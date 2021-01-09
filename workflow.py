import logging
import os

import aiofiles
import virtool_core.samples.db
import virtool_core.utils
import virtool_workflow.analysis.cache
from virtool_workflow import startup, step, cleanup

import utils

AODP_MAX_HOMOLOGY = 0
AODP_OLIGO_SIZE = 8

logger = logging.getLogger(__name__)


@startup
async def make_analysis_dir(params, run_in_executor):
    """
    Make a directory for the analysis in the sample/analysis directory.

    """
    await run_in_executor(os.mkdir, params["temp_analysis_path"])
    await run_in_executor(os.mkdir, params["raw_path"])
    await run_in_executor(os.mkdir, params["reads_path"])


@step
async def join_reads(params, proc, run_subprocess, results):
    """
    Join overlapping paired reads into single reads.

    """
    max_overlap = round(0.65 * params["sample_read_length"])

    command = [
        "flash",
        "--max-overlap", str(max_overlap),
        "-d", params["temp_analysis_path"],
        "-o", "flash",
        "-t", str(proc - 1),
        *params["read_paths"]
    ]

    await run_subprocess(command)

    joined_path = params["temp_analysis_path"] / "flash.extendedFrags.fastq"
    remainder_path = params["temp_analysis_path"] / "flash.notCombined_1.fastq"
    hist_path = params["temp_analysis_path"] / "flash.hist"

    results.update({
        "join_histogram": await utils.parse_flash_histogram(hist_path),
        "joined_pair_count": await virtool_core.utils.file_length(joined_path) / 4,
        "remainder_pair_count": await virtool_core.utils.file_length(remainder_path) / 4
    })


@step
async def deduplicate_reads(params, run_in_executor):
    """
    Remove duplicate reads. Store the counts for unique reads.

    """
    joined_path = params["temp_analysis_path"] / "flash.extendedFrags.fast"
    output_path = params["temp_analysis_path"] / "unique.fa"

    counts = await run_in_executor(
        utils.run_deduplication,
        joined_path,
        output_path
    )

    params["intermediate"]["sequence_counts"] = counts


@step
async def aodp(params, proc, run_subprocess, results):
    """
    Run AODP, parse the ouput, and update the .

    """
    cwd = params["temp_analysis_path"]

    aodp_output_path = params["aodp_output_path"]
    base_name = params["temp_analysis_path"] / "aodp"
    target_path = params["temp_analysis_path"] / "unique.fa"

    command = [
        "aodp",
        f"--basename={base_name}",
        f"--threads={proc}",
        f"--oligo-size={AODP_OLIGO_SIZE}",
        f"--match={target_path}",
        f"--match-output={aodp_output_path}",
        f"--max-homolo={AODP_MAX_HOMOLOGY}",
        params["temp_index_path"]
    ]

    await run_subprocess(command, cwd=cwd)

    parsed = list()

    async with aiofiles.open(params["aodp_output_path"], "r") as f:
        async for line in f:
            split = line.rstrip().split("\t")
            assert len(split) == 7

            sequence_id = split[1]

            if sequence_id == "-":
                continue

            identity = split[2]

            if identity[0] == "<":
                continue
            else:
                identity = float(identity.replace("%", ""))

            read_id = split[0]

            sequence_id = split[1]

            otu_id = params["sequence_otu_map"][sequence_id]
            otu_version = params["manifest"][otu_id]

            parsed.append({
                "id": read_id,
                "sequence_id": sequence_id,
                "identity": identity,
                "matched_length": int(split[3]),
                "read_length": int(split[4]),
                "min_cluster": int(split[5]),
                "max_cluster": int(split[6]),
                "count": params["intermediate"]["sequence_counts"][read_id],
                "otu": {
                    "version": otu_version,
                    "id": otu_id
                }
            })

    results["results"] = parsed


@step
async def import_results(db, params, results):
    """
    Import the analysis results to the database.

    """
    analysis_id = params["analysis_id"]
    sample_id = params["sample_id"]

    # Update the database document with the small data.
    await db.analyses.update_one({"_id": analysis_id}, {
        "$set": {
            **results,
            "ready": True
        }
    })

    await virtool_core.samples.db.recalculate_workflow_tags(db, sample_id)


@cleanup
async def delete_analysis(params):
    await virtool_workflow.analysis.cache.delete_analysis(
        params["analysis_id"],
        params["analysis_path"],
        params["sample_id"]
    )
