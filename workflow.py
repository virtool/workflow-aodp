import logging
from pathlib import Path
from typing import List

import aiofiles
import virtool_core.samples.db
import virtool_core.utils
from virtool_workflow import step
from virtool_workflow.analysis.indexes import Index

import utils

AODP_MAX_HOMOLOGY = 0
AODP_OLIGO_SIZE = 8

logger = logging.getLogger(__name__)


@step
async def join_reads(
        joined_path: Path,
        proc: int,
        run_subprocess,
        results: dict,
        sample,
        work_path: Path
):
    """
    Join overlapping paired reads into single reads.

    """
    max_overlap = round(0.65 * sample.read_length)

    command = [
        "flash",
        "--max-overlap", str(max_overlap),
        "-d", work_path,
        "-o", "flash",
        "-t", str(proc - 1),
        sample.read_paths
    ]

    await run_subprocess(command)

    hist_path = work_path / "flash.hist"
    remainder_path = work_path / "flash.notCombined_1.fastq"

    results.update({
        "join_histogram": await utils.parse_flash_histogram(hist_path),
        "joined_pair_count": await virtool_core.utils.file_length(joined_path) / 4,
        "remainder_pair_count": await virtool_core.utils.file_length(remainder_path) / 4
    })


@step
async def deduplicate_reads(
        run_in_executor,
        joined_path: Path,
        results: dict,
        unique_path: Path,
        work_path: Path
):
    """Remove duplicate reads. Store the counts for unique reads."""
    counts = await run_in_executor(
        utils.run_deduplication,
        joined_path,
        unique_path
    )

    results["sequence_counts"] = counts


@step
async def aodp(
        proc: int,
        indexes: List[Index],
        run_subprocess,
        results: dict,
        unique_path: Path,
        work_path: Path
):
    """
    Run AODP, parse the output, and update the.

    TODO: Upload data and finalize analysis

    """
    aodp_path = work_path / "aodp.out"
    base_name = work_path / "aodp"
    index_path = indexes[0].fasta_path

    command = [
        "aodp",
        f"--basename={base_name}",
        f"--threads={proc}",
        f"--oligo-size={AODP_OLIGO_SIZE}",
        f"--match={unique_path}",
        f"--match-output={aodp_path}",
        f"--max-homolo={AODP_MAX_HOMOLOGY}",
        index_path
    ]

    await run_subprocess(command, cwd=work_path)

    parsed = list()

    async with aiofiles.open(aodp_path, "r") as f:
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

            otu_id = indexes[0].get_otu_id_by_sequence_id(sequence_id)
            otu_version = indexes[0].get_otu_version_by_sequence_id(sequence_id)

            parsed.append({
                "id": read_id,
                "sequence_id": sequence_id,
                "identity": identity,
                "matched_length": int(split[3]),
                "read_length": int(split[4]),
                "min_cluster": int(split[5]),
                "max_cluster": int(split[6]),
                "count": results["sequence_counts"][read_id],
                "otu": {
                    "version": otu_version,
                    "id": otu_id
                }
            })

    results["hits"] = parsed
