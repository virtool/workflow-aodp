import collections
import os
import pathlib
import shutil
import typing

import aiofiles
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import virtool_core.utils
from virtool_workflow import startup, step, cleanup

AODP_MAX_HOMOLOGY = 0
AODP_OLIGO_SIZE = 8


@startup
async def prepare_index(params, run_in_executor):
    """
    Copy reference index from its location in the Virtool application data path to a local temporary path.

    """

    await run_in_executor(
        os.makedirs,
        pathlib.Path(params["temp_index_path"]).parent
    )

    await run_in_executor(
        shutil.copy,
        params["index_path"],
        params["temp_index_path"]
    )


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
        "join_histogram": await parse_flash_histogram(hist_path),
        "joined_pair_count": await virtool_core.utils.file_length(joined_path) / 4,
        "remainder_pair_count": await virtool_core.utils.file_length(remainder_path) / 4
    })


@step
async def deduplicate_reads():
    pass


@step
async def aodp():
    pass


@step
async def import_results():
    pass


async def parse_flash_histogram():
    pass


async def parse_joined_fastq():
    pass


async def run_deduplication():
    pass


@cleanup
async def delete_analysis():
    pass


@cleanup
async def delete_analysis():
    pass
