import collections
import pathlib
from typing import Generator

import aiofiles
import virtool_core.caches
import virtool_core.db
import virtool_core.history
import virtool_core.samples
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


async def get_index_info(db, settings, index_id) -> dict:
    document = await db.indexes.find_one(index_id, ["manifest", "sequence_otu_map"])

    try:
        sequence_otu_map = document["sequence_otu_map"]
    except KeyError:
        sequence_otu_map = await get_sequence_otu_map(db, settings, document["manifest"])

    return {
        "manifest": document["manifest"],
        "sequence_otu_map": sequence_otu_map
    }


async def get_sequence_otu_map(db, settings, manifest) -> dict:
    app_dict = {
        "db": db,
        "settings": settings
    }

    sequence_otu_map = dict()

    for otu_id, otu_version in manifest.items():
        _, patched, _ = await virtool_core.history.db.patch_to_version(
            app_dict,
            otu_id,
            otu_version
        )

        for isolate in patched["isolates"]:
            for sequence in isolate["sequences"]:
                sequence_id = sequence["_id"]
                sequence_otu_map[sequence_id] = patched["_id"]

        return sequence_otu_map


async def parse_flash_histogram(path: pathlib.Path):
    """
    Parse the histogram output file from FLASH.

    :param path: the path to the histogram file
    :return: a list-bsed representation of the histogram data.
    """
    histogram = list()

    async with aiofiles.open(path, "r") as f:
        async for line in f:
            histogram.append([int(i) for i in line.rstrip().split()])

    return histogram


async def parse_joined_fastq(path: pathlib.Path, counts: collections.defaultdict) -> Generator[SeqRecord, None, None]:
    """
    Parse the joined FASTQ file at `path` and yield Biopython `SeqRecord` objects. Does not yield duplicate reads,
    de-duplicating the input.

    Updates `counts` with observed duplicate count for each sequence.

    :param path: the path to the input FASTQ file
    :param counts: a dict to track the duplication count for each read
    """
    sequence_id_map = dict()

    for record in SeqIO.parse(path, format="fastq"):
        try:
            sequence_id = sequence_id_map[str(record.seq)]
        except KeyError:
            sequence_id = f"read_len_{len(sequence_id_map) + 1}"
            sequence_id_map[str(record.seq)] = sequence_id

            yield SeqRecord(record.seq, id=sequence_id)

        counts[sequence_id] += 1


async def run_deduplication(joined_path: pathlib.Path, output_path: pathlib.Path):
    """
    Deduplicate the reads at `joined_path` and output at `output_path`. Tis function is computationally intensive and
    should be executed in a separate process.

    :param joined_path: the path to the file containing the reads to be joined
    :param output_path: the path to a file to write the deduplicated reads to
    :return: the sequence-wise duplicate counts
    """
    counts = collections.defaultdict(int)

    with open(output_path, "w") as f:
        for record in parse_joined_fastq(joined_path, counts):
            SeqIO.write(record, f, format="fasta")

    return dict(counts)
