from collections import defaultdict
from pathlib import Path
from typing import Generator, List, Mapping

import aiofiles
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


async def parse_flash_histogram(path: Path) -> List[List[int]]:
    """
    Parse the histogram output file from FLASH.

    :param path: the path to the histogram file
    :return: a list-based representation of the histogram data.

    """
    histogram = list()

    async with aiofiles.open(path, "r") as f:
        async for line in f:
            histogram.append([int(i) for i in line.rstrip().split()])

    return histogram


def parse_joined_fastq(path: Path, counts: Mapping[str, int]) -> Generator[SeqRecord, None, None]:
    """
    Parse the joined FASTQ file at `path` and yield Biopython `SeqRecord` objects. Does not yield
    duplicate reads, de-duplicating the input.

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


async def run_deduplication(joined_path: Path, output_path: Path):
    """
    Deduplicate the reads at `joined_path` and output at `output_path`. Tis function is computationally intensive and
    should be executed in a separate process.

    :param joined_path: the path to the file containing the reads to be joined
    :param output_path: the path to a file to write the deduplicated reads to
    :return: the sequence-wise duplicate counts

    """
    counts = defaultdict(int)

    with open(output_path, "w") as f:
        for record in parse_joined_fastq(joined_path, counts):
            SeqIO.write(record, f, format="fasta")

    return dict(counts)
