import logging
from pathlib import Path

from virtool_workflow import fixture

logger = logging.getLogger(__name__)


@fixture
def joined_path(work_path: Path) -> Path:
    """
    The path to an intermediate file containing joined output from FLASH.

    """
    return work_path / "flash.extendedFrags.fastq"


@fixture
def unique_path(work_path: Path) -> Path:
    """
    The path to an intermediate file containing all unique, joined reads.

    """
    return work_path / "unique.fa"
