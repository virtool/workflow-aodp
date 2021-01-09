import logging

from virtool_workflow import fixture
from virtool_workflow_runtime.config.configuration import db_name, db_connection_string
from virtool_workflow_runtime.db import VirtoolDatabase

from utils import get_index_info

logger = logging.getLogger(__name__)


@fixture
def db():
    return VirtoolDatabase(db_name(), db_connection_string())


@fixture
async def params(db, job_args, analysis_path, index_path, sample_path, temp_analysis_path, data_path, reads_path,
                 raw_path, temp_cache_path, temp_path_str):
    logger.debug("Retrieving job parameters")

    # The document for the sample being analyzed. Assigned after database connection is made.
    sample = await db.samples.find_one(job_args["sample_id"])

    analysis = await db.analyses.find_one(job_args["analysis_id"], ["subtraction"])

    params = dict(job_args)

    params.update({
        # The Path to the directory where all analysis result files will be written.
        "analysis_path": analysis_path,
        "index_path": index_path,
        "sample_path": sample_path,
        "paired": sample["paired"],
        # The number of reads in the sample library. Assigned after database connection is made.
        "read_count": int(sample["quality"]["count"]),
        "sample_read_length": int(sample["quality"]["length"][1]),
        "library_type": sample["library_type"],
        "reads_path": reads_path,
        "subtraction_path":
            data_path /
            "subtractions" /
            analysis["subtraction"]["id"].replace(" ", "_").lower() /
            "reference",
        "raw_path": raw_path,
        "temp_cache_path": temp_cache_path,
        "temp_analysis_path": temp_analysis_path
    })

    index_info = await get_index_info(
        db,
        params["settings"],
        params["index_id"]
    )

    params.update(index_info)

    read_paths = [params["reads_path"] / "reads_1.fq.gz"]

    if params["paired"]:
        read_paths.append(params["reads_path"] / "reads_2.fq.gz")

    params["read_paths"] = read_paths
    params["temp_index_path"] = temp_path_str / "reference" / "reference.fa"
    params["aodp_output_path"] = params["temp_analysis_path"] / "aodp.out"

    return params
