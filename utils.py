from copy import deepcopy
from virtool_core.history.db import patch_to_version


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
        _, patched, _ = await patch_to_version(
            app_dict,
            otu_id,
            otu_version
        )

        for isolate in patched["isolates"]:
            for sequence in isolate["sequences"]:
                sequence_id = sequence["_id"]
                sequence_otu_map[sequence_id] = patched["_id"]

        return sequence_otu_map
