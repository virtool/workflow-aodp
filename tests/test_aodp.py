import pytest
from virtool_workflow.fixtures.scope import WorkflowFixtureScope

import workflow


@pytest.mark.asyncio
async def test_make_analysis_dir(tmpdir):
    with WorkflowFixtureScope as scope:
        scope["params"] = {
            "temp_analysis_path": f"{tmpdir}/foo1",
            "raw_path": f"{tmpdir}/foo2",
            "reads_path": f"{tmpdir}/foo3"
        }

        bound_function = await scope.bind(workflow.make_analysis_dir)
        await bound_function()

        assert scope["params"]["temp_analysis_path"].exists()
        assert scope["params"]["raw_path"].exists()
        assert scope["params"]["reads_path"].exists()
