from pathlib import Path
from setuptools import setup, find_packages

AUTHORS = ["Jake Alexander"]

CLASSIFIERS = [
    "Topic :: Software Development :: Libraries",
    "Programming Language:: Python:: 3.9"
]

PACKAGES = find_packages(exclude="tests")

INSTALL_REQUIRES = [
    "aiofiles==0.6.0",
    "biopython==1.78",
    "virtool-core==0.1.1",
    "virtool-workflow==0.4.0"
]

setup(
    name="vt-workflow-aodp",
    version="0.1.0",
    description="A workflow for using amplicon analysis in Virtool using AODP.",
    long_description=Path("README.md").read_text(),
    long_description_context_type="text/markdown",
    url="https://github.com/virtool/workflow-aodp",
    author=", ".join(AUTHORS),
    license="MIT",
    platforms="linux",
    packages=PACKAGES,
    install_requires=INSTALL_REQUIRES,
    python_requires=">=3.9",
)