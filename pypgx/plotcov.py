import os
from subprocess import run, PIPE
from .common import logging

logger = logging.getLogger(__name__)

def plotcov(sdf: str, out: str) -> None:
    """
    Plot coverage data.

    Args:
        sdf (str): SDF file.
        out (str): PDF file.
    """
    command = [
        "Rscript",
        f"{os.path.dirname(__file__)}/resources/r/plotcov.R",
        sdf, 
        out,
    ]
    output = run(command, stdout=PIPE, stderr=PIPE)

    logger.info("plotcov:")
    for line in output.stderr.decode("utf-8").split("\n"):
        logger.info(f"    {line}")
