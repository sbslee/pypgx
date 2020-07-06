import os
from subprocess import run, PIPE
from .common import logging

logger = logging.getLogger(__name__)

def cpa(rdata: str) -> str:
    """
    Run change point analysis in copy number.

    Returns:
        str: Result file.

    Args:
        rdata (str): RData file.
    """

    command = [
        "Rscript", f"{os.path.dirname(__file__)}/resources/r/cpa.R", rdata]
    output = run(command, stdout=PIPE, stderr=PIPE)

    logger.info("changepoint:")
    for line in output.stderr.decode("utf-8").split("\n"):
        logger.info(f"    {line}")

    result = output.stdout.decode("utf-8")

    return result
