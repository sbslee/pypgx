import statistics

from typing import Optional, List, TextIO

def sdf2gdf(
        fn: str,
        id: List[str],
        f: Optional[TextIO] = None, 
    ) -> str:
    """
    Create GDF file from SDF file.

    Returns:
        str: GDF file.

    Args:
        fn (str): SDF file.
        id (list[str]): Sample ID(s).
        f (TextIO, optional): SDF file.
    """

    if fn:
        f = open(fn)

    # Get the header.
    result = "Locus\tTotal_Depth\tAverage_Depth_sample"
    for x in id:
        result += f"\tDepth_for_{x}"
    result += "\n"

    # Check the sample count with the first line.
    fields1 = next(f).strip().split("\t")
    if len(fields1) - 2 != len(id):
        raise ValueError("incorrect sample count")
    locus = f"{fields1[0]}:{fields1[1]}"
    depth = [int(x) for x in fields1[2:]]
    avg = round(statistics.mean(depth), 2)
    fields2 = [locus, sum(depth), avg] + depth
    result += "\t".join([str(x) for x in fields2]) + "\n"

    # Read rest of the lines.
    for line in f:
        fields1 = line.strip().split("\t")
        locus = f"{fields1[0]}:{fields1[1]}"
        depth = [int(x) for x in fields1[2:]]
        avg = round(statistics.mean(depth), 2)
        fields2 = [locus, sum(depth), avg] + depth
        result += "\t".join([str(x) for x in fields2]) + "\n"

    if fn:
        f.close()

    return result
