import statistics
import sys

from .common import is_namespace, get_parser, is_file

def _sdf2gdf(sdf, id):
    # Check the sample count.
    fields = next(sdf).strip().split("\t")
    if len(fields) != len(id) + 2:
        raise ValueError("incorrect sample count") 

    # Get the header.
    result = "Locus\tTotal_Depth\tAverage_Depth_sample"
    for x in id:
        result += f"\tDepth_for_{x}"
    result += "\n"

    # Get the data.
    for line in sdf:
        fields1 = line.strip().split("\t")
        locus = f"{fields1[0]}:{fields1[1]}"
        depth = [int(x) for x in fields1[2:]]
        avg = round(statistics.mean(depth), 2)
        fields2 = [locus, sum(depth), avg] + depth
        result += "\t".join([str(x) for x in fields2]) + "\n"

    return result

def sdf2gdf(*args):
    is_console = is_namespace(args[0])

    if is_console:
        sdf = open(args[0].sdf)
        id = args[0].id
    elif is_file(args[0]):
        sdf = args[0]
        id = args[1]
    else:
        parser = get_parser()
        args = parser.parse_args(args)
        sdf = open(args.sdf)
        id = args.id

    result = _sdf2gdf(sdf, id)
    
    sdf.close()

    if is_console:
        if args[0].o:
            with open(args[0].o, "w") as f:
                f.write(result)
        else:
            sys.stdout.write(result)
    else:
        return result
