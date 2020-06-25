from .bam2sdf import bam2sdf
from .sdf2gdf import sdf2gdf
from .common import sm_tag, str2file, is_namespace, get_parser, stdout

def _bam2gdf(tg, cg, bam):
    sdf = bam2sdf("bam2sdf", tg, cg, *bam)
    sm = [sm_tag(x) for x in bam]
    result = sdf2gdf(str2file(sdf), sm)
    return result

def bam2gdf(*args):
    is_console = is_namespace(args[0])

    if is_console:
        result = _bam2gdf(args[0].tg, args[0].cg, args[0].bam)
    else:
        parser = get_parser()
        args = parser.parse_args(args)
        result = _bam2gdf(args.tg, args.cg, args.bam)

    if is_console:
        if args[0].o:
            with open(args[0].o, "w") as f:
                f.write(result)
        else:
            stdout(result)
    else:
        return result