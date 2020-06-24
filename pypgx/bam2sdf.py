import pysam
import sys

from .common import is_namespace, get_parser, is_file, read_gene_table, logging, parse_region, sort_regions

def _bam2sdf(tg, cg, bam):
    logger = logging.getLogger(__name__)

    genes = read_gene_table()
    targets = [k for k, v in genes.items() if v["type"] == "target"]
    controls = [k for k, v in genes.items() if v["control"] == "yes"]
    if tg not in targets:
        raise ValueError(f"'{tg}' is not among target genes: {targets}")
    if cg not in controls:
        raise ValueError(f"'{cg}' is not among control genes: {controls}")

    # Get sample and sequence names from BAM headers.
    sm = []
    sn = []
    for x in bam:
        has_sm = False
        result = pysam.view("-H", x).strip().split("\n")
        for line in result:
            fields = line.split("\t")
            if "@RG" == fields[0] and not has_sm:
                for field in fields:
                    if "SM:" in field:
                        y = field.replace("SM:", "")
                        sm.append(y)
                        has_sm = True
            if "@SQ" == fields[0]:
                for field in fields:
                    if "SN:" in field:
                        y = field.replace("SN:", "")
                        if y not in sn:
                            sn.append(y)
        if not has_sm:
            raise ValueError(f"SM tag not found from BAM file: {x}")

    logger.info(f"Sample IDs: {sm}")
    logger.info(f"Contigs: {sn}")

    # Determine whether the "chr" string should be used.
    if any(["chr" in x for x in sn]):
        chr = "chr"
    else:
        chr = ""

    tr = genes[tg]["hg19_region"].replace("chr", "")
    cr = genes[cg]["hg19_region"].replace("chr", "")

    regions = sort_regions([tr, cr])

    result = ""

    for region in regions:
        temp = pysam.depth("-a", "-Q", "1", "-r", f"{chr}{region}", *bam)
        result += temp

    return result

def bam2sdf(*args):
    is_console = is_namespace(args[0])

    if is_console:
        tg = args[0].tg
        cg = args[0].cg
        bam = args[0].bam
    else:
        parser = get_parser()
        args = parser.parse_args(args)
        tg = args.tg
        cg = args.cg
        bam = args.bam

    result = _bam2sdf(tg, cg, bam)

    if is_console:
        if args[0].o:
            with open(args[0].o, "w") as f:
                f.write(result)
        else:
            sys.stdout.write(result)
    else:
        return result