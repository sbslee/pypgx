import logging
import io
import pkgutil
import gzip
from typing import Optional, TextIO, List

import pysam

LINE_BREAK1 = "-" * 70
LINE_BREAK2 = "*" * 70

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def str2file(x):
    return io.StringIO(x)

def sort_star_names(names):
    def f(x):
        cn = 1
        if '*' not in x or x == '*DEL':
            n = 999
        else:
            n = int(''.join([y for y in x.split('+')[0].split('x')[0] if y.isdigit()]))
            if 'x' in x.split('+')[0]:
                cn = int(x.split('+')[0].split('x')[1])
        return (n, cn, len(x))

    return sorted(names, key = f)

def read_gene_table():
    result = {}
    text = pkgutil.get_data(__name__, "resources/sg/gene_table.txt").decode()
    for line in text.strip().split("\n"):
        fields = line.split("\t")
        gene = fields[1]
        if gene == "name":
            header = fields
            continue
        result[gene] = dict(zip(header, fields))
    return result

def read_pt_table():
    result = {}
    text = pkgutil.get_data(__name__, "resources/sg/pt_table.txt").decode()
    for line in text.strip().split("\n"):
        fields = line.split("\t")
        gene = fields[0]
        name = fields[2]
        rules = fields[3]
        if gene == "gene":
            header = fields
            continue
        if gene not in result:
            result[gene] = {}
        result[gene][name] = dict(zip(header, fields))
    return result

def parse_region(region: str) -> List[str]:
    """
    Parse region.

    Returns:
        list[str]: Parsed region [chr, start, end].

    Args:
        region (str): Region.
    """

    return (
        region.split(":")[0],
        int(region.split(":")[1].split("-")[0]),
        int(region.split(":")[1].split("-")[1]),
    )

def sort_regions(regions: List[str]) -> List[str]:
    """
    Sort regions.

    Returns:
        list[str]: Sorted regions.

    Args:
        regions (list[str]): Regions.
    """

    def f(x):
        r = parse_region(x)
        if "X" in r[0]:
            chr = 23
        elif "Y" in r[0]:
            chr = 24
        else:
            chr = int(r[0].replace("chr", ""))
        return (chr, r[1], r[2])
    return sorted(regions, key = f)

def sm_tag(bam: str) -> str:
    """
    Extract SM tag from BAM file.

    Returns:
        str: SM tag.

    Args:
        bam (str): BAM file.
    """

    header = pysam.view("-H", bam).strip().split("\n")

    l = []

    for line in header:
        fields = line.split("\t")
        if "@RG" == fields[0]:
            for field in fields:
                if "SM:" in field:
                    l.append(field.replace("SM:", ""))

    l = list(set(l))

    if not l:
        raise ValueError(f"SM tag not found: {bam}")

    if len(l) > 1:
        logger.warning(
            f"Multiple SM tags found (will return the first one): {bam}")
        result = l[0]
    else:
        result = l[0]

    return result




def is_chr(bam: str) -> bool:
    """
    Check whether SN tags in BAM file contain "chr" string.

    Returns:
        bool: True if found.

    Args:
        bam (str): BAM file.
    """

    header = pysam.view("-H", bam).strip().split("\n")

    l = []

    for line in header:
        fields = line.split("\t")
        if "@SQ" == fields[0]:
            for field in fields:
                if "SN:" in field:
                    l.append(field.replace("SN:", ""))

    return any(["chr" in x for x in l])

class VCFFile:
    """
    VCF file ojbect.

    Attributes:
       f (TextIO): VCF file.
       meta (list[str]): Meta data.
       header (list[str]): Header.
       data (list[str]): Genotype data.

    Examples:

        >>> vcf = VCFFile("in.vcf") # also works with "in.vcf.gz"
        >>> vcf.read("chr10:96519437-96615962") # read CYP2C19 region only
        >>> vcf.unphase()
        >>> result = vcf.to_str()
        >>> vcf.to_file("out.vcf")
        >>> vcf.close()
    """

    def __init__(self, fn: str, f: Optional[TextIO] = None) -> None:
        """
        Initialize VCF file object.

        Args:
            fn (str): VCF file.
            f (TextIO, optional): VCF file.
        """

        if fn:
            if ".gz" in fn:
                self.f = gzip.open(fn, "rt")
            else:
                self.f = open(fn)
        else:
            self.f = f
        self.meta = []
        self.header = []
        self.data = []

    def read(self, region: Optional[str] = None) -> None:
        """
        Read VCF file.

        Args:
            region (str, optional): Target region.
        """

        if region:
            r = parse_region(region)
            for line in self.f:
                if "##" in line:
                    self.meta.append(line)
                    continue
                fields = line.strip().split("\t")
                if fields[0] == "#CHROM":
                    self.header = fields
                    continue
                chr = fields[0]
                pos = int(fields[1])
                if chr != r[0] or pos < r[1]:
                    continue
                if pos > r[2]:
                    break
                self.data.append(fields)
        else:
            for line in self.f:
                if "##" in line:
                    self.meta.append(line)
                    continue
                fields = line.strip().split("\t")
                if fields[0] == "#CHROM":
                    self.header = fields
                    continue
                self.data.append(fields)

    def to_str(self) -> str:
        """
        Return VCF file.

        Returns:
            str: VCF file.
        """

        string = ""
        for line in self.meta:
            string += line
        string += "\t".join(self.header) + "\n"
        for fields in self.data:
            string += "\t".join(fields) + "\n"
        return string

    def to_file(self, fn: str) -> None:
        """
        Write VCF file.

        Args:
            fn (str): VCF file.
        """

        string = self.to_str()
        with open(fn, "w") as f:
            f.write(string)

    def unphase(self) -> None:
        """
        Change genotype separator from '|' to '/'.
        """

        for i in range(len(self.data)):
            self.data[i][9:] = [x.replace("|", "/") for x in self.data[i][9:]]

    def phase(self) -> None:
        """
        Change genotype separator from '/' to '|'.

        .. warning::
            This is not statistcal phasing.
        """

        for i in range(len(self.data)):
            self.data[i][9:] = [x.replace("/", "|") for x in self.data[i][9:]]

    def close(self):
        """
        Close VCF file.
        """

        self.f.close()
        self.f = None
