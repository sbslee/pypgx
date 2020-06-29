import logging
import io
import pkgutil
import gzip
from typing import Optional, TextIO

import pysam

LINE_BREAK1 = "-" * 70
LINE_BREAK2 = "*" * 70

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def str2file(x):
    return io.StringIO(x)

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

def parse_region(x):
    return (
        x.split(":")[0],
        int(x.split(":")[1].split("-")[0]),
        int(x.split(":")[1].split("-")[1]),
    )

def sort_regions(x):
    def f(x):
        r = parse_region(x)
        if "X" in r[0]:
            chr = 23
        elif "Y" in r[0]:
            chr = 24
        else:
            chr = int(r[0].replace("chr", ""))
        return (chr, r[1], r[2])
    return sorted(x, key = f)

def sm_tag(x):
    l = []
    header = pysam.view("-H", x).strip().split("\n")
    for line in header:
        fields = line.split("\t")
        if "@RG" == fields[0]:
            for field in fields:
                if "SM:" in field:
                    y = field.replace("SM:", "")
                    l.append(y)
    l = list(set(l))
    if len(l) == 0:
        raise ValueError(f"SM tag not found: {x}")
    elif len(l) > 1:
        logger.warning(f"Multiple SM tags found (will use the 1st one): {x}")
        return l[0]
    else:
        return l[0]

class VCFFile:
    def __init__(
        self,
        fn: str,
        f: Optional[TextIO] = None,
    ) -> None:
        """
        Initialize VCF file.

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

    def read(self, region = None):
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
        string = ""
        for line in self.meta:
            string += line
        string += "\t".join(self.header) + "\n"
        for fields in self.data:
            string += "\t".join(fields) + "\n"
        return string

    def to_file(self, fn: str) -> None:
        string = self.to_str()
        with open(fn, "w") as f:
            f.write(string)

    def unphase(self) -> None:
        for i in range(len(self.data)):
            self.data[i][9:] = [x.replace("|", "/") for x in self.data[i][9:]]

    def phase(self) -> None:
        for i in range(len(self.data)):
            self.data[i][9:] = [x.replace("/", "|") for x in self.data[i][9:]]

    def close(self):
        self.f.close()
        self.f = None
