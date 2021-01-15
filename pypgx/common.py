import os
import logging
import random
import string
import configparser
from typing import Dict, List, Optional
from tempfile import TemporaryDirectory
import pysam
from functools import wraps

import copy
import gzip
import statistics

from .sglib import (
    read_gene_table,
    read_snp_table,
    build_snpdb,
    read_star_table,
    build_stardb,
    StarAllele,
)

LINE_BREAK1 = "-" * 70
LINE_BREAK2 = "*" * 70






class Record:
    """Stores data for single SNP record.

    Args:
        fields: SNP data in columns.

    Attributes:
        fields (list of str): SNP data in columns.
    """

    def __init__(self, fields: List[str]):
        """Inits a Record."""
        self.fields = fields

    @property
    def chrom(self):
        """str: Chromosome."""
        return self.fields[0]

    @chrom.setter
    def chrom(self, x: str):
        self.fields[0] = x

    @property
    def pos(self):
        """int: Position."""
        return int(self.fields[1])

    @pos.setter
    def pos(self, x: int):
        self.fields[1] = str(x)

    @property
    def id(self):
        """str: Identifier."""
        return self.fields[2]

    @id.setter
    def id(self, x: str):
        self.fields[2] = x

    @property
    def ref(self):
        """str: Reference allele."""
        return self.fields[3]

    @ref.setter
    def ref(self, x: str):
        self.fields[3] = x

    @property
    def alt(self):
        """list of str: Alternate allele(s)."""
        return self.fields[4].split(",")

    @alt.setter
    def alt(self, x: List[str]):
        self.fields[4] = ",".join(x)

    @property
    def qual(self):
        """int: Quality. 0 is equivalent to the missing value."""
        s = self.fields[5]
        if s == ".":
            return 0
        else:
            return int(s)

    @qual.setter
    def qual(self, x: int):
        self.fields[5] = str(x) if x else "."

    @property
    def filter(self):
        """list of str: Filter status. An empty list is equivalent to the
        missing value.
        """
        s = self.fields[6]
        if s == ".":
            return []
        else:
            return s.split(";")

    @filter.setter
    def filter(self, x: List[str]):
        self.fields[6] = ";".join(x) if x else "."

    @property
    def info(self):
        """dict of str to str: Additional information. An empty dictionary is
        equivalent to the missing value.
        """
        s = self.fields[7]
        if s == ".":
            return {}
        d = {}
        for x in s.split(";"):
            k = x.split("=")[0]
            v = x.split("=")[1]
            d[k] = v
        return d

    @info.setter
    def info(self, x: Dict[str, str]):
        if x:
            l = []
            for k, v in x.items():
                l.append(f"{k}={v}")
            self.fields[7] = ";".join(l)
        else:
            self.fields[7] = "."

    @property
    def format(self):
        """list of str: Types of genotype data."""
        return self.fields[8].split(":")

    @format.setter
    def format(self, x: List[str]):
        self.fields[8] = ":".join(x)

    @property
    def data(self):
        """list of str: Genotype data."""
        return self.fields[9:]

    @data.setter
    def data(self, x: List[str]):
        self.fields[9:] = x

    def has_indel(self) -> bool:
        """Returns true if the REF or ALT allele is an indel."""
        return len(self.ref) > 1 or any([len(x) > 1 for x in self.alt])





class VCFFile:
    """Versatile object for reading, writing and manipulating a VCF file.

    Attributes:
        meta (list of str): Meta information lines.
        header (list of str): Column names.
        data (list of Record): SNP records.
    """

    def __init__(self, *args, **kwargs):
        """Inits a VCFFile."""
        self.args   = args
        self.kwargs = kwargs
        self.meta   = []
        self.header = []
        self.data   = []

    @property
    def filepath(self):
        """str: VCF file path."""
        if "filepath" in self.kwargs:
            return self.kwargs["filepath"]
        elif self.args:
            return self.args[0]
        else:
            return ""

    def clear(self) -> None:
        """Clears any loaded VCF data."""
        self.meta.clear()
        self.header.clear()
        self.data.clear()

    def read(self,
             region: Optional[str] = None,
             tidy: bool = False) -> None:
        """Reads the provided VCF file.

        Args:
            region: If provided, only read VCF data in the region.
            tidy: If true, remove ``chr`` string from ``CHROM`` column.
        """

        self.clear()

        if region:
            chr = region.split(":")[0].replace("chr", "")
            start = int(region.split(":")[1].split("-")[0])
            end = int(region.split(":")[1].split("-")[1])

        if self.filepath.endswith(".gz"):
            f = gzip.open(self.filepath, "rt")
        else:
            f = open(self.filepath)

        for line in f:
            if "##" in line:
                self.meta.append(line)
                continue

            fields = line.strip().split("\t")

            if fields[0] == "#CHROM":
                self.header = fields
                continue

            record = Record(fields)

            if region:
                if (chr == record.chrom.replace("chr", "") and
                    start <= record.pos <= end):
                    pass
                elif (chr == record.chrom.replace("chr", "") and
                      end < record.pos):
                      break
                else:
                    continue

            if tidy:
                record.fields[0] = record.fields[0].replace("chr", "")

            self.data.append(record)

        f.close()

    def missing_filter(self,
                       threshold: float = 0.0) -> List[Record]:
        """Filters out records with missing genotype ``./.``.

        Returns:
            list of Record: Filtered records.

        Args:
            threshold: A record will be removed if the fraction of samples
                with missing genotype is greater than this threshold.
        """

        filtered = []

        for i in reversed(range(len(self.data))):
            record = self.data[i]
            missing = ["." in x.split(":")[0] for x in record.data]
            missing = missing.count(True) / len(missing)

            if missing > threshold:
                if not record.filter or "PASS" in record.filter:
                    record.filter = ["MissingFilter"]
                else:
                    record.filter += ["MissingFilter"]

                filtered.append(self.data.pop(i))

        return filtered

    def multiallelic_filter(self) -> List[Record]:
        """Filters out records with more than two ALT alleles.

        Returns:
            list of Record: Filtered records.
        """

        filtered = []

        for i in reversed(range(len(self.data))):
            record = self.data[i]

            if len(record.alt) > 1:
                if not record.filter or "PASS" in record.filter:
                    record.filter = ["MultiallelicFilter"]
                else:
                    record.filter += ["MultiallelicFilter"]

                filtered.append(self.data.pop(i))

        return filtered

    def invalid_allele_filter(self) -> List[Record]:
        """Filters out records with invalid alleles (e.g. ``I``, ``D``).

        Returns:
            list of Record: Filtered records.
        """

        filtered = []

        for i in reversed(range(len(self.data))):
            record = self.data[i]

            if (record.ref in ["I", "D"] or
                "." in record.alt or "D" in record.alt):

                if not record.filter or "PASS" in record.filter:
                    record.filter = ["InvalidAlleleFilter"]
                else:
                    record.filter += ["InvalidAlleleFilter"]

                filtered.append(self.data.pop(i))

        return filtered

    def allelic_imbalance_filter(self,
                                 threshold: float = 0.8,
                                 count: int = 3) -> List[Record]:
        """Filters out records with high allelic imbalance.

        Returns:
            list of Record: Filtered records.

        Args:
            threshold: A record is said to have high allelic imbalance if the
                median of allele fractions is greather than this threshold.
            count: Records will not be filtered out if they are indel and
                have the number of supporting samples less than this count.
        """

        def get_af(x: str, record: Record) -> float:
            i = record.format.index("GT")
            gt = x.split(":")[i]

            if "|" in gt:
                gt = gt.split("|")
            else:
                gt = gt.split("/")

            i = record.format.index("AD")
            ad = x.split(":")[i]
            ad = [int(x) for x in ad.split(",")]

            if "." in gt or gt[0] == gt[1] or sum(ad) == 0:
                af = 0.0
            else:
                af = max(ad) / sum(ad)

            return af

        filtered = []

        for i in reversed(range(len(self.data))):
            record = self.data[i]

            if "AD" not in record.format:
                continue

            afs = [get_af(x, record) for x in record.data]
            afs = [x for x in afs if x]

            if not afs:
                continue

            median = statistics.median(afs)

            if median <= threshold:
                continue

            if record.has_indel() and len(afs) < count:
                continue

            if not record.filter or "PASS" in record.filter:
                record.filter = ["AllelicImbalanceFilter"]
            else:
                record.filter += ["AllelicImbalanceFilter"]

            filtered.append(self.data.pop(i))

        return filtered

    def search_meta(self, key: str) -> str:
        """Returns the value of the meta-information key, if present."""
        for line in self.meta:
            if key in line:
                return line.strip().replace(f"##{key}=", "")
        return ""

    def copy(self,
             meta: Optional[List[str]] = None,
             header: Optional[List[str]] = None,
             data: Optional[List[Record]] = None) -> "VCFFile":
        """Returns a copy of the VCFFile."""
        new = VCFFile()
        new.meta = copy.deepcopy(self.meta) if meta is None else meta
        new.header = copy.deepcopy(self.header) if header is None else header
        new.data = copy.deepcopy(self.data) if data is None else data
        return new

    def to_str(self) -> str:
        """Returns a string representation of the VCFFile."""
        s = ""
        for line in self.meta:
            s += line
        s += "\t".join(self.header) + "\n"
        for record in self.data:
            s += "\t".join(record.fields) + "\n"
        return s

    def to_file(self, filepath: str) -> None:
        """Writes the VCFFile to a file."""
        with open(filepath, "w") as f:
            f.write(self.to_str())

    def phase(self) -> None:
        """Changes genotype separator from ``/`` to ``|``."""
        for record in self.data:
            record.fields[9:] = [x.replace("/", "|") for x
                                 in record.fields[9:]]

    def unphase(self) -> None:
        """Changes genotype separator from ``|`` to ``/``."""
        for record in self.data:
            record.fields[9:] = [x.replace("|", "/") for x
                                 in record.fields[9:]]

    def sort(self) -> None:
        """Sorts SNP records by chromosome and then by position."""
        self.data.sort(key = lambda x: (x.chrom, x.pos))




















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

def get_stardb(tg: str, gb: str) -> Dict[str, StarAllele]:
    """
    Get StarAllele database.

    Returns:
        dict[str, StarAllele]: StarAllele objects.

    Args:
        tg (str): Target gene.
        gb (str): Genome build.

    Examples:

        >>> stardb = get_stardb("cyp2d6", "hg19")
        >>> print(stardb["*2"].score)
        1.0
    """

    p = os.path.dirname(__file__)
    gene_table = read_gene_table(f"{p}/resources/sg/gene_table.txt")
    snp_table = read_snp_table(f"{p}/resources/sg/snp_table.txt", gene_table)
    snpdb = build_snpdb(tg, gb, snp_table)
    star_table = read_star_table(f"{p}/resources/sg/star_table.txt")

    return build_stardb(tg, gb, star_table, snpdb)

def get_gene_table() -> Dict[str, Dict[str, str]]:
    """
    Get gene table object.

    Returns:
        dict[str, dict[str, str]]: Gene table object.
    """

    p = os.path.dirname(__file__)
    return read_gene_table(f"{p}/resources/sg/gene_table.txt")

def get_target_genes() -> List[str]:
    """Get the list of target gene names.

    Returns:
        list[str]: A list of gene names.
    """
    gene_table = get_gene_table()
    return [k for k, v in gene_table.items() if v["type"] == "target"]

def get_target_region(tg: str, gb: str) -> str:
    """Get the genomic region for the target gene.

    Returns:
        str: Genomic region.

    Args:
        tg (str): Target gene.
        gb (str): Genome build (hg19, hg38).
    """
    gene_table = get_gene_table()
    target_genes = [k for k, v in gene_table.items() if v["type"] == "target"]

    if tg not in target_genes:
        raise ValueError(f"'{tg}' is not among target genes: {target_genes}")

    return gene_table[tg][f"{gb}_region"]

def get_file_list(
        td: str,
        fe: Optional[str] = None
    ) -> List[str]:
    """ Get the list of files from the target directory.

    Returns:
        list[str]: List of files.

    Args:
        td (str): Target directory.
        fe (str, optional): File extension.
    """
    result = []
    for r, d, f in os.walk(td):
        for fn in f:
            if fe and not fn.endswith(fe):
                continue
            result.append(os.path.join(r, fn))
    return result

def read_file_list(fl: str) -> List[str]:
    """ Get the list of files from the list file.

        Returns:
            list[str]: List of files.

        Args:
            fl (str): List file.
    """
    result = []
    with open(fl) as f:
        for line in f:
            if not line.strip():
                continue
            result.append(line.strip())
    return result

def randstr(
    chars: str = string.ascii_uppercase + string.digits,
    n: int = 5
) -> str:
    """Generate a random string of length n."""
    first = random.choice(string.ascii_lowercase)
    return first + "".join(random.choice(chars) for _ in range(n))

def temp_env(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        if "temp_dir" in kwargs and kwargs["temp_dir"]:
            temp_obj = None
            temp_path = os.path.realpath(kwargs["temp_dir"])
        else:
            temp_obj = TemporaryDirectory()
            temp_path = temp_obj.name

        func(*args, **kwargs, temp_path=temp_path)

        if temp_obj:
            temp_obj.cleanup()

    return wrapper

def bam_getter(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        if kwargs["bam_file"]:
            input_files = kwargs["bam_file"]
        elif kwargs["bam_dir"]:
            bam_path = os.path.realpath(kwargs["bam_dir"])
            input_files = []
            for r, d, f in os.walk(bam_path):
                for x in f:
                    if x.endswith("bam"):
                        input_files.append(f"{bam_path}/{x}")
        elif kwargs["bam_list"]:
            input_files = []
            with open(kwargs["bam_list"]) as f:
                for line in f:
                    input_files.append(line.strip())
        else:
            input_files = []

        if not input_files:
            raise ValueError("No input BAM files found")

        result = func(*args, **kwargs, input_files=input_files)

        return result

    return wrapper

def get_logger():
    logging.basicConfig(
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    return logging.getLogger()

def conf_env(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        logger = get_logger()
        logger.info("Configuration:")

        with open(kwargs["conf_file"]) as f:
            for line in f:
                logger.info("    " + line.strip())

        config = configparser.ConfigParser()
        config.read(kwargs["conf_file"])

        func(*args, **kwargs, config=config)

    return wrapper
