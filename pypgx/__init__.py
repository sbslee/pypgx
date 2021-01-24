import pkgutil
import pandas as pd
from io import BytesIO
from pypgx.sdk import Locus

from ._compare_stargazer_calls import compare_stargazer_calls
from ._calculate_read_depth import calculate_read_depth

from._phenotyper import phenotyper

gene_df = pd.read_table(BytesIO(
    pkgutil.get_data(__name__, "resources/gene_table.tsv")))

snp_df = pd.read_table(BytesIO(
    pkgutil.get_data(__name__, "resources/snp_table.tsv")))

star_df = pd.read_table(BytesIO(
    pkgutil.get_data(__name__, "resources/star_table.tsv")))

target_genes = gene_df[gene_df["type"] == "target"]["name"]
control_genes = gene_df[gene_df["control"] == "yes"]["name"]

commands = {
    "compare-stargazer-calls": compare_stargazer_calls,
    "calculate-read-depth": calculate_read_depth,
}

__all__ = [
    "compare_stargazer_calls",
    "calculate_read_depth",
    "phenotyper",
]
