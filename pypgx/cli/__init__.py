from ._compare_stargazer_calls import compare_stargazer_calls
from ._calculate_read_depth import calculate_read_depth
from ._call_variants_gatk_sge import call_variants_gatk_sge

commands = {
    "compare-stargazer-calls": compare_stargazer_calls,
    "calculate-read-depth": calculate_read_depth,
    "call-variants-gatk-sge": call_variants_gatk_sge,
}
