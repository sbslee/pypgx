import configparser
import os

from .common import logging, LINE_BREAK1, is_chr
from .sglib import read_gene_table

logger = logging.getLogger(__name__)

def sges(conf: str) -> None:
    """
    Run per-sample genotyping with Stargazer.

    Args:
        conf (str): Configuration file.

    Examples:
        .. code-block:: python

            # File: example_conf.txt
            # Do not make any changes to this section.
            [DEFAULT]
            mapping_quality = 1
            output_prefix = pypgx
            target_genes = ALL
            genome_build = hg19
            control_gene = NONE
            vcf_only = FALSE
            qsub_options = NONE

            # Make any necessary changes to this section.
            [USER]
            fasta_file = reference.fa
            project_path = path/to/project/
            control_gene = vdr
            data_type = wgs
            dbsnp_file = dbsnp.vcf
            stargazer_tool = stargazer.py
            gatk_tool = gatk.jar
            bam_file = in.bam
    """

    gene_table = read_gene_table(
        f"{os.path.dirname(__file__)}/resources/sg/gene_table.txt")

    # Log the configuration data.
    logger.info(LINE_BREAK1)
    logger.info("Configureation:")
    with open(conf) as f:
        for line in f:
            logger.info("    " + line.strip())
    logger.info(LINE_BREAK1)

    # Read the configuration file.
    config = configparser.ConfigParser()
    config.read(conf)

    # Parse the configuration data.
    target_genes = config["USER"]["target_genes"]
    project_path = config["USER"]["project_path"]
    bam_file = config["USER"]["bam_file"]
    gatk_tool = config["USER"]["gatk_tool"]
    fasta_file = config["USER"]["fasta_file"]
    dbsnp_file = config["USER"]["dbsnp_file"]
    output_prefix = config["USER"]["output_prefix"]
    control_gene = config["USER"]["control_gene"]
    mapping_quality = config["USER"]["mapping_quality"]
    stargazer_tool = config["USER"]["stargazer_tool"]
    data_type = config["USER"]["data_type"]
    genome_build = config["USER"]["genome_build"]
    vcf_only = config["USER"].getboolean("vcf_only")
    qsub_options = config["USER"]["qsub_options"]

    t = [k for k, v in gene_table.items() if v["type"] == "target"]
    
    if target_genes == "ALL":
        select_genes = t
    else:
        select_genes = []
        for gene in target_genes.split(","):
            select_genes.append(gene.strip().lower())
        for gene in select_genes:
            if gene not in t:
                raise ValueError(f"Unrecognized target gene found: {gene}")

    # Make the project directories.
    os.mkdir(project_path)
    os.mkdir(f"{project_path}/log")
    os.mkdir(f"{project_path}/gene")

    for select_gene in select_genes:
        os.mkdir(f"{project_path}/gene/{select_gene}")
        os.mkdir(f"{project_path}/gene/{select_gene}/temp")
        os.mkdir(f"{project_path}/gene/{select_gene}/shell")
        os.mkdir(f"{project_path}/gene/{select_gene}/log")

    if is_chr(bam_file):
        chr_str = "chr"
    else:
        chr_str = ""

# -- Shell script for each gene ----------------------------------------------

    for select_gene in select_genes:
        target_region = gene_table[select_gene][f"{genome_build}_region"].replace("chr", "")

        s = (
            f"p={project_path}\n"
            f"tg={select_gene}\n"
            f"gatk={gatk_tool}\n"
            f"fasta={fasta_file}\n"
            f"dbsnp={dbsnp_file}\n"
            f"bam={bam_file}\n"
            f"region={chr_str}{target_region}\n"
            f"gvcf=$p/gene/$tg/temp/{output_prefix}.g.vcf\n"
            "\n"
            "java -jar $gatk -T HaplotypeCaller \\\n"
            "  -R $fasta \\\n"
            "  -D $dbsnp \\\n"
            "  -I $bam \\\n"
            "  -o $gvcf \\\n"
            "  -L $region \\\n"
            "  --emitRefConfidence GVCF \\\n"
            "  -U ALLOW_SEQ_DICT_INCOMPATIBILITY \\\n"
            "\n"
            f"vcf1=$p/gene/$tg/temp/{output_prefix}.joint.vcf\n"
            "\n"
            "java -jar $gatk -T GenotypeGVCFs \\\n"
            "  -R $fasta \\\n"
            "  -D $dbsnp \\\n"
            "  -o $vcf1 \\\n"
            "  -L $region \\\n"
            "  --variant $gvcf \\\n"
            "\n"
            f"vcf2=$p/gene/$tg/{output_prefix}.joint.filtered.vcf\n"
            "\n"
            "java -jar $gatk -T VariantFiltration \\\n"
            "  -R $fasta \\\n"
            "  --filterExpression 'QUAL <= 50.0' \\\n"
            "  --filterName QUALFilter \\\n"
            "  --variant $vcf1 \\\n"
            "  -o $vcf2 \\\n"
            "  -L $region \\\n"
            "\n"
        )

        if not vcf_only:
            s += (
                f"pypgx bam2gdf $tg {control_gene} $bam \\\n"
                "  -o $p/gene/$tg/$tg.gdf \\\n"
                "\n"
            )

        s += (
            f"stargazer={stargazer_tool}\n"
            "\n"
            f"python3 $stargazer \\\n"
            f"  {data_type} \\\n"
            f"  {genome_build} \\\n"
            "  $tg \\\n"
            "  $vcf2 \\\n"
            "  $p/gene/$tg/stargazer\n"
        )

        if not vcf_only:
            s += (
                f"  --cg {control_gene} \\\n"
                "  --gdf $p/gene/$tg/$tg.gdf \\\n"
            )

        with open(
            f"{project_path}/gene/{select_gene}/shell/run.sh", "w"
        ) as f:
            f.write(s)

# -- Shell script for reporting ----------------------------------------------

    s = (
        f"p={project_path}\n"
        "\n"
        f'''head -n1 $p/gene/{select_genes[0]}/stargazer/genotype.txt | awk '{{print "gene",$0}}' OFS="\\t" > $p/merged.sg-genotype.txt\n'''
        "\n"
        f"for gene in {' '.join(select_genes)}\n"
        "do\n"
        f'''  tail -n+2 $p/gene/$gene/stargazer/genotype.txt | awk -v var="$gene" '{{print var,$0}}' OFS="\\t" >> $p/merged.sg-genotype.txt\n'''
        "done\n"
        "\n"
        "pypgx report $p/merged.sg-genotype.txt -o $p/report.html"
    )

    with open(f"{project_path}/report.sh", "w") as f:
        f.write(s)

# -- Shell script for qsub ---------------------------------------------------

    if qsub_options == "NONE":
        q = "qsub"
    else:
        q = f"qsub {qsub_options}"

    s = (
        f"p={project_path}\n"
        "\n"
    )

    for select_gene in select_genes:
        s += f"{q} -o $p/gene/{select_gene}/log -e $p/gene/{select_gene}/log -N run $p/gene/{select_gene}/shell/run.sh\n"

    s += (
        f"{q} -o $p/log -e $p/log -hold_jid run $p/report.sh"
    )

    with open(f"{project_path}/example-qsub.sh", "w") as f:
        f.write(s)
