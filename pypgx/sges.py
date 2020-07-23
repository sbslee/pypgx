import configparser
from os import mkdir
from os.path import realpath

from .common import logging, LINE_BREAK1, is_chr, get_gene_table

logger = logging.getLogger(__name__)

def sges(conf: str) -> None:
    """Run per-sample genotyping for multiple genes with SGE.

    This command runs the per-sample genotyping pipeline by submitting 
    jobs to the Sun Grid Engine (SGE) cluster. After genotype analysis by 
    Stargazer, it will generate a HTML report using the ``report`` tool.

    Args:
        conf (str): Configuration file.

    .. note::

        BCFtools, SGE and Stargazer must be pre-installed.

    This is what a typical configuration file for ``sges`` looks like:

        .. code-block:: python

            # File: example_conf.txt
            # Do not make any changes to this section.
            [DEFAULT]
            mapping_quality = 1
            output_prefix = pypgx
            target_genes = ALL
            control_gene = NONE
            qsub_options = NONE

            # Make any necessary changes to this section.
            [USER]
            fasta_file = reference.fa
            project_path = /path/to/project/
            genome_build = hg19
            data_type = wgs
            bam_file = in.bam
            qsub_options = -l mem_requested=2G
            target_genes = cyp2b6, cyp2d6, dpyd
            control_gene = vdr

    This table summarizes the configuration parameters specific to ``sges``:

        .. list-table::
           :widths: 25 75
           :header-rows: 1

           * - Parameter
             - Summary
           * - bam_file
             - BAM file.
           * - control_gene
             - Control gene or region.
           * - data_type
             - Data type (wgs, ts, chip).
           * - fasta_file
             - Reference FASTA file.
           * - genome_build
             - Genome build (hg19, hg38).
           * - mapping_quality
             - Minimum mapping quality for read detph.
           * - output_prefix
             - Output prefix.
           * - project_path
             - Output project directory.
           * - qsub_options
             - Options for qsub command.
           * - target_genes
             - Target genes.
    """

    gene_table = get_gene_table()

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
    project_path = realpath(config["USER"]["project_path"])
    bam_file = realpath(config["USER"]["bam_file"])
    fasta_file = realpath(config["USER"]["fasta_file"])
    target_genes = config["USER"]["target_genes"]
    output_prefix = config["USER"]["output_prefix"]
    control_gene = config["USER"]["control_gene"]
    mapping_quality = config["USER"]["mapping_quality"]
    data_type = config["USER"]["data_type"]
    genome_build = config["USER"]["genome_build"]
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
    mkdir(project_path)
    mkdir(f"{project_path}/gene")

    for select_gene in select_genes:
        mkdir(f"{project_path}/gene/{select_gene}")

    if is_chr(bam_file):
        chr_str = "chr"
    else:
        chr_str = ""

    # Write the shell script for genotyping pipeline.
    for select_gene in select_genes:
        target_region = gene_table[select_gene][f"{genome_build}_region"].replace("chr", "")

        s = (
            "pypgx genotype \\\n"
            f"  {fasta_file} \\\n"
            f"  {data_type} \\\n"
            f"  {genome_build} \\\n"
            f"  {select_gene} \\\n"
            f"  {project_path}/gene/{select_gene}/stargazer \\\n"
            f"  {bam_file} \\\n"
        )

        if control_gene != "NONE":
            s += f"  --cg {control_gene}\n"

        with open(
            f"{project_path}/gene/{select_gene}/run.sh", "w"
        ) as f:
            f.write(s)

    # Write the shell script for report.
    s = (
        f"p={project_path}\n"
        "\n"
        f'''head -n1 $p/gene/{select_genes[0]}/stargazer/genotype.txt > $p/genotype.merged.txt\n'''
        "\n"
        f"for gene in {' '.join(select_genes)}\n"
        "do\n"
        f'''  tail -n+2 $p/gene/$gene/stargazer/genotype.txt >> $p/genotype.merged.txt\n'''
        "done\n"
        "\n"
        "pypgx report $p/genotype.merged.txt -o $p/report.html"
    )

    with open(f"{project_path}/report.sh", "w") as f:
        f.write(s)

    # Write the shell script for qsub.
    if qsub_options == "NONE":
        q = "qsub"
    else:
        q = f"qsub {qsub_options}"

    s = (
        f"p={project_path}\n"
        "\n"
    )

    for select_gene in select_genes:
        s += f"{q} -o $p/gene/{select_gene} -e $p/gene/{select_gene} -N run $p/gene/{select_gene}/run.sh\n"

    s += (
        f"{q} -o $p -e $p -hold_jid run $p/report.sh"
    )

    with open(f"{project_path}/example-qsub.sh", "w") as f:
        f.write(s)
