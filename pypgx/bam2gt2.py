import configparser
from os import mkdir
from os.path import realpath
from .bam2vcf2 import bam2vcf2
from .common import logging, get_target_genes, sm_tag, LINE_BREAK1, randstr

logger = logging.getLogger(__name__)

def _write_bam2gdf_shell(
        genome_build,
        target_gene,
        control_gene,
        bam_files,
        gdf_file,
        shell_file
    ):
    s = (
        "pypgx bam2gdf \\\n"
        f"  {genome_build} \\\n"
        f"  {target_gene} \\\n"
        f"  {control_gene} \\\n"
        f"  {gdf_file} \\\n"
    )

    for name in bam_files:
        s += f"  {bam_files[name]} \\\n"

    with open(shell_file, "w") as f:
        f.write(s)

def _write_bam2vcf_shell(
        fasta_file,
        target_gene,
        project_path,
        genome_build,
        bam_files
    ):
    s = (
        "pypgx bam2vcf \\\n"
        f"  bcftools \\\n"
        f"  {fasta_file} \\\n"
        f"  {target_gene} \\\n"
        f"  {project_path}/pypgx.vcf \\\n"
        f"  {genome_build} \\\n"
    )

    for name in bam_files:
        s += f"  {bam_files[name]} \\\n"

    with open(f"{project_path}/shell/bam2vcf.sh", "w") as f:
        f.write(s)

def _write_bam2vcf2_shell(
        fasta_file,
        bam_list,
        project_path,
        target_gene,
        genome_build,
        qsub_options,
        java_options,
        dbsnp_file
    ):
    s = (
        "# Do not make any changes to this section.\n"
        "[DEFAULT]\n"
        "qsub_options = NONE\n"
        "java_options = NONE\n"
        "dbsnp_file = NONE\n"
        "\n"
        "# Make any necessary changes to this section.\n"
        "[USER]\n"
        f"fasta_file = {fasta_file}\n"
        f"bam_list = {bam_list}\n"
        f"project_path = {project_path}/bam2vcf2\n"
        f"target_gene = {target_gene}\n"
        f"genome_build = {genome_build}\n"
        f"qsub_options = {qsub_options}\n"
        f"java_options = {java_options}\n"
        f"dbsnp_file = {dbsnp_file}\n"
    )

    with open(f"{project_path}/conf.txt", "w") as f:
        f.write(s)

    bam2vcf2(f"{project_path}/conf.txt")

def _write_stargazer_shell(
        snp_caller,
        data_type,
        genome_build,
        target_gene,
        project_path,
        control_gene
    ):

    if snp_caller == "bcftools":
        vcf_file = "$p/pypgx.vcf"
    else:
        vcf_file = "$p/bam2vcf2/bam2vcf2.vcf"

    s = (
        f"p={project_path}\n"
        "\n"
        "stargazer \\\n"
        f"  {data_type} \\\n"
        f"  {genome_build} \\\n"
        f"  {target_gene} \\\n"
        f"  {vcf_file} \\\n"
        "  $p/stargazer \\\n"
    )

    if control_gene != "NONE":
        s += (
            f"  --cg {control_gene} \\\n"
            f"  --gdf $p/pypgx.gdf \\\n"
        )

    with open(f"{project_path}/shell/stargazer.sh", "w") as f:
        f.write(s)

def _write_qsub_shell(
        snp_caller,
        qsub_options,
        control_gene,
        project_path
    ):
    q = "qsub -e $p/log -o $p/log"

    if qsub_options != "NONE":
        q += f" {qsub_options}"

    s = (
        "#!/bin/bash\n"
        "\n"
        f"p={project_path}\n"
        f"j={randstr()}\n"
        "\n"
    )

    if control_gene != "NONE":
        s += f"{q} -N $j-bam2gdf $p/shell/bam2gdf.sh\n"

    if snp_caller == "bcftools":
        s += f"{q} -N $j-bam2vcf $p/shell/bam2vcf.sh\n"

        if control_gene == "NONE":
            s += f"{q} -hold_jid $j-bam2vcf -N $j-stargazer $p/shell/stargazer.sh\n"
        else:
            s += f"{q} -hold_jid $j-bam2gdf,$j-bam2vcf -N $j-stargazer $p/shell/stargazer.sh\n"

    else:
        with open(f"{project_path}/bam2vcf2/example-qsub.sh") as f:
            for line in f:
                if line.startswith("qsub"):
                    s += line.replace("$p", "$p/bam2vcf2")

        if control_gene == "NONE":
            s += f"{q} -hold_jid $j-post-hc -N $j-stargazer $p/shell/stargazer.sh\n"
        else:
            s += f"{q} -hold_jid $j-bam2gdf,$j-post-hc -N $j-stargazer $p/shell/stargazer.sh\n"

    with open(f"{project_path}/example-qsub.sh", "w") as f:
        f.write(s)

def bam2gt2(conf: str) -> None:
    """Convert BAM files to genotype files [SGE].

    This command runs the entire genotyping pipeline for BAM files 
    with the Sun Grid Engine (SGE) cluster. By default, it will genotype 
    all genes currently targeted by the Stargazer program (you can specify 
    select genes too). For each gene, the command runs under the hood 
    ``bam2vcf`` with ``bcftools`` caller (i.e. BCFtools) or ``bam2vcf2`` 
    (i.e. GATK) to create the input VCF file. The input GDF file is 
    created with ``bam2gdf``.

    Args:
        conf (str): Configuration file.

    .. warning::

        SGE, Stargazer and BCFtools/GATK must be pre-installed.

    This is what a typical configuration file for ``bam2gt2`` looks like:

        .. code-block:: python

            # File: example_conf.txt
            # To execute:
            #   $ pypgx bam2gt2 example_conf.txt
            #   $ sh ./myproject/example-qsub.sh

            # Do not make any changes to this section.
            [DEFAULT]
            control_gene = NONE
            dbsnp_file = NONE
            java_options = NONE
            qsub_options = NONE
            target_genes = ALL

            # Make any necessary changes to this section.
            [USER]
            bam_list = bam-list.txt
            control_gene = vdr
            data_type = wgs
            fasta_file = hs37d5.fa
            genome_build = hg19
            project_path = ./myproject
            qsub_options = -l mem_requested=2G
            snp_caller = gatk
            target_genes = cyp2b6, cyp2d6

    This table summarizes the configuration parameters specific to 
    ``bam2gt2``:

        .. list-table::
           :widths: 25 75
           :header-rows: 1

           * - Parameter
             - Summary
           * - bam_list
             - List of input BAM files, one file per line.
           * - control_gene
             - Control gene or region.
           * - data_type
             - Data type ('wgs' or 'ts').
           * - dbsnp_file
             - dbSNP VCF file.
           * - fasta_file
             - Reference FASTA file.
           * - genome_build
             - Genome build ('hg19' or 'hg38').
           * - java_options
             - Java-specific arguments for GATK (e.g. ‘-Xmx4G’).
           * - project_path
             - Output project directory.
           * - qsub_options
             - Options for qsub command (e.g. '-l mem_requested=2G').
           * - snp_caller
             - SNP caller (‘gatk’ or ‘bcftools’).
           * - target_genes
             - Names of target genes (e.g. 'cyp2d6').
    """
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
    fasta_file = realpath(config["USER"]["fasta_file"])
    genome_build = config["USER"]["genome_build"]
    target_genes = config["USER"]["target_genes"]
    control_gene = config["USER"]["control_gene"]
    data_type = config["USER"]["data_type"]
    qsub_options = config["USER"]["qsub_options"]
    snp_caller = config["USER"]["snp_caller"]
    java_options = config["USER"]["java_options"]
    dbsnp_file = config["USER"]["dbsnp_file"]
    bam_list = realpath(config["USER"]["bam_list"])

    bam_files = {}

    with open(bam_list) as f:
        for line in f:
            bam = line.strip()
            name = sm_tag(bam)
            bam_files[name] = bam

    all_genes = get_target_genes()

    if target_genes == "ALL":
        select_genes = all_genes
    else:
        select_genes = []

        for gene in target_genes.split(","):
            select_genes.append(gene.strip().lower())

        for gene in select_genes:
            if gene not in all_genes:
                raise ValueError(f"Unrecognized target gene found: {gene}")

    # Sort the samples by name since GATK does this.
    bam_files = {k: v for k, v in sorted(bam_files.items(), key=lambda x: x[0])}

    # Log the number of samples.
    logger.info(f"Number of samples: {len(bam_files)}")

    # Make the project directories.
    mkdir(project_path)
    mkdir(f"{project_path}/gene")

    s = (
        "#!/bin/bash\n"
        "\n"
    )

    for select_gene in select_genes:
        s += f"sh {project_path}/gene/{select_gene}/example-qsub.sh\n"

    with open(f"{project_path}/example-qsub.sh", "w") as f:
        f.write(s)

    for select_gene in select_genes:
        gene_path = f"{project_path}/gene/{select_gene}"

        mkdir(gene_path)
        mkdir(f"{gene_path}/shell")
        mkdir(f"{gene_path}/log")

        if control_gene != "NONE":
            _write_bam2gdf_shell(
                genome_build,
                select_gene,
                control_gene,
                bam_files,
                f"{gene_path}/pypgx.gdf",
                f"{gene_path}/shell/bam2gdf.sh"
            )

        if snp_caller == "bcftools":
            _write_bam2vcf_shell(
                fasta_file,
                select_gene,
                gene_path,
                genome_build,
                bam_files
            )

        elif snp_caller == "gatk":
            _write_bam2vcf2_shell(
                fasta_file,
                bam_list,
                gene_path,
                select_gene,
                genome_build,
                qsub_options,
                java_options,
                dbsnp_file
            )

        else:
            raise ValueError(f"Incorrect SNP caller: {snp_caller}")

        _write_stargazer_shell(
            snp_caller,
            data_type,
            genome_build,
            select_gene,
            gene_path,
            control_gene
        )

        _write_qsub_shell(
            snp_caller,
            qsub_options,
            control_gene,
            gene_path
        )
