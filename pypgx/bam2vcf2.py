import configparser
from os import mkdir
from os.path import realpath
from .common import logging, LINE_BREAK1, is_chr, get_gene_table, randstr

logger = logging.getLogger(__name__)

def bam2vcf2(conf: str) -> None:
    """Convert BAM files to a VCF file.

    This command outputs a single- or multi-sample VCF file from one or 
    more input BAM files. The output VCF file will only contain variants 
    within the target gene or region. This command is essentially a 
    wrapper with pre-specified parameters for the Genome Analysis Toolkit 
    (GATK). It also uses Sun Grid Engine (SGE) for parallelism to make 
    GATK run faster.

    Args:
        conf (str): Configuration file.

    .. warning::
        GATK and SGE must be pre-installed.

    This is what a typical configuration file for ``bam2vcf2`` looks like:

        .. code-block:: python

            # File: example_conf.txt
            # To execute:
            #   $ pypgx bam2vcf2 example_conf.txt
            #   $ sh ./myproject/example-qsub.sh

            # Do not make any changes to this section.
            [DEFAULT]
            qsub_options = NONE
            java_options = NONE
            dbsnp_file = NONE

            # Make any necessary changes to this section.
            [USER]
            fasta_file = reference.fa
            manifest_file = manifest.txt
            project_path = ./myproject
            target_gene = cyp2d6
            genome_build = hg19
            qsub_options = -l mem_requested=4G
            java_options = -Xmx4G
            dbsnp_file = dbsnp.vcf

    This table summarizes the configuration parameters specific to 
    ``bam2vcf2``:

        .. list-table::
           :widths: 25 75
           :header-rows: 1

           * - Parameter
             - Summary
           * - dbsnp_file
             - dbSNP VCF file.
           * - fasta_file
             - Reference FASTA file.
           * - genome_build
             - Genome build (hg19, hg38).
           * - java_options
             - Java-specific arguments for GATK (e.g. ‘-Xmx4G’).
           * - manifest_file
             - Manifest file.
           * - project_path
             - Output project directory.
           * - qsub_options
             - Options for qsub command.
           * - target_gene
             - Target gene.
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
    manifest_file = realpath(config["USER"]["manifest_file"])
    fasta_file = realpath(config["USER"]["fasta_file"])
    target_gene = config["USER"]["target_gene"]
    genome_build = config["USER"]["genome_build"]
    qsub_options = config["USER"]["qsub_options"]
    java_options = config["USER"]["java_options"]
    dbsnp_file = config["USER"]["dbsnp_file"]

    # Read the manifest file.
    bam_files = {}
    with open(manifest_file) as f:
        header = next(f).strip().split("\t")
        i1 = header.index("sample_id")
        i2 = header.index("bam")
        for line in f:
            fields = line.strip().split("\t")
            sample_id = fields[i1]
            bam = fields[i2]
            bam_files[sample_id] = bam

    # Make the project directories.
    mkdir(project_path)
    mkdir(f"{project_path}/shell")
    mkdir(f"{project_path}/log")
    mkdir(f"{project_path}/temp")

    t = [is_chr(v) for k, v in bam_files.items()]
    if all(t):
        chr_str = "chr"
    elif not any(t):
        chr_str = ""
    else:
        raise ValueError("Mixed types of SN tags found.")

    gene_table = get_gene_table()
    target_region = gene_table[target_gene][f"{genome_build}_region"].replace("chr", "")

    # Write the shell script for HaplotypeCaller.
    for i, sample_id in enumerate(bam_files):
        bam_file = bam_files[sample_id]
        s = (
            "gatk HaplotypeCaller \\\n"
            f"  -R {fasta_file} \\\n"
            f"  --emit-ref-confidence GVCF \\\n"
            f"  -I {bam_file} \\\n"
            f"  -O {project_path}/temp/{i}.g.vcf \\\n"
            f"  -L {target_region} \\\n"
            "  --QUIET \\\n"
        )

        if java_options != "NONE":
            s += f"  --java-options {java_options} \\\n"

        with open(f"{project_path}/shell/haplotypecaller-{i}.sh", "w") as f:
            f.write(s)

    # Write the shell script for post-HaplotypeCaller.
    s = (
        "#!/bin/bash\n"
        "\n"
        f"p={project_path}\n"
        "\n"
    )

    s += (
        "gatk GenomicsDBImport \\\n"
        f"  --intervals {target_region} \\\n"
        f"  --genomicsdb-workspace-path $p/temp/datastore \\\n"
        "  --QUIET \\\n"
    )

    for i, sample_id in enumerate(bam_files):
        s += f"  -V $p/temp/{i}.g.vcf \\\n"

    s += (
        "\n"
        "gatk GenotypeGVCFs \\\n"
        f"  -R {fasta_file} \\\n"
        f"  -V gendb://$p/temp/datastore \\\n"
        f"  -O $p/temp/pypgx.joint.vcf \\\n"
        "  --QUIET \\\n"
    )

    if java_options != "NONE":
        s += f"  --java-options {java_options} \\\n"

    if dbsnp_file != "NONE":
        s += f"  -D {dbsnp_file} \\\n"

    s += (
        "\n"
        "gatk VariantFiltration \\\n"
        f"  -R {fasta_file} \\\n"
        f"  -L {target_region} \\\n"
        f"  -O $p/bam2vcf2.vcf \\\n"
        f"  --variant $p/temp/pypgx.joint.vcf \\\n"
        "  --filter-expression 'QUAL <= 50.0' \\\n"
        "  --filter-name QUALFilter \\\n"
        "  --QUIET\n"
    )

    with open(f"{project_path}/shell/post-haplotypecaller.sh", "w") as f:
        f.write(s)

    # Write the shell script for qsub.
    q = "qsub -e $p/log -o $p/log"

    if qsub_options != "NONE":
        q += f" {qsub_options}"

    s = (
        "#!/bin/bash\n"
        "\n"
        f"p={project_path}\n"
        "\n"
    )

    jid = randstr()

    for i, sample_id in enumerate(bam_files):
        s += f"{q} -N {jid}-haplotypecaller $p/shell/haplotypecaller-{i}.sh\n"

    s += f"{q} -hold_jid {jid}-haplotypecaller -N {jid}-post-haplotypecaller $p/shell/post-haplotypecaller.sh\n"

    with open(f"{project_path}/example-qsub.sh", "w") as f:
        f.write(s)
