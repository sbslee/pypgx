import configparser
from os import mkdir
from os.path import realpath
from .common import is_chr, randstr, conf_env, get_target_region

@conf_env
def bam2vcf2(conf_file: str, **kwargs) -> None:
    """Convert BAM files to a VCF file [SGE].

    This command outputs a single- or multi-sample VCF file from one or 
    more input BAM files. The output VCF file will only contain variants 
    within the target gene or region. This command is essentially a 
    wrapper with pre-specified parameters for the Genome Analysis Toolkit 
    (GATK). It also uses Sun Grid Engine (SGE) for parallelism to make 
    GATK run faster.

    Args:
        conf_file (str): Configuration file.

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
            conda_env = NONE
            dbsnp_file = NONE
            java_options = NONE
            qsub_options = NONE

            # Make any necessary changes to this section.
            [USER]
            bam_list = bam-list.txt
            conda_env = env_name
            dbsnp_file = dbsnp.vcf
            fasta_file = reference.fa
            genome_build = hg19
            java_options = -Xmx4G
            project_path = ./myproject
            qsub_options = -l mem_requested=4G
            target_gene = cyp2d6

    This table summarizes the configuration parameters specific to 
    ``bam2vcf2``:

        .. list-table::
           :widths: 25 75
           :header-rows: 1

           * - Parameter
             - Summary
           * - bam_list
             - List of input BAM files, one file per line.
           * - conda_env
             - Name of conda environment to be activated.
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
           * - target_gene
             - Name of target gene (e.g. 'cyp2d6'). 
               Also accepts a BED file.
    """
    config = kwargs["config"]

    # Parse the configuration data.
    project_path = realpath(config["USER"]["project_path"])
    bam_list = realpath(config["USER"]["bam_list"])
    fasta_file = realpath(config["USER"]["fasta_file"])
    target_gene = config["USER"]["target_gene"]
    genome_build = config["USER"]["genome_build"]
    qsub_options = config["USER"]["qsub_options"]
    java_options = config["USER"]["java_options"]
    dbsnp_file = config["USER"]["dbsnp_file"]
    conda_env = config["USER"]["conda_env"]

    bam_files = []

    with open(bam_list) as f:
        for line in f:
            bam_files.append(line.strip())

    # Make the project directories.
    mkdir(project_path)
    mkdir(f"{project_path}/shell")
    mkdir(f"{project_path}/log")
    mkdir(f"{project_path}/temp")

    t = [is_chr(x) for x in bam_files]
    if all(t):
        chr_str = "chr"
    elif not any(t):
        chr_str = ""
    else:
        raise ValueError("Mixed types of SN tags found.")

    if target_gene.endswith(".bed"):
        target_region = realpath(target_gene)
    else:
        target_region = chr_str + get_target_region(target_gene, genome_build).replace("chr", "")

    # Write the shell script for HaplotypeCaller.
    for i, bam_file in enumerate(bam_files):
        s = "#!/bin/bash\n"

        if conda_env != "NONE":
            s += (
                "\n"
                f"conda activate {conda_env}\n"
            )

        s += (
            "\n"
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
    s = "#!/bin/bash\n"

    if conda_env != "NONE":
        s += (
            "\n"
            f"conda activate {conda_env}\n"
        )

    s += (
        "\n"
        f"p={project_path}\n"
        "\n"
        "gatk GenomicsDBImport \\\n"
        f"  --intervals {target_region} \\\n"
        f"  --genomicsdb-workspace-path $p/temp/datastore \\\n"
        "  --merge-input-intervals \\\n"
        "  --QUIET \\\n"
    )

    if java_options != "NONE":
        s += f"  --java-options {java_options} \\\n"

    for i, bam_file in enumerate(bam_files):
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
        f"  -O $p/pypgx.vcf \\\n"
        f"  --variant $p/temp/pypgx.joint.vcf \\\n"
        "  --filter-expression 'QUAL <= 50.0' \\\n"
        "  --filter-name QUALFilter \\\n"
        "  --QUIET\n"
    )

    if java_options != "NONE":
        s += f"  --java-options {java_options} \\\n"

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
        f"j={randstr()}\n"
        "\n"
    )

    for i, bam_file in enumerate(bam_files):
        s += f"{q} -N $j-hc $p/shell/haplotypecaller-{i}.sh\n"

    s += f"{q} -hold_jid $j-hc -N $j-post-hc $p/shell/post-haplotypecaller.sh\n"

    with open(f"{project_path}/example-qsub.sh", "w") as f:
        f.write(s)
