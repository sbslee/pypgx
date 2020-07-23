import configparser
from os import mkdir
from os.path import realpath

from .common import logging, LINE_BREAK1, is_chr, get_gene_table

logger = logging.getLogger(__name__)

def sgea(conf: str) -> None:
    """Run per-project genotyping for single gene with SGE (2).

    This command runs the per-project genotyping pipeline by submitting 
    jobs to the Sun Grid Engine (SGE) cluster. The main difference between
    ``sgea`` and ``sgep`` is that the former uses Genome Analysis Tool 
    Kit (GATK) v4 while the latter uses BCFtools.

    Args:
        conf (str): Configuration file.

    .. note::

        SGE, Stargazer, and GATK must be pre-installed.

    This is what a typical configuration file for ``sgea`` looks like:

        .. code-block:: python

            # File: example_conf.txt
            # Do not make any changes to this section.
            [DEFAULT]
            mapping_quality = 1
            output_prefix = pypgx
            control_gene = NONE
            qsub_options = NONE

            # Make any necessary changes to this section.
            [USER]
            fasta_file = reference.fa
            manifest_file = manifest.txt
            project_path = /path/to/project/
            target_gene = cyp2d6
            genome_build = hg19
            data_type = wgs
            dbsnp_file = dbsnp.vcf
            control_gene = vdr

    This table summarizes the configuration parameters specific to ``sgea``:

        .. list-table::
           :widths: 25 75
           :header-rows: 1

           * - Parameter
             - Summary
           * - control_gene
             - Control gene or region.
           * - data_type
             - Data type (wgs, ts, chip).
           * - dbsnp_file
             - dbSNP VCF file.
           * - fasta_file
             - Reference FASTA file.
           * - genome_build
             - Genome build (hg19, hg38).
           * - manifest_file
             - Manifest file.
           * - mapping_quality
             - Minimum mapping quality for read depth.
           * - output_prefix
             - Output prefix.
           * - project_path
             - Output project directory.
           * - qsub_options
             - Options for qsub command.
           * - target_gene
             - Target gene.
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
    manifest_file = realpath(config["USER"]["manifest_file"])
    dbsnp_file = realpath(config["USER"]["dbsnp_file"])
    fasta_file = realpath(config["USER"]["fasta_file"])
    mapping_quality = config["USER"]["mapping_quality"]
    output_prefix = config["USER"]["output_prefix"]
    target_gene = config["USER"]["target_gene"]
    control_gene = config["USER"]["control_gene"]
    data_type = config["USER"]["data_type"]
    genome_build = config["USER"]["genome_build"]
    qsub_options = config["USER"]["qsub_options"]

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

    # Log the number of samples.
    logger.info(f"Number of samples: {len(bam_files)}")

    # Make the project directories.
    project_path = f"{project_path}"
    mkdir(project_path)
    mkdir(f"{project_path}/shell")
    mkdir(f"{project_path}/log")
    mkdir(f"{project_path}/temp")

    if control_gene != "NONE":
        # Write the shell script for bam2gdf.
        s = f"pypgx bam2gdf {genome_build} {target_gene} {control_gene} \\\n"

        for sample_id in bam_files:
            s += f"  {bam_files[sample_id]} \\\n"

        s += f"  -o {project_path}/{output_prefix}.gdf\n"

        with open(f"{project_path}/shell/doc.sh", "w") as f:
            f.write(s)

    # Write the shell script for HaplotypeCaller.
    target_region = gene_table[target_gene][f"{genome_build}_region"].replace("chr", "")

    t = [is_chr(v) for k, v in bam_files.items()]
    if all(t):
        chr_str = "chr"
    elif not any(t):
        chr_str = ""
    else:
        raise ValueError("Mixed types of SN tags found.")

    for sample_id in bam_files:
        s = (
            f"project={project_path}\n"
            f"bam={bam_files[sample_id]}\n"
            f"fasta={fasta_file}\n"
            f"dbsnp={dbsnp_file}\n"
            f"gvcf=$project/temp/{sample_id}.g.vcf\n"
            f"region={chr_str}{target_region}\n"
            "\n"
            "gatk HaplotypeCaller \\\n"
            "  -R $fasta \\\n"
            "  -D $dbsnp \\\n"
            "  --emit-ref-confidence GVCF \\\n"
            "  -I $bam \\\n"
            "  -O $gvcf \\\n"
            "  -L $region\n"
        )

        with open(f"{project_path}/shell/hc-{sample_id}.sh", "w") as f:
            f.write(s)

    # Write the schell script for post-HaplotypeCaller.
    s = (
        f"project={project_path}\n"
        f"fasta={fasta_file}\n"
        f"dbsnp={dbsnp_file}\n"
        f"gvcf=$project/temp/{output_prefix}.g.vcf\n"
        f"region={chr_str}{target_region}\n"
        "\n"
        "gatk CombineGVCFs \\\n"
        "  -R $fasta \\\n"
        "  -D $dbsnp \\\n"
        "  -L $region \\\n"
    )

    for sample_id in bam_files:
        s += f"  --variant $project/temp/{sample_id}.g.vcf \\\n"

    s += (
        "  -O $gvcf\n"
        "\n"
        f"vcf1=$project/temp/{output_prefix}.joint.vcf\n"
        "\n"
        "gatk GenotypeGVCFs \\\n"
        "  -R $fasta \\\n"
        "  -D $dbsnp \\\n"
        "  -O $vcf1 \\\n"
        "  -L $region \\\n"
        "  --variant $gvcf\n"
        "\n"
        f"vcf2=$project/{output_prefix}.joint.filtered.vcf\n"
        "\n"
        "gatk VariantFiltration \\\n"
        "  -R $fasta \\\n"
        "  -O $vcf2 \\\n"
        "  -L $region \\\n"
        "  --filter-expression 'QUAL <= 50.0' \\\n"
        "  --filter-name QUALFilter \\\n"
        "  --variant $vcf1\n"
    )

    with open(f"{project_path}/shell/post.sh", "w") as f:
        f.write(s)

    # Write the shell script for Stargazer.
    s = (
        f"project={project_path}\n"
        f"vcf=$project/{output_prefix}.joint.filtered.vcf\n"
        "\n"
        "stargazer \\\n"
        f"  {data_type} \\\n"
        f"  {genome_build} \\\n"
        f"  {target_gene} \\\n"
        "  $vcf \\\n"
        "  $project/stargazer \\\n"
    )

    if control_gene != "NONE":
        s += (
            f"  --cg {control_gene} \\\n"
            f"  --gdf $project/{output_prefix}.gdf \\\n"
        )

    with open(f"{project_path}/shell/rs.sh", "w") as f:
        f.write(s)

    # Write the shell script for qsub.
    q = "qsub -e $p/log -o $p/log"

    if qsub_options != "NONE":
        q += f" {qsub_options}"

    s = (
        f"p={project_path}\n"
        "\n"
    )

    if control_gene != "NONE":
        s += f"{q} -N doc $p/shell/doc.sh\n"

    for sample_id in bam_files:
        s += f"{q} -N hc $p/shell/hc-{sample_id}.sh\n"

    s += f"{q} -hold_jid hc -N post $p/shell/post.sh\n"

    if control_gene == "NONE":
        s += f"{q} -hold_jid post -N rs $p/shell/rs.sh\n"

    else:
        s += f"{q} -hold_jid doc,post -N rs $p/shell/rs.sh\n"

    with open(f"{project_path}/example-qsub.sh", "w") as f:
        f.write(s)
