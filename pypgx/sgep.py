import configparser
from os import mkdir
from os.path import realpath

from .common import logging, LINE_BREAK1, is_chr, get_gene_table

logger = logging.getLogger(__name__)

def sgep(conf: str) -> None:
    """Run per-project genotyping for single gene with SGE (1).

    This command runs the per-project genotyping pipeline by submitting 
    jobs to the Sun Grid Engine (SGE) cluster.

    Args:
        conf (str): Configuration file.

    .. note::

        BCFtools, SGE and Stargazer must be pre-installed.

    This is what a typical configuration file for ``sgep`` looks like:

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
            control_gene = vdr
            qsub_options = -V -l mem_requested=10G

    This table summarizes the configuration parameters specific to ``sgep``:

        .. list-table::
           :widths: 25 75
           :header-rows: 1

           * - Parameter
             - Summary
           * - control_gene
             - Control gene or region.
           * - data_type
             - Data type (wgs, ts, chip).
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
    fasta_file = realpath(config["USER"]["fasta_file"])
    mapping_quality = config["USER"]["mapping_quality"]
    output_prefix = config["USER"]["output_prefix"]
    genome_build = config["USER"]["genome_build"]
    target_gene = config["USER"]["target_gene"]
    control_gene = config["USER"]["control_gene"]
    data_type = config["USER"]["data_type"]
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
    mkdir(project_path)
    mkdir(f"{project_path}/shell")
    mkdir(f"{project_path}/log")

    t = [is_chr(v) for k, v in bam_files.items()]
    if all(t):
        chr_str = "chr"
    elif not any(t):
        chr_str = ""
    else:
        raise ValueError("Mixed types of SN tags found.")

    target_region = gene_table[target_gene][f"{genome_build}_region"].replace("chr", "")

    if control_gene != "NONE":
        # Write the shell script for bam2gdf.
        s = f"pypgx bam2gdf {genome_build} {target_gene} {control_gene} \\\n"

        for sample_id in bam_files:
            s += f"  {bam_files[sample_id]} \\\n"

        s += f"  -o {project_path}/{output_prefix}.gdf\n"

        with open(f"{project_path}/shell/bam2gdf.sh", "w") as f:
            f.write(s)

    # Write the shell script for bam2vcf.
    s = f"pypgx bam2vcf {genome_build} {target_gene} {fasta_file} \\\n"

    for sample_id in bam_files:
        s += f"  {bam_files[sample_id]} \\\n"    

    s += f"  -o {project_path}/{output_prefix}.vcf\n"

    with open(f"{project_path}/shell/bam2vcf.sh", "w") as f:
        f.write(s)

    # Write the shell script for Stargazer.
    s = (
        f"project={project_path}\n"
        f"vcf=$project/{output_prefix}.vcf\n"
        "\n"
        f"stargazer \\\n"
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

    with open(f"{project_path}/shell/stargazer.sh", "w") as f:
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
        s += f"{q} -N bam2gdf $p/shell/bam2gdf.sh\n"

    s += f"{q} -N bam2vcf $p/shell/bam2vcf.sh\n"

    if control_gene == "NONE":
        s += f"{q} -hold_jid bam2vcf -N stargazer $p/shell/stargazer.sh\n"

    else:
        s += f"{q} -hold_jid bam2gdf,bam2vcf -N stargazer $p/shell/stargazer.sh\n"

    with open(f"{project_path}/example-qsub.sh", "w") as f:
        f.write(s)
