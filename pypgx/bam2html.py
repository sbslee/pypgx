import configparser
from os import mkdir
from os.path import realpath
from .common import conf_env, get_target_genes, randstr

@conf_env
def bam2html(conf_file: str, **kwargs) -> None:
    """Run per-sample genotyping for multiple genes with SGE.

    This command runs the per-sample genotyping pipeline by submitting 
    jobs to the Sun Grid Engine (SGE) cluster. This essentially deploys 
    the ``genotype`` command to multiple genes in parallel. After genotype 
    analysis is complete, it will merge the genotype results and then 
    generate a HTML report using the ``report`` command.

    Args:
        conf (str): Configuration file.

    .. note::

        BCFtools, SGE and Stargazer must be pre-installed.

    This is what a typical configuration file for ``sges`` looks like:

        .. code-block:: python

            # File: example_conf.txt
            # To execute:
            #   $ pypgx sges example_conf.txt
            #   $ sh ./myproject/example-qsub.sh

            # Do not make any changes to this section.
            [DEFAULT]
            target_genes = ALL
            control_gene = NONE
            qsub_options = NONE

            # Make any necessary changes to this section.
            [USER]
            snp_caller = gatk
            fasta_file = reference.fa
            project_path = ./myproject
            genome_build = hg19
            data_type = wgs
            bam_file = in.bam
            qsub_options = -l mem_requested=2G
            target_genes = cyp2b6, cyp2d6
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
             - Name or region of control gene 
               (e.g. 'vdr', 'chr12:48232319-48301814').
           * - data_type
             - Data type ('wgs' or 'ts').
           * - fasta_file
             - Reference FASTA file.
           * - genome_build
             - Genome build ('hg19' or 'hg38').
           * - project_path
             - Output project directory.
           * - qsub_options
             - Options for qsub command (e.g. '-l mem_requested=2G').
           * - target_genes
             - Names of target genes (e.g. 'cyp2d6')..
    """
    config = kwargs["config"]

    # Parse the configuration data.
    snp_caller = config["USER"]["snp_caller"]
    project_path = realpath(config["USER"]["project_path"])
    bam_file = realpath(config["USER"]["bam_file"])
    fasta_file = realpath(config["USER"]["fasta_file"])
    target_genes = config["USER"]["target_genes"]
    control_gene = config["USER"]["control_gene"]
    data_type = config["USER"]["data_type"]
    genome_build = config["USER"]["genome_build"]
    qsub_options = config["USER"]["qsub_options"]

    t = get_target_genes()
    
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

    # Write the shell script for genotyping pipeline.
    for select_gene in select_genes:
        s = (
            "pypgx bam2gt \\\n"
            f"  {snp_caller} \\\n"
            f"  {fasta_file} \\\n"
            f"  {select_gene} \\\n"
            f"  {genome_build} \\\n"
            f"  {data_type} \\\n"
            f"  {project_path}/gene/{select_gene}/stargazer \\\n"
            f"  {bam_file} \\\n"
        )

        if control_gene != "NONE":
            s += f"  --control_gene {control_gene}\n"

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
        "pypgx gt2html $p/genotype.merged.txt -o $p/report.html"
    )

    with open(f"{project_path}/report.sh", "w") as f:
        f.write(s)

    # Write the shell script for qsub.
    if qsub_options == "NONE":
        q = "qsub"
    else:
        q = f"qsub {qsub_options}"

    s = (
        "#!/bin/bash\n"
        "\n"
        f"p={project_path}\n"
        "\n"
    )

    jid = randstr()

    for select_gene in select_genes:
        s += f"{q} -o $p/gene/{select_gene} -e $p/gene/{select_gene} -N {jid} $p/gene/{select_gene}/run.sh\n"

    s += (
        f"{q} -o $p -e $p -hold_jid {jid} $p/report.sh"
    )

    with open(f"{project_path}/example-qsub.sh", "w") as f:
        f.write(s)
