import os
import shutil
import sys
from pathlib import Path

from pypgx import get_sn_tags, Locus, randstr

def call_variants_gatk_sge(target_gene, bam_path, fasta_file, output_dir,
                           genome_build="hg19", dbsnp_file=None,
                           java_options=None, qsub_options=None,
                           conda_env=None):
    """Create a VCF (Variant Call Format) file for Stargazer from BAM files
    by calling SNVs and indels.

    Parameters
    ----------
    target_gene : str
        Name of the target gene. Choices: {'abcb1', 'cacna1s',
        'cftr', 'cyp1a1', 'cyp1a2', 'cyp1b1', 'cyp2a6',
        'cyp2a13', 'cyp2b6', 'cyp2c8', 'cyp2c9', 'cyp2c19',
        'cyp2d6', 'cyp2e1', 'cyp2f1', 'cyp2j2', 'cyp2r1',
        'cyp2s1', 'cyp2w1', 'cyp3a4', 'cyp3a5', 'cyp3a7',
        'cyp3a43', 'cyp4a11', 'cyp4a22', 'cyp4b1', 'cyp4f2',
        'cyp17a1', 'cyp19a1', 'cyp26a1', 'dpyd', 'g6pd',
        'gstm1', 'gstp1', 'gstt1', 'ifnl3', 'nat1', 'nat2',
        'nudt15', 'por', 'ptgis', 'ryr1', 'slc15a2',
        'slc22a2', 'slco1b1', 'slco1b3', 'slco2b1', 'sult1a1',
        'tbxas1', 'tpmt', 'ugt1a1', 'ugt1a4', 'ugt2b7',
        'ugt2b15', 'ugt2b17', 'vkorc1', 'xpc'}.    bam_path : str
    bam_path : str
        Read BAM files from PATH, one file path per line.
    fasta_file : str
        Path to a reference FASTA file.
    output_dir : str
        Path to the output directory.
    genome_build : str, default: 'hg19'
        Build of the reference genome assembly. Choices:
        {'hg19', 'hg38'}.
    dbsnp_file : str, optional
        Path to a dbSNP file (.vcf or .vcf.gz). Used to assign rs ID to
        observed variants.
    java_options : str, optional
        Options passed to Java to run GATK.
    qsub_options : str, optional
        Options passed to SGE.
    conda_env : str, optional
        Name of the conda environment to be activated when the jobs
        are submitted to SGE.

    """
    # Get absolute path for the input files.
    _fasta_file = Path(fasta_file).resolve()
    _output_dir = Path(output_dir).resolve()
    if dbsnp_file is None:
        _dbsnp_file = None
    else:
        _dbsnp_file = Path(dbsnp_file).resolve()

    # Get the input BAM files.
    bam_files = []
    with open(bam_path) as f:
        for line in f:
            bam_files.append(line.strip())

    # Make the output directoires.
    try:
        shutil.rmtree(_output_dir)
    except OSError:
        pass
    os.mkdir(_output_dir)
    os.mkdir(f"{_output_dir}/shell")
    os.mkdir(f"{_output_dir}/log")
    os.mkdir(f"{_output_dir}/temp")

    # Record the command line.
    with open(f"{_output_dir}/command-line.txt", 'w') as f:
        f.write(' '.join(sys.argv))

    # Determine whether to include the 'chr' string in chromosome name.
    sn_tags = []
    for bam_file in bam_files:
        sn_tags += get_sn_tags(bam_file)
    if any(["chr" in x for x in list(set(sn_tags))]):
        chr = "chr"
    else:
        chr = ""

    # Parse the target locus.
    target_locus = Locus.from_input(target_gene, genome_build)

    # Define the conda environment line.
    conda_env_line = f"conda activate {conda_env}\n"
    if conda_env is None:
        conda_env_line = f"# {conda_env_line}"

    # Define the Java options line.
    java_options_line = f"--java-options {java_options}\n"
    if java_options is None:
        java_options_line = f"# {java_options_line}"

    # Define the dbSNP files line.
    dbsnp_file_line = f"-D {_dbsnp_file}\n"
    if _dbsnp_file is None:
        dbsnp_file_line = f"# {dbsnp_file_line}"

    # Define the gVCF file lines.
    gvcf_file_lines = ''
    for i, bam_file in enumerate(bam_files):
        gvcf_file_lines += f"-V $p/temp/{i}.g.vcf\n"

    # Write the shell script for HaplotypeCaller.
    for i, bam_file in enumerate(bam_files):
        with open(f"{_output_dir}/shell/hc-{i}.sh", 'w') as f:
            f.write(("#!/bin/bash\n"
                     '\n'
                     f"{conda_env_line}"
                     '\n'
                     "arguments=(\n"
                     f"-R {_fasta_file}\n"
                     f"--emit-ref-confidence GVCF\n"
                     f"-I {bam_file}\n"
                     f"-O {_output_dir}/temp/{i}.g.vcf\n"
                     f"-L {chr}{target_locus.region}\n"
                     "--QUIET\n"
                     f"{java_options_line}"
                     ")\n"
                     '\n'
                     "gatk HaplotypeCaller ${arguments[@]}\n"))

    # Write the shell script for post-HaplotypeCaller.
    with open(f"{_output_dir}/shell/post-hc.sh", "w") as f:
        f.write(("#!/bin/bash\n"
                 '\n'
                 f"{conda_env_line}"
                 '\n'
                 f"p={_output_dir}\n"
                 '\n'
                 "arguments=(\n"
                 f"--intervals {chr}{target_locus.region}\n"
                 f"--genomicsdb-workspace-path $p/temp/datastore\n"
                 "--merge-input-intervals\n"
                 "--QUIET\n"
                 f"{java_options_line}"
                 f"{gvcf_file_lines}"
                 ")\n"
                 '\n'
                 "gatk GenomicsDBImport ${arguments[@]}\n"
                 '\n'
                 "arguments=(\n"
                 f"-R {_fasta_file}\n"
                 f"-V gendb://$p/temp/datastore\n"
                 f"-O $p/temp/pypgx.joint.vcf\n"
                 "--QUIET\n"
                 f"{java_options_line}"
                 f"{dbsnp_file_line}"
                 ")\n"
                 '\n'
                 "gatk GenotypeGVCFs ${arguments[@]}\n"
                 '\n'
                 "arguments=(\n"
                 f"-R {_fasta_file}\n"
                 f"-L {chr}{target_locus.region}\n"
                 f"-O $p/pypgx.vcf \\\n"
                 f"--variant $p/temp/pypgx.joint.vcf\n"
                 "--filter-expression 'QUAL <= 50.0'\n"
                 "--filter-name QUALFilter\n"
                 "--QUIET\n"
                 f"{java_options_line}"
                 ")\n"
                 '\n'
                 "gatk VariantFiltration ${arguments[@]}\n"))

    # Write the shell script for qsub.
    q = "qsub -e $p/log -o $p/log"

    if qsub_options is not None:
        q += f" {qsub_options}"

    s = (
        "#!/bin/bash\n"
        '\n'
        f"p={_output_dir}\n"
        f"j={randstr()}\n"
        '\n'
    )

    for i, bam_file in enumerate(bam_files):
        s += f"{q} -N $j-hc $p/shell/hc-{i}.sh\n"

    s += f"{q} -hold_jid $j-hc -N $j-post-hc $p/shell/post-hc.sh\n"

    with open(f"{_output_dir}/example-qsub.sh", "w") as f:
        f.write(s)
