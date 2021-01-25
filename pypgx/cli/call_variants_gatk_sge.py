import os
import shutil
import sys

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
    bam_files = []
    with open(bam_path) as f:
        for line in f:
            bam_files.append(line.strip())

    try:
        shutil.rmtree(output_dir)
    except OSError:
        pass

    os.mkdir(output_dir)
    os.mkdir(f"{output_dir}/shell")
    os.mkdir(f"{output_dir}/log")
    os.mkdir(f"{output_dir}/temp")

    with open(f"{output_dir}/command-line.txt", 'w') as f:
        f.write(" ".join(sys.argv))

    sn_tags = []

    for bam_file in bam_files:
        sn_tags += get_sn_tags(bam_file)

    if any(["chr" in x for x in list(set(sn_tags))]):
        chr = "chr"
    else:
        chr = ""

    target_locus = Locus.from_input(target_gene, genome_build)

    # Write the shell script for HaplotypeCaller.
    for i, bam_file in enumerate(bam_files):
        s = "#!/bin/bash\n"

        if conda_env is not None:
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
            f"  -O {output_dir}/temp/{i}.g.vcf \\\n"
            f"  -L {target_locus.region} \\\n"
            "  --QUIET \\\n"
        )

        if java_options is not None:
            s += f"  --java-options {java_options} \\\n"

        with open(f"{output_dir}/shell/haplotypecaller-{i}.sh", "w") as f:
            f.write(s)


    # Write the shell script for post-HaplotypeCaller.
    s = "#!/bin/bash\n"

    if conda_env is not None:
        s += (
            "\n"
            f"conda activate {conda_env}\n"
        )

    s += (
        "\n"
        f"p={output_dir}\n"
        "\n"
        "gatk GenomicsDBImport \\\n"
        f"  --intervals {target_locus.region} \\\n"
        f"  --genomicsdb-workspace-path $p/temp/datastore \\\n"
        "  --merge-input-intervals \\\n"
        "  --QUIET \\\n"
    )

    if java_options is not None:
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

    if java_options is not None:
        s += f"  --java-options {java_options} \\\n"

    if dbsnp_file is not None:
        s += f"  -D {dbsnp_file} \\\n"

    s += (
        "\n"
        "gatk VariantFiltration \\\n"
        f"  -R {fasta_file} \\\n"
        f"  -L {target_locus.region} \\\n"
        f"  -O $p/pypgx.vcf \\\n"
        f"  --variant $p/temp/pypgx.joint.vcf \\\n"
        "  --filter-expression 'QUAL <= 50.0' \\\n"
        "  --filter-name QUALFilter \\\n"
        "  --QUIET\n"
    )

    if java_options is not None:
        s += f"  --java-options {java_options} \\\n"

    with open(f"{output_dir}/shell/post-haplotypecaller.sh", "w") as f:
        f.write(s)

    # Write the shell script for qsub.
    q = "qsub -e $p/log -o $p/log"

    if qsub_options is not None:
        q += f" {qsub_options}"

    s = (
        "#!/bin/bash\n"
        "\n"
        f"p={output_dir}\n"
        f"j={randstr()}\n"
        "\n"
    )

    for i, bam_file in enumerate(bam_files):
        s += f"{q} -N $j-hc $p/shell/haplotypecaller-{i}.sh\n"

    s += f"{q} -hold_jid $j-hc -N $j-post-hc $p/shell/post-haplotypecaller.sh\n"

    with open(f"{output_dir}/example-qsub.sh", "w") as f:
        f.write(s)
