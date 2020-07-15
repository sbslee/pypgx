import configparser
from os import mkdir
from os.path import realpath

from .common import logging, LINE_BREAK1

logger = logging.getLogger(__name__)

def remap(conf: str) -> None:
    """
    Remap BAM file(s) to different reference.

    Args:
        conf (str): Configuration file.

    This is what a typical configuration file for ``remap`` looks like:

        .. code-block:: python

            # File: example_conf.txt
            # Do not make any changes to this section.
            [DEFAULT]
            platform = illumina
            threads = 1
            java_heap = -Xmx2g
            qsub_options1 = NONE
            qsub_options2 = NONE

            # Make any necessary changes to this section.
            [USER]
            fasta_file = reference.fa
            manifest_file = manifest.txt
            project_path = /path/to/project/
            vcf_files = in1.vcf, in2.vcf, in3.vcf
            library = awesome_experiment
            gatk_tool = GenomeAnalysisTK.jar
            picard_tool = picard.jar
            qsub_options1 = -q nick-grad.q -l mem_requested=2G -pe serial 1
            qsub_options2 = -q nick-grad.q -l mem_requested=2G

    This table summarizes the configuration parameters specific to ``remap``:

        .. list-table::
           :widths: 25 75
           :header-rows: 1

           * - Parameter
             - Summary
           * - fasta_file
             - Reference FASTA file.
           * - gatk_tool
             - GATK program.
           * - java_heap
             - Java heap size.
           * - library
             - Sequencing library name.
           * - manifest_file
             - Manifest file.
           * - picard_tool
             - Picard program.
           * - platform
             - Sequencing platform.
           * - project_path
             - Output project directory.
           * - qsub_options1
             - Options for the first qsub command. Recommended to set a parallel environment.
           * - qsub_options2
             - Options for the second qsub command.
           * - threads
             - Number of threads.
           * - vcf_files
             - Reference VCF files used for base quality score recalibration.
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
    gatk_tool = realpath(config["USER"]["gatk_tool"])
    picard_tool = realpath(config["USER"]["picard_tool"])
    threads = config["USER"]["threads"]
    platform = config["USER"]["platform"]
    java_heap = config["USER"]["java_heap"]
    library = config["USER"]["library"]
    qsub_options1 = config["USER"]["qsub_options1"]
    qsub_options2 = config["USER"]["qsub_options2"]

    vcf_files = []
    for vcf_file in config["USER"]["vcf_files"].split(","):
        vcf_files.append(realpath(vcf_file.strip()))

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
    mkdir(f"{project_path}/bam")
    mkdir(f"{project_path}/log")
    mkdir(f"{project_path}/temp")
    mkdir(f"{project_path}/fastq")

    # Write the first qsub script.
    for id in bam_files:
        s = (
            "#!/bin/bash\n"
            "\n"
            f"name={id}\n"
            f"project={project_path}\n"
            f"threads={threads}\n"
            f"bam1={bam_files[id]}\n"
            f"bam2=$project/temp/$name.collated.bam\n"
            f"bam3=$project/temp/$name.sorted.bam\n"
            f"fastq=$project/fastq/$name.fq\n"
            f"fasta={fasta_file}\n"
            "\n"
            "# Collate the input BAM file by read name.\n"
            "samtools collate -@ $threads $bam1 -o $bam2\n"
            "\n"
            "# Convert the new BAM file to a FASTQ file.\n"
            "samtools fastq -0 /dev/null $bam2 > $fastq\n"
            "\n"
            "# Get the read group.\n"
            "read_group1=`samtools view -H $bam1 | grep -m 1 '^@RG'`\n"
            "id_field=`echo $read_group1 | awk '{for (i=1; i<=NF; i++) "
                "{if ($i ~ /ID/) {print $i}}}' | sed 's/ID://g'`\n"
            "pu_field=`echo $read_group1 | awk '{for (i=1; i<=NF; i++) "
                "{if ($i ~ /PU/) {print $i}}}' | sed 's/PU://g'`\n"
            f"platform={platform}\n"
            f"library={library}\n"
            'read_group2="@RG\\tID:$id_field\\tPU:$pu_field\\tSM:$name'
                '\\tPL:$platform\\tLB:$library"\n'
            "\n"
           "# Align the sequence reads.\n"
           "bwa mem -M -t $threads -R $read_group2 -p $fasta $fastq | "
                "samtools sort -@ $threads -o $bam3 -\n"
            )

        with open(f"{project_path}/shell/run-{id}-1.sh", "w") as f:
            f.write(s)


    # Write the second qsub script.
    for id in bam_files:
        s = (
            "#!/bin/bash\n"
            "\n"
            f"name={id}\n"
            f"project={project_path}\n"
            "\n"
            "# Mark duplicate reads.\n"
            f"bam1=$project/temp/$name.sorted.bam\n"
            f"bam2=$project/temp/$name.sorted.markeddups.bam\n"
            f"metrics=$project/temp/$name.metrics\n"
            f"picard={picard_tool}\n"
            f"java {java_heap} -jar $picard MarkDuplicates \\\n"
            "  I=$bam1 \\\n"
            "  M=$metrics \\\n"
            "  O=$bam2 \\\n"
            "  ASSUME_SORTED=true\n"
            "\n"
            "# Index the resulting BAM file.\n"
            "samtools index $bam2\n"
            "\n"
            "# Build the BQSR model.\n"
        )

        for i, vcf_file in enumerate(vcf_files, 1):
            s += f"vcf{i}={vcf_file}\n"

        s += (
            f"fasta={fasta_file}\n"
            f"bqsr=$project/temp/$name.table\n"
            f"gatk={gatk_tool}\n"
            f"java {java_heap} -jar $gatk -T BaseRecalibrator \\\n"
            "  -I $bam2 \\\n"
            "  -R $fasta \\\n"
        )

        for i in range(len(vcf_files)):
            s += f"  --knownSites $vcf{i + 1} \\\n"

        s += (
            "  -o $bqsr\n"
            "\n"
            "# Apply the BQSR model.\n"
            f"bam3=$project/bam/$name.sorted.markeddups.recal.bam\n"
            f"java {java_heap} -jar $gatk -T PrintReads \\\n"
            "  -R $fasta \\\n"
            "  -I $bam2 \\\n"
            "  -o $bam3 \\\n"
            "  -BQSR $bqsr\n"
        )

        with open(f"{project_path}/shell/run-{id}-2.sh", "w") as f:
            f.write(s)

    # Write the shell script for qsub.
    s = f"p={project_path}\n"

    for sample_id in bam_files:
        q1 = f"qsub -e $p/log -o $p/log"

        if qsub_options1 != "NONE":
            q1 += f" {qsub_options1}"

        q2 = f"qsub -e $p/log -o $p/log"

        if qsub_options2 != "NONE":
            q2 += f" {qsub_options2}"

        s += (
            "\n"
            f"n1=run-{sample_id}-1\n"
            f"n2=run-{sample_id}-2\n"
            f"{q1} -N $n1 $p/shell/$n1.sh\n"
            f"{q2} -N $n2 -hold_jid $n1 $p/shell/$n2.sh\n"
        )

    with open(f"{project_path}/example-qsub.sh", "w") as f:
        f.write(s)
