import configparser
import os

from .common import logging, LINE_BREAK1

logger = logging.getLogger(__name__)

def remap(conf: str) -> None:
    """
    Remap BAM file(s) to different reference.

    Args:
        conf (str): Configuration file.

    Examples:
        .. code-block:: python

            # example_conf.txt
            # Do not make any changes to this section.
            [DEFAULT]
            output_prefix = pypgx
            output_directory = CURRENT
            platform = illumina
            thread = 1
            heap = -Xmx8g
            resource = mem_requested=2G

            # Make any necessary changes to this section.
            [USER]
            fasta = reference.fa
            manifest = manifest.txt
            gatk = gatk.jar
            picard = picard.jar
            vcf = in1.vcf, in2.vcf, in3.vcf
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

    manifest = config["USER"]["manifest"]
    output_prefix = config["USER"].get("output_prefix")
    
    if config["USER"].get("output_directory") == "CURRENT":
        output_directory = os.getcwd()
    else:
        output_directory = config["USER"].get("output_directory")

    thread = config["USER"]["thread"]
    fasta = config["USER"]["fasta"]
    platform = config["USER"]["platform"]
    picard = config["USER"]["picard"]
    heap = config["USER"]["heap"]
    gatk = config["USER"]["gatk"]
    resource = config["USER"]["resource"]

    vcf_files = []
    for vcf_file in config["USER"]["vcf"].split(","):
        vcf_files.append(vcf_file.strip())

    # Read the manifest file.
    bams = {}
    with open(manifest) as f:
        header = next(f).strip().split("\t")
        i1 = header.index("sample_id")
        i2 = header.index("bam")
        for line in f:
            fields = line.strip().split("\t")
            id = fields[i1]
            bam = fields[i2]
            bams[id] = bam

    # Log the number of samples.
    logger.info(f"Number of samples: {len(bams)}")

    # Make the project directories.
    project_directory = f"{output_directory}/{output_prefix}-project"
    os.mkdir(project_directory)
    os.mkdir(f"{project_directory}/shell")
    os.mkdir(f"{project_directory}/bam")
    os.mkdir(f"{project_directory}/log")
    os.mkdir(f"{project_directory}/temporary")
    os.mkdir(f"{project_directory}/fastq")

    # Write the first qsub script.
    for id in bams:
        with open(f"{project_directory}/shell/run-{id}-1.sh", "w") as f:
            f.write((
                "#!/bin/bash\n"
                "\n"
                f"name={id}\n"
                f"project_dir={project_directory}\n"
                f"thread={thread}\n"
                f"bam1={bams[id]}\n"
                f"bam2=$project_dir/temp/$name.collated.bam\n"
                f"bam3=$project_dir/temp/$name.sorted.bam\n"
                f"fastq=$project_dir/fastq/$name.fq\n"
                f"fasta={fasta}\n"
                "\n"
                "# Collate the input BAM file by read name.\n"
                "samtools collate -@ $thread $bam1 -o $bam2\n"
                "\n"
                "# Convert the new BAM file to a FASTQ file.\n"
                "samtools fastq -0 /dev/null $bam2 > $fastq\n"
                "\n"
                "# Get the read group.\n"
                "read_group1=`samtools view -H $bam1 | grep -m 1 '^@RG'`\n"
                "id_field=`echo $read_group1 | awk '{for (i=1; i<=NF; i++)"
                    "{if ($i ~ /ID/) {print $i}}}' | sed 's/ID://g'`\n"
                "pu_field=`echo $read_group1 | awk '{for (i=1; i<=NF; i++) "
                    "{if ($i ~ /PU/) {print $i}}}' | sed 's/PU://g'`\n"
                f"platform={platform}\n"
                f"library={output_prefix}\n"
                "read_group2='@RG\\tID:$id_field\\tPU:$pu_field\\tSM:$name"
                    "\\tPL:$platform\\tLB:$library'\n"
                "\n"
                "# Align the sequence reads.\n"
                "bwa mem -M -t $thread -R $read_group2 -p $fasta $fastq | "
                    "samtools sort -@ $thread -o $bam3 -\n"
            ))

    # Write the second qsub script.
    for id in bams:
        temp = (
            "#!/bin/bash\n"
            "\n"
            f"name={id}\n"
            f"project_dir={project_directory}\n"
            "\n"
            "# Mark duplicate reads.\n"
            f"bam1=$project_dir/temp/$name.sorted.bam\n"
            f"bam2=$project_dir/temp/$name.sorted.markeddups.bam\n"
            f"metrics=$project_dir/temp/$name.metrics\n"
            f"picard=\n{picard}"
            f"heap={heap}\n"
            "java $heap -jar $picard MarkDuplicates \\\n"
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
            temp += f"vcf{i}={vcf_file}\n"

        temp += (
            f"fasta={fasta}\n"
            f"bqsr=$project_dir/temp/$name.table\n"
            f"gatk={gatk}\n"
            f"java $heap -jar $gatk -T BaseRecalibrator \\\n"
            "  -I $bam2 \\\n"
            "  -R $fasta \\\n"
        )

        for i in range(len(vcf_files)):
            temp += f"  --knownSites $vcf{i + 1} \\\n"

        temp += (
            "  -o $bqsr\n"
            "\n"
            "# Apply the BQSR model.\n"
            f"bam3=$project_dir/bam/$name.sorted.markeddups.recal.bam\n"
            "java $heap -jar $gatk -T PrintReads \\\n"
            "  -R $fasta \\\n"
            "  -I $bam2 \\\n"
            "  -o $bam3 \\\n"
            "  -BQSR $bqsr\n"
        )

        with open(f"{project_directory}/shell/run-{id}-2.sh", "w") as f:
            f.write(temp)

    with open(f"{project_directory}/example-qsub.sh", "w") as f:
        f.write(f"log_directory={project_directory}/log\n")
        f.write(f"shell_directory={project_directory}/shell\n")
        f.write(f"resource={resource}\n")
        f.write("\n")
        for id in bams:
            f.write((
                "qsub -q nick-grad.q -e $log_directory -o $log_directory "
                f"-N run-{id}-1 -l $resource -pe serial {thread} "
                f"$shell_directory/run-{id}-1.sh\n"
            ))
            f.write((
                "qsub -q nick-grad.q -e $log_directory -o $log_directory "
                f"-N run-{id}-2 -l $resource -hold_jid run-{id}-1 "
                f"$shell_directory/run-{id}-2.sh\n"
                "\n"
            ))
