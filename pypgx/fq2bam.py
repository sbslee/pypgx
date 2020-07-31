from os import mkdir
from os.path import realpath
from .common import conf_env

@conf_env
def fq2bam(conf_file: str, **kwargs) -> None:
    """
    Create BAM file(s) from FASTQ file(s).

    Args:
        conf_file (str): Configuration file.

    This is what a typical configuration file for ``fq2bam`` looks like:

        .. code-block:: python

            # File: example_conf.txt
            # Do not make any changes to this section.
            [DEFAULT]
            platform = illumina
            threads = 1
            read_length = 150
            qsub_options1 = NONE
            qsub_options2 = NONE

            # Make any necessary changes to this section.
            [USER]
            fasta_file = reference.fa
            manifest_file = manifest.txt
            project_path = /path/to/project/
            vcf_files = in1.vcf, in2.vcf, in3.vcf
            library = awesome_experiment
            bed_file = in.bed
            threads = 15
            qsub_options1 = -V -q biall.q -S /bin/bash -pe pePAC 15
            qsub_options2 = -V -q biall.q -S /bin/bash

    This table summarizes the configuration parameters specific to ``fq2bam``:

        .. list-table::
           :widths: 25 75
           :header-rows: 1

           * - Parameter
             - Summary
           * - bed_file
             - BED file.
           * - fasta_file
             - Reference FASTA file.
           * - library
             - Sequencing library name.
           * - manifest_file
             - Manifest file.
           * - platform
             - Sequencing platform.
           * - project_path
             - Output project directory.
           * - qsub_options1
             - Options for the first qsub command. Recommended to set a parallel environment.
           * - qsub_options2
             - Options for the second qsub command.
           * - read_length
             - Sequence read length.
           * - threads
             - Number of threads.
           * - vcf_files
             - Reference VCF files used for base quality score recalibration.
    """
    # Parse the configuration data.
    config = kwargs["config"]
    manifest_file = realpath(config["USER"]["manifest_file"])
    project_path = realpath(config["USER"]["project_path"])
    fasta_file = realpath(config["USER"]["fasta_file"])
    bed_file = realpath(config["USER"]["bed_file"])
    platform = config["USER"]["platform"]
    library = config["USER"]["library"]
    threads = config["USER"]["threads"]
    read_length = config["USER"]["read_length"]
    qsub_options1 = config["USER"]["qsub_options1"]
    qsub_options2 = config["USER"]["qsub_options2"]

    vcf_files = []

    for vcf_file in config["USER"]["vcf_files"].split(","):
        vcf_files.append(realpath(vcf_file.strip()))

    # Read the manifest file.
    fastq_files = {}

    with open(manifest_file) as f:
        header = next(f).strip().split("\t")
        i1 = header.index("sample_id")
        i2 = header.index("fastq")
        for line in f:
            fields = line.strip().split("\t")
            sample_id = fields[i1]
            fastq = fields[i2]
            if sample_id not in fastq_files:
                fastq_files[sample_id] = []
            fastq_files[sample_id].append(fastq)
            if len(fastq_files[sample_id]) > 2:
                raise ValueError(
                    f"More than two FASTQ files found: {sample_id}")

    # Log the number of samples.
    logger.info(f"Number of samples: {len(fastq_files)}")

    # Make the project directories.
    mkdir(project_path)
    mkdir(f"{project_path}/shell")
    mkdir(f"{project_path}/bam")
    mkdir(f"{project_path}/log")
    mkdir(f"{project_path}/temp")
    mkdir(f"{project_path}/stat")

    # Write the first shell script.
    for sample_id in fastq_files:
        s = (
            "#!/bin/bash\n"
            "\n"
            f"name={sample_id}\n"
            f"read1={fastq_files[sample_id][0]}\n"
            f"read2={fastq_files[sample_id][1]}\n"
            f"threads={threads}\n"
            f"fasta={fasta_file}\n"
            f"platform={platform}\n"
            f"library={library}\n"
            f"bam1={project_path}/temp/$name.sorted.bam\n"
            "\n"
            "# Get the read group information.\n"
            "first=`zcat $read1 | head -1`\n"
            '''flowcell=`echo "$first" | awk -F " " '{print $1}' | '''
                '''awk -F ":" '{print $3}'`\n'''
            '''barcode=`echo "$first" | awk -F " " '{print $2}' | '''
                '''awk -F ":" '{print $4}'`\n'''
            f'''group="@RG\\tID:$flowcell\\tPU:$flowcell.$barcode'''
                '''\\tSM:$name\\tPL:$platform\\tLB:$library"\n'''
            "\n"
            "# Align and sort the seuqnece reads. Note that this step \n"
            "# will also assign a read group to the mapped reads.\n"
            "bwa mem -M -R $group -t $threads $fasta $read1 $read2 | "
                "samtools sort -@ $threads -o $bam1 -\n"
        )

        with open(f"{project_path}/shell/run-{sample_id}-1.sh", "w") as f:
            f.write(s)

    # Write the second shell script.
    for sample_id in fastq_files:
        s = (
            "#!/bin/bash\n"
            "\n"
            f"name={sample_id}\n"
            f"length={read_length}\n"
            f"fasta={fasta_file}\n"
            f"bed={bed_file}\n"
            f"project={project_path}\n"
            "bam1=$project/temp/$name.sorted.bam\n"
            "bam2=$project/temp/$name.sorted.markeddups.bam\n"
            "bam3=$project/bam/$name.sorted.markeddups.recal.bam\n"
            "table=$project/temp/$name.table\n"
            "metrics=$project/temp/$name.metrics\n"
            "stat=$project/stat/$name.txt\n"
        )

        for i, vcf_file in enumerate(vcf_files, 1):
            s += f"vcf{i}={vcf_file}\n"

        s += (
            "\n"
            "# Mark the duplicate reads.\n"
            f"gatk MarkDuplicates -I $bam1 -M $metrics \\\n"
            "  --ASSUME_SORT_ORDER coordinate \\\n"
            "  -O $bam2\n"
            "\n"
            "# Index the BAM file.\n"
            "samtools index $bam2\n"
            "\n"
            "# Build the BQSR model.\n"
            f"gatk BaseRecalibrator -I $bam2 \\\n"
            "  --intervals $bed \\\n"
        )

        for i, vcf_file in enumerate(vcf_files, 1):
            s += f"  --known-sites $vcf{i} \\\n"

        s += (
            "  -R $fasta \\\n"
            "  -O $table\n"
            "\n"
            "# Apply the BQSR model. Note that this step will remove any\n"
            "# sequence reads that were mapped outside of targeted regions.\n"
            f"gatk ApplyBQSR -I $bam2 -O $bam3 \\\n"
            "  --intervals $bed \\\n"
            "  -bqsr $table\n"
            "\n"
            "# Get the stats.\n"
            "flagstat1=`samtools flagstat $bam2`\n"
            "flagstat2=`samtools flagstat $bam3`\n"
            '''total_n=`echo "$flagstat1" | awk 'NR == 1 {print $1}'`\n'''
            '''mapped_n=`echo "$flagstat1" | awk 'NR == 5 {print $1}'`\n'''
            '''mapped_p=`echo "$flagstat1" | awk 'NR == 5 {print $5}' | '''
                '''sed "s/(//g" | sed "s/%//g"`\n'''
            '''unique_n=`samtools view -q 1 -c $bam2`\n'''
            '''unique_p=`bc <<< "scale = 2; $unique_n * 100 / $total_n"`\n'''
            '''dup_n=`echo "$flagstat1" | awk 'NR == 4 {print $1}'`\n'''
            '''dup_p=`bc <<< "scale = 2; $dup_n * 100 / $total_n"`\n'''
            '''total_b=`echo "$(($total_n * $length))"`\n'''
            '''panel_b=`awk -F '\\t' 'BEGIN {sum = 0} {sum += $3 - $2} '''
                '''END {print sum}' $bed`\n'''
            '''total_x=$(($total_b / $panel_b))\n'''
            '''target_n=`echo "$flagstat2" | awk 'NR == 1 {print $1}'`\n'''
            '''target_p=`bc <<< "scale = 2; $target_n * 100 / $total_n"`\n'''
            '''target_b=`echo "$(($target_n * $length))"`\n'''
            '''target_x=`echo "$(($target_b / $panel_b))"`\n'''
            '''mean_x=`samtools depth $bam3 -b $bed | awk '{sum += $3} '''
                '''END {print sum / NR}' | sed 's/\.[^\.]*$//'`\n'''
            '''x0001_p=`printf "%.2f\\n" <<< echo "$(samtools depth '''
                '''-b $bed $bam3 | awk '{if ($3 >= 1) c++} '''
                '''END {print c / NR * 100}')"`\n'''
            '''x0010_p=`printf "%.2f\\n" <<< echo "$(samtools depth '''
                '''-b $bed $bam3 | awk '{if ($3 >= 10) c++} '''
                '''END {print c / NR * 100}')"`\n'''
            '''x0020_p=`printf "%.2f\\n" <<< echo "$(samtools depth '''
                '''-b $bed $bam3 | awk '{if ($3 >= 20) c++} '''
                '''END {print c / NR * 100}')"`\n'''
            '''x0050_p=`printf "%.2f\\n" <<< echo "$(samtools depth '''
                '''-b $bed $bam3 | awk '{if ($3 >= 50) c++} '''
                '''END {print c / NR * 100}')"`\n'''
            '''x0100_p=`printf "%.2f\\n" <<< echo "$(samtools depth '''
                '''-b $bed $bam3 | awk '{if ($3 >= 100) c++} '''
                '''END {print c / NR * 100}')"`\n'''
            '''x0200_p=`printf "%.2f\\n" <<< echo "$(samtools depth '''
                '''-b $bed $bam3 | awk '{if ($3 >= 200) c++} '''
                '''END {print c / NR * 100}')"`\n'''
            '''x0300_p=`printf "%.2f\\n" <<< echo "$(samtools depth '''
                '''-b $bed $bam3 | awk '{if ($3 >= 300) c++} '''
                '''END {print c / NR * 100}')"`\n'''
            '''x0400_p=`printf "%.2f\\n" <<< echo "$(samtools depth '''
                '''-b $bed $bam3 | awk '{if ($3 >= 400) c++} '''
                '''END {print c / NR * 100}')"`\n'''
            '''x0500_p=`printf "%.2f\\n" <<< echo "$(samtools depth '''
                '''-b $bed $bam3 | awk '{if ($3 >= 500) c++} '''
                '''END {print c / NR * 100}')"`\n'''
            '''x1000_p=`printf "%.2f\\n" <<< echo "$(samtools depth '''
                '''-b $bed $bam3 | awk '{if ($3 >= 1000) c++} '''
                '''END {print c / NR * 100}')"`\n'''
            "\n"
            "# Report the stats.\n"
            '''echo "total_n $total_n" > $stat\n'''
            '''echo "mapped_n $mapped_n" >> $stat\n'''
            '''echo "mapped_p $mapped_p" >> $stat\n'''
            '''echo "unique_n $unique_n" >> $stat\n'''
            '''echo "unique_p $unique_p" >> $stat\n'''
            '''echo "dup_n $dup_n" >> $stat\n'''
            '''echo "dup_p $dup_p" >> $stat\n'''
            '''echo "total_b $total_b" >> $stat\n'''
            '''echo "panel_b $panel_b" >> $stat\n'''
            '''echo "total_x $total_x" >> $stat\n'''
            '''echo "target_n $target_n" >> $stat\n'''
            '''echo "target_p $target_p" >> $stat\n'''
            '''echo "target_b $target_b" >> $stat\n'''
            '''echo "target_x $target_x" >> $stat\n'''
            '''echo "mean_x $mean_x" >> $stat\n'''
            '''echo "x0001_p $x0001_p" >> $stat\n'''
            '''echo "x0010_p $x0010_p" >> $stat\n'''
            '''echo "x0020_p $x0020_p" >> $stat\n'''
            '''echo "x0050_p $x0020_p" >> $stat\n'''
            '''echo "x0100_p $x0100_p" >> $stat\n'''
            '''echo "x0200_p $x0200_p" >> $stat\n'''
            '''echo "x0300_p $x0300_p" >> $stat\n'''
            '''echo "x0400_p $x0400_p" >> $stat\n'''
            '''echo "x0500_p $x0500_p" >> $stat\n'''
            '''echo "x1000_p $x1000_p" >> $stat\n'''
        )

        with open(f"{project_path}/shell/run-{sample_id}-2.sh", "w") as f:
            f.write(s)

    # Write the qsub script.
    s = f"p={project_path}\n"

    for sample_id in fastq_files:
        q1 = "qsub -e $p/log -o $p/log"

        if qsub_options1 != "NONE":
            q1 += f" {qsub_options1}"

        q2 = "qsub -e $p/log -o $p/log"

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
