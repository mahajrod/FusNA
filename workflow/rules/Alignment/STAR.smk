

rule STAR:
    input:
        #fastq_dir=rules.create_fastq_links.output,
        forward_fastq=rules.cutadapt.output.forward_fastq,
        reverse_fastq=rules.cutadapt.output.reverse_fastq,
        genome_dir=lambda wildcards: reference_dict[wildcards.reference]["STAR_index"]
    output:
        unsorted_bam=out_dir_path/ "alignment/STAR/{reference}/{sample}/{sample}.unsorted.bam",
    log:
        std=log_dir_path / "{sample}/STAR.{sample}.{reference}.log",
        cluster_log=cluster_log_dir_path / "{sample}/STAR.{sample}.{reference}.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/STAR.{sample}.{reference}.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/STAR.{sample}.{reference}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["threads"]["STAR"],
        time=config["time"]["STAR"],
        mem=config["memory_mb"]["STAR"],
        io=1
    threads:
        config["threads"]["STAR"]
    shell:
        " OUT_DIR=`dirname {output.unsorted_bam}`; "
        " LOG=`realpath {log.std}`; "
        " GENOME_DIR=`realpath {input.genome_dir}`; "
        " FORWARD_FASTQ=`realpath {input.forward_fastq}`; "
        " REVERSE_FASTQ=`realpath {input.reverse_fastq}`; "
        " BAM=`basename {output.unsorted_bam}`; "
        " cd ${{OUT_DIR}}; "
        " STAR --runThreadN {threads} --genomeDir ${{GENOME_DIR}} --genomeLoad NoSharedMemory "
        " --readFilesIn ${{FORWARD_FASTQ}} ${{REVERSE_FASTQ}} --readFilesCommand zcat "
        " --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 "
        " --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 "
        " --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 10 --chimOutType WithinBAM HardClip "
        " --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 "
        " --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50  > ${{BAM}} 2>${{LOG}} "