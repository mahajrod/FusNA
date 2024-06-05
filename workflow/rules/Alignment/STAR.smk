

rule STAR:
    input:
        #fastq_dir=rules.create_fastq_links.output,
        forward_fastqs=rules.cutadapt.output.forward_fastq,
        reverse_fastqs=rules.cutadapt.output.reverse_fastq,
        genome_dir=lambda wildcards: reference_dict[wildcards.reference]["STAR_index"]
    output:
        unsorted_bam=out_dir_path/ "alignment/STAR/{reference}{sample}/{sample}.unsorted.bam",
    log:
        std=log_dir_path / "{sample}/STAR.{sample}.{reference}.log",
        cluster_log=cluster_log_dir_path / "{sample}/STAR.{sample}.{reference}.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/STAR.{sample}.{reference}.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/STAR.{sample}.{reference},benchmark.txt"
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
        " STAR --runThreadN {threads} --genomeDir {input.genome_dir} --genomeLoad NoSharedMemory "
        " --readFilesIn {input.forward_fastqs} {input.reverse_fastqs} --readFilesCommand zcat "
        " --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 "
        " --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 "
        " --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 10 --chimOutType WithinBAM HardClip "
        " --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 "
        " --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50  > {output.unsorted_bam} 2>{log.std} "