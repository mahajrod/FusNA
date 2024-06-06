rule index_bam:
    input:
        bam="{bam_prefix}.bam"
    output:
        bai="{bam_prefix}.bam.bai"
    log:
        std="{bam_prefix}.index_bam.log",
        cluster_log="{bam_prefix}.index_bam.cluster.log",
        cluster_err="{bam_prefix}.index_bam.cluster.err"
    benchmark:
        "{bam_prefix}.index_bam.benchmark.txt"
    conda:
        "../../../%s" % config["arriba_config"]
    resources:
        cpus=config["threads"]["samtools_index"] ,
        time=config["time"]["samtools_index"],
        mem=config["memory_mb"]["samtools_index"]
    threads: config["threads"]["samtools_index"]
    shell:
        " samtools index -@ {threads} {input} > {log.std} 2>&1; "

rule sort_bam:
    input:
        bam="{bam_prefix}.bam"
    output:
        bam="{bam_prefix}.sorted.bam"
    params:
        per_thread_sort_mem=config["memory_mb"]["sort_bam_per_thread"]
    log:
        std="{bam_prefix}.sort_bam.log",
        cluster_log="{bam_prefix}.sort_bam.cluster.log",
        cluster_err="{bam_prefix}.sort_bam.cluster.err"
    benchmark:
        "{bam_prefix}.sort_bam.benchmark.txt"
    conda:
        "../../../%s" % config["arriba_config"]
    resources:
        cpus=config["threads"]["sort_bam"] ,
        time=config["time"]["sort_bam"],
        mem=config["threads"]["sort_bam"] * config["memory_mb"]["sort_bam_per_thread"] + 500
    threads: config["threads"]["sort_bam"]
    shell:
        " TMP_PREFIX=`basename {output.bam}`'.sort_tmp'; "
        " samtools sort -T ${{TMP_PREFIX}} -@ {threads} -m {params.per_thread_sort_mem} "
        " -o {output.bam} {input.bam} > {log.std} 2>&1; "