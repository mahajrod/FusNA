rule umi_tools_dedup:
    input:
        sorted_bam=out_dir_path/ "alignment/STAR/{reference}/{sample}/{sample}.sorted.bam",
        sorted_bam_bai=out_dir_path/ "alignment/STAR/{reference}/{sample}/{sample}.sorted.bam.bai",
    priority: 1000
    output:
        rmdup_bam=out_dir_path/ "alignment/STAR/{reference}/{sample}/{sample}.rmdup.bam",
        #rmdup_stats=out_dir_path/ "alignment/STAR/{reference}/{sample}/{sample}.rmdup.stats"
    log:
        std=log_dir_path / "{sample}/umi_tools_dedup.{sample}.{reference}.log",
        err_log=log_dir_path / "{sample}/umi_tools_dedup.{sample}.{reference}.err.log",
        cluster_log=cluster_log_dir_path / "{sample}/umi_tools_dedup.{sample}.{reference}.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/umi_tools_dedup.{sample}.{reference}.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/umi_tools_dedup.{sample}.{reference}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["threads"]["umi_tools_dedup"],
        time=config["time"]["umi_tools_dedup"],
        mem=config["memory_mb"]["umi_tools_dedup"]
    threads: config["threads"]["umi_tools_dedup"]
    shell:
        " umi_tools dedup --paired -I {input.sorted_bam} -S {output.rmdup_bam} -L {log.std} > {log.err_log} 2>&1; "
