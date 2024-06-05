

rule arriba:
    input:
        #fastq_dir=rules.create_fastq_links.output,
        bam=out_dir_path/ "alignment/{aligner}/{reference}/{sample}/{sample}.unsorted.bam",
        reference_fasta=lambda wildcards: reference_dict[wildcards.reference]["fasta"],
        reference_annotation=lambda wildcards: reference_dict[wildcards.reference]["annotation"],
        blacklist=lambda wildcards: reference_dict[wildcards.reference]["blacklist"],
        known_fusions=lambda wildcards: reference_dict[wildcards.reference]["known_fusions"],
        protein_domains=lambda wildcards: reference_dict[wildcards.reference]["protein_domains"],
    output:
        fusions=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.fusions.tsv",
        fusions_discarded=out_dir_path/ "alignment/{aligner}..arriba/{reference}/{sample}/{sample}.fusions.discarded.tsv"
    log:
        std=log_dir_path / "{sample}/arriba.{sample}.{aligner}.{reference}.log",
        cluster_log=cluster_log_dir_path / "{sample}/arriba.{sample}.{aligner}.{reference}.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/arriba.{sample}.{aligner}.{reference}.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/arriba.{sample}.{aligner}.{reference}.benchmark.txt"
    conda:
        "../../../%s" % config["arriba_config"]
    resources:
        cpus=config["threads"]["arriba"],
        time=config["time"]["arriba"],
        mem=config["memory_mb"]["arriba"],
        io=1
    threads:
        config["threads"]["arriba"]
    shell:
        " arriba -x {input.bam} -o {output.fusions} -O {output.fusions_discarded} "
        " -a {input.reference_fasta} -g {input.reference_annotation} "
        " -b {input.blacklist} -k {input.known_fusions} -t {input.known_fusions} "
        " -p {input.protein_domains} >{log.std} 2>&1; "