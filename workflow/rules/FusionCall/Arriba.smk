

rule arriba:
    input:
        #fastq_dir=rules.create_fastq_links.output,
        bam=out_dir_path/ "alignment/{alignment_tool}/{reference}/{sample}/{sample}.unsorted.bam",
        reference_fasta=,
        reference_annotation=
    output:
        fusions=out_dir_path/ "fusion_call/{alignment_tool}..arriba/{reference}/{sample}/{sample}.fusions.tsv",
        fusions_discarded=out_dir_path/ "alignment/{alignment_tool}..arriba/{reference}/{sample}/{sample}.fusions.discarded.tsv"
    log:
        std=log_dir_path / "{sample}/arriba.{sample}.{alignment_tool}.{reference}.log",
        cluster_log=cluster_log_dir_path / "{sample}/arriba.{sample}.{alignment_tool}.{reference}.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/arriba.{sample}.{alignment_tool}.{reference}.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/arriba.{sample}.{alignment_tool}.{reference}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
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
        " -b $BLACKLIST_TSV -k $KNOWN_FUSIONS_TSV -t $TAGS_TSV -p $PROTEIN_DOMAINS_GFF3  2>{log.std} "