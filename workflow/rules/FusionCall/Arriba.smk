localrules: add_sample_id

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
        " -p {input.protein_domains} > {log.std} 2>&1; "

rule add_sample_id:
    input:
        fusions=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.fusions.tsv",
    output:
        fusions=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.fusions.labeled.tsv",

    log:
        std=log_dir_path / "{sample}/add_sample_id.{sample}.{aligner}.{reference}.log",
        cluster_log=cluster_log_dir_path / "{sample}/add_sample_id.{sample}.{aligner}.{reference}.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/add_sample_id.{sample}.{aligner}.{reference}.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/add_sample_id.{sample}.{aligner}.{reference}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["threads"]["add_sample_id"],
        time=config["time"]["add_sample_id"],
        mem=config["memory_mb"]["add_sample_id"],
        io=1
    threads:
        config["threads"]["add_sample_id"]
    shell:
        " workflow/scripts/add_sample_and_fusion_ids.py -i {input.fusions} -s {wildcards.sample}"
        "  -o {output.fusions} > {log.std} 2>&1; "

rule combine_arriba_fusion_files:
    input:
        fusions=expand(out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.fusions.labeled.tsv",
                       sample=sample_list),
    output:
        fusions=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/all_samples.fusions.labeled.tsv",
    log:
        std=log_dir_path / "combine_arriba_fusion_files.{aligner}.{reference}.log",
        cluster_log=cluster_log_dir_path / "combine_arriba_fusion_files.{aligner}.{reference}.cluster.log",
        cluster_err=cluster_log_dir_path / "combine_arriba_fusion_files.{aligner}.{reference}.cluster.err"
    benchmark:
        benchmark_dir_path / "combine_arriba_fusion_files.{aligner}.{reference}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["threads"]["combine_arriba_fusion_files"],
        time=config["time"]["combine_arriba_fusion_files"],
        mem=config["memory_mb"]["combine_arriba_fusion_files"],
        io=1
    threads:
        config["threads"]["combine_arriba_fusion_files"]
    shell:
        " FUSION_FILE_ARRAY=({input.fusions}); "
        " head -n 1 ${{FUSION_FILE_ARRAY[0]}} > {output.fusions} 2>{log.std}; "
        " for FILE in ${{FUSION_FILE_ARRAY[@]}}; "
        "   do "
        "   tail -n +2 ${{FILE}} >> {output.fusions} 2>>{log.std}; "
        "   done "