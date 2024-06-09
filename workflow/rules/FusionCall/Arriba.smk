localrules: add_sample_id

def arriba_duplicate_handling_mode(wildcards):
    if wildcards.stage == "sorted":
        # at sorted stage no deduplication were done, so internal arriba's deduplication (enabled in arriba by default) is used
        return ""
    if config["panel_parameters"][config["panel"]]["UMI"] and config["panel_parameters"][config["panel"]]["use_UMI"]:
        return " -u "
    return ""

rule arriba:
    input:
        bam=out_dir_path/ "alignment/{aligner}/{reference}/{sample}/{sample}.{stage}.bam",
        reference_fasta=lambda wildcards: reference_dict[wildcards.reference]["fasta"],
        reference_annotation=lambda wildcards: reference_dict[wildcards.reference]["annotation"],
        blacklist=lambda wildcards: reference_dict[wildcards.reference]["blacklist"],
        known_fusions=lambda wildcards: reference_dict[wildcards.reference]["known_fusions"],
        protein_domains=lambda wildcards: reference_dict[wildcards.reference]["protein_domains"],
    params:
        use_external_duplicate_flag=arriba_duplicate_handling_mode
    output:
        fusions=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.{stage}.fusions.tsv",
        fusions_discarded=out_dir_path/ "alignment/{aligner}..arriba/{reference}/{sample}/{sample}.{stage}.fusions.discarded.tsv"
    log:
        std=log_dir_path / "{sample}/arriba.{sample}.{aligner}.{reference}.{stage}.log",
        cluster_log=cluster_log_dir_path / "{sample}/arriba.{sample}.{aligner}.{reference}.{stage}.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/arriba.{sample}.{aligner}.{reference}.{stage}.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/arriba.{sample}.{aligner}.{reference}.{stage}.benchmark.txt"
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
        " arriba {params.use_external_duplicate_flag} "
        " -x {input.bam} -o {output.fusions} -O {output.fusions_discarded} "
        " -a {input.reference_fasta} -g {input.reference_annotation} "
        " -b {input.blacklist} -k {input.known_fusions} -t {input.known_fusions} "
        " -p {input.protein_domains} > {log.std} 2>&1; "

rule filter_arriba_output:
    input:
        fusions=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.{stage}.fusions.tsv"
    output:
        filtered_fusions=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.{stage}.fusions.filtered.tsv",
        filtered_out_fusions=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.{stage}.fusions.filtered_out.tsv",
        all_with_filters_fusions=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.{stage}.fusions.all_with_filters.tsv",
    params:
        min_coverage=config["filter_arriba_output"]["min_coverage"],
        min_non_zero_sup_cov=config["filter_arriba_output"]["min_non_zero_sup_cov"],
        coverage_ratio_threshold=config["filter_arriba_output"]["coverage_ratio_threshold"],
        min_supporting_ratios_above_threshold=config["filter_arriba_output"]["min_supporting_ratios_above_threshold"]
    log:
        std=log_dir_path / "{sample}/filter_arriba_output.{sample}.{aligner}.{reference}.{stage}.log",
        cluster_log=cluster_log_dir_path / "{sample}/filter_arriba_output.{sample}.{aligner}.{reference}.{stage}.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/filter_arriba_output.{sample}.{aligner}.{reference}.{stage}.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/filter_arriba_output.{sample}.{aligner}.{reference}.{stage}.benchmark.txt"
    conda:
        "../../../%s" % config["arriba_config"]
    resources:
        cpus=config["threads"]["filter_arriba_output"],
        time=config["time"]["filter_arriba_output"],
        mem=config["memory_mb"]["filter_arriba_output"],
        io=1
    threads:
        config["threads"]["arriba"]
    shell:
        " workflow/scripts/filter_arriba_fusion.py -x {params.min_coverage} -n {params.min_non_zero_sup_cov} "
        " -a {params.coverage_ratio_threshold} -m {params.min_supporting_ratios_above_threshold} "
        " -i {input.fusions} -r {output.filtered_out_fusions} "
        " -w {output.all_with_filters_fusions} -f {output.filtered_fusions} > {log.std} 2>&1; "

rule classify_genes:
    input:
        fusions=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.{stage}.fusions.{filtering}.tsv"
    output:
        control_genes=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.{stage}.fusions.{filtering}.control.tsv",
        target_genes=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.{stage}.fusions.{filtering}.target.tsv",
        offtarget_genes=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.{stage}.fusions.{filtering}.offtarget.tsv",
        all_genes_classified=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.{stage}.fusions.{filtering}.all_classified.tsv",
    params:
        control_gene_file=panel_dir_path / config["panel"] / "control_genes.ids",
        target_gene_file=panel_dir_path / config["panel"] / "target_genes.ids",
    log:
        std=log_dir_path / "{sample}/classify_genes.{sample}.{aligner}.{reference}.{filtering}.{stage}.log",
        cluster_log=cluster_log_dir_path / "{sample}/classify_genes.{sample}.{aligner}.{reference}.{filtering}.{stage}.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/classify_genes.{sample}.{aligner}.{reference}.{filtering}.{stage}.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/classify_genes.{sample}.{aligner}.{reference}.{filtering}.{stage}.benchmark.txt"
    conda:
        "../../../%s" % config["arriba_config"]
    resources:
        cpus=config["threads"]["classify_genes"],
        time=config["time"]["classify_genes"],
        mem=config["memory_mb"]["classify_genes"],
        io=1
    threads:
        config["threads"]["arriba"]
    shell:
        " workflow/scripts/filter_by_panel_genes.py -i {input.fusions} "
        " -c {params.control_gene_file} -t {params.target_gene_file} "
        " --offtarget_output {output.offtarget_genes} --control_output {output.control_genes} "
        "--target_output {output.target_genes} --all_output {output.all_genes_classified} > {log.std} 2>&1; "

rule visualization:
    input:
        fusions=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.{stage}.fusions.{filter}.{gene_type}.tsv",
        sorted_bam=out_dir_path/ "alignment/{aligner}/{reference}/{sample}/{sample}.{stage}.bam",
        sorted_bam_index=out_dir_path/ "alignment/{aligner}/{reference}/{sample}/{sample}.{stage}.bam.bai",
        cytobands=lambda wildcards: reference_dict[wildcards.reference]["cytobands"],
        reference_annotation=lambda wildcards: reference_dict[wildcards.reference]["annotation"],
        protein_domains=lambda wildcards: reference_dict[wildcards.reference]["protein_domains"],
    output:
        fusions_pdf=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.{stage}.fusions.{filter}.{gene_type}.pdf",
    log:
        std=log_dir_path / "{sample}/visualization.{sample}.{aligner}.{reference}.{stage}.{filter}.{gene_type}.log",
        cluster_log=cluster_log_dir_path / "{sample}/avisualization.{sample}.{aligner}.{reference}.{stage}.{filter}.{gene_type}.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/visualization.{sample}.{aligner}.{reference}.{stage}.{filter}.{gene_type}.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/visualization.{sample}.{aligner}.{reference}.{stage}.{filter}.{gene_type}.benchmark.txt"
    conda:
        "../../../%s" % config["arriba_config"]
    resources:
        cpus=config["threads"]["visualization"],
        time=config["time"]["visualization"],
        mem=config["memory_mb"]["visualization"],
        io=1
    threads:
        config["threads"]["visualization"]
    shell:
        " draw_fusions.R --fusions={input.fusions} --alignments={input.sorted_bam} "
        " --annotation={input.reference_annotation} "
        " --cytobands={input.cytobands} "
        " --proteinDomains={input.protein_domains} "
        " --output={output.fusions_pdf} > {log.std} 2>&1; "

rule add_sample_id:
    input:
        fusions=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.{stage}.fusions.{filter}.{gene_type}.tsv",
    output:
        fusions=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.{stage}.fusions.{filter}.{gene_type}.labeled.tsv",
    log:
        std=log_dir_path / "{sample}/add_sample_id.{sample}.{aligner}.{reference}.{stage}.{filter}.{gene_type}.log",
        cluster_log=cluster_log_dir_path / "{sample}/add_sample_id.{sample}.{aligner}.{reference}.{stage}.{filter}.{gene_type}.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/add_sample_id.{sample}.{aligner}.{reference}.{stage}.{filter}.{gene_type}.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/add_sample_id.{sample}.{aligner}.{reference}.{stage}.{filter}.{gene_type}.benchmark.txt"
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
        fusions=expand(out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/{sample}/{sample}.{stage}.fusions.{filter}.{gene_type}.labeled.tsv",
                       sample=sample_list, allow_missing=True),
    output:
        fusions=out_dir_path/ "fusion_call/{aligner}..arriba/{reference}/all_samples.{stage}.fusions.{filter}.{gene_type}.labeled.tsv",
    log:
        std=log_dir_path / "combine_arriba_fusion_files.{aligner}.{reference}.{stage}.{filter}.{gene_type}.log",
        cluster_log=cluster_log_dir_path / "combine_arriba_fusion_files.{aligner}.{reference}.{stage}.{filter}.{gene_type}.cluster.log",
        cluster_err=cluster_log_dir_path / "combine_arriba_fusion_files.{aligner}.{reference}.{stage}.{filter}.{gene_type}.cluster.err"
    benchmark:
        benchmark_dir_path / "combine_arriba_fusion_files.{aligner}.{reference}.{stage}.{filter}.{gene_type}.benchmark.txt"
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