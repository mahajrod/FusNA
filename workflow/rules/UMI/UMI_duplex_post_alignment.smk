
rule merge_alignment_with_umi_bam:
    input:
        mapped_unsorted_bam=out_dir_path/ "alignment/{aligner}/{reference}/{sample}/{sample}.bam",
        unmapped_bam=rules.extract_umis_from_bam.output.bam,
        reference_fasta=lambda wildcards: reference_dict[wildcards.reference]["fasta"],
    priority: 1000
    output:
        sorted_umi_bam=out_dir_path/ "alignment/{aligner}/{reference}/{sample}/{sample}.with_umi.sorted.bam"
    params:
        fixmate_threads=config["threads"]["fixmate"],
        sort_threads=config["threads"]["sort"],
        bwa_threads=config["threads"]["bwa"],
        picard_threads=config["threads"]["picard"],
        picard_mem=config["memory_mb"]["picard"],
        per_thread_sort_mem="%sG" % config["per_thread_sort_mem"],
    log:
        samtofastq=log_dir_path / "{sample}/bwa_dna_map.samtofastq.{sample}.{aligner}.{reference}.log",
        mergebamalignments=log_dir_path / "{sample}/bwa_dna_map.mergebamalignments.{sample}.{aligner}.{reference}.log",
        fixmate=log_dir_path / "{sample}/bwa_dna_map.fixmate.{sample}.{aligner}.{reference}.log",
        bwa=log_dir_path / "{sample}/bwa_dna_map.bwa.{sample}.{aligner}.{reference}.log",
        sort=log_dir_path / "{sample}/bwa_dna_map.sort.{sample}.{aligner}.{reference}.log",
        cluster_log=cluster_log_dir_path / "{sample}/bwa_dna_map.{sample}.{aligner}.{reference}.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/bwa_dna_map.{sample}.{aligner}.{reference}.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/bwa_dna_map.benchmark.{sample}.{aligner}.{reference}.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["threads"]["sort"] + config["threads"]["fixmate"] + 2 * config["threads"]["picard"],
        time=config["time"]["merge_alignment_with_umi_bam"],
        mem=config["per_thread_sort_mem"] * config["threads"]["sort"] * 1024 +  config["memory_mb"]["fixmate"] + 2 * config["memory_mb"]["picard"]
    threads: config["threads"]["bwa"]
    shell:
        " TMP_PREFIX={output.sorted_umi_bam}.tmp; "
        " picard -Xmx{params.picard_mem}m MergeBamAlignment UNMAPPED={input.unmapped_bam} "
        " ALIGNED={input.mapped_unsorted_bam} R={input.reference_fasta} SO=unsorted "
        " ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR "
        " VALIDATION_STRINGENCY=STRICT CREATE_INDEX=false O=/dev/stdout 2>{log.mergebamalignments} | "
        " samtools fixmate -@ {params.fixmate_threads} -m - -  2>{log.fixmate} | "
        " samtools sort -T ${{TMP_PREFIX}} -@ {params.sort_threads} "
        " -m {params.per_thread_sort_mem} -o {output.sorted_umi_bam} 2>{log.sort} "



rule group_reads_by_umi:
    input:
        bam=alignment_dir_path / "{library_id}/{library_id}.{cut_type}.{mol_type}.raw.sorted.bam"
    output:
        bam=temp(alignment_dir_path / "{library_id}/{library_id}.{cut_type}.{mol_type}.umi.grouped.bam"),
        read_per_umi_hist=alignment_dir_path / "{library_id}/{library_id}.{cut_type}.{mol_type}.umi.grouped.hist",
        stats=alignment_dir_path / "{library_id}/{library_id}.{cut_type}.{mol_type}.group_umi.stats",
        hist_stats=alignment_dir_path / "{library_id}/{library_id}.{cut_type}.{mol_type}.group_umi.hist.stats"
    log:
        std=log_dir_path / "{library_id}/{cut_type}.{mol_type}.group_reads_by_umi.log",
        stats=log_dir_path / "{library_id}/{cut_type}.{mol_type}.group_reads_by_umi.stats.log",
        cluster_log=cluster_log_dir_path / "{library_id}.{cut_type}.{mol_type}.group_reads_by_umi.cluster.log",
        cluster_err=cluster_log_dir_path / "{library_id}.{cut_type}.{mol_type}.group_reads_by_umi.cluster.err"
    params:
        max_umi_mismatch=config["fgbio"]["UMI_duplex"]["GroupReadsByUmi"]["max_umi_mismatch"],
        min_mapping_quality=config["fgbio"]["UMI_duplex"]["GroupReadsByUmi"]["min_mapping_quality"],
        strategy=config["fgbio"]["UMI_duplex"]["GroupReadsByUmi"]["strategy"]
    benchmark:
        benchmark_dir_path / "{library_id}/{cut_type}.{mol_type}.group_reads_by_umi.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["threads"]["group_reads_by_umi"],
        time=config["time"]["group_reads_by_umi"],
        mem=config["memory_mb"]["group_reads_by_umi"],
    threads: config["threads"]["group_reads_by_umi"]
    shell:
        " fgbio -Xmx{resources.mem}m GroupReadsByUmi -i {input} -o {output.bam} -s {params.strategy} "
        " --family-size-histogram {output.read_per_umi_hist} "
        " -e {params.max_umi_mismatch} -m {params.min_mapping_quality} > {log.std} 2>&1;"
        " ONE=`head -n 2 {output.read_per_umi_hist} | tail -n 1 | cut -f 2`; "
        " ONE_FR=`head -n 2 {output.read_per_umi_hist} | tail -n 1 | cut -f 3`; "
        " TWO=`head -n 3 {output.read_per_umi_hist} | tail -n 1 | cut -f 2`; "
        " TWO_FR=`head -n 3 {output.read_per_umi_hist} | tail -n 1 | cut -f 3`; "
        " THREE_PLUS=`awk 'BEGIN {{SUM=0}}; NR>3 {{SUM+=$2}}; END {{print SUM}}' {output.read_per_umi_hist}`; "
        " THREE_PLUS_FR=`awk 'BEGIN {{SUM=0}}; NR>3 {{SUM+=$3}}; END {{print SUM}}' {output.read_per_umi_hist}`; "
        " ACCEPTED=`grep 'Accepted' {log.std} | sed 's/.*Accepted \([0-9,]\+\).*/\\1/;s/,//g'`; "
        " FILTERED_OUT_MAPPINGS=`grep 'due to mapping issues' {log.std} | sed 's/.*out \([0-9,]\+\).*/\\1/;s/,//g'`; "
        " FILTERED_OUT_AMBIGUOUS_UMI=`grep 'that contained one or more Ns in their UMIs' {log.std} | sed 's/.*out \([0-9,]\+\).*/\\1/;s/,//g'`; "
        " echo 'library_id	accepted_reads(grouping)	filtered_out_mapping_issues(grouping)	filtered_out_ambiguous_umi(grouping)' > {output.stats} 2>{log.stats}; "
        " echo \"{wildcards.library_id}	${{ACCEPTED}}	${{FILTERED_OUT_MAPPINGS}}	${{FILTERED_OUT_AMBIGUOUS_UMI}}\" >> {output.stats} 2>>{log.stats};"
        " echo 'library_id	1	1(fraction)	2	2(fraction)	3+	3+(fraction)' > {output.hist_stats} 2>> {log.stats}; " 
        " echo \"{wildcards.library_id}	${{ONE}}	${{ONE_FR}}	${{TWO}}	${{TWO_FR}}	${{THREE_PLUS}}	${{THREE_PLUS_FR}}\" >> {output.hist_stats} 2>>{log.stats};"
"""
rule create_duplex_consensus_reads:
    input:
        rules.group_reads_by_umi.output.bam
    output:
        bam=alignment_dir_path / "{library_id}/{library_id}.{cut_type}.{mol_type}.consensus.unmaped.bam",#TODO: add temp after testing
        stats=alignment_dir_path / "{library_id}/{library_id}.{cut_type}.{mol_type}.consensus.stats"
    log:
        std=log_dir_path / "{library_id}/{cut_type}.{mol_type}.create_duplex_consensus_reads.log",
        stats=log_dir_path / "{library_id}/{cut_type}.{mol_type}.create_duplex_consensus_reads.log",
        cluster_log=cluster_log_dir_path / "{library_id}.{cut_type}.{mol_type}.create_duplex_consensus_reads.cluster.log",
        cluster_err=cluster_log_dir_path / "{library_id}.{cut_type}.{mol_type}.create_duplex_consensus_reads.cluster.err"
    params:
        min_reads=config["fgbio"]["UMI_duplex"]["CallDuplexConsensusReads"]["min_reads"],
        error_rate_post_umi=config["fgbio"]["UMI_duplex"]["CallDuplexConsensusReads"]["error_rate_post_umi"]
    benchmark:
        benchmark_dir_path / "{library_id}/{cut_type}.{mol_type}.create_duplex_consensus_reads.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["threads"]["create_consensus_reads"],
        time=config["time"]["create_consensus_reads"],
        mem=config["memory_mb"]["create_consensus_reads"],
    threads: config["threads"]["create_consensus_reads"]
    shell:
        " fgbio -Xmx{resources.mem}m CallDuplexConsensusReads -i {input} -o {output.bam} "
        " --error-rate-post-umi {params.error_rate_post_umi} --min-reads {params.min_reads} > {log.std} 2>&1;"
        " COUNTS=`grep 'Consensus reads emitted' {log.std} | sed 's/.*emitted: \([0-9,]\+\).*/\\1/;s/,//g'`; "
        " echo 'library_id	consensus_reads' > {output.stats} 2>{log.stats}; "
        " echo \"{wildcards.library_id}	${{COUNTS}}\" >> {output.stats} 2>>{log.stats} "

"""