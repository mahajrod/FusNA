

rule merge_files:
    input:
        #fastq_dir=rules.create_fastq_links.output,
        forward_fastqs=lambda wildcards: sample_dict[wildcards.sample]["input_files"][::2],
        reverse_fastqs=lambda wildcards: sample_dict[wildcards.sample]["input_files"][1::2]
    output:
        forward_fastq=out_dir_path/ "data/merged_raw/{sample}%s%s" % (config["data_type_description"]["fastq"]["output"]["suffix_list"]["forward"],
                                                                      config["data_type_description"]["fastq"]["output"]["extension"]),
        reverse_fastq=out_dir_path/ "data/merged_raw/{sample}%s%s" % (config["data_type_description"]["fastq"]["output"]["suffix_list"]["reverse"],
                                                                      config["data_type_description"]["fastq"]["output"]["extension"]),
        #stats=merged_raw_fastqc_dir_path / "{library_id}/{library_id}.raw.fastqc.stats"
    log:
        std=out_dir_path/ "log/merge_files.{sample}.log",
        cluster_log=out_dir_path/ "cluster_log/merge_files.{sample}.cluster.log",
        cluster_err=out_dir_path/ "cluster_err/merge_files.{sample}.cluster.err"
    benchmark:
        out_dir_path/ "benchmark/merge_files.{sample}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=config["threads"]["merge_files"],
        time=config["time"]["merge_files"],
        mem=config["memory_mb"]["merge_files"],
    threads:
        config["threads"]["merge_files"]
    shell:
        " cat {input.forward_fastqs} > {output.forward_fastq} 2>{log.std}; "
        " cat {input.reverse_fastqs} > {output.reverse_fastq} 2>{log.std}; "
