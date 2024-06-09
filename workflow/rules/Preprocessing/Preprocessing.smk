
if config["panel_parameters"][config["panel"]]["UMI"] and (not config["panel_parameters"][config["panel"]]["use_UMI"]):
    # trim UMIs if they will not be used
    rule trim_UMIs:
        input:
            #fastq_dir=rules.create_fastq_links.output,
            forward_fastqs=lambda wildcards: sample_dict[wildcards.sample]["input_files"][::2],
            reverse_fastqs=lambda wildcards: sample_dict[wildcards.sample]["input_files"][1::2]
        output:
            forward_fastq=out_dir_path/ ("data/merged_raw/{sample}/{sample}%s%s" % (config["data_type_description"]["fastq"]["output"]["suffix_list"]["forward"],
                                                                          config["data_type_description"]["fastq"]["output"]["extension"])),
            reverse_fastq=out_dir_path/ ("data/merged_raw/{sample}/{sample}%s%s" % (config["data_type_description"]["fastq"]["output"]["suffix_list"]["reverse"],
                                                                          config["data_type_description"]["fastq"]["output"]["extension"])),
            #stats=merged_raw_fastqc_dir_path / "{library_id}/{library_id}.raw.fastqc.stats"
        params:
            forward_read_trim=config["panel_parameters"][config["panel"]]["forward_read_trim"],
            reverse_read_trim=config["panel_parameters"][config["panel"]]["reverse_read_trim"],
        log:
            forward_log=log_dir_path / "{sample}/trim_UMIs.forward.{sample}.log",
            reverse_log=log_dir_path / "{sample}/trim_UMIs.reverse.{sample}.log",
            cluster_log=cluster_log_dir_path / "{sample}/trim_UMIs.{sample}.cluster.log",
            cluster_err=cluster_log_dir_path / "{sample}/trim_UMIs.{sample}.cluster.err"
        benchmark:
            benchmark_dir_path / "{sample}/trim_UMIs.{sample}.benchmark.txt"
        conda:
            "../../../%s" % config["conda_config"]
        resources:
            cpus=config["threads"]["trim_UMIs"],
            time=config["time"]["trim_UMIs"],
            mem=config["memory_mb"]["trim_UMIs"],
            io=1
        threads:
            config["threads"]["trim_UMIs"]
        shell:
            " zcat {input.forward_fastqs} | fastx_trimmer -f {params.forward_read_trim} | "
            " pigz -p {threads} > {output.forward_fastq} 2>{log.forward_log}; "
            " zcat {input.reverse_fastqs} | fastx_trimmer -f {params.reverse_read_trim} | "
            " pigz -p {threads} > {output.reverse_fastq} 2>{log.reverse_log}; "

else:
    rule merge_files:
        input:
            #fastq_dir=rules.create_fastq_links.output,
            forward_fastqs=lambda wildcards: sample_dict[wildcards.sample]["input_files"][::2],
            reverse_fastqs=lambda wildcards: sample_dict[wildcards.sample]["input_files"][1::2]
        output:
            forward_fastq=out_dir_path/ ("data/trimmed/{sample}/{sample}%s%s" % (config["data_type_description"]["fastq"]["output"]["suffix_list"]["forward"],
                                                                          config["data_type_description"]["fastq"]["output"]["extension"])),
            reverse_fastq=out_dir_path/ ("data/trimmed/{sample}/{sample}%s%s" % (config["data_type_description"]["fastq"]["output"]["suffix_list"]["reverse"],
                                                                          config["data_type_description"]["fastq"]["output"]["extension"])),
            #stats=merged_raw_fastqc_dir_path / "{library_id}/{library_id}.raw.fastqc.stats"
        log:
            std=log_dir_path / "{sample}/merge_files.{sample}.log",
            cluster_log=cluster_log_dir_path / "{sample}/merge_files.{sample}.cluster.log",
            cluster_err=cluster_log_dir_path / "{sample}/merge_files.{sample}.cluster.err"
        benchmark:
            benchmark_dir_path / "{sample}/merge_files.{sample}.benchmark.txt"
        conda:
            "../../../%s" % config["conda_config"]
        resources:
            cpus=config["threads"]["merge_files"],
            time=config["time"]["merge_files"],
            mem=config["memory_mb"]["merge_files"],
            io=1
        threads:
            config["threads"]["merge_files"]
        shell:
            " cat {input.forward_fastqs} > {output.forward_fastq} 2>{log.std}; "
            " cat {input.reverse_fastqs} > {output.reverse_fastq} 2>{log.std}; "
