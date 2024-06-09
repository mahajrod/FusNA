rule umi_tools_extract_umi:
    input:
        forward_fastq=out_dir_path/ ("data/merged_raw/{sample}/{sample}%s%s" % (config["data_type_description"]["fastq"]["output"]["suffix_list"]["forward"],
                                                                          config["data_type_description"]["fastq"]["output"]["extension"])),
        reverse_fastq=out_dir_path/ ("data/merged_raw/{sample}/{sample}%s%s" % (config["data_type_description"]["fastq"]["output"]["suffix_list"]["reverse"],
                                                                          config["data_type_description"]["fastq"]["output"]["extension"])),
    priority: 1000
    params:
        forward_read_regex=lambda wildcards: "'(?P<umi_1>.{%s})(?P<discard_1>.{%s}).*'" % (config["panel_parameters"][config["panel"]]["forward_barcode_length"],
                                                                             config["panel_parameters"][config["panel"]]["forward_linker_length"],),
        reverse_read_regex=lambda wildcards: "'(?P<umi_2>.{%s})(?P<discard_2>.{%s}).*'" % (config["panel_parameters"][config["panel"]]["reverse_barcode_length"],
                                                                             config["panel_parameters"][config["panel"]]["reverse_linker_length"],),
    output:
        forward_fastq=out_dir_path/ ("data/trimmed/{sample}/{sample}%s%s" % (config["data_type_description"]["fastq"]["output"]["suffix_list"]["forward"],
                                                                          config["data_type_description"]["fastq"]["output"]["extension"])),
        reverse_fastq=out_dir_path/ ("data/trimmed/{sample}/{sample}%s%s" % (config["data_type_description"]["fastq"]["output"]["suffix_list"]["reverse"],
                                                                          config["data_type_description"]["fastq"]["output"]["extension"])),
    log:
        std=log_dir_path / "{sample}/umi_tools_extract_umi.{sample}.log",
        err_log=log_dir_path / "{sample}/umi_tools_extract_umi.{sample}.err.log",
        cluster_log=cluster_log_dir_path / "{sample}/umi_tools_extract_umi{sample}.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/umi_tools_extract_umi.{sample}.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/umi_tools_extract_umi.{sample}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["threads"]["umi_tools_extract_umi"],
        time=config["time"]["umi_tools_extract_umi"],
        mem=config["memory_mb"]["umi_tools_extract_umi"]
    threads: config["threads"]["umi_tools_extract_umi"]
    shell:
        " OUTPUT_FORWARD_FASTQ={output.forward_fastq}; "
        " OUTPUT_FORWARD_FASTQ=${{OUTPUT_FORWARD_FASTQ%.gz}}; "
        " OUTPUT_REVERSE_FASTQ={output.reverse_fastq}; "
        " OUTPUT_REVERSE_FASTQ=${{OUTPUT_REVERSE_FASTQ%.gz}}; "
        " umi_tools extract -I {input.forward_fastq} --read2-in {input.reverse_fastq} "
        " --extract-method=regex --bc-pattern={params.forward_read_regex} "
        " --bc-pattern2={params.forward_read_regex} "
        " --read2-out ${{OUTPUT_REVERSE_FASTQ}} -L {log.std}  > ${{OUTPUT_FORWARD_FASTQ}} 2>{log.err_log}; "
        " pigz -p {threads} ${{OUTPUT_FORWARD_FASTQ}}; "
        " pigz -p {threads} ${{OUTPUT_REVERSE_FASTQ}}; "

