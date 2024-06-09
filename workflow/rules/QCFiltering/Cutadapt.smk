
rule cutadapt:
    input:
        forward_fastq=out_dir_path/ ("data/trimmed/{sample}/{sample}%s%s" % (config["data_type_description"]["fastq"]["output"]["suffix_list"]["forward"],
                                                                          config["data_type_description"]["fastq"]["output"]["extension"])),
        reverse_fastq=out_dir_path/ ("data/trimmed/{sample}/{sample}%s%s" % (config["data_type_description"]["fastq"]["output"]["suffix_list"]["reverse"],
                                                                          config["data_type_description"]["fastq"]["output"]["extension"])),
    output:
        forward_fastq=out_dir_path/ ("data/filtered/{sample}/{sample}%s%s" % (config["data_type_description"]["fastq"]["output"]["suffix_list"]["forward"],
                                                                              config["data_type_description"]["fastq"]["output"]["extension"])),
        reverse_fastq=out_dir_path/ ("data/filtered/{sample}/{sample}%s%s" % (config["data_type_description"]["fastq"]["output"]["suffix_list"]["reverse"],
                                                                              config["data_type_description"]["fastq"]["output"]["extension"])),
        stats=out_dir_path/ "data/filtered/{sample}/{sample}.stats"
    params:
        adapters=config["cutadapt"]["adapter_seq"],
        minlength=config["cutadapt"]["minlength"]
    log:
        std=log_dir_path / "{sample}/{sample}.cutadapt.log",
        stats=log_dir_path / "{sample}/{sample}.cutadapt.stats.log",
        cluster_log=cluster_log_dir_path / "{sample}.cutadapt.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.cutadapt.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/cutadapt.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["threads"]["cutadapt"],
        time=config["time"]["cutadapt"],
        mem=config["memory_mb"]["cutadapt"],
    threads:
        config["threads"]["cutadapt"]
    shell:
         "cutadapt -j {threads} -m {params.minlength} -a {params.adapters} -A {params.adapters} "
         "-o {output.forward_fastq} -p {output.reverse_fastq} {input.forward_fastq} {input.reverse_fastq} > {log.std} 2>&1; "
         " echo 'library_id	input read pairs	read1 with adapter	read2 with adapter	too short pairs	pairs passed filters' > {output.stats} 2>{log.stats}; "
         " COUNTS=`grep -A 6  '=== Summary ===' {log.std} | sed -ne  's/.*: \+//g;s/ .*//g;s/,//g;3,$p' | tr '\\n' '\t' | sed 's/\t$//'`; "
         " echo \"{wildcards.sample}	${{COUNTS}}\" >> {output.stats} 2>>{log.stats}"
