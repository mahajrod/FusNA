

rule fastqc:
    input:
        fastq=out_dir_path/ ("data/{stage}/{fastq_prefix}%s" % config["data_type_description"]["fastq"]["output"]["extension"]),
    output:
        fastqc=out_dir_path/ "qc/fastqc/{stage}/{fastq_prefix}_fastqc.zip"
    params:
        kmer=config["fastqc"]["kmer_length"],
    log:
        std=out_dir_path/ "log/fastqc.{stage}.{fastq_prefix}.log",
        cluster_log=out_dir_path/ "cluster_log/fastqc.{stage}.{fastq_prefix}.log",
        cluster_err=out_dir_path/ "cluster_err/fastqc.{stage}.{fastq_prefix}.err"
    benchmark:
        out_dir_path/ "benchmark/fastqc.{stage}.{fastq_prefix}.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["threads"]["fastqc"],
        time=config["time"]["fastqc"],
        mem=config["memory_mb"]["fastqc"],
        io=1
    threads:
        config["threads"]["fastqc"]
    shell:
        #" mkdir -p {output.dir}; "
        " OUTDIR=`dirname {output.fastqc}`; "
        " fastqc --nogroup -k {params.kmer} -t {threads} -o ${{OUTDIR}} {input}.fastq 1>{log.std} 2>&1; "
        #" workflow/scripts/convert_fastqc_output.py -f {output.forward_fastqc} -r {output.reverse_fastqc} "
        #" -s {wildcards.library_id} -o {output.stats} 1>{log.stats} 2>&1 "
