

rule fastqc:
    input:
        fastq=out_dir_path/ ("data/{qc_stage}/{sample}/{fastq_prefix}%s" % config["data_type_description"]["fastq"]["output"]["extension"]),
    output:
        fastqc=out_dir_path/ "qc/fastqc/{qc_stage}/{sample}/{fastq_prefix}_fastqc.zip",
        fastqc_dir=directory(out_dir_path/ "qc/fastqc/{qc_stage}/{sample}/{fastq_prefix}_fastqc")
    params:
        kmer=config["fastqc"]["kmer_length"],
    log:
        std=log_dir_path / "{sample}/fastqc.{qc_stage}.{sample}.{fastq_prefix}.log",
        unzip_log=log_dir_path / "{sample}/fastqc.{qc_stage}.{sample}.{fastq_prefix}.unzip.log",
        cluster_log=cluster_log_dir_path / "{sample}/fastqc.{qc_stage}.{sample}.{fastq_prefix}.log",
        cluster_err=cluster_log_dir_path / "{sample}/fastqc.{qc_stage}.{sample}.{fastq_prefix}.err"
    benchmark:
        benchmark_dir_path / "{sample}/fastqc.{qc_stage}.{sample}.{fastq_prefix}.benchmark.txt"
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
        " fastqc --nogroup -k {params.kmer} -t {threads} -o ${{OUTDIR}} {input.fastq} 1>{log.std} 2>&1; "
        " unzip {output.fastqc} -d ${{OUTDIR}} > {log.unzip_log} 2>&1; "
        #" workflow/scripts/convert_fastqc_output.py -f {output.forward_fastqc} -r {output.reverse_fastqc} "
        #" -s {wildcards.library_id} -o {output.stats} 1>{log.stats} 2>&1 "

