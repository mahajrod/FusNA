import glob
from collections import OrderedDict
from copy import deepcopy
import pandas as pd
from pathlib import Path

input_rna_dir = Path(config["input_rna_dir"])

sample_list = [element.name for element in input_rna_dir.glob("*")]
print(sample_list)

"""
def check_dna_index_presence(dna_reference_path):
    for index_extension in bwa_index_extension_list:
        if not (dna_reference_path.parent / (dna_reference_path.name + index_extension)).exists():
            return False
    return True

def check_rna_index_presence(rna_reference_path):
    for index_extension in bowtie2_index_extension_list:
        if not (rna_reference_path.parent / (rna_reference_path.name + index_extension)).exists():
            return False
    return True

#-- Initialization of path variables from config file --

#---- Initialization of path variables for input----
seq_run_dir_path = Path(config["seq_run_dir"])
input_fastq_dir_path = Path(config["input_fastq_dir"]) # used if fastq is starting point

sample_sheet_path = Path(config["sample_sheet"])
gene_so_retain_file_path = Path(config["gene_so_retain_file"])
ignore_so_terms_file_path = Path(config["ignore_so_terms_file"])
input_replicate_file_path = Path(config["replicate_file"])
dna_target_regions_dir_path = Path(config["target_dna_regions_dir"])
rna_target_regions_dir_path = Path(config["target_rna_regions_dir"])
rna_target_fusions_dir_path = Path(config["target_rna_fusions_dir"])
dna_target_regions_raw_reads_dir_path = Path(config["target_dna_raw_reads_regions_dir"])
ontarget_dna_regions_dir_path = Path(config["ontarget_dna_regions_dir"])
ontarget_rna_regions_dir_path = Path(config["ontarget_rna_regions_dir"])
mutations_of_interest_path = Path(config["mutations_of_interest"])
mutations_of_interest_short_ids_path = Path(config["mutations_of_interest_short_ids"])
fusions_of_interest_short_ids_path = Path(config["fusions_of_interest_short_ids"])
mRNA_of_interest_short_ids_path = Path(config["mRNA_of_interest_short_ids"])
mRNA_syn_file_path = Path(config["mRNA_syn_file"])

#---- Initialization of path variables for test mode ----
test_fastq_dir_path = Path(config["test_fastq_dir"]) / config["umi_type"] # used in test mode
test_gene_so_retain_file_path = Path(config["test_gene_so_retain_file"])
test_ignore_so_terms_file_path = Path(config["test_ignore_so_terms_file"])
test_dna_target_regions_dir_path = Path(config["test_target_dna_regions_dir"])
test_rna_target_regions_dir_path = Path(config["test_target_rna_regions_dir"])
test_rna_target_fusions_dir_path = Path(config["test_target_rna_fusions_dir"])
test_dna_target_regions_raw_reads_dir_path = Path(config["test_target_dna_raw_reads_regions_dir"])
test_ontarget_dna_regions_dir_path = Path(config["test_ontarget_dna_regions_dir"])
test_ontarget_rna_regions_dir_path = Path(config["test_ontarget_rna_regions_dir"])
test_mutations_of_interest_path = Path(config["test_mutations_of_interest"])
test_mutations_of_interest_short_ids_path = Path(config["test_mutations_of_interest_short_ids"])
test_fusions_of_interest_short_ids_path = Path(config["test_fusions_of_interest_short_ids"])
test_mRNA_of_interest_short_ids_path = Path(config["test_mRNA_of_interest_short_ids"])
test_mRNA_syn_file_path = Path(config["test_mRNA_syn_file"])
#---- Initialization of path variables for output directories ----
out_dir_path = Path(config["out_dir"])

log_dir_path = out_dir_path / config["log_dir"]
error_dir_path  = out_dir_path / config["error_dir"]
benchmark_dir_path =  out_dir_path / config["benchmark_dir"]
cluster_log_dir_path = out_dir_path / config["cluster_log_dir"]

preprocessing_dir_path = out_dir_path / config["preprocessing_dir"]
adjusted_sample_sheet_path = preprocessing_dir_path / 'SampleSheet.adjusted.csv' # sample_sheet_path.parent / (sample_sheet_path.stem + '.adjusted.csv')

rna_target_regions_gff_path = preprocessing_dir_path / "RNA.target.gff"
dna_target_regions_gff_path = preprocessing_dir_path  / "DNA.target.gff"
rna_ontarget_regions_gff_path = preprocessing_dir_path / "RNA.ontarget.gff"
dna_ontarget_regions_gff_path = preprocessing_dir_path  / "DNA.ontarget.gff"

basecall_dir_path = out_dir_path / config["basecall_dir"]
basecall_stat_folder = basecall_dir_path / "Stats"
basecall_report_dir = basecall_dir_path / "Reports"

per_lane_raw_fastq_dir_path = out_dir_path / config["fastq_dir"] / config["per_lane_raw_fastq_dir"]
merged_raw_fastq_dir_path = out_dir_path / config["fastq_dir"] / config["merged_raw_fastq_dir"]
filtered_fastq_dir_path = out_dir_path / config["fastq_dir"] / config["filtered_fastq_dir"]
consensus_fastq_dir_path = out_dir_path / config["fastq_dir"] / config["consensus_fastq_dir"]
per_lane_raw_fastqc_dir_path = out_dir_path / config["fastqc_dir"] / config["per_lane_raw_fastqc_dir"]
merged_raw_fastqc_dir_path = out_dir_path / config["fastqc_dir"] / config["merged_raw_fastqc_dir"]
filtered_fastqc_dir_path = out_dir_path / config["fastqc_dir"] / config["filtered_fastqc_dir"]
alignment_dir_path = out_dir_path / config["alignment_dir"]
variantcall_dir_path = out_dir_path / config["variantcall_dir"]
final_dir_path = out_dir_path / config["final_dir"]
final_per_sample_dir_path = final_dir_path / "per_sample"
#------- Initialization of path variables for output variables --------
basecall_stat_json_path = basecall_stat_folder / "Stats.json"
basecall_parsed_stats = final_dir_path / "basecalling.read.stats"

dna_target_regions_path = preprocessing_dir_path / "DNA.consensus.target.bed"
rna_target_regions_path = preprocessing_dir_path / "RNA.target.bed"
rna_target_fusions_path = preprocessing_dir_path / "RNA.fusions.bed"
dna_ontarget_regions_path = preprocessing_dir_path / "DNA.ontarget.bed"
rna_ontarget_regions_path = preprocessing_dir_path / "RNA.ontarget.bed"
dna_target_regions_raw_reads_path = preprocessing_dir_path / "DNA.raw.target.bed"

replicate_file_path = final_dir_path / "replicates.tab"

#------------ Initialization of path variable for panel of normals creation mode ------------
panel_of_normals_dir_path = out_dir_path / config["panel_of_normals_dir"]

#---- Initialization path variables for resources ----

#-------- Reference path variables --------
reference_dir_path = Path(config["reference_dir"])
rna_reference_dir_path = reference_dir_path / "RNA"
dna_reference_dir_path = reference_dir_path / "DNA"
rna_reference_path = Path(config["rna_reference_fasta"]) #reference_splited_dir_path / "reference.RNA.fasta"
dna_reference_path = Path(config["dna_reference_fasta"]) # reference_splited_dir_path / "reference.DNA.fasta"

rna_reference_fai_path = rna_reference_path.parent / (rna_reference_path.name + '.fai')
dna_reference_fai_path = dna_reference_path.parent / (dna_reference_path.name + '.fai')
rna_reference_dict_path = rna_reference_path.parent / (rna_reference_path.stem + '.dict')
dna_reference_dict_path = dna_reference_path.parent / (dna_reference_path.stem + '.dict')

bwa_index_dna_prefix = dna_reference_path
bwa_index_dna_ref_bwt_path = dna_reference_path.parent / (dna_reference_path.name + ".bwt")

bowtie2_index_rna_prefix = rna_reference_path
bowtie2_index_rna_ref_bt2_path = rna_reference_path.parent / (rna_reference_path.name + ".rev.1.bt2")

dna_reference_genome_size_file_path = reference_dir_path / "GenomeSize.xml"

# gff File describing lengths of transcripts (Required for whole transcriptome read count)

rna_total_gff_path = preprocessing_dir_path  /  "RNA.total.gff"
# gff File describing coordinates of mutations of interest (Required for read count on this positions)
mutations_of_interest_gff_path = preprocessing_dir_path / "mutations_of_interest.gff"
mutations_of_interest_counts_path = final_dir_path / "mutation_of_interest.read.counts"

rna_target_fusions_gff_path = preprocessing_dir_path / "RNA.fusions.gff"
rna_target_fusions_counts_path = final_dir_path / "rna.fusions.read.counts"

#-------- Variant filtration path variables ---------
gnomad_all_path = Path(config["germline_resource"]["gnomad_all"])
gnomad_all_index_path = gnomad_all_path.parent / (gnomad_all_path.name + ".tbi")
gnomad_target_path = preprocessing_dir_path / (gnomad_all_path.stem + ".target_regions.vcf")

cosmic_coding_all_path = Path(config["somatic_resource"]["cosmic_coding_all"])
cosmic_coding_all_index_path = cosmic_coding_all_path.parent / (cosmic_coding_all_path.name + ".tbi")
cosmic_target_path = preprocessing_dir_path / (cosmic_coding_all_path.stem + ".target_regions.vcf")

panel_of_normals_WGS = Path(config["mutect2"]["panel_of_normals"]["WGS"])
panel_of_normals_exome = Path(config["mutect2"]["panel_of_normals"]["exome"])
panel_of_normals_hybrid_capture = Path(config["mutect2"]["panel_of_normals"]["hybrid_capture"])

#-------- Setting desired panel of normals --------
if config["panel_of_normals"] in ["WGS", "exome", "hybrid_capture"]:
    panel_of_normals = Path(config["mutect2"]["panel_of_normals"][config["panel_of_normals"]])
else:
    panel_of_normals = Path(config["panel_of_normals"])
    if not panel_of_normals.is_file():
        raise ValueError("ERROR!!! Panel of normals is not available! Tried {}".format(str(panel_of_normals)))

#-------- Annotation path variables --------
cannonical_transcripts_path = Path(config["ensembl_cannonical_transcripts"])
#----

#-------- Verify index availability --------
dna_index_presence = check_dna_index_presence(dna_reference_path)
dna_genome_size_presence = dna_reference_genome_size_file_path.exists()
rna_index_presence = check_rna_index_presence(rna_reference_path)
rna_total_gff_presence = rna_total_gff_path.exists()

if not (dna_genome_size_presence and dna_index_presence):

    old_dna_reference_dir = dna_reference_dir_path
    old_dna_reference_path = dna_reference_path
    old_reference_dir_path = reference_dir_path

    reference_dir_path  = preprocessing_dir_path
    dna_reference_dir_path = preprocessing_dir_path / "DNA"
    dna_reference_path = dna_reference_dir_path / "reference.DNA.fasta" # reference_splited_dir_path / "reference.DNA.fasta"
    dna_reference_fai_path = dna_reference_path.parent / (dna_reference_path.name + '.fai')
    dna_reference_dict_path = dna_reference_path.parent / (dna_reference_path.stem + '.dict')
    bwa_index_dna_prefix = dna_reference_path
    bwa_index_dna_ref_bwt_path = dna_reference_path.parent / (dna_reference_path.name + ".bwt")
    dna_reference_genome_size_file_path = reference_dir_path / "GenomeSize.xml"

if not rna_index_presence:
    old_rna_reference_dir = rna_reference_dir_path
    old_rna_reference_path = rna_reference_path

    rna_reference_dir_path = preprocessing_dir_path / "RNA"
    rna_reference_path = rna_reference_dir_path / "reference.RNA.fasta" #reference_splited_dir_path / "reference.RNA.fasta"
    rna_reference_fai_path = rna_reference_path.parent / (rna_reference_path.name + '.fai')
    rna_reference_dict_path = rna_reference_path.parent / (rna_reference_path.stem + '.dict')
    bowtie2_index_rna_prefix = rna_reference_path
    bowtie2_index_rna_ref_bt2_path = rna_reference_path.parent / (rna_reference_path.name + ".rev.1.bt2")
"""
#---- Setting mode of pipeline ----


# ----
#-------------------------------------------

#cut_list = config["variant_calling"][config["umi_type"]]["three_prime_cut"]

#---- Functions ----
#----
localrules: all
results_list = []
if config["pipeline_mode"] in ["qc"]:

    results_list += []

#---- Final rule ----
rule all:
    input:
        results_list


"""
if pipeline_mode in ["index", "index_rna", "index_dna"]:
    include: "workflow/rules/Preprocessing/Reference.smk"
else:
    include: "workflow/rules/Preprocessing/Target.smk"
    include: "workflow/rules/Preprocessing/AnnotationSources.smk"
    include: "workflow/rules/BaseCall/BaseCall.smk"
    include: "workflow/rules/QCFiltering/Cutadapt.smk"
    include: "workflow/rules/QCFiltering/FastQC.smk"
    if config["umi_type"] == "UMI_UDI":
        include: "workflow/rules/Alignment/ConsensusCall/UMI_UDI.smk"
    elif config["umi_type"] == "UMI_duplex":
        include: "workflow/rules/Alignment/ConsensusCall/UMI_duplex.smk"
    include: "workflow/rules/Alignment/ConsensusCall/Common.smk"
    include: "workflow/rules/Alignment/RawReads.smk"
    include: "workflow/rules/VariantCall/Pisces.smk"
    include: "workflow/rules/VariantCall/Mutect2.smk"
    include: "workflow/rules/VariantCall/PanelOfNormals.smk"
    include: "workflow/rules/Alignment/CountReads.smk"
    include: "workflow/rules/Annotation/Filter.smk"
    include: "workflow/rules/Annotation/SnpEff.smk"
    include:  "workflow/rules/Stats/CollectStats.smk"
"""

