import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")
configfile: "config.json"

GLOBAL_REF_PATH = config["globalResources"]
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

if not 'aligner' in config:
  config['aligner'] = "bowtie2" # BWA, Bowtie2
if not 'max_len_frags' in config:
  config['max_len_frags'] = "120"
if not 'bowtie2_sens' in config:
  config['bowtie2_sens'] = "very"
if not 'dovetailing' in config:
  config['dovetailing'] = True

# Reference processing
#
# if config["lib_ROI"] != "wgs":
#     # setting reference from lib_ROI
#     f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","lib_ROI.json"))
#     lib_ROI_dict = json.load(f)
#     f.close()
#     config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]


# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference2.json"),)
reference_dict = json.load(f)
f.close()
config["species_name"] = [organism_name for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
config["organism"] = config["species_name"].split(" (")[0].lower().replace(" ","_")
if len(config["species_name"].split(" (")) > 1:
    config["species"] = config["species_name"].split(" (")[1].replace(")","")

if not "min_qual" in config:
    config['min_qual'] = "20"

##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])
fastq_dir = "cleaned_fastq" if (config["preprocess"]!="none") else "raw_fastq"

# Samples
#
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

if not config["is_paired"]:
    read_pair_tags = [""]
    pair_tag = ["SE"]
else:
    read_pair_tags = ["_R1","_R2"]
    pair_tag = ["R1","R2"]

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),
    read_pair_tag = "(_R.)?"


##### Target rules #####
rule all:
    input:  expand("mapped/{sample}.bam",sample = sample_tab.sample_name),
            expand("mapped/{sample}.bam.bai", sample = sample_tab.sample_name),
            "qc_reports/all_samples/alignment_chip_multiqc/multiqc.html"

##### Modules #####

include: "rules/alignment_ChIP.smk"
# include: "rules/prepare_reference.smk"
