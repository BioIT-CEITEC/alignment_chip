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
if not "min_qual" in config:
    config['min_qual'] = "30"

##### BioRoot utilities - basics #####
module BR:
    snakefile: github("BioIT-CEITEC/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*

fastq_dir = "processed_fastq" if (config["preprocess"]!="none") else "raw_fastq"

# Samples
sample_tab = BR.load_sample()

config = BR.load_organism()

spikein_ref = BR.spikein_reference()

if not config["is_paired"]:
    read_pair_tags = [""]
    pair_tag = ["SE"]
else:
    read_pair_tags = ["_R1","_R2"]
    pair_tag = ["R1","R2"]

wildcard_constraints:
    sample = ("|".join(sample_tab.sample_name)+'|'+'|'.join([name+".spike" for name in sample_tab.sample_name])),
    read_pair_tag = "(_R.)?"


##### Target rules #####
rule all:
    input:  expand("mapped/{sample}.bam",sample = sample_tab.sample_name),
            expand("mapped/{sample}.bam.bai", sample = sample_tab.sample_name),
            "qc_reports/all_samples/alignment_chip_multiqc/multiqc.html"

##### Modules #####

include: "rules/alignment_ChIP.smk"

##### BioRoot utilities - prepare reference #####
module PR:
    snakefile: github("BioIT-CEITEC/bioroots_utilities", path="prepare_reference.smk",branch="master")
    config: config

use rule * from PR as other_*
