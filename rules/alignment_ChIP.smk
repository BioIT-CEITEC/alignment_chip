
rule alignment_bwa:
    input:  fastq = expand("{dir}/{{sample}}{rt}.fastq.gz", dir=fastq_dir, rt=read_pair_tags),
            ref = config['organism_bwa'],
    output: bam = "mapped/{sample}.not_markDups.BWA.bam",
    log:    "logs/{sample}/alignment_bwa.log"
    params: entity_name=config["entity_name"],
    threads: 40
    conda: "../wrappers/alignment_bwa/env.yaml"
    script: "../wrappers/alignment_bwa/script.py"
    

rule alignment_bowtie2:
    input:  fastq = expand("{dir}/{{sample}}{rt}.fastq.gz", dir=fastq_dir, rt=read_pair_tags),
            ref = config['organism_bowtie2'],
    output: bam = "mapped/{sample}.not_markDups.bowtie2.bam",
    log:    "logs/{sample}/alignment_bowtie2.log"
    params: entity_name = config["entity_name"],
            protocol = config['protocol'],
            dovetailing = config['dovetailing'],
            sensitivity = config['bowtie2_sens'],
            unmapped = "mapped/{sample}.unmapped",
    threads: 40
    conda: "../wrappers/alignment_bowtie2/env.yaml"
    script: "../wrappers/alignment_bowtie2/script.py"
    

rule alignment_spike_bowtie2:
    input:  fastq = expand("{dir}/{{sample}}{rt}.fastq.gz", dir=fastq_dir, rt=read_pair_tags),
            ref = spikein_ref,
    output: bam = "mapped/{sample}.spike.not_markDups.bowtie2.bam",
    log:    "logs/{sample}/alignment_spike_bowtie2.log"
    params: entity_name = config["entity_name"],
            protocol = "spike",
            dovetailing = False,
            sensitivity = 'very',
            unmapped = "mapped/{sample}.spike.unmapped",
    threads: 40
    conda: "../wrappers/alignment_bowtie2/env.yaml"
    script: "../wrappers/alignment_bowtie2/script.py"
    
    
def index_bam_before_deduplication_input(wc):
    bam = "mapped/{sample}.not_markDups.BWA.bam"
    if config['aligner'] == "bowtie2":
      bam = "mapped/{sample}.not_markDups.bowtie2.bam"
    return bam

rule index_bam_before_deduplication:
    input:  bam = index_bam_before_deduplication_input,
    output: bam = "mapped/{sample}.not_markDups.bam",
            bai = "mapped/{sample}.not_markDups.bam.bai"
    log:    "logs/{sample}/index_bam_before_deduplication.log"
    threads: 4
    conda: "../wrappers/index_bam_before_deduplication/env.yaml"
    script: "../wrappers/index_bam_before_deduplication/script.py"
    

rule mark_duplicates:
    input:  bam = "mapped/{sample}.not_markDups.bam",
            bai = "mapped/{sample}.not_markDups.bam.bai",
    output: bam = "mapped/{sample}.markDups.bam",
            mtx = "qc_reports/{sample}/MarkDuplicates/{sample}.markDups_metrics.txt"
    log:    "logs/{sample}/mark_duplicates.log"
    resources: mem=10
    params: mark_duplicates=config["mark_duplicates"],
            UMI=config["UMI"],
            umi_usage=config["umi_usage"],
            keep_not_markDups_bam=config["keep_not_markDups_bam"],
            tmpd = GLOBAL_TMPD_PATH
    conda: "../wrappers/mark_duplicates/env.yaml"
    script: "../wrappers/mark_duplicates/script.py"


rule max_length_filter:
    input:  bam = "mapped/{sample}.markDups.bam",
    output: bam = "mapped/{sample}.markDups.maxLen.bam"
    log:    "logs/{sample}/max_length_filter.log"
    threads:  4
    params: max_len_frags = config['max_len_frags'],
    conda:  "../wrappers/max_length_filter/env.yaml"
    script: "../wrappers/max_length_filter/script.py"


rule index_and_stats:
    input:  "mapped/{sample}.markDups.maxLen.bam",
    output: bam = "mapped/{sample}.bam",
            bai = "mapped/{sample}.bam.bai",
            idxstats = "qc_reports/{sample}/index_and_stats/{sample}.idxstats.tsv",
            flagstats = "qc_reports/{sample}/index_and_stats/{sample}.flagstat.tsv",
            stats = "qc_reports/{sample}/index_and_stats/{sample}.stats.txt",
    log:    "logs/{sample}/index_and_stats.log"
    threads:    8
    conda: "../wrappers/index_and_stats/env.yaml"
    script: "../wrappers/index_and_stats/script.py"


def alignment_chip_multiqc_inputs(wc):
    inputs = {
      'bam' : expand("mapped/{sample}.bam",sample = sample_tab.sample_name),
      'idxstats' : expand("qc_reports/{sample}/index_and_stats/{sample}.idxstats.tsv",sample = sample_tab.sample_name),
      'flagstats' : expand("qc_reports/{sample}/index_and_stats/{sample}.flagstat.tsv",sample = sample_tab.sample_name),
      'stats' : expand("qc_reports/{sample}/index_and_stats/{sample}.stats.txt",sample = sample_tab.sample_name)
    }
    if config["spikein"]:
      inputs['spikeins'] = expand("mapped/{sample}.spike.bam",sample = sample_tab.sample_name)
    if config["mark_duplicates"]:
      inputs['dup_log'] = expand("qc_reports/{sample}/MarkDuplicates/{sample}.markDups_metrics.txt",sample = sample_tab.sample_name)
    return inputs

rule alignment_chip_multiqc:
    input:  unpack(alignment_chip_multiqc_inputs)
    output: html= "qc_reports/all_samples/alignment_chip_multiqc/multiqc.html"
    log:    "logs/all_samples/alignment_chip_multiqc.log"
    params: mark_duplicates = config["mark_duplicates"],
            umi_usage = config["umi_usage"],
            config = workflow.basedir+"/wrappers/alignment_chip_multiqc/multiqc_config.yaml",
    conda: "../wrappers/alignment_chip_multiqc/env.yaml"
    script: "../wrappers/alignment_chip_multiqc/script.py"

