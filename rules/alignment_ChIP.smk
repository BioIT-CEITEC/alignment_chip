

# def trim_adapters_input(wildcards):
#     if read_pair_tags == [""]:
#         return "raw_fastq/{sample}.fastq.gz"
#     else:
#         return ["raw_fastq/{sample}_R1.fastq.gz","raw_fastq/{sample}_R2.fastq.gz"]

rule automatic_preprocess:
    input:  fastq = expand("raw_fastq/{{sample}}{rt}.fastq.gz",rt=read_pair_tags),
    output: fastq = expand("cleaned_fastq/{{sample}}{rt}.auto.fastq.gz",rt = read_pair_tags),
            trim_stats = "qc_reports/{{sample}}/trim_galore/trim_stats.log"
    log:    "logs/{sample}/automatic_preprocess.log",
    threads: 8,
    params: paired = config["is_paired"],
            outdir = "cleaned_fastq",
    conda: "../wrappers/automatic_preprocess/env.yaml"
    script: "../wrappers/automatic_preprocess/script.py"
    

rule advanced_preprocess:
    input:  raw = expand("raw_fastq/{{sample}}{rt}.fastq.gz",rt=read_pair_tags),
    output: cleaned = expand("cleaned_fastq/{{sample}}{rt}.advanced.fastq.gz",rt=read_pair_tags),
            trim_stats = "qc_reports/{sample}/trimmomatic/trim_stats.log"
    log:    "logs/{sample}/advanced_preprocess.log",
    threads: 10,
    resources:  mem = 10,
    params: adaptors = config["trim_adapters"],
            r1u = "cleaned_fastq/trimmed/{sample}_R1.discarded.fastq.gz",
            r2u = "cleaned_fastq/trimmed/{sample}_R2.discarded.fastq.gz",
            trim_left1 = config["trim_left1"], # Applied only if trim left is true, trimming from R1 (different for classic:0, quant:10, sense:9)
            trim_right1 = config["trim_right1"], # Applied only if trim right is true, trimming from R1; you should allow this if you want to trim the last extra base and TRIM_LE is true as RD_LENGTH is not effective
            trim_left2 = config["trim_left2"], # Applied only if trim left is true, trimming from R2 (different for classic:0, quant:?, sense:7)
            trim_right2 = config["trim_right2"], # Applied only if trim right is true, trimming from R2; you should allow this if you want to trim the last extra base and TRIM_LE is true as RD_LENGTH is not effective
            phred = "-phred33",
            leading = 3,
            trailing = 3,
            crop = 250,
            minlen = config["min_length"],
            slid_w_1 = 4,
            slid_w_2 = 5,
    conda:  "../wrappers/advanced_preprocess/env.yaml"
    script: "../wrappers/advanced_preprocess/script.py"


# def cleaned_fastq_qc_input(wildcards):
#     preprocessed = "cleaned_fastq"
#     if read_pair_tags == ["SE"]:
#         return os.path.join(preprocessed,"{sample}.fastq.gz")
#     else:
#         return expand("cleaned_fastq/{{sample}}{rt}.fastq.gz", rt=read_pair_tags)

rule cleaned_fastq_qc:
    input:  fastq = expand("cleaned_fastq/{{sample}}{rt}.{type}.fastq.gz", rt=read_pair_tags, type="advanced" if config["preprocess"]=="trimmomatic" else "auto"),
    output: html = "qc_reports/{sample}/cleaned_fastqc/{read_pair_tags}_fastqc.html",
            fastq = expand("cleaned_fastq/{{sample}}{rt}.fastq.gz", rt=read_pair_tags)
    log:    "logs/{sample}/cleaned_fastqc_{read_pair_tags}.log"
    params: extra = "--noextract --format fastq --nogroup",
    threads:  1
    conda:  "../wrappers/cleaned_fastq_qc/env.yaml"
    script: "../wrappers/cleaned_fastq_qc/script.py"


def alignment_chip_input(wc):
    inputs = {
      'ref': expand("{ref_dir}/index/BWA/{ref}.bwt",ref_dir=reference_directory,ref = config["reference"])[0]
    }
    preprocessed = "cleaned_fastq" if config["preprocess"] != "None" else "raw_fastq"
    inputs['fastqs'] = expand("{dir}/{{sample}}{rt}.fastq.gz". dir=preprocessed, rt=read_pair_tags)
    # inputs['fastqc'] = expand("qc_reports/{{sample}}/{dir}c/{rt}_fastqc.html". dir=preprocessed, rt=pair_tag)
    return inputs


rule alignment_chip:
    input:  unpack(alignment_chip_input),
    # input:  fastqs = alignment_chip_input,
    #         ref = expand("{ref_dir}/index/BWA/{ref}.bwt",ref_dir=reference_directory,ref = config["reference"])[0],
    output: bam = "mapped/{sample}.not_markDups.bam",
            bai = "mapped/{sample}.not_markDups.bam.bai"
    log:    "logs/{sample}/alignment_chip.log"
    params: entity_name=config["entity_name"]
    threads: 40
    conda: "../wrappers/alignment_chip/env.yaml"
    script: "../wrappers/alignment_chip/script.py"


rule mark_duplicates:
    input:  bam = "mapped/{sample}.not_markDups.bam",
            bai = "mapped/{sample}.not_markDups.bam.bai",
    output: bam = "mapped/{sample}.markDups.bam",
            mtx = "qc_reports/{sample}/MarkDuplicates/{sample}.markDups_metrics.txt"
    log:    "logs/{sample}/mark_duplicates.log"
    resources: mem=10
    params: mark_duplicates=config["mark_duplicates"],
            rmDup=config["remove_duplicates"],
            UMI=config["UMI"],
            umi_usage=config["umi_usage"],
            keep_not_markDups_bam=config["keep_not_markDups_bam"],
            tmpd = GLOBAL_TMPD_PATH
    conda: "../wrappers/mark_duplicates/env.yaml"
    script: "../wrappers/mark_duplicates/script.py"


# rule umi_concensus:
#     input:  bam = "mapped/{sample}.not_markDups.bam",
#             bai = "mapped/{sample}.not_markDups.bam.bai",
#             ref = expand("{ref_dir}/index/BWA/{ref}.bwt",ref_dir=reference_directory,ref=config["reference"])[0],
#             lib_ROI = expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
#             fa = expand("{ref_dir}/seq/{ref}.fa",ref_dir=reference_directory,ref=config["reference"])[0],
#     output: bam = "mapped/{sample}.concensus.bam",
#             html = "qc_reports/{sample}/umi_concensus/umi_concensus.html",
#             json = "qc_reports/{sample}/umi_concensus/umi_concensus.json",
#     log: "logs/{sample}/umi_concensus.log"
#     params: umi_consensus_min_support = config["umi_consensus_min_support"],
#             keep_not_markDups_bam=config["keep_not_markDups_bam"],
#             tmpd = GLOBAL_TMPD_PATH
#     conda: "../wrappers/umi_concensus/env.yaml"
#     script: "../wrappers/umi_concensus/script.py"


# def index_and_stats_input(wildcards):
#     if not config["umi_usage"] == "umi_concensus":
#         return "mapped/{sample}.markDups.bam"
#     else:
#         return "mapped/{sample}.concensus.bam"

rule index_and_stats:
    # input:  index_and_stats_input,
    input:  "mapped/{sample}.markDups.bam",
    output: bam = "mapped/{sample}.bam",
            bai = "mapped/{sample}.bam.bai",
            idxstats = "qc_reports/{sample}/index_and_stats/{sample}.idxstats.tsv",
            flagstats = "qc_reports/{sample}/index_and_stats/{sample}.flagstat.tsv",
    log:    "logs/{sample}/index_and_stats.log"
    threads:    8
    conda: "../wrappers/index_and_stats/env.yaml"
    script: "../wrappers/index_and_stats/script.py"


rule alignment_chip_multiqc:
    input:  bam = expand("mapped/{sample}.bam",sample = sample_tab.sample_name),
            idxstats = expand("qc_reports/{sample}/index_and_stats/{sample}.idxstats.tsv",sample = sample_tab.sample_name),
            fastqc = expand("qc_reports/{sample}/{dir}/{rt}_fastqc.html". dir="cleaned_fastqc" if config['preprocess']!="none" else "raw_fastqc", rt=pair_tag, sample=sample_tab.sample_name)
    output: html= "qc_reports/all_samples/alignment_chip_multiqc/multiqc.html"
    log:    "logs/all_samples/alignment_chip_multiqc.log"
    params: trim_adapters=config["trim_adapters"],
            mark_duplicates=config["mark_duplicates"],
            umi_usage = config["umi_usage"]
    conda: "../wrappers/alignment_chip_multiqc/env.yaml"
    script: "../wrappers/alignment_chip_multiqc/script.py"

