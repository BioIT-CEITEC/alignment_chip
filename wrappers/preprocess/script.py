######################################
# wrapper for rule: preprocess
######################################
import subprocess
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: preprocess \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

paired_flag = ""
r1 = snakemake.output.fastq[0]
if snakemake.params.paired:
    paired_flag = " --paired"
    r2 = snakemake.output.fastq[1]
    
adapters = ""
if snakemake.params.adaptors != "":
    adaptors = snakemake.params.adaptors.split(",")
    for i, adaptor in enumerate(adaptors):
        if adaptor.split("-")[0] == "1":
            adapters += " -a "+adaptor.split("-")[1]
        elif adaptor.split("-")[0] == "2":
            adapters += " -a2 "+adaptor.split("-")[1]
        else:
            adapters += " -a "+adaptor
          
trim_ends = ""
trim_ends+= " --clip_R1 "+str(snakemake.params.trim_left1) if int(snakemake.params.trim_left1) > 0 else ""
trim_ends+= " --three_prime_clip_R1 "+str(snakemake.params.trim_right1) if int(snakemake.params.trim_right1) > 0 else ""
max_length = " --max_length "+str(snakemake.params.crop)
if snakemake.params.paired:
    trim_ends+= " --clip_R2 "+str(snakemake.params.trim_left2) if int(snakemake.params.trim_left2) > 0 else ""
    trim_ends+= " --three_prime_clip_R2 "+str(snakemake.params.trim_right2) if int(snakemake.params.trim_right2) > 0 else ""
    #max_length = ""
    
# ABout cores usage see https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
command = "$(which time) trim_galore --gzip" +\
          paired_flag +\
          " --"+snakemake.params.phred +\
          " --quality "+str(snakemake.params.qual)+\
          adapters+\
          max_length+\
          " --length "+str(snakemake.params.minlen)+\
          " --trim-n"+\
          trim_ends+\
          " --retain_unpaired"+\
          " -j 4"+\
          " "+str(snakemake.input.fastq) +\
          " -o "+snakemake.params.outdir+\
          " --basename "+snakemake.wildcards.sample+\
          " >> "+log_filename+" 2>&1"
with open(log_filename, 'at') as f:
        f.write("## COMMAND: " + command + "\n")
shell(command)

prefix = os.path.join(snakemake.params.outdir, snakemake.wildcards.sample)
if snakemake.params.paired:
    command = "mv " + prefix + "_val_1.fq.gz" + " " + r1 + " >> "+log_filename+" 2>&1"
    with open(log_filename, 'at') as f:
            f.write("## COMMAND: " + command + "\n")
    shell(command)

    command = "mv " + prefix + "_val_2.fq.gz" + " " + r2 + " >> "+log_filename+" 2>&1"
    with open(log_filename, 'at') as f:
            f.write("## COMMAND: " + command + "\n")
    shell(command)

    # report for multiQC
    command = "cp " + prefix + "_R1.fastq.gz_trimming_report.txt" + " " + snakemake.output.trim_stats[0] + " >> "+log_filename+" 2>&1"
    with open(log_filename, 'at') as f:
            f.write("## COMMAND: " + command + "\n")
    shell(command)

    command = "cp " + prefix + "_R2.fastq.gz_trimming_report.txt" + " " + snakemake.output.trim_stats[1] + " >> "+log_filename+" 2>&1"
    with open(log_filename, 'at') as f:
            f.write("## COMMAND: " + command + "\n")
    shell(command)
    
    if os.path.isfile(prefix+"_R1_unpaired_1.fq.gz"):
        command = "mv " + prefix + "_R1_unpaired_1.fq.gz" + " " + snakemake.params.r1u + " >> "+log_filename+" 2>&1"
        with open(log_filename, 'at') as f:
            f.write("## COMMAND: " + command + "\n")
        shell(command)
        
    if os.path.isfile(prefix+"_R2_unpaired_2.fq.gz"):
        command = "mv " + prefix + "_R2_unpaired_2.fq.gz" + " " + snakemake.params.r2u + " >> "+log_filename+" 2>&1"
        with open(log_filename, 'at') as f:
            f.write("## COMMAND: " + command + "\n")
        shell(command)

else:
    command = "mv " + prefix + "_trimmed.fq.gz" + " " + r1 + " >> "+log_filename+" 2>&1"
    with open(log_filename, 'at') as f:
        f.write("## COMMAND: " + command + "\n")
    shell(command)

    # report for multiQC
    command = "cp " + prefix + ".fastq.gz_trimming_report.txt" + " " + snakemake.output.trim_stats[0] + " >> "+log_filename+" 2>&1"
    with open(log_filename, 'at') as f:
        f.write("## COMMAND: " + command + "\n")
    shell(command)

