######################################
# wrapper for rule: alignment_bowtie2
######################################
import re
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: alignment_bowtie2 \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

index_prefix = re.sub(".1.bt2$","",snakemake.input.ref)

input_reads = ""
if len(snakemake.input.fastq) == 2:
  input_reads =  " -1 "+snakemake.input.fastq[0]
  input_reads += " -2 "+snakemake.input.fastq[1]
  unmapped_reads = " --un-conc-gz "+snakemake.params.unmapped+"%.fastq.gz"
else:
  input_reads = " -U "+snakemake.input.fastq[0]
  unmapped_reads = " --un-gz "+snakemake.params.unmapped+".fastq.gz"
  
special_arguments = " --no-mixed --no-discordant --phred33 -I 10 -X 1000"
if snakemake.params.protocol == "car":
  special_arguments += " --end-to-end"
else:
  special_arguments += " --local"

if snakemake.params.sensitivity == "very":
  special_arguments += " --very-sensitive"
else:
  special_arguments += " --sensitive"
  
if snakemake.params.dovetailing:
  special_arguments += " --dovetail"

command = "$(which time) bowtie2 -t -p "+str(snakemake.threads)+\
          " -x " + index_prefix +\
          " " + input_reads + unmapped_reads + special_arguments +\
          " 2>> " + log_filename + \
          " | $(which time) samtools sort -@ " +str(snakemake.threads)+" -o "+snakemake.output.bam+" /dev/stdin 2>> "+log_filename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

