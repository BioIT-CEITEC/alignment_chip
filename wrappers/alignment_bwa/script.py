######################################
# wrapper for rule: alignment_chip
######################################
import re
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: alignment_chip \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

bwa_ref_prefix = re.sub(".bwt$","",snakemake.input.ref)

command = "$(which time) bwa mem -t "+str(snakemake.threads)+\
          " -R '@RG\\tID:"+snakemake.params.entity_name+"_"+snakemake.wildcards.sample+"\\tSM:"+snakemake.wildcards.sample+"\\tPL:illumina' -v 1" +\
          " "+bwa_ref_prefix + " " + " ".join(snakemake.input.fastq) + " 2>> " + log_filename + \
          " | $(which time) samtools sort -@ " +str(snakemake.threads)+" -o "+snakemake.output.bam+" /dev/stdin 2>> "+log_filename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# command = "$(which time) samtools index -@ " +str(snakemake.threads)+" -b "+snakemake.output.bam+" >> "+log_filename+ " 2>&1"
# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
