######################################
# wrapper for rule: cleaned_fastq_qc
######################################
import subprocess
from os.path import dirname
from os.path import basename
from snakemake.shell import shell
import re
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: cleaned_fastq_qc \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

command = "mkdir -p " + dirname(snakemake.output.html) + " >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

tags = ""
if re.search("_R1.fastq",str(snakemake.input.fastq)):
    tags = "_R1"
elif re.search("_R2.fastq",str(snakemake.input.fastq)):
    tags = "_R2"

tag_input_fastq = dirname(snakemake.input.fastq[0]) + "/" + snakemake.wildcards.sample + tags + ".fastq.gz"

f = open(log_filename, 'at')
f.write("## tag_input_fastq: "+tag_input_fastq+"\n")
f.close()

command = "$(which time) fastqc -o " + dirname(snakemake.output.html) + " "+snakemake.params.extra+\
          " --threads "+str(snakemake.threads)+" "+tag_input_fastq+" >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv " + dirname(snakemake.output.html) + "/" + basename(tag_input_fastq).replace(".fastq.gz","_fastqc.html") +\
          " " + snakemake.output.html+" >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv " + dirname(snakemake.output.html) + "/" + basename(tag_input_fastq).replace(".fastq.gz","_fastqc.zip") + \
          " " + snakemake.output.arch+" >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
