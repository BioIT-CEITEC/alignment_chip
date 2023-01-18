######################################
# wrapper for rule: index_and_stats
######################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: index_and_stats \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()


command = "$(which time) samtools sort -@ " + str(snakemake.threads) + \
          " -o " +snakemake.output.bam + " " + str(snakemake.input)+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) samtools index -@ " + str(snakemake.threads) + " " + snakemake.output.bam+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


command = "$(which time) samtools idxstats "+snakemake.output.bam+" > "+snakemake.output.idxstats+" 2>> "+log_filename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command ="$(which time) samtools flagstat "+snakemake.output.bam+" > "+snakemake.output.flagstats+" 2>> "+log_filename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


command = "rm " + str(snakemake.input)+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)
