######################################
# wrapper for rule: max_length_filter
######################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: max_length_filter \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()


command = "$(which time) samtools view -hb -e 'tlen < "+str(snakemake.params.max_len_frags)+" && tlen > -"+str(snakemake.params.max_len_frags)+"' -@ " + str(snakemake.threads) + \
          " -o " +snakemake.output.bam + " " + snakemake.input.bam+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
