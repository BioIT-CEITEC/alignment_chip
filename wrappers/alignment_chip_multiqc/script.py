######################################
# wrapper for rule: alignment_chip_multiqc
######################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: alignment_chip_multiqc \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

multiqc_search_paths = "./mapped/*" + " ./qc_reports/*/index_and_stats/*" + " ./qc_reports/*/"+snakemake.params.fastqc_dir+"/*"

if snakemake.params.trim_adapters:
    multiqc_search_paths += " ./qc_reports/*/trim_galore/*"
if snakemake.params.mark_duplicates:
    multiqc_search_paths += " ./qc_reports/*/MarkDuplicates/*"


command = "multiqc -f -n " + snakemake.output.html +\
          " --config "+snakemake.params.config+\
          " " + multiqc_search_paths + \
          " >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
