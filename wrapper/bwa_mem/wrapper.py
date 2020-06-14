from os import path
from snakemake.shell import shell


# Extract arguments.
para_bwa = snakemake.params.get("para_bwa", "")
para_sambambasort = snakemake.params.get("para_sambambasort", "")
para_samblaster = snakemake.params.get("para_samblaster", "")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Check inputs/arguments.
if not isinstance(snakemake.input.reads, str) and len(snakemake.input.reads) not in {1,2}:
    raise ValueError("Input must have 1 (SE) or 2 (PE) input fastq files.")

shell(
    "module load bwa/0.7.17 samblaster/0.1.22 sambamba"
    "(bwa mem"
    " -t {snakemake.threads}"
    " {para_bwa}"
    " {snakemake.params.index}"
    " {snakemake.input.reads}"
    " | samblaster"
    " {para_samblaster}"
    " | sambamba view -S -f bam /dev/stdin"
    " -t {snakemake.threads}"
    " | sambamba sort /dev/stdin"
    " -t {snakemake.threads}"
    " -o {snakemake.output.bam}"
    " {para_sambambasort}"
    ") {log}"
)