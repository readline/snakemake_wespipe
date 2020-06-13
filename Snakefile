from os.path import join
import pandas as pd

GENOME = 'genome/HWB.chromosome.fa'
GTF = 'genes/HWB.gene.models.gtf'


(SAMPLES,) = glob_wildcards('pairedDIR/{sample}_1P.fq.gz')
PATTERN_R1 = join('pairedDIR', '{sample}_1P.fq.gz')
PATTERN_R2 = join('pairedDIR', '{sample}_2P.fq.gz')



#sbcmd="sbatch --cpus-per-task={threads} --mem={cluster.mem}"
#sbcmd+=" --time={cluster.time} --partition={cluster.partition}"
#sbcmd+=" --out={cluster.out} {cluster.extra}"


rule all:
    input:
        expand("01.CleanData/{run}/{run}.1.cln.fq.gz", run=RUN),
        expand("01.CleanData/{run}/{run}.2.cln.fq.gz", run=RUN),
        expand("02.Alignment/Level3/{sample}/{sample}.BQSR.bam", sample=SAMPLE),
        expand("03.Germline/{sample}/{sample}.gvcf.gz", sample=SAMPLES),
        expand("03.Germline.samtools/{sample}/{sample}.gsam.vcf.gz", sample=RUN),
        "03.Germline/Merge.flt.vqsr.vcf.gz",
        "03.Germline.samtools/Merge.gsam.flt.snv.vcf.gz",
        "03.Germline.samtools/Merge.gsam.flt.indel.vcf.gz",
        "03.Germline/Merge.flt.vqsr.vcf.anno/Merge.Anno.matrix.gz",
        "03.Germline.samtools/Merge.gsam.flt.snv.vcf.anno/Merge.Anno.matrix.gz"
        
rule QC:
    input:
        reads1="00.RawData/{run}/{run}.R1.fq.gz",
        reads2="00.RawData/{run}/{run}.R2.fq.gz"
    output:
        reads1out="01.CleanData/{run}/{run}.R1.cln.fq.gz",
        reads2out="01.CleanData/{run}/{run}.R2.cln.fq.gz",
        htmlout="01.CleanData/{run}/{run}.QC.html",
        jsonout="01.CleanData/{run}/{run}.QC.json"
    log: "logs/QC.{run}.log"
    threads: 8
    resources:
        mem  = 16,
        time = '2-00:00:00'
    message: "Executing fastq QC with {threads} threads on the following files {input}."
    shell:
        """
        /data/yuk5/script/fastp \
            -i {input.reads1} \
            -I {input.reads2} \
            -o {output.reads1out} \
            -O {output.reads2out} \
            -h {output.htmlout} \
            -j {output.jsonout} \
            -w {threads}
        """

rule bwa_mem:
    input:
        reads=["01.CleanData/{run}/{run}.R1.cln.fq.gz", "01.CleanData/{run}/{run}.R2.cln.fq.gz"]
    output:
        runbam=temp("02.Alignment/Level1/{run}/{run}.sort.bam"),
        runindex=temp("02.Alignment/Level1/{run}/{run}.sort.bam.bai")
    log:
        "logs/bwa_mem.{run}.log"
    params:
        index="/data/yuk5/pipeline/wgs_germline/ref/gatk_bundle_hg38_v0/Homo_sapiens_assembly38.fasta",
        para_bwa=' -K 100000000 -v 3 -R "@RG\tID:{run}\tLB:{lib}\tPL:illumina\tPU:{run}\tSM:{sample}" ',
        para_samblaster=""
        para_sambambasort=" --tmpdir /lscratch/$SLURM_JOB_ID "
    threads: 16
    resources:
        mem  = 32,
        time = '2-00:00:00'
        gres = 'lscratch:40'
    wrapper:
        "wrapper/bwa_mem"
        
rule merge_bam:
    input:
        ["mapped/A.bam", "mapped/B.bam"]
    output:
        "merged.bam"
    params:
        "" # optional additional parameters as string
    threads:  8
    resources:
        mem  = 32,
        time = '2-00:00:00'
        gres = 'lscratch:40'
    wrapper:
        "0.60.0/bio/samtools/merge"


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample}"
        " -O bam {input} > {output}"
        
        
SBT=["wt1","wt2","epcr1","epcr2"]

rule all:
    input:
        expand("02_clean/{nico}_1.paired.fq.gz", nico=SBT),
        expand("02_clean/{nico}_2.paired.fq.gz", nico=SBT),
        expand("03_align/{nico}.sam", nico=SBT),
        expand("04_exp/{nico}_count.txt", nico=SBT),
        expand("05_ft/{nico}_gene.gtf", nico=SBT),
        expand("05_ft/{nico}_transcript.gtf", nico=SBT)

rule trim:
    input:
        "01_raw/{nico}_1.fastq",
        "01_raw/{nico}_2.fastq"
    output:
        "02_clean/{nico}_1.paired.fq.gz",
        "02_clean/{nico}_1.unpaired.fq.gz",
        "02_clean/{nico}_2.paired.fq.gz",
        "02_clean/{nico}_2.unpaired.fq.gz",
    threads: 20
    shell:
        "java -jar /software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads {threads} {input[0]} {input[1]} {output[0]} {output[1]} {output[2]} {output[3]} ILLUMINACLIP:/software/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 "

rule map:
    input:
        "02_clean/{nico}_1.paired.fq.gz",
        "02_clean/{nico}_2.paired.fq.gz"
    output:
        "03_align/{nico}.sam"
    log:
        "logs/map/{nico}.log"
    threads: 20
    shell:
        "hisat2 -p {threads} --dta -x /root/s/r/p/A_th/WT-Al_VS_WT-CK/index/tair10 -1 {input[0]} -2 {input[1]} -S {output} >{log} 2>&1 "
        
rule sort2bam:
    input:
        "03_align/{nico}.sam"
    output:
        "03_align/{nico}.bam"
    threads: 8
    shell:
        "samtools sort -@ {threads} -m 2G -o {output} {input}"

rule count:
    input:
        "03_align/{nico}.bam"
    output:
        "04_exp/{nico}_count.txt"
    threads: 20
    shell:
        "featureCounts -T {threads} -p -t exon -g gene_id -a /root/s/r/p/A_th/WT-Al_VS_WT-CK/genome/tair10.gtf -o {output} {input}"

rule fpkm:
    input:
        "03_align/{nico}.bam"
    output:
        "05_ft/{nico}_gene.gtf",
        "05_ft/{nico}_transcript.gtf"
    threads: 20
    shell:
        "stringtie -e -p {threads} -G /root/s/r/p/A_th/WT-Al_VS_WT-CK/genome/tair10.gtf -A {output[0]} -o {output[1]} {input}"