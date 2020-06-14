from os.path import join
import pandas as pd
from scripts.Load import samplesheet

configfile: 'config.yml'
samplesheetpath: 'samplesheet.tsv'

sampledic, libdic, rundic = samplesheet(samplesheetpath)

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
    log: "logs/QC.{run}.snakelog"
    threads: 8
    cluster:
        mem  = 16,
        time = '2-00:00:00'
        partition = 'norm'
        out = 'logs/QC.{run}.log'
        error = 'logs/QC.{run}.err'
        extra = ' --gres=lscratch:40 '
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
        bam=temp("02.Alignment/Level1/{run}/{run}.sort.bam"),
        bai=temp("02.Alignment/Level1/{run}/{run}.sort.bam.bai")
    log:
        "logs/bwa_mem.{run}.log"
    params:
        index="/data/yuk5/pipeline/wgs_germline/ref/gatk_bundle_hg38_v0/Homo_sapiens_assembly38.fasta",
        para_bwa=' -K 100000000 -v 3 -R "@RG\tID:{run}\tLB:{lib}\tPL:illumina\tPU:{run}\tSM:{sample}" ',
        para_samblaster=""
        para_sambambasort=" --tmpdir /lscratch/$SLURM_JOB_ID "
    message: "Executing fastq QC with {threads} threads on the following files {input}."
    threads: 16
    resources:
        mem  = 32,
        time = '2-00:00:00'
        partition = 'norm'
        out = 'logs/QC.{run}.log'
        error = 'logs/QC.{run}.err'
        extra = ' --gres=lscratch:40 '
    shell:
        """
        module load {config[modules][bwa]} {config[modules][samblaster]} {config[modules][sambamba]}
        bwa mem -t {threads} \
            {para_bwa} \
            config[references][bwaidx] \
            {input.reads} \
        |samblaster \
            {params.para_samblaster} \
        |sambamba view -S -f bam /dev/stdin \
            -t {threads} \
        |sambamba sort /dev/stdin \
            -t {threads} \
            -o {output.bam} \
            {params.para_sambambasort}
        """
        

rule merge_level3:
    input:
        expand("02.Alignment/Level1/{lib}/{lib}.sort.bam", lib=sampledic[{sample}])
    inputs:
        " ".join(input)
    output:
        bam="02.Alignment/Level3/{sample}/{sample}.sort.md.bam",
        bai="02.Alignment/Level3/{sample}/{sample}.sort.md.bam.bai",
    threads:  4
    resources:
        mem  = 16,
        time = '2-00:00:00'
        partition = 'norm'
        out = 'logs/B2.mergelib.{lib}.log'
        error = 'logs/B2.mergelib.{lib}.err'
        extra = ' --gres=lscratch:40 '
    shell:
        if len(sampledic[{sample}]) > 1:
            """
            module load {config[modules][sambamba]}
            sambamba merge \
                -t {threads} \
                {output.bam} \
                {inputs}
            sambamba index \
                -t {threads} \
                {output.bam}
            """
        else:
            """
            mv {inputs[0]} {output.bam}
            sambamba index \
                -t {threads} \
                {output.bam}
            """


rule bqsr:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.sort.md.bam"
    output:
        metrics="02.Alignment/Level3/{sample}/{sample}.BQSR.metrics"
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
        bai="02.Alignment/Level3/{sample}/{sample}.BQSR.bai",
    threads:  4
    resources:
        mem  = 16,
        time = '2-00:00:00'
        partition = 'norm'
        out = 'logs/B2.mergelib.{lib}.log'
        error = 'logs/B2.mergelib.{lib}.err'
        extra = ' --gres=lscratch:40 '
    shell:
        """
        module load {config[modules][gatk]}
        gatk --java-options "-Xmx12000m -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" BaseRecalibrator \
            -R config[references][fasta] \
            -I {input.bam} \
            -O {output.metrics} \
            --use-original-qualities \
            --known-sites config[references][gatk_dbsnp] \
            --known-sites config[references][gatk_1000g] \
            --known-sites config[references][gatk_indel] \
            --intervals   config[references][flankbed]
        gatk --java-options "-Xmx12000m -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" ApplyBQSR \
            --add-output-sam-program-record \
            -R config[references][fasta] \
            -I {input.bam} \
            -O {output.bam} \
            --use-original-qualities \
            -bqsr {output.metrics} \
            -L config[references][flankbed]
        """            
            

