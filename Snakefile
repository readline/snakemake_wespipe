from os.path import join
import pandas as pd
from scripts.Load import samplesheet

configfile: 'config.yaml'
print(config)
sampledic, libdic, rundic = samplesheet(config['samplesheet'])
workdir: config['workdir']
#sbcmd="sbatch --cpus-per-task={threads} --mem={cluster.mem}"
#sbcmd+=" --time={cluster.time} --partition={cluster.partition}"
#sbcmd+=" --out={cluster.out} {cluster.extra}"

"""
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
        """
rule all:
    input:
        expand("02.Alignment/Level3/{sample}/{sample}.BQSR.bam", sample=list(sampledic.keys())),
        
rule QC:
    input:
        reads1=lambda wildcards: rundic[wildcards.run]['Read1'],
        reads2=lambda wildcards: rundic[wildcards.run]['Read2'],
    output:
        reads1out="01.CleanData/{run}/{run}.R1.cln.fq.gz",
        reads2out="01.CleanData/{run}/{run}.R2.cln.fq.gz",
        htmlout="01.CleanData/{run}/{run}.QC.html",
        jsonout="01.CleanData/{run}/{run}.QC.json"
    log: "logs/QC.{run}.snakelog"
    threads: 8
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
        para_bwa=lambda wildcards:' -K 100000000 -v 3 -R "@RG\tID:{}\tLB:{}\tPL:illumina\tPU:{}\tSM:{}" '.format(
                                wildcards.run,  rundic[wildcards.run]['Lib'],  wildcards.run,  rundic[wildcards.run]['Sample']),
        para_samblaster="",
        para_sambambasort=" --tmpdir /lscratch/$SLURM_JOB_ID "
    threads: 16
    message: "Executing bwa_mem with {threads} threads on the run {wildcards.run} with fastq files: {input}."
    resources: mem_mb = 32*1024
    shell:
        """
        module load {config[modules][bwa]} {config[modules][samblaster]} {config[modules][sambamba]}
        bwa mem -t {threads} \
            {params.para_bwa} \
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


rule markdup:
    input:
        lambda wildcards: \
            ["02.Alignment/Level1/{}/{}.sort.bam".format(run, run) \
                for run in libdic[wildcards.lib] ]
    output:
        bam="02.Alignment/Level2/{lib}/{lib}.sort.md.bam",
        bai="02.Alignment/Level2/{lib}/{lib}.sort.md.bam.bai",
        metrics="02.Alignment/Level2/{lib}/{lib}.sort.md.metrics",
    threads:  4
    resources:
        mem_mb  = 16*1024,
    message: "Executing markdup with {threads} threads on the lib {wildcards.lib} bam files: {input}."
    run:
        inputs = " ".join("INPUT={}".format(f) for f in snakemake.input)
        shell(
            """
            module load {config[modules][picard]} {config[modules][sambamba]}
            java -Xmx15G  -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID  -jar $PICARD_JAR MarkDuplicates \
                {inputs} \
                OUTPUT={output.bam} \
                METRICS_FILE={output.metrics} \
                MAX_RECORDS_IN_RAM=2000000 \
                VALIDATION_STRINGENCY=SILENT \
                OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
                ASSUME_SORT_ORDER=coordinate \
                CLEAR_DT=false \
                ADD_PG_TAG_TO_READS=false
            sambamba index -t 4 {output.bam}
            """)

rule merge_level3:
    input:
        lambda wildcards: \
            ["02.Alignment/Level2/{0}/{0}.sort.md.bam".format(lib, lib) \
                for lib in sampledic[wildcards.sample] ]
    output:
        bam="02.Alignment/Level3/{sample}/{sample}.sort.md.bam",
        bai="02.Alignment/Level3/{sample}/{sample}.sort.md.bam.bai",
    threads:  4
    resources:
        mem_mb  = 16*1024,
    message: "Executing merge bam with {threads} threads on the sample {wildcards.sample} bam files: {input}."
    run:
        if len(sampledic[{sample}]) > 1:
            shell("""
            module load {config[modules][sambamba]}
            sambamba merge \
                -t {threads} \
                {output.bam} \
                {inputs}
            sambamba index \
                -t {threads} \
                {output.bam}
            """)
        else:
            shell("""
            mv {inputs[0]} {output.bam}
            sambamba index \
                -t {threads} \
                {output.bam}
            """)


rule bqsr:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.sort.md.bam"
    output:
        metrics="02.Alignment/Level3/{sample}/{sample}.BQSR.metrics",
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
        bai="02.Alignment/Level3/{sample}/{sample}.BQSR.bai",
    threads:  4
    resources:
        mem  = 16,
    shell:
        """
        module load {config[modules][gatk]}
        gatk --java-options "-Xmx12000m -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" BaseRecalibrator \
            -R {config[references][fasta]} \
            -I {input.bam} \
            -O {output.metrics} \
            --use-original-qualities \
            --known-sites {config[references][gatk_dbsnp]} \
            --known-sites {config[references][gatk_1000g]} \
            --known-sites {config[references][gatk_indel]} \
            --intervals   {config[references][flankbed]}
        gatk --java-options "-Xmx12000m -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" ApplyBQSR \
            --add-output-sam-program-record \
            -R {config[references][fasta]} \
            -I {input.bam} \
            -O {output.bam} \
            --use-original-qualities \
            -bqsr {output.metrics} \
            -L {config[references][flankbed]}
        """            

        
rule stat_bqsr:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
    output:
        stat1="02.Alignment/Level3/{sample}/stat/{sample}.stat1",
        stat2="02.Alignment/Level3/{sample}/stat/{sample}.stat2",
        stat3="02.Alignment/Level3/{sample}/stat/{sample}.stat3",
        flags="02.Alignment/Level3/{sample}/{sample}.BQSR.bam.flagstat",
        metrics="02.Alignment/Level3/{sample}/{sample}.BQSR.metrics",
    threads:  4
    resources:
        mem  = 16,
    shell:
        """
        module load {config[modules][picard]} {config[modules][sambamba]} {config[modules][samtools]}
        java -Xmx8000m -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $PICARD_JAR CollectMultipleMetrics \
            INPUT={input.bam} \
            REFERENCE_SEQUENCE={config[references][fasta]} \
            OUTPUT={output.stat1} \
            ASSUME_SORTED=true \
            PROGRAM="null" \
            PROGRAM="CollectAlignmentSummaryMetrics" \
            PROGRAM="CollectInsertSizeMetrics" \
            PROGRAM="CollectSequencingArtifactMetrics" \
            PROGRAM="CollectGcBiasMetrics" \
            PROGRAM="QualityScoreDistribution" \
            METRIC_ACCUMULATION_LEVEL="null" \
            METRIC_ACCUMULATION_LEVEL="SAMPLE" \
            METRIC_ACCUMULATION_LEVEL="LIBRARY"
        java -Xmx8000m -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $PICARD_JAR CollectMultipleMetrics \
            INPUT={input.bam} \
            REFERENCE_SEQUENCE={config[references][fasta]} \
            OUTPUT={output.stat2} \
            ASSUME_SORTED=true \
            PROGRAM="null" \
            PROGRAM="CollectAlignmentSummaryMetrics" \
            PROGRAM="CollectGcBiasMetrics" \
            METRIC_ACCUMULATION_LEVEL="null" \
            METRIC_ACCUMULATION_LEVEL="READ_GROUP" &&
        sambamba flagstat -t {threads} {input.bam} > {output.flags}
        /data/yuk5/script/BamStat/bam_cov_estimate.py \
            {input.bam} \
            wxs {output.stat3} \
            {threads} 
        """        
        
rule HaplotypeCaller:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
    output:
        gvcf="03.Germline/{sample}/{sample}.gvcf.gz",
        stat="03.Germline/{sample}/stat/{sample}.germline"
    threads:  4
    resources:
        mem  = 16,
    shell:
        """
        module load {config[modules][gatk]} {config[modules][picard]}
        gatk --java-options "-Xmx8000m -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" HaplotypeCaller \
            -R {config[references][fasta]} \
            -O {output.gvcf} \
            -I {input.bam} \
            --emit-ref-confidence GVCF \
            --max-alternate-alleles 3 \
            --read-filter OverclippedReadFilter \
            --dbsnp {config[references][gatk_dbsnp]} \
            -L {config[references][flankbed]}
        gatk --java-options "-Xmx8000m -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" ValidateVariants \
            -V {output.gvcf} \
            -R {config[references][fasta]} \
            -gvcf \
            --validation-type-to-exclude ALLELES \
            --dbsnp {config[references][gatk_dbsnp]} \
            -L {config[references][flankbed]}
        java -Xms2000m -Xmx8000m -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $PICARD_JAR \
            CollectVariantCallingMetrics \
            INPUT={output.gvcf} \
            OUTPUT={output.stat} \
            DBSNP={config[references][gatk_dbsnp]} \
            SEQUENCE_DICTIONARY={config[references][fadict]} \
            TARGET_INTERVALS={config[references][itvlist]} \
            GVCF_INPUT=true 
        """
        
        
rule GenomicsDBImport:
    input:
        bed="/data/yuk5/pipeline/wxs_phase1/ref/IDT_xGen_Exome_Research_Panel/hg38/intervals_withflank150/hg38.xGen.{itv}.bed",
        gvcf=lambda wildcards: ["03.Germline/{}/{}.gvcf.gz".format(sample, sample) for sample in sampledic ]
    output:
        itvvcf=temp("03.Germline/VQSR/Merge.itv_{itv}.vcf.gz"),
        itvmf=temp("03.Germline/VQSR/Merge.itv_{itv}.mf.vcf.gz"),
        itvflt="03.Germline/VQSR/Merge.itv_{itv}.flt.vcf.gz"
    threads:  4
    resources:
        mem  = 16
    run:
        inputs = " ".join("-V {}".format(f) for f in snakemake.input.gvcf)
        shell(
        """
        module load {config[modules][gatk]}
        gatk --java-options "-Xmx4g -Xms4g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            GenomicsDBImport \
            {inputs} \
            -L /data/yuk5/pipeline/wxs_phase1/ref/IDT_xGen_Exome_Research_Panel/hg38/intervals_withflank150/hg38.xGen.{itv}.bed \
            --genomicsdb-workspace-path /lscratch/$SLURM_JOB_ID/gdb.itv_{itv} &&
        gatk --java-options "-Xmx5g -Xms5g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            GenotypeGVCFs \
            -R {config[references][fasta]} \
            -O {output.itvvcf} \
            -D {config[references][gatk_dbsnp]} \
            -G StandardAnnotation \
            --only-output-calls-starting-in-intervals \
            --use-new-qual-calculator \
            -V gendb:///lscratch/$SLURM_JOB_ID/gdb.itv_3 \
            -L /data/yuk5/pipeline/wxs_phase1/ref/IDT_xGen_Exome_Research_Panel/hg38/intervals_withflank150/hg38.xGen.{itv}.bed &&
        gatk --java-options "-Xmx3g -Xms3g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            VariantFiltration \
            --filter-expression "ExcessHet>54.69" \
            --filter-name ExcessHet \
            -V {output.itvvcf} \
            -O {output.itvmf}
        gatk --java-options "-Xmx3g -Xms3g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            SelectVariants \
            --exclude-filtered \
            -V {output.itvmf} \
            -O {output.itvflt}
        """
        )
            
            
            
            
rule genotyping:
    input:
        itvs=lambda wildcards: \
             ["03.Germline/VQSR/Merge.itv_{}.flt.vcf.gz".format(itv) for itv in range(1, config['references']['interval']+1)]
    output:
        mvcf1="03.Germline/VQSR/Merge.flt.vcf.gz",
        sidrecal="03.Germline/VQSR/Merge.flt.sid.recal",
        snvrecal="03.Germline/VQSR/Merge.flt.snv.recal",
        sidtranch="03.Germline/VQSR/Merge.flt.sid.tranches",
        snvtranch="03.Germline/VQSR/Merge.flt.snv.tranches",
        indrecal="03.Germline/VQSR/Merge.flt.tmp.indel.recalibrated.vcf.gz",
        vqsr="03.Germline/Merge.flt.vqsr.vcf.gz"
    threads:  4
    resources:
        mem  = 16,
    run:
        inputs = " ".join("--INPUT {}".format(i) for i in snakemake.input.itvs)
        shell(
        """
        module load {config[modules][gatk]} {config[modules][samtools]}
        
        gatk --java-options "-Xmx3g -Xms3g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            GatherVcfs \
            {inputs}
            --OUTPUT {output.mvcf1}
            
        tabix -p vcf {output.mvcf1}
        
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            VariantRecalibrator \
            -R {config[references][fasta]} \
            -V {output.mvcf1} \
            -O {output.sidrecal} \
            --tranches-file {output.sidtranch} \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 \
            -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 \
            -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 \
            -tranche 91.0 -tranche 90.0 \
            -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
            -mode INDEL \
            --max-gaussians 4 \
            -resource mills,known=false,training=true,truth=true,prior=12:{config[references][gatk_1000g]} \
            -resource axiomPoly,known=false,training=true,truth=false,prior=10:{config[references][gatk_axiom]} \
            -resource dbsnp,known=true,training=false,truth=false,prior=2:{config[references][gatk_dbsnp]}

        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            VariantRecalibrator \
            -V {output.mvcf1} \
            -O {output.snvrecal} \
            --tranches-file {output.snvtranch} \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 \
            -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
            -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
            -mode SNP \
            --max-gaussians 6 \
            -resource hapmap,known=false,training=true,truth=true,prior=15:{config[references][gatk_hapmap]} \
            -resource omni,known=false,training=true,truth=true,prior=12:{config[references][gatk_omni]} \
            -resource 1000G,known=false,training=true,truth=false,prior=10:{config[references][gatk_1000hc]} \
            -resource dbsnp,known=true,training=false,truth=false,prior=7:{config[references][gatk_dbsnp]}
            
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            ApplyVQSR \
            -V {output.mvcf1} \
            -O {output.indrecal} \
            --recal-file {output.sidrecal} \
            --tranches-file {output.sidtranch} \
            --truth-sensitivity-filter-level 99.7 \
            --create-output-variant-index true \
            -mode INDEL
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            ApplyVQSR \
            -V {output.indrecal} \
            -O {output.vqsr} \
            --recal-file {output.snvrecal} \
            --tranches-file {output.snvtranch} \
            --truth-sensitivity-filter-level 99.7 \
            --create-output-variant-index true \
            -mode SNP
        """
        )

            
rule annotation:
    input:
        "03.Germline/Merge.flt.vqsr.vcf.gz"
    output:
        "03.Germline/Merge.flt.vqsr.vcf.anno"
    threads:  48
    resources:
        mem  = 16,
    shell:
        """
        {config[bins][vcfanno]} \
            {input} \
            {output} \
            {threads} n 
        """
        