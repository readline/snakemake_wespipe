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


rule all:
    input:
        #lambda wildcards: \
        #    ["01.CleanData/{}/{}.1.cln.fq.gz".format(run, run) \
        #        for run in rundic ],
        #lambda wildcards: \
        #    ["01.CleanData/{}/{}.2.cln.fq.gz".format(run, run) \
        #        for run in rundic ],
        #lambda wildcards: \
        #    ["02.Alignment/Level3/{}/{}.BQSR.bam".format(sample, sample) \
        #        for sample in sampledic ],
        #lambda wildcards: \
        #    ["03.Germline/{}/{}.gvcf.gz".format(sample, sample) \
        #        for sample in sampledic ],
        #"03.Germline/Merge.flt.vqsr.vcf.gz",
        #"03.Germline.samtools/Merge.gsam.flt.snv.vcf.gz",
        #"03.Germline.samtools/Merge.gsam.flt.indel.vcf.gz",
        "03.Germline/Merge.flt.vqsr.vcf.anno/Merge.Anno.matrix.gz",
        lambda wildcards: ["02.Alignment/Level3/{}/{}.BQSR.bam.flagstat".format(sample, sample) for sample in sampledic],
        lambda wildcards: ["02.Alignment/Callable/{}/{}.bed".format(sample, sample) for sample in sampledic],
        lambda wildcards: ["03.Germline.Lofreq/{}/{}.lofreq.vcf.gz.anno/Merge.Anno.matrix.gz".format(sample, sample) for sample in sampledic],
        lambda wildcards: ["03.Germline.freebayes/{}/{}.vcf.gz.anno/Merge.Anno.matrix.gz".format(sample, sample) for sample in sampledic],
        lambda wildcards: ["03.Germline.chrM/{}/{}.vcf".format(sample, sample) for sample in sampledic],

rule QC:
    input:
        reads1=lambda wildcards: rundic[wildcards.run]['Read1'],
        reads2=lambda wildcards: rundic[wildcards.run]['Read2'],
    output:
        reads1out=temp("01.CleanData/{run}/{run}.R1.cln.fq.gz"),
        reads2out=temp("01.CleanData/{run}/{run}.R2.cln.fq.gz"),
        htmlout="01.CleanData/{run}/{run}.QC.html",
        jsonout="01.CleanData/{run}/{run}.QC.json"
    log: 
        out = "logs/A1.QC/{run}.o",
        err = "logs/A1.QC/{run}.e"
    threads: 8
    resources: 
        mem = '16g',
        extra = ' --gres=lscratch:10 '
    shell:
        "/data/yuk5/script/fastp "
        " -i {input.reads1}" 
        " -I {input.reads2}" 
        " -o {output.reads1out}" 
        " -O {output.reads2out}" 
        " -h {output.htmlout}" 
        " -j {output.jsonout}" 
        " -w {threads} > {log.out} 2> {log.err}"

rule bwa_mem:
    input:
        reads=["01.CleanData/{run}/{run}.R1.cln.fq.gz", 
               "01.CleanData/{run}/{run}.R2.cln.fq.gz"]
    output:
        #bam="02.Alignment/Level1/{run}/{run}.sort.bam",
        #bai="02.Alignment/Level1/{run}/{run}.sort.bam.bai"
        bam=temp("02.Alignment/Level1/{run}/{run}.sort.bam"),
        bai=temp("02.Alignment/Level1/{run}/{run}.sort.bam.bai"),
    log:
        out = "logs/B1.bwa/{run}.o",
        err = "logs/B1.bwa/{run}.e"
    params:
        bwa=lambda wildcards:' -K 100000000 -v 3 -R "@RG\\tID:{}\\tLB:{}\\tPL:illumina\\tPU:{}\\tSM:{}" '.format(
                                wildcards.run,  rundic[wildcards.run]['Lib'],  wildcards.run,  rundic[wildcards.run]['Sample']),
        samblaster="",
        sambambasort=" --tmpdir /lscratch/$SLURM_JOB_ID "
    threads: 16
    resources: 
        mem = '16g',
        extra = ' --gres=lscratch:20 '
    shell:
        "module load {config[modules][bwa]} {config[modules][samblaster]} {config[modules][sambamba]}\n"
        "bwa mem -t {threads} "
        " {params.bwa} "
        " {config[references][bwaidx]} "
        " {input.reads} 2> {log.err} "
        " |samblaster "
        " {params.samblaster} "
        " |sambamba view -S -f bam /dev/stdin "
        " -t {threads} "
        " |sambamba sort /dev/stdin "
        " -t {threads} "
        " -o {output.bam} "
        " {params.sambambasort} >> {log.out} 2>> {log.err}"

rule markdup:
    input:
        bam= lambda wildcards: ["02.Alignment/Level1/{}/{}.sort.bam".format(run, run) for run in libdic[wildcards.lib]],
        bai= lambda wildcards: ["02.Alignment/Level1/{}/{}.sort.bam.bai".format(run, run) for run in libdic[wildcards.lib]],
    output:
        bam    =temp("02.Alignment/Level2/{lib}/{lib}.sort.md.bam"),
        bai    =temp("02.Alignment/Level2/{lib}/{lib}.sort.md.bam.bai"),
        metrics=temp("02.Alignment/Level2/{lib}/{lib}.sort.md.metrics"),
    log:
        out = "logs/B2.mkdup/{lib}.o",
        err = "logs/B2.mkdup/{lib}.e",
    threads:  4
    resources:
        mem = '16g',
        extra = ' --gres=lscratch:40 ',
    run:
        #inputs = " ".join(["INPUT={}".format(f) for f in {input}]),
        inputs = " ".join("INPUT={}".format(in_) for in_ in input.bam),
        shell(
            """
            module load {config[modules][picard]} {config[modules][sambamba]}
            java -Xmx{resources.mem}  -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID  -jar $PICARD_JAR MarkDuplicates \
                {inputs} \
                OUTPUT={output.bam} \
                METRICS_FILE={output.metrics} \
                MAX_RECORDS_IN_RAM=2000000 \
                VALIDATION_STRINGENCY=SILENT \
                OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
                ASSUME_SORT_ORDER=coordinate \
                CLEAR_DT=false \
                ADD_PG_TAG_TO_READS=false >> {log.out} 2>> {log.err}
            sambamba index -t 4 {output.bam} >> {log.out} 2>> {log.err}
            """)

rule merge_level3:
    input:
        bam = lambda wildcards: ["02.Alignment/Level2/{}/{}.sort.md.bam".format(lib, lib) for lib in sampledic[wildcards.sample]],
        bai = lambda wildcards: ["02.Alignment/Level2/{}/{}.sort.md.bam.bai".format(lib, lib) for lib in sampledic[wildcards.sample]],
    output:
        bam=temp("02.Alignment/Level3/{sample}/{sample}.sort.md.bam"),
        bai=temp("02.Alignment/Level3/{sample}/{sample}.sort.md.bam.bai"),
    log:
        out = "logs/B3.merge_level3/{sample}.o",
        err = "logs/B3.merge_level3/{sample}.e",
    threads:  8
    resources:
        mem = '16g',
        extra = ' --gres=lscratch:10 ',
    run:
        if len(sampledic[wildcards.sample]) > 1:
            inputs = " ".join(input.bam),
            shell("""
            module load {config[modules][sambamba]}
            sambamba merge \
                -t {threads} \
                {output.bam} \
                {inputs} >> {log.out} 2>> {log.err}
            """)
        else:
            shell("""
            mv {snakemake.input[0]} {output.bam}
            sambamba index \
                -t {threads} \
                {output.bam} 2>> {log.err}
            """)
            
rule chrMbam:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.sort.md.bam",
        bai="02.Alignment/Level3/{sample}/{sample}.sort.md.bam.bai",
    output:
        mbam="02.Alignment/chrM/{sample}/{sample}.bam",
    log:
        out = "logs/B4.chrMbam/{sample}.o",
        err = "logs/B4.chrMbam/{sample}.e",
    threads:  8
    resources:
        mem = '16g',
        extra = ' --gres=lscratch:10 ',
    shell:
        """
        module load samtools
        samtools view -Sb -h {input.bam} chrM > {output.mbam} 2>> {log.err}
        samtools index {output.mbam}
        """
        
rule bqsr:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.sort.md.bam",
        bai="02.Alignment/Level3/{sample}/{sample}.sort.md.bam.bai",
    output:
        metrics="02.Alignment/Level3/{sample}/{sample}.BQSR.metrics",
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
        bai="02.Alignment/Level3/{sample}/{sample}.BQSR.bai",
    log:
        out = "logs/B5.BQSR/{sample}.o",
        err = "logs/B5.BQSR/{sample}.e",
    threads:  4
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:40 ',
    shell:
        """
        module load {config[modules][gatk]}
        gatk --java-options "-Xmx{resources.mem} -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" BaseRecalibrator \
            -R {config[references][fasta]} \
            -I {input.bam} \
            -O {output.metrics} \
            --use-original-qualities \
            --known-sites {config[references][gatk_dbsnp]} \
            --known-sites {config[references][gatk_1000g]} \
            --known-sites {config[references][gatk_indel]} \
            --intervals   {config[references][flankbed]} >> {log.out} 2>> {log.err}
        gatk --java-options "-Xmx12000m -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" ApplyBQSR \
            --add-output-sam-program-record \
            -R {config[references][fasta]} \
            -I {input.bam} \
            -O {output.bam} \
            --use-original-qualities \
            -bqsr {output.metrics} \
            -L {config[references][flankbed]} >> {log.out} 2>> {log.err}
        """            

        
rule stat_bqsr:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
    params:
        stat1="02.Alignment/Level3/{sample}/stat/{sample}.stat1",
        stat2="02.Alignment/Level3/{sample}/stat/{sample}.stat2",
        stat3="02.Alignment/Level3/{sample}/stat/{sample}.stat3",
    output:
        statdir=directory("02.Alignment/Level3/{sample}/stat"),
        flags="02.Alignment/Level3/{sample}/{sample}.BQSR.bam.flagstat",
    log:
        out = "logs/B6.BQSRstat/{sample}.o",
        err = "logs/B6.BQSRstat/{sample}.e",
    threads:  4
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:10 ',
    shell:
        """
        module load {config[modules][picard]} {config[modules][sambamba]} {config[modules][samtools]}
        mkdir -p {output.statdir}
        java -Xmx8000m -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $PICARD_JAR CollectMultipleMetrics \
            INPUT={input.bam} \
            REFERENCE_SEQUENCE={config[references][fasta]} \
            OUTPUT={params.stat1} \
            ASSUME_SORTED=true \
            PROGRAM="null" \
            PROGRAM="CollectAlignmentSummaryMetrics" \
            PROGRAM="CollectInsertSizeMetrics" \
            PROGRAM="CollectSequencingArtifactMetrics" \
            PROGRAM="CollectGcBiasMetrics" \
            PROGRAM="QualityScoreDistribution" \
            METRIC_ACCUMULATION_LEVEL="null" \
            METRIC_ACCUMULATION_LEVEL="SAMPLE" \
            METRIC_ACCUMULATION_LEVEL="LIBRARY" >> {log.out} 2>> {log.err}
        java -Xmx8000m -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $PICARD_JAR CollectMultipleMetrics \
            INPUT={input.bam} \
            REFERENCE_SEQUENCE={config[references][fasta]} \
            OUTPUT={params.stat2} \
            ASSUME_SORTED=true \
            PROGRAM="null" \
            PROGRAM="CollectAlignmentSummaryMetrics" \
            PROGRAM="CollectGcBiasMetrics" \
            METRIC_ACCUMULATION_LEVEL="null" \
            METRIC_ACCUMULATION_LEVEL="READ_GROUP"  >> {log.out} 2>> {log.err}
        sambamba flagstat -t {threads} {input.bam} > {output.flags}
        module load {config[modules][python37]} 
        ./scripts/ExomeBamCovStat.py hg38idt \
            {input.bam} {params.stat3} >> {log.out} 2>> {log.err}
        """        
        
rule lofreqindelbam:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
    output:
        bam=temp("02.Alignment/Lofreq/{sample}/{sample}.li.bam"),
        bai=temp("02.Alignment/Lofreq/{sample}/{sample}.li.bam.bai"),
    log:
        out = "logs/B7.lofreqindelbam/{sample}.o",
        err = "logs/B7.lofreqindelbam/{sample}.e",
    threads:  4
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:10 ',
    shell:
        """
        module load {config[modules][lofreq]} {config[modules][sambamba]}
        lofreq indelqual \
          -f {config[references][fasta]} \
          --dindel \
          -o {output.bam} \
          {input.bam} >> {log.out} 2>> {log.err}
        sambamba index -t {threads} {output.bam}  >> {log.out} 2>> {log.err}
        """        

rule CallableLoci:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
    output:
        bed="02.Alignment/Callable/{sample}/{sample}.bed",
        summary="02.Alignment/Callable/{sample}/{sample}.summary",
    log:
        out = "logs/B8.Callable/{sample}.o",
        err = "logs/B8.Callable/{sample}.e",
    threads:  4
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:10 ',
    shell:
        """module load GATK/3.8-1
        GATK CallableLoci \
          -R {config[references][fasta]} \
          -L {config[references][flankitv]} \
          -I {input.bam} \
          --maxDepth 1000 \
          --minBaseQuality 10 \
          --minMappingQuality 10 \
          --minDepth 10 \
          --minDepthForLowMAPQ 20 \
          --summary {output.summary} \
          -o {output.bed}  >> {log.out} 2>> {log.err}
          """

rule HaplotypeCaller:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
    output:
        gvcf="03.Germline/{sample}/{sample}.gvcf.gz",
        statdir=directory("03.Germline/{sample}/stat")
    params:
        stat="03.Germline/{sample}/stat/{sample}.germline"
    threads:  8
    log:
        out = "logs/C1.HCaller/{sample}.o",
        err = "logs/C1.HCaller/{sample}.e",
    resources:
        mem  = '32g',
        extra = ' --gres=lscratch:40 ',
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
            -L {config[references][flankbed]} >> {log.out} 2>> {log.err}
        gatk --java-options "-Xmx8000m -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" ValidateVariants \
            -V {output.gvcf} \
            -R {config[references][fasta]} \
            -gvcf \
            --validation-type-to-exclude ALLELES \
            --dbsnp {config[references][gatk_dbsnp]} \
            -L {config[references][flankbed]} >> {log.out} 2>> {log.err}
        mkdir -p {output.statdir}
        java -Xms2000m -Xmx8000m -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $PICARD_JAR \
            CollectVariantCallingMetrics \
            INPUT={output.gvcf} \
            OUTPUT={params.stat} \
            DBSNP={config[references][gatk_dbsnp]} \
            SEQUENCE_DICTIONARY={config[references][fadict]} \
            TARGET_INTERVALS={config[references][itvlist]} \
            GVCF_INPUT=true  >> {log.out} 2>> {log.err}
        """
        
        
rule GenomicsDBImport:
    input:
        gvcf= lambda wildcards: ["03.Germline/{}/{}.gvcf.gz".format(sample, sample) for sample in sampledic],
    params:
        bed = "/data/yuk5/pipeline/wgs_germline/ref/hg38_chr_intervals/itv_{itv}.bed",
    output:
        itvvcf=temp("03.Germline/VQSR/Merge.itv_{itv}.vcf.gz"),
        itvmf=temp("03.Germline/VQSR/Merge.itv_{itv}.mf.vcf.gz"),
        itvso="03.Germline/VQSR/Merge.itv_{itv}.siteonly.vcf.gz"
    log:
        out = "logs/C2.GDBimport/{itv}.o",
        err = "logs/C2.GDBimport/{itv}.e",
    threads:  4
    resources:
        mem  = '16g',
        extra = ' --gres=lscratch:40 ',
    run:
        inputs = " ".join("-V {}".format(f) for f in input.gvcf)
        shell(
        """
        module load {config[modules][gatk]}
        gatk --java-options "-Xmx4g -Xms4g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            GenomicsDBImport \
            {inputs} \
            -L {params.bed} \
            --genomicsdb-workspace-path /lscratch/$SLURM_JOB_ID/gdb.itv_{wildcards.itv} \
            --merge-input-intervals \
            --consolidate >> {log.out} 2>> {log.err}
        gatk --java-options "-Xmx5g -Xms5g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            GenotypeGVCFs \
            -R {config[references][fasta]} \
            -O {output.itvvcf} \
            -D {config[references][gatk_dbsnp]} \
            -G StandardAnnotation -G AS_StandardAnnotation \
            --allow-old-rms-mapping-quality-annotation-data \
            --merge-input-intervals \
            -V gendb:///lscratch/$SLURM_JOB_ID/gdb.itv_{wildcards.itv} \
            -L {params.bed}   >> {log.out} 2>> {log.err}
            
        gatk --java-options "-Xmx16g -Xms16g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            VariantFiltration \
            --filter-expression "ExcessHet>54.69" \
            --filter-name ExcessHet \
            -V {output.itvvcf} \
            -O {output.itvmf} >> {log.out} 2>> {log.err}
        gatk --java-options "-Xmx16g -Xms16g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            MakeSitesOnlyVcf \
             -I {output.itvmf} \
             -O {output.itvso} >> {log.out} 2>> {log.err}
        """
        )
            
            
            
            
rule genotyping:
    input:
        mfitvs=lambda wildcards: \
             ["03.Germline/VQSR/Merge.itv_{}.mf.vcf.gz".format('%.4d'%itv) for itv in range(1, config['references']['interval']+1)],
        soitvs=lambda wildcards: \
             ["03.Germline/VQSR/Merge.itv_{}.siteonly.vcf.gz".format('%.4d'%itv) for itv in range(1, config['references']['interval']+1)],
    output:
        sovcf  ="03.Germline/VQSR/Merge.sites_only.vcf.gz",
        mfvcf  ="03.Germline/VQSR/Merge.mf.vcf.gz",
        mir = "03.Germline/VQSR/Merge.indels.recal",
        mit = "03.Germline/VQSR/Merge.indels.tranches",
        msr = "03.Germline/VQSR/Merge.snps.recal",
        mst = "03.Germline/VQSR/Merge.snps.tranches",
        msmr = "03.Germline/VQSR/Merge.snps.model.report",
        mirv = "03.Germline/VQSR/tmp.indel.recalibrated.vcf",
        vqsr="03.Germline/Merge.flt.vqsr.vcf.gz",
    log:
        out = "logs/C3.Genotype/Genotype.o",
        err = "logs/C3.Genotype/Genotype.e",
    threads:  8
    resources:
        mem  = '64g',
        extra = ' --gres=lscratch:100 ',
    run:
        mfinputs = " ".join("--input {}".format(i) for i in input.mfitvs)
        soinputs = " ".join("--input {}".format(i) for i in input.soitvs)
        shell(
        """
        module load {config[modules][gatk]} {config[modules][samtools]}
        
        gatk --java-options "-Xmx16g -Xms16g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          GatherVcfsCloud \
          --ignore-safety-checks \
          --gather-type BLOCK \
          --output {output.sovcf} \
          {soinputs}  >> {log.out} 2>> {log.err}
        gatk --java-options "-Xmx16g -Xms16g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          GatherVcfsCloud \
          --ignore-safety-checks \
          --gather-type BLOCK \
          --output {output.mfvcf} \
          {mfinputs}  >> {log.out} 2>> {log.err}
        tabix -p vcf {output.sovcf}
        tabix -p vcf {output.mfvcf}
        gatk --java-options -Xms24g \
          VariantRecalibrator \
          -V {output.sovcf} \
          -O {output.mir} \
          --tranches-file {output.mit} \
          --trust-all-polymorphic \
          -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
          -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
          -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP\
          --use-allele-specific-annotations \
          -mode INDEL \
          --max-gaussians 4 \
          --resource:mills,known=false,training=true,truth=true,prior=12 {config[references][gatk_1000g]} \
          --resource:axiomPoly,known=false,training=true,truth=false,prior=10 {config[references][gatk_axiom]} \
          --resource:dbsnp,known=true,training=false,truth=false,prior=2 {config[references][gatk_dbsnp]} >> {log.out} 2>> {log.err}
        gatk --java-options -Xms50g \
          VariantRecalibrator \
          -V {output.sovcf} \
          -O {output.msr} \
          --tranches-file {output.mst} \
          --trust-all-polymorphic \
          -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 \
          -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
          -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
          --use-allele-specific-annotations \
          -mode SNP \
          --sample-every-Nth-variant 10 \
          --output-model {output.msmr} \
          --max-gaussians 6 \
          -resource:hapmap,known=false,training=true,truth=true,prior=15 {config[references][gatk_hapmap]} \
          -resource:omni,known=false,training=true,truth=true,prior=12 {config[references][gatk_omni]} \
          -resource:1000G,known=false,training=true,truth=false,prior=10 {config[references][gatk_1000hc]} \
          -resource:dbsnp,known=true,training=false,truth=false,prior=7 {config[references][gatk_dbsnp]} >> {log.out} 2>> {log.err}
        gatk --java-options -Xms80g \
          VariantRecalibrator \
          -V {output.sovcf} \
          -O {output.msr} \
          --tranches-file {output.mst} \
          --trust-all-polymorphic \
          -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 \
          -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
          -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
          --use-allele-specific-annotations \
          -mode SNP \
          --input-model {output.msmr} \
          --max-gaussians 6 \
          -resource:hapmap,known=false,training=true,truth=true,prior=15 {config[references][gatk_hapmap]} \
          -resource:omni,known=false,training=true,truth=true,prior=12 {config[references][gatk_omni]} \
          -resource:1000G,known=false,training=true,truth=false,prior=10 {config[references][gatk_1000hc]} \
          -resource:dbsnp,known=true,training=false,truth=false,prior=7 {config[references][gatk_dbsnp]} >> {log.out} 2>> {log.err}
        gatk --java-options -Xms5g \
          ApplyVQSR \
          -O {output.mirv} \
          -V {output.mfvcf} \
          --recal-file {output.mir} \
          --use-allele-specific-annotations \
          --tranches-file {output.mit} \
          --truth-sensitivity-filter-level 95.0 \
          --create-output-variant-index true \
          -mode INDEL >> {log.out} 2>> {log.err}
        gatk --java-options -Xms5g \
          ApplyVQSR \
          -O {output.vqsr} \
          -V {output.mirv} \
          --recal-file {output.msr} \
          --use-allele-specific-annotations \
          --tranches-file {output.mst} \
          --truth-sensitivity-filter-level 99.7 \
          --create-output-variant-index true \
          -mode SNP >> {log.out} 2>> {log.err}
        """
        )

rule lofreq:
    input:
        bam="02.Alignment/Lofreq/{sample}/{sample}.li.bam",
        bai="02.Alignment/Lofreq/{sample}/{sample}.li.bam.bai",
    params:
        vcf="03.Germline.Lofreq/{sample}/{sample}.lofreq.vcf",
    output:
        vgz="03.Germline.Lofreq/{sample}/{sample}.lofreq.vcf.gz",
        tbi="03.Germline.Lofreq/{sample}/{sample}.lofreq.vcf.gz.tbi",
    log:
        out = "logs/C05.lofreq/{sample}.o",
        err = "logs/C05.lofreq/{sample}.e",
    threads:  16
    resources:
        mem  = '32g',
        extra = ' --gres=lscratch:10 ',
    shell:
        """
        module load {config[modules][lofreq]}
        lofreq call-parallel \
          --pp-threads {threads} \
          -f {config[references][fasta]} \
          -o {params.vcf} \
          --call-indels \
          {input.bam}  >> {log.out} 2>> {log.err}
        module load {config[modules][samtools]}
        bgzip {params.vcf}
        tabix -p vcf {output.vgz}  >> {log.out} 2>> {log.err}
        """
        
rule mtoolbox:
    input:
        bam="02.Alignment/chrM/{sample}/{sample}.bam",
    output:
        config="03.Germline.chrM/{sample}/config",
        path=directory("03.Germline.chrM/{sample}"),
        vcf="03.Germline.chrM/{sample}/{sample}.vcf",
    log:
        out = "logs/C6.mtoolbox/{sample}.o",
        err = "logs/C6.mtoolbox/{sample}.e",
    threads:  4
    resources:
        mem  = '32g',
        extra = ' --gres=lscratch:10 ',
    shell:
        """
        module load {config[modules][python27]}
        mkdir -p {output.path}/OUT_{wildcards.sample}
        set +eu
        source /home/yuk5/yuk5/app/anaconda3/etc/profile.d/conda.sh >> {log.out} 2>> {log.err}
        conda activate mtoolbox >> {log.out} 2>> {log.err}
        source /home/yuk5/yuk5/app/anaconda3/envs/mtoolbox/MToolBox-1.2.1/setenv.sh >> {log.out} 2>> {log.err}
        set -eu
        ./scripts/write_mtoolbox_config.py {input.bam} {output.path} {wildcards.sample} > {output.config} 2>> {log.err}
        which python
        MToolBox.sh -i {output.config} >> {log.out} 2>> {log.err}
        """
        
rule freebayes:
    input:
        bam="02.Alignment/Level3/{sample}/{sample}.BQSR.bam",
    output:
        vgz="03.Germline.freebayes/{sample}/{sample}.vcf.gz",
    log:
        out = "logs/C7.freebayes/{sample}.o",
        err = "logs/C7.freebayes/{sample}.e",
    threads:  16
    resources:
        mem  = '64g',
        extra = ' --gres=lscratch:10 ',
    shell:
        """
        module load {config[modules][freebayes]} {config[modules][samtools]}
        freebayes-parallel \
          <(fasta_generate_regions.py {config[references][fasta]}.fai \
          100000) \
          {threads} \
          -f {config[references][fasta]} \
          {input.bam} 2>> {log.err} |\
        bgzip > {output.vgz}
        tabix -p vcf {output.vgz}
        """
        
            
rule anno_gatk:
    input:
        "03.Germline/Merge.flt.vqsr.vcf.gz",
    output:
        folder=directory("03.Germline/Merge.flt.vqsr.vcf.anno"),
        result="03.Germline/Merge.flt.vqsr.vcf.anno/Merge.Anno.matrix.gz"
    log:
        out = "logs/D1.anno_gatk/Anno_GATK.o",
        err = "logs/D1.anno_gatk/Anno_GATK.e",
    threads:  24
    resources:
        mem  = '48g',
        extra = ' --gres=lscratch:10 ',
    shell:
        """
        {config[bins][vcfanno]} \
            {input} \
            {output.folder} \
            {threads} n 
        """

rule anno_lofreq:
    input:
        "03.Germline.Lofreq/{sample}/{sample}.lofreq.vcf.gz",
    output:
        folder=directory("03.Germline.Lofreq/{sample}/{sample}.lofreq.vcf.gz.anno"),
        result="03.Germline.Lofreq/{sample}/{sample}.lofreq.vcf.gz.anno/Merge.Anno.matrix.gz"
    log:
        out = "logs/D2.anno_lofreq/{sample}.o",
        err = "logs/D2.anno_lofreq/{sample}.e",
    threads:  16
    resources:
        mem  = '64g',
        extra = ' --gres=lscratch:10 ',
    shell:
        """
        {config[bins][vcfanno]} \
            {input} \
            {output.folder} \
            {threads} n 
        """
        
rule anno_freebayes:
    input:
        "03.Germline.freebayes/{sample}/{sample}.vcf.gz",
    output:
        folder=directory("03.Germline.freebayes/{sample}/{sample}.vcf.gz.anno"),
        result="03.Germline.freebayes/{sample}/{sample}.vcf.gz.anno/Merge.Anno.matrix.gz",
    log:
        out = "logs/D3.anno_freebayes/{sample}.o",
        err = "logs/D3.anno_freebayes/{sample}.e",
    threads:  16
    resources:
        mem  = '64g',
        extra = ' --gres=lscratch:10 ',
    shell:
        """
        {config[bins][vcfanno]} \
            {input} \
            {output.folder} \
            {threads} n 
        """
        
        