#!/usr/bin/env python3
# =============================================================================
# Filename: vcf_annotation.v1.wxs.py
# Version: 
# Author: Kai Yu - finno@live.cn
# https://github.com/readline
# Last modified: 2018-10-29 09:46
# Description: 
# 
# =============================================================================
import os,sys,gzip
import pandas
import math
import multiprocessing

def pen(task_name, line):
    task_path = 'Tasks/%s.sh'%task_name
    with open(task_path, 'w') as savefile:
        savefile.write('#!/bin/bash\necho Job init: `date` at `hostname` &&\ncd %s &&\n'%projdir)
        tmp = line.rstrip('\n').split('\n')
        for line in tmp:
            if line[-2:] == ' \\':
                savefile.write(line + '\n')
            else:
                savefile.write(line + ' &&\n')
        savefile.write('echo Job end: `date` at `hostname`')
        

def split_vcf(vcfpath):
    vcfwc = 0
    if split == 'y':
        if vcfpath[-3:] == 'vcf':
            vcfwc = int(os.popen('grep -v "^#" %s|wc -l|cut -d " " -f1'%vcfpath).readline())
            base = 'cat %s|grep -v "^#"|cut -f1-7|'%vcfpath
        elif vcfpath[-3:] == '.gz':
            vcfwc = int(os.popen('zcat %s|grep -v "^#"|wc -l|cut -d " " -f1'%vcfpath).readline())
            base = 'zcat %s|grep -v "^#"|cut -f1-7|'%vcfpath
    elif split == 'n':
        vcfwc = int(os.popen('/data/yuk5/app/anaconda3/bin/vcf_parser %s --split|grep -v "^#"|wc -l|cut -d " " -f1'%vcfpath).readline())
        base = '/data/yuk5/app/anaconda3/bin/vcf_parser %s --split|grep -v "^#"|cut -f1-7|'%vcfpath
    
    splitsize = int(math.ceil(vcfwc*1.0/20))
    if vcfpath[-3:] == 'vcf':
        infile = open(vcfpath)
    elif vcfpath[-3:] == '.gz':
        infile = gzip.open(vcfpath)
    cmd = base + 'split -l %d - tmp/itv_'%splitsize
    tmp = os.popen(cmd)
    tmp.readline()
    os.system('pwd')
    rawnames = [i.rstrip() for i in os.popen('ls tmp/itv_??').readlines()]
    itvlist = []
    for n in range(len(rawnames)):
        os.rename(rawnames[n],'tmp/itv_%d.base'%(n+1))
        itvlist.append('tmp/itv_%d'%(n+1))
        
    print('VCF splited')
    return itvlist

       

def base_anno1(itv):
    print('%s, SnpEff'%itv)
    cmd = 'java -Xmx2g -jar %s/snpEff.jar hg38 %s.base > %s.ann'%(bindir, itv, itv)
    tmp = os.popen(cmd)
    return tmp.readline()

def base_anno2(itv):
    print('%s, SnpSift dbSNP'%itv)
    cmd = 'java -Xmx2g -jar %s/SnpSift.jar annotate /fdb/dbSNP/organisms/human_9606_b150_GRCh38p7/00-All.vcf.gz %s.base > %s.dbsnp'%(bindir, itv, itv)
    tmp = os.popen(cmd)
    return tmp.readline()

def base_anno3(itv):
    print('itv_%s, SnpSift 1000g'%itv)
    cmd = 'java -Xmx2g -jar %s/SnpSift.jar annotate /data/RUNX1/data/1000genome_GRCh38.genotype/ALLchr.GRCh38_sites.20170504.vcf.gz %s.base > %s.tgp'%(bindir, itv, itv)
    tmp = os.popen(cmd)
    return tmp.readline()

def base_anno4(itv):
    print('%s, SnpSift ExAC'%itv)
    cmd = 'java -Xmx2g -jar %s/SnpSift.jar annotate /data/RUNX1/data/ExAC/ExAC_nonTCGA.r1.hg38.sites.vep.vcf.gz %s.base > %s.exac'%(bindir, itv, itv)
    tmp = os.popen(cmd)
    return tmp.readline()

def base_anno5(itv):
    print('%s, SnpSift Clinvar'%itv)
    cmd = 'java -Xmx2g -jar %s/SnpSift.jar annotate /fdb/clinvar/vcf_GRCh38/clinvar_20180805.vcf.gz %s.base > %s.clinvar'%(bindir, itv, itv)
    tmp = os.popen(cmd)
    return tmp.readline()

def base_anno6(itv):
    print('%s, SnpSift dbNSFP'%itv)
    cmd = 'java -Xmx2g -jar %s/SnpSift.jar dbnsfp -v -db /fdb/dbNSFP2/3.5a/dbNSFP3.5a.txt.gz %s.base > %s.dbnsfp'%(bindir, itv, itv)
    tmp = os.popen(cmd)
    return tmp.readline()

def base_anno7(itv):
    print('%s, SnpSift COSMIC'%itv)
    cmd = 'java -Xmx2g -jar %s/SnpSift.jar annotate /data/RUNX1/data/COSMIC/hg38/CosmicCodingMuts.vcf.gz %s.base > %s.cosmic'%(bindir, itv, itv)
    tmp = os.popen(cmd)
    return tmp.readline()



def main():
    global vcfpath,bindir,split
    try:
        vcfpath = sys.argv[1]
        outdir  = sys.argv[2]
        threads = sys.argv[3]
        split   = sys.argv[4]
    except:
        sys.exit(sys.argv[0] + ' [vcf path] [output dir] [threads|takes at least 2*t gb memory] [split: y|n]')

    if vcfpath[0] != '/':
        vcfpath = os.getcwd()+'/'+vcfpath
    bindir = '/usr/local/apps/snpEff/5.1d'
    intervals = ['/data/RUNX1/pipeline/wxs_phase1/ref/IDT_xGen_Exome_Research_Panel/hg38/intervals_withflank150/hg38.xGen.%d.bed'%n for n in range(1,21)]
    threads = int(threads)

    os.system('mkdir -p %s/tmp'%outdir)
    os.chdir(outdir)
    
    if threads //2 >0:
        halfthread = threads //2
    else:
        halfthread = 1
    itvlist = split_vcf(vcfpath)
    pool = multiprocessing.Pool(processes=halfthread)
    results = []
    for itv in itvlist:
        results.append(pool.apply_async(base_anno1,(itv,)))
        results.append(pool.apply_async(base_anno2,(itv,)))
        results.append(pool.apply_async(base_anno3,(itv,)))
        results.append(pool.apply_async(base_anno4,(itv,)))
        results.append(pool.apply_async(base_anno5,(itv,)))
        results.append(pool.apply_async(base_anno6,(itv,)))
        results.append(pool.apply_async(base_anno7,(itv,)))
    pool.close()
    pool.join()

    filetypes = ['ann','clinvar','dbsnp','exac','tgp','dbnsfp','cosmic','base']
    for i in filetypes:
        os.system('cat %s.%s > Merge.%s.vcf'%(itvlist[0], i, i))
        for itv in itvlist[1:]:
            os.system('cat %s.%s|grep -v "^#" >> Merge.%s.vcf'%(itv, i, i))
    if vcfpath[-3:] == '.gz':
        cmd = '''zcat %s|head -n 10000|grep '^#'|tail -1|sed 's/^#//g' > tmp.1 &&\n'''%(vcfpath)
    else:
        cmd = '''cat %s|head -n 10000|grep '^#'|tail -1|sed 's/^#//g' > tmp.1 &&\n'''%(vcfpath)
    if split == 'y':
        if vcfpath[-3:] == '.gz':
            cmd += 'zcat %s|grep -v "^#" >> tmp.1 &&\n'%(vcfpath)
        else:
            cmd += 'cat %s|grep -v "^#" >> tmp.1 &&\n'%(vcfpath)
    elif split == 'n':
        cmd += '/data/yuk5/app/anaconda3/bin/vcf_parser %s --split|grep -v "^#" >> tmp.1 &&\n'%(vcfpath)
        
    cmd += '''echo SNPEFF > tmp.7 &&
cat Merge.ann.vcf|grep -v '^#'  |cut -f8 >> tmp.7 &&
echo dbSNP > tmp.2 &&
cat Merge.dbsnp.vcf|grep -v '^#'  |cut -f3,8|sed 's/\t/^/g'|sed 's/\.^\././g' >> tmp.2 &&
echo TGP > tmp.3 &&
cat Merge.tgp.vcf|grep -v '^#'    |cut -f8 >> tmp.3 &&
echo ExAC > tmp.4 &&
cat Merge.exac.vcf|grep -v '^#'   |cut -f8 >> tmp.4 &&
echo dbNSFP > tmp.5 &&
cat Merge.dbnsfp.vcf|grep -v '^#' |cut -f8 >> tmp.5 &&
echo Clinvar > tmp.6 &&
cat Merge.clinvar.vcf|grep -v '^#'|cut -f8 >> tmp.6 &&
echo COSMIC > tmp.8 &&
cat Merge.cosmic.vcf|grep -v '^#'|cut -f3,8|sed 's/\t/^/g'|sed 's/\.^\././g' >> tmp.8 &&
paste tmp.1 tmp.7 tmp.2 tmp.3 tmp.4 tmp.5 tmp.6 tmp.8 > Merge.Anno.matrix &&
gzip Merge.Anno.matrix &&
rm tmp.1 tmp.7 tmp.2 tmp.3 tmp.4 tmp.5 tmp.6 tmp.8 Merge.ann.vcf Merge.dbsnp.vcf Merge.tgp.vcf Merge.exac.vcf Merge.dbnsfp.vcf Merge.clinvar.vcf Merge.cosmic.vcf Merge.base.vcf &&
rm -rf tmp'''
    os.system(cmd)
if __name__ == '__main__':
    main()

