#!/usr/bin/env python3
import os,sys
import pandas as pd
import numpy as np
from io import StringIO
import multiprocessing
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def font(size):
    return matplotlib.font_manager.FontProperties(size=size,fname='/data/yuk5/app/font/arial.ttf')

def calc_depth(bedpath, bampath):
    cmd = '%s depth -b %s -q 20 %s'%(samtools, bedpath, bampath)
    df = pd.read_csv(StringIO(os.popen(cmd).read()), sep="\t", header=None)
    df.columns = ['chr','pos','depth']
    return df

def get_sample(bampath):
    cmd = '%s view -H %s'%(samtools, bampath)
    heads = os.popen(cmd).readlines()
    for line in heads:
        c = line.strip().split('\t')
        if c[0] != '@RG':
            continue
        for i in c:
            if i[:3] == 'SM:':
                return i[3:]
            else:
                continue
    return 'Sample'

def main():
    bedfiles = {}
    bedfiles['hg38idt'] = '/data/yuk5/pipeline/wxs_phase1/ref/IDT_xGen_Exome_Research_Panel/hg38/xGen_Exome_Research_Panel.targets.hg38.10pct.bed'
    bedfiles['mm10SeqCap'] = '/data/yuk5/pipeline/wxs_phase1_mm10/ref/bedfile/SeqCapEZ.target.mm10.10pct.bed'
    bedfiles['hg38TrueSeq'] = '/data/yuk5/pipeline/wxs_phase1_truseqexome/TruSeq_Exome_Kit/truseq-dna-exome-targeted-regions-manifest-v1-2.hg38.10pct.bed'
    global samtools
    samtools= '/usr/local/apps/samtools/1.6/bin/samtools'
    try:
        ref     = sys.argv[1]
        bampath = sys.argv[2]
        prefix  = sys.argv[3]
    except:
        sys.exit(sys.argv[0] + ' [Ref: hg38idt|mm10SeqCap|hg38TrueSeq] [bam file] [prefix]')
    if ref not in 'hg38idt|mm10SeqCap|hg38TrueSeq'.split('|'):
        sys.exit(sys.argv[0] + ' [Ref: hg38idt|mm10SeqCap|hg38TrueSeq] [bam file] [prefix]')
    sample = get_sample(bampath)
    ############## Calc depth ##############
    df = calc_depth(bedfiles[ref], bampath)
    chrlist = ['chr1','chr2','chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
        'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 
        'chr18', 'chr19','chr20', 'chr21', 'chr22','chrX', 'chrY']
    for c in set(df.chr):
        if c not in chrlist:
            df = df.drop(df.loc[df.chr==c].index)
    
    ############## Load Bed ##############
    beddf = pd.read_csv(bedfiles[ref], sep='\t', header=None)
    beddf.columns = ['chr','start','end']
    beddf['range'] = beddf.end-beddf.start
    beddf = pd.Series(data=[beddf.loc[beddf.chr==i,'range'].sum() for i in chrlist], index=chrlist)
    
    ############## Calc Chr ##############
    chrdpdf = pd.Series(data=[df.loc[df.chr==i,'depth'].sum()/beddf[i] for i in chrlist], index=chrlist)
    
    ############## Calc Autosome chr ##############
    chr22df = df.loc[(df.chr!='chrX')&(df.chr!='chrY')]
    regsize = beddf[chrlist[:22]].sum()
    chr22d = pd.Series(data=list(chr22df.depth) + [0]*(regsize-chr22df.shape[0])).sort_values()
    fracdf = pd.Series(data=[chr22d.loc[chr22d==i].shape[0]/chr22d.sum() for i in range(chr22d.max())])
    covdf = pd.Series(data=[chr22d[chr22d>i].shape[0]/regsize for i in range(chr22d.max())])
    tmpqt = list(chr22d.quantile([.25, .5, .75]))
    tmpdq = tmpqt[2]-tmpqt[0]
    tmpuo = tmpqt[2]+1.5*tmpdq
    tmpdo = tmpqt[0]-1.5*tmpdq
    meandepth = chr22d.sum()/regsize
    adjmeandepth = np.mean([n for n in df.depth if n > tmpdo and n < tmpuo])
    tmplinex,tmpliney = [],[]
    tmp50,tmp100,tmp150,tmp200 = None,None,None,None
    for x,y in zip(covdf.index,covdf*100):
        if not tmp50 and x >= 50:
            tmp50 = True
            tmplinex.append(x)
            tmpliney.append(y)
        if not tmp100 and x >= 100:
            tmp100 = True
            tmplinex.append(x)
            tmpliney.append(y)
        if not tmp150 and x >= 150:
            tmp150 = True
            tmplinex.append(x)
            tmpliney.append(y)
        if not tmp200 and x >= 200:
            tmp200 = True
            tmplinex.append(x)
            tmpliney.append(y)
            
    ############## Plot ##############
    fig = plt.figure(figsize=(9,6), dpi=100)
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(212)

    ax1.plot(fracdf.index,fracdf*100,color='grey')
    ax1.axvline(x=adjmeandepth,color='red',linestyle=':',alpha=0.6)
    ax1.text(adjmeandepth,0,' Adj.Depth. %.1fX'%adjmeandepth, va='bottom',ha='left', font_properties=font(13))
    ax1.set_title('Base Depth Distribution', font_properties=font(13))
    ax1.set_xlim(0,int(adjmeandepth*4))
    ax1.set_ylabel('Fraction of Bases (%)', font_properties=font(13))
    [i.set_fontproperties(font(11)) for i in ax1.get_xticklabels()]
    [i.set_fontproperties(font(11)) for i in ax1.get_yticklabels()]


    ax2.plot(covdf.index,covdf*100,color='grey')
    #ax2.set_ylim(0,100)
    for x,y in zip(tmplinex,tmpliney):
        ax2.axvline(x=x,ymax=y*1.0/max(covdf*100),color='red',linestyle=':',alpha=0.6)
        ax2.axhline(y=y,xmax=x*1.0/(adjmeandepth*4),color='red',linestyle=':',alpha=0.6)
        ax2.text(x,y,'      %.1f%% (%dX)'%(y,int(x)),va='top',ha='left', font_properties=font(11))
    ax2.set_title('Coverage', font_properties=font(13))
    ax2.set_xlabel('Depth', font_properties=font(13))
    ax2.set_ylabel('Covered (%)', font_properties=font(13))
    ax2.set_xlim(0,int(adjmeandepth*4))
    ax2.set_ylim(0,100)
    [i.set_fontproperties(font(11)) for i in ax2.get_xticklabels()]
    [i.set_fontproperties(font(11)) for i in ax2.get_yticklabels()]

    if chrdpdf['chrY']< 0.2*chrdpdf[['chr%s'%str(i) for i in list(range(1,23))]].mean():
        color='#ff7c7c'
    else:
        color='#88e1f2'
    ax3.bar(range(chrdpdf.shape[0]), chrdpdf, color=color)
    ax3.set_xticks(range(chrdpdf.shape[0]))
    ax3.set_xticklabels([i.replace('chr','') for i in chrdpdf.index], font_properties=font(13))
    ax3.set_ylabel(sample, font_properties=font(14))
    ax3.axhline(y=chrdpdf[['chr%s'%str(i) for i in list(range(1,23))]].mean(),
                    linestyle=':', color='black')
    [i.set_fontproperties(font(11)) for i in ax3.get_yticklabels()]
    ax3.set_xlim(-1,chrdpdf.shape[0])

    fig.tight_layout()
    plt.savefig(prefix+'.xcov.pdf')
    
    ############## Log ##############
    with open(prefix+'.xcov.log','w') as savefile:
        savefile.write('#Sample\t%s\n'%sample)
        savefile.write('A1. Target region effective bases\t%d\n'%(df.depth.sum()*10))
        savefile.write('A2. Target region average depth\t%.2f\n'%(meandepth))
        savefile.write('A3. Target region QTadj depth\t%.2f\n'%(adjmeandepth))
        for n in range(12):
            cutoff = [10,20,30,40,50,60,70,80,90,100,150,200][n]
            savefile.write('B%d. Target region coverage >=%dx\t%.2f%%\n'%(n+1,cutoff, chr22df.loc[chr22df.depth>=cutoff].shape[0]/regsize*100))
        for c in chrdpdf.index:
            savefile.write('C.%s\t%.2f\n'%(c, chrdpdf[c]))
            
if __name__ == '__main__':
    main()
