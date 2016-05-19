#!/usr/bin/env python
#file name: get_snp.py
from __future__ import division
import time
#f_in=open("dqSRR1197554.pileup","r")
date=time.strftime("%Y%m%d",time.localtime())
f_in=open("./output/reads_1.pileup","r")
#f_in=f_in.readline()
depth=3
f_out=open("./output/reads_1.vcf","w")
print("vcf file output to folder: output")
f_out.write('##fileformat=VCFv4.3\n')
f_out.write('##source=LambdaPhageV01\n')
f_out.write('##reference=./sample/lambda_virus.fa\n')
f_out.write('##fileDate='+str(date)+'\n')
f_out.write('##INFO=<ID=DEPTH,Number=1,Type=Integer,Description="sequencing at this postition">\n'
        '##INFO=<ID=SNP,Number=0,Type=Flag,Description="classified as SNP">\n'
        '##INFO=<ID=HOMO,Number=0,Type=Flag,Description="classified as HOMO">\n'
        '##INFO=<ID=HETE,Number=0,Type=Flag,Description="classified as HETE">\n'
        '##INFO=<ID=VF,Number=1,Type=Float,Description="varant frequency">\n')

f_out.write('#CHROM\t'+'POS\t'+'ID\t'+'REF\t'+'ALT\t'+'QUAL\t'+'FILTER'+'\t'+'INFO\n')
j=0
for line in f_in:
    line=line.rstrip('\n')
    ls=line.split('\t')
    #print(int(ls[3]))
    if int(ls[3])>depth:
        j+=1
        line_str=str(ls[4])
        ref_count=0
        A=0
        T=0
        G=0
        C=0
        N=0
        r_count=0
        alt=''
        info=''
        info="DEPTH="+str(len(line_str))+';'
        for char in line_str:
            char=char.upper()
            ref=str(ls[2]).upper()
            if char=="," or char==".":
                r_count+=1
                #ref_count+=1
            else:
                alt=alt+char
                if char=="A":
                    A+=1
                if char=="T":
                    T+=1
                if char=="G":
                    G+=1
                if char=="C":
                    C+=1
                if char=="N":
                    N+=1
        total=A+T+C+G+N+r_count
        var=A+T+C+G
        if total==0:
            print(line)
            continue
        frac=var/total
        if frac>0.8:
            #if frac>0.8:
            info=info+'SNP;'
            #elif 0.2<frac <0.8:
            #    info=info+'HETE;'
            #else:
            #    info=info+'HETE;'
            info=info+"FRACTION="+str(frac)
    
            
            f_out.write(str(ls[0])+'\t'+str(ls[1])+'\t'+'.\t'+str(ls[2])+'\t'+alt+'\t.\tPASS\t'+ info+'\n')
        
#print("j is "+str(j))




