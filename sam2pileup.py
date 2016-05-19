#!/usr/bin/env python
#filename: sam_pileup.py

sam_f=open("./output/reads_1.sam","r")
sam_sort=open('./output/reads_1.sorted.sam','w')
pileup=open('./output/reads_1.pileup','w')
ref=open('./sample/lambda_virus.fa','r')
lnum=0
dic={}
for line in sam_f:
    #print(line)
    lnum+=1
    if lnum>3:
        line=line.strip('\n')
        line_ls=line.split('\t')
        if int(line_ls[3])!=0:
            dic[int(line_ls[3])]=line
#print(lnum)
dic_sort=sorted(dic.iteritems(),key=lambda asd:asd[0],reverse=False)
for e in dic_sort:
    sam_sort.write(str(e[1])+'\n')
sam_sort.close()
sam_sort=open('./output/reads_1.sorted.sam','r')
refs=''
for line in ref:
    if str(line[0])!=">":
        line=line.strip('\n')
        #print(line)
        refs=refs+line
print("length of reference genome: "+str(len(refs)))
print("sorted sam file outputs to folder: output")
print("Pileup file outputs to folder: output")
chrname='gi|9626243|ref|NC_001416.1|'
dic_base={}
dic_qual={}
marks=''
for line in sam_sort:
    line=line.strip('\n')
    line_ls=line.split('\t')
    flag=line_ls[1]
    #print(flag)
    pos=int(line_ls[3])
    #print(pos)
    #print(line_ls[4])
    ci=str(line_ls[5])
    if ci=='*':
        ci='0*'
    else:
        ci=str(len(line))+'*'
    #print(ci)
    ciga=int(ci[0:-1])
    base=line_ls[9]
    bq=line_ls[10]
    if flag=='0':
        mark='.'
    elif flag=='16':
        mark=','
    else:
        mark=''
    loc=0
    marks=''
    for b in base:
        coverage=1
        #if pos>48502:
        #    pos=48502
        if b==refs[pos-1]:
            if dic_base.has_key(pos):
                a=str(dic_base[pos])
                a=a.split('\t')
                coverage=int(a[2])+1
                marks=str(a[3])
                baseq=str(a[4])+bq[loc]
                marks=marks+mark
                dic_base[pos]=str(pos)+'\t'+b+'\t'+str(coverage)+'\t'+marks+'\t'+str(baseq)
            else:
                #print('pos '+str(pos))
                dic_base[pos]=str(pos)+'\t'+b+'\t'+str(coverage)+'\t'+mark+'\t'+str(bq[loc])
        if b!=refs[pos-1]:
        #else:
            if dic_base.has_key(pos):
                a=str(dic_base[pos])
                a=a.split('\t')
                coverage=int(a[2])+1
                marks=str(a[3])
                baseq=str(a[4])+bq[loc]
                if flag=='16':
                    b=b.lower()
                marks=marks+b
                dic_base[pos]=str(pos)+'\t'+refs[pos-1]+'\t'+str(coverage)+'\t'+marks+'\t'+str(baseq)
            else:
                #print('pos: '+str(pos))
                dic_base[pos]=str(pos)+'\t'+refs[pos-1]+'\t'+str(coverage)+'\t'+b+'\t'+str(bq[loc])
        loc+=1
        pos+=1

dic_base_sort=sorted(dic_base.iteritems(),key=lambda asd:asd[0],reverse=False)
for e in dic_base_sort:
    pile=str(e[1])
    #print(pile)
    pileup.write(chrname+'\t'+pile+'\n')

    
    

    




