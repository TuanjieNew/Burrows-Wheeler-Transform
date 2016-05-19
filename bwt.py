#!/usr/bin/env python
#fn:bwt.py
#python2
import re

file_i=open("./sample/lambda_virus.fa",'r')#open reference genome
#file_i=open("lam_2.fa",'r')
file_i.readline()
list1=[]

for line in file_i:
    if line[0]=='>':
        pattern=re.compile(">*[ ]")
        rname=pattern.search(str(line)).groups()
        #print(rname)
        break

    line=line.strip('\n')
    for char in line:
         list1.append(char)
list1.append('$')
length=len(list1)
stri="".join(list1)
dic={}
l=20
s=20 #fixed length
j=1
i=0
dic[stri[:s]]=['$']
while i < length-1:
    if(l>1):
        stri3=stri[length-j:length]+stri[0:l-1]
    else:
        stri3=stri[length-j:length-j+s+1]

    dic[stri3]=stri[length-j-1]
    l-=1
    j+=1
    i+=1
ls_dic=sorted(dic.iteritems(),key=lambda asd:asd[0],reverse=False)
result=[]
for e in ls_dic:
    result.append(e[1])

#print 'The first five: ',result[:10] 
#print 'The last five: ',result[length-5:] 
file_i.close()
###############################################################suffix index
si=open("./output/suffix_index.txt","w")
print("Suffix index file outputs to foleder: output")
dic2={}
dic3={}
l=20
s=20
j=1
i=0
k=48502
while i < length:
    if(l>1):
        stri3=stri[length-j:length]+stri[0:l-1]
    else:
        stri3=stri[length-j:length-j+s+1]
    dic2[k]=stri3
    dic3[k]=stri[length-j-1]
    l-=1
    j+=1
    i+=1
    k-=1

ls_dic2=sorted(dic2.iteritems(),key=lambda asd:asd[1],reverse=False)
#print(ls_dic2)
#print(dic2[1])
sf=[]
for e in ls_dic2:
    sf.append(e[0])
    si.write(str(e[0])+"\n")
#print('sflength: '+str(len(sf)))

#print(sf[48500:])
#print(j)
#print(len(result))
################################################################Full Tally
#suffix_index=open("suffix_cyr.txt","r")
suffix_index=open("./output/suffix_index.txt","r")
#suffix_index.readline()
suffix_ls=[]
suf_dic={}
a=0
for line in suffix_index:
    a+=1
    line=line.strip("\n")
    '''
    if str(line)=='0':
        print(a)
    if a==0:
        print("line: "+str(line)+" "+str(a))
    '''
    suffix_ls.append(a)
    suf_dic[a]=line
    #a+=1

#print("a: "+str(a))
#print(suffix_ls[0])


aa=[]
cc=[]
gg=[]
tt=[]
ai=0
ci=0
gi=0
ti=0
for n in result:
    if n=="A":
        ai+=1
        aa.append(ai)
    else:
        aa.append(ai)
    if n=="C":
        ci+=1
        cc.append(ci)
    else:
        cc.append(ci)
    if n=="G":
        gi+=1
        gg.append(gi)
    else:
        gg.append(gi)
    if n=="T":
        ti+=1
        tt.append(ti)
    else:
        tt.append(ti)
#print(tt[48500:])
out_list=open("./output/out_list.txt","w")
print("Out list file outputs to folder: output")
o=0
while o< 48503:
    out_list.write(str(sf[o])+'\t'+str(aa[o])+'\t'+str(cc[o])+'\t'+str(gg[o])+'\t'+str(tt[o])+'\n')

    o+=1
##########################################################backtrack
def getpos(query):
        sa=0
        a=12334 + sa #12334
        c=11362+a  #23696
        g=12820+c  #36516
        t=11986+g  #48502
        #query="CCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTT"
        #query="TTTTCGCTATTTATGAAAATTTTCCG" #19
        #query="ATCCGGCGCGTGAGTTCACCATGA"#7071
        #q_ls=list(query)
        #q_ls.reverse()
        #query="".join(q_ls)
        #query="GCTGTA"
        #print(query)
        que_len=len(query)-1
        if query[que_len]=="A":
            for_start=sa
            for_end=a
        if query[que_len]=="C":
            for_start=a
            for_end=c
        if query[que_len]=="G":
            for_start=c
            for_end=g
        if query[que_len]=="T":
            for_start=g
            for_end=t
        if query[que_len]=="N":
            for_start=sa
            for_end=a
        r=0
        que_len=que_len+1
        #print(que_len)
        stop=0
        while que_len>1:
            r+=1
            #if for_end>48502:
            #    for_end=48502
            #if for_start>48502:
            #    for_start=48502
            #print(str(r)+" "+str(for_start))
            #print(str(r)+" "+str(for_end))
            que_len -=1
            q=query[que_len]
            #print(q)
            if q=="A":
                start=for_start
                end=for_end
            if q=="C":
                start=for_start
                end  =for_end
            if q=="G":
                start=for_start
                end  =for_end
            if q=="T":
                start=for_start
                end  =for_end
            if q=="N":
                start=for_start
                end  =for_end

            q2=query[que_len-1]
            #stop=0
            if q2=="N":
                q2=="A"
            if q2=="A":
                num=aa[end]-aa[start]
                #print('A'+' '+str(num))
                if num>0:
                    for_start=sa+aa[start]
                    for_end  =sa+aa[start]+num
                else:
                    q2="C"
                    stop+=1
            if q2=="C":
                num=cc[end]-cc[start]
                #print('C'+' '+str(num))
                if num>0:
                    for_start=a+cc[start]
                    for_end  =a+cc[start]+num
                else:
                    q2="G"
                    stop+=1
            if q2=="G":
                num=gg[end]-gg[start]
                #print('G'+' '+str(num))
                if num>0:
                    for_start=c+gg[start]
                    for_end  =c+gg[start]+num
                else:
                    q2="T"
                    stop+=1
            if q2=="T":
                num=tt[end]-tt[start]
                #print('T'+' '+str(num))
                if num>0:
                    for_start=g+tt[start]
                    for_end  = g+tt[start]+num
                else:
                    q2="A"
                    stop+=2
            
            if stop>1 and stop < 5:
                if q2=="A":
                    num=aa[end]-aa[start]
                    #print('A'+' '+str(num))
                    if num>0:
                        for_start=sa+aa[start]
                        for_end  =sa+aa[start]+num
                    else:
                        q2="C"
                        stop+=1
                if q2=="C":
                    num=cc[end]-cc[start]
                    #print('C'+' '+str(num))
                    if num>0:
                        for_start=a+cc[start]
                        for_end  =a+cc[start]+num
                    else:
                        q2="G"
                        stop+=1
                if q2=="G":
                    num=gg[end]-gg[start]
                    #print('G'+' '+str(num))
                    if num>0:
                        for_start=c+gg[start]
                        for_end  =c+gg[start]+num
                    else:
                        q2="T"
                        stop+=1
            #print('stop: '+str(stop))
            if stop>5:
                pos=1
                return pos
            #query[0]
            #print(r)
            #print(for_start)
            #print(for_end)
            #print(len(suffix_ls))
            #print(sf[for_end]+1)
            pos=sf[for_end]+1
        return pos
def getposn(query):
    que_len=len(query)-1
    
def revl(seq):
    seq=seq[::-1]
    pseq=''
    for char in seq:
        if char=="A":
            pseq=pseq+'T'
        if char=="T":
            pseq=pseq+'A'
        if char=='C':
            pseq=pseq+'G'
        if char=='G':
            pseq=pseq+'C'
    return pseq
#print(revl(seq))
#print(getpos(revl(seq)))
############################################################test part
'''
seq='TTTTCGCTATTTATGAAAATTTTCCG' #test sequence
seq='GATCCGGCGCGTGAGTTCACCATGANTCAGTCAGCACCGCT' #test sequence
#seq=revl(seq)
pos=getpos(seq)
#print(seq)
if pos==1:
    seq=revl(seq)
    pos=getpos(seq)
else:
    pos=pos
if pos==1:
    print('over')
'''

#print(int(suf_dic[for_start-1])+1)
#########################################################fastq files
reads=open("./sample/reads_1.fq","r")
sam=open("./output/reads_1.sam","w")
rname="gi|9626243|ref|NC_001416.1|"
sam.write("@HD"+"\t"+"VN:1.0"+"\t"+"SO:unsorted"+"\n")
sam.write("@SQ"+"\t"+"SN:"+rname+"\t"+"LN:48502"+"\n")
sam.write("@PG"+"\t"+"ID:nt"+"\t"+"PN:nt"+"\t"+"VN:1.0.0"+"\n")
ln=0
pos=0
seq=""
cigar=''
flag=0

for line in reads:
    ln+=1
    line=line.strip("\n")
    if ln%4==1:
        qname=str(line[1:])
    if ln%4==2:
        seq=str(line)
        #if "N" not in line:
        flag=0
        rname="gi|9626243|ref|NC_001416.1|"
        #print("seq: "+str(line))
        pos=getpos(seq)
    #print("ln: "+str(ln))
    #pos=sf[for_end-1]+1
    #cigar=str(len(seq))+"M"
    cigar=str(len(seq))+"M"
    #print(cigar)
    #flag=0
    mapq=0
   # cigar=str(len(line))+"M"
    if ln%4==0:
        qual=line
    rnext='*'
    pnext='0'
    tlen='0'
    #seq='0'
    if(ln%4==0):
        sam.write(str(qname)+'\t'+str(flag)+'\t'+str(rname)+'\t'+str(pos)+'\t'+str(mapq)+'\t'+str(cigar)+'\t'+str(rnext)+'\t'+str(pnext)+'\t'+str(tlen)+'\t'+str(seq)+'\t'+str(qual)+'\n')
sam.close()
#sfile=open("./output/suffix_ls.txt","w")
print("sam file outputs to folder: output")
#sfile.write("@HD"+"\t"+"VN:1.0"+"\t"+"SO:unsorted"+"\n")
#sfile.write("SQ"+"\t"+"SN:"+rname+"\t"+"LN:48502"+"\n")
#sfile.write("PG"+"\t"+"ID:nt"+"\t"+"PN:nt"+"\t"+"VN:1.0.0"+"\n")
#sfile.write(str(qname)+'\t'+str(flag)+'\t'+str(rname)+'\t'+str(pos)+'\t'+str(mapq)+'\t'+str(cigar)+'\t'+str(rnext)+'\t'+str(pnext)+'\t'+str(tlen)+'\t'+str(seq)+'\t'+str(qual)+'\n')
'''
for i in sf:
    sfile.write(str(i))
    sfile.write('\n')
'''
sam.close()
reads.close()
