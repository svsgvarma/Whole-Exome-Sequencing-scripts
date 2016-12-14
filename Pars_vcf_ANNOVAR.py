#!/usr/bin/python

#Script to Pars vcf...
#####------Inputs-------
# script.py INPUT1=sampily_num INPUT2=sampleIDs INPUT3=orderedsampleIDs_withcomma(,) INPUT4=Filteringtype(1-2-3)

#./Pars_vcf_ANNOVAR.py 123456789101112 FSCD1FSCD2FSCD3FSCD4FSCD5FSCD6FSCD7FSCD8FSCD9FSCD10FSCD11FSCD12 FSCD1,FSCD2,FSCD3,FSCD4,FSCD5,FSCD6,FSCD7,FSCD8,FSCD9,FSCD10,FSCD11,FSCD12
#./Pars_vcf_ANNOVAR.py 123456 FSCD1FSCD2FSCD3FSCD4FSCD5FSCD6 FSCD1,FSCD2,FSCD3,FSCD4,FSCD5,FSCD6
#./Pars_vcf_ANNOVAR.py 789101112 FSCD7FSCD8FSCD9FSCD10FSCD11FSCD12 FSCD7,FSCD8,FSCD9,FSCD10,FSCD11,FSCD12


import sys
import re
import os
import tempfile
import commands
import subprocess
import subprocess32
from subprocess import *
from subprocess import call

c1 = str(sys.argv[1])
b1 = str(sys.argv[2])
b2 = str(sys.argv[3])

WDIR="/home/varma/proj_exome/work_WES/"
TMP="/home/varma/proj_exome/work_WES/TMP/"
AN_ANNOVAR="/home/varma/proj_exome/DVT2_WES_ANNOTATION_ANNOVAR/"
AN_CADD="/home/varma/proj_exome/DVT2_WES_ANNOTATION_ANNOVAR_CADD/"

bb= list(map(str, b2.split(',')))
def magic(numList):
    s = ''.join(map(str, numList))
    return str(s)
c2 = str(magic(bb))

lst1=[];
for xx in range(len(bb)):
	lst1.append("Sample_"+str(bb[xx])+"_genotype"+str("	")+"Sample_"+str(bb[xx])+"_AD:DP:GQ:PL:TP")
	#lst2.append("(re.compile("+str("[|/]")+").split(el1["+str(xx)+"])[0] != re.compile("+str("[|/]")+").split(el1["+str(xx)+"])[1])")
	#lst3.append("elm2["+str(xx)+"][0:3]"+str("	")+"elm2["+str(xx)+"][4:]")
Joinall_gvcfs1 = '\t'.join(lst1)
#Joinall_gvcfs2 = ' and '.join(lst2)
#Joinall_gvcfs3 = '\t'.join(lst3)

#samp7_2153339_variants_GATK_step1_norm.vcf
#samp7_2153339_variants_GATK_step2_annovar.avinput
#samp7_2153339_variants_GATK_step2_annovar.hg19_multianno.vcf

infile1  = open(AN_ANNOVAR+"samp"+c1+"_"+b1+"_variants_GATK_step1_norm.vcf",'r').readlines()[141:] #[140:]  #[144:]
infile2  = open(AN_ANNOVAR+"samp"+c1+"_"+b1+"_variants_GATK_step2_annovar.hg19_multianno.vcf",'r').readlines()[164:] #[151:]  #[155:]  
infile3  = open(AN_ANNOVAR+"samp"+c1+"_"+b1+"_variants_GATK_step2_annovar.avinput",'r').readlines()
owrite = open(AN_ANNOVAR+"samp"+c1+"_"+c2+"_Variants_GATK_ANNOVAR.txt",'w')

owrite.write(str("annovar_chr	annovar_pos_1	annovar_ref	annovar_alt	allele_freq	nb_homoz	nb_heteroz	quality	position	ref	alt	rsID	DP	QD	ReadPosRankSum	SOR	")+Joinall_gvcfs1+str("	annovar_loc	annovar_gene_name	annovar_effect(ExonicFunc)	annovar_exonic_effects(AAChange)	Clinvar	ExAC_ALL	1000genome	ESP(exome sequencing project) 	SIFT(mean a variant with score<0.05 is predicted as deleterious)	 PolyPhen2:HVAR	 PolyPhen2:HDIV	GERP++_annotation(Higher the score, the more conserved the site)")+"\n")

out = open(TMP+"xa", "w")
totlns= len(infile3)
count1 = 0; count2 = 0;count3 = 0; dit1={}

while (count1 < totlns):
	#-----file1------
	hom = 0; het = 0;
	#if infile1[count1][0]!="#":
	elm1 = infile1[count1].strip()
	#[9:14] is for 5 individuals , [9:15] is for 6 individuals,[9:21] is for 12 individuals. change here if num of individuals increases.  
	elm2 = elm1.split("\t")[9:9+len(bb)] #or [9:21]
	#elm3 = elm2[0][0:3],elm2[0][4:] ,elm2[1][0:3],elm2[1][4:], elm2[2][0:3],elm2[2][4:], elm2[3][0:3],elm2[3][4:], elm2[4][0:3],elm2[4][4:]
	elm4 = elm1.split("\t")[0:5]
	#elm2_1= elm2[0] #elm2_5= elm2[0]
	#elm2_2= elm2[1] #elm2_3= elm2[1]
	#elm2_3= elm2[2] #elm2_1= elm2[2]
	#elm2_4= elm2[3] #elm2_4= elm2[3]
	#elm2_5= elm2[4] #elm2_2= elm2[4]
	el1 =[]
	for zyg in elm2:
		el = zyg[0:3]
		el1.append(el)
		#----count hom and het
		el2 = re.compile("[|/]").split(el)
		if el2[0] == "0" and el2[1] == "0":
			False
		elif el2[0] == el2[1]:
			hom+=1
		elif el2[0] != el2[1]:
			het+=1
	#----check het in all--------------------
	#ss1 = re.compile("[|/]").split(el1[0]) #ss5 (V6O)1 #B00GV6K=5=1 -het 4
	#ss2 = re.compile("[|/]").split(el1[1]) #ss3 (V6Q)2 #B00GV6N=3=2 -hom 2
	#ss3 = re.compile("[|/]").split(el1[2]) #ss1 (V6N)3 #B00GV6O=1=3 -hom 0	
	#ss4 = re.compile("[|/]").split(el1[3]) #ss4 (V6P)4 #B00GV6P=4=4 -het 3
	#ss5 = re.compile("[|/]").split(el1[4]) #ss2 (V6K)5 #B00GV6Q=2=5 -het 1
	lst2=[];lst3=[];lst4=[]
	for xx in range(len(bb)):
		lst3.append(elm2[xx][0:3]+str("	")+elm2[xx][4:])
	Joinall_gvcfs3 = '\t'.join(lst3)

	#-----MAF-----------------------------------
	col=elm1.split()
	#print  col[7].split(";")[1]
	genos=[]
	#-----file2-------
	#if infile2[count1][0]!="#":
	flm1 = infile2[count1]
	flm2 = flm1.split("\t")[7].split(";")
	#print elm2[1], elm2[3], elm2[5], elm2[15:18],elm2[21:30]
	for fl in flm2:
		fl = fl.split("=")
		if len(fl) == 1: #fl[0] == "NEGATIVE_TRAIN_SITE" or fl[0] == "ALLELE_END":
			False
		else:
			dit1[fl[0]] = fl[1]
	if "ReadPosRankSum" not in dit1:
		dit1['ReadPosRankSum'] = "NA"
		print "NA"
	#print dit1['ReadPosRankSum']
	#
	#-----file3-------
	glm1= infile3[count1]
	glm2 = glm1.split("\t")
	glm3 = '\t'.join(glm2[9:10])
	glm4 = '\t'.join(glm2[11:13])
	owrite.write(str(glm2[0])+"\t"+str(glm2[2])+"\t"+str(glm2[3])+"\t"+str(glm2[4])+"\t"+str(glm2[5])+"\t"+str(hom)+"\t"+str(het)+"\t"+str(glm2[6])+"\t"+str(glm3)+"\t"+str(glm4)+"\t"+dit1['snp138']+"\t"+str(dit1['DP'])+"\t"+str(dit1['QD'])+"\t"+str(dit1['ReadPosRankSum'])+"\t"+str(dit1['SOR'])+"\t"+Joinall_gvcfs3.rstrip('\n')+"\t"+dit1['Func.refGene']+"\t"+dit1['Gene.refGene']+"\t"+dit1['ExonicFunc.refGene']+"\t"+dit1['AAChange.refGene']+"\t"+dit1['clinvar_20150629']+"\t"+dit1['ExAC_ALL']+"\t"+dit1['1000g2014oct_all']+"\t"+dit1['esp6500siv2_all']+"\t"+dit1['ljb23_sift']+"\t"+dit1['ljb23_pp2hvar']+"\t"+dit1['ljb23_pp2hdiv']+"\t"+dit1['ljb23_gerp++']+"\n") #+"\t"+str(srchcadd)) #+"\n")#
	count1 = count1 + 1

print "Total_Number of lines","\t",count1
owrite.close()



