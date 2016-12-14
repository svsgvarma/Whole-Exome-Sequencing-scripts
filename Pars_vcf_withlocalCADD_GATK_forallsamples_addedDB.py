#!/usr/bin/python

#Script to Pars vcf...
#####------Inputs-------
# script.py INPUT1=sampily_num INPUT2=sampleIDs INPUT3=orderedsampleIDs_withcomma(,) INPUT4=Filteringtype(1-2-3)


#./Pars_vcf_withlocalCADD_GATK_forallsamples_addedDB.py 1234 DBRIDJERDMARDVIN DBRI,DJER,DMAR,DVIN 1

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
b3 = str(sys.argv[4])

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
owrite = open(TMP+"samp"+c1+"_"+c2+"_variants_Parser_ANNOVAR_OUT.txt",'w')

srchdb={}
def createCADD(CADD):
	######-----CADD--database--annotation
	#cmd = CADD > "+TMP+"xa
	#for header 
	cmd1 = "head -141 "+AN_ANNOVAR+"samp"+c1+"_"+b1+"_variants_GATK_step1_norm.vcf > "+TMP+"CADD_"+c2+"_step1_vcfhdr" #-144 "+WDIR+"samp"+c1+"_"+c2+"_variants_step1_norm.vcf
	#---add header to the file
	cmd2 = "cat "+TMP+"CADD_"+c2+"_step1_vcfhdr "+TMP+"xa > "+TMP+"xa_in"
	#---rm ''$TMP''xa$i$i.gz
	cmd3 = "gzip "+TMP+"xa_in"
	#---Run CADD-----
	cmd4 = "/home/varma/softwares/CADD_v1.2/bin/score_anno.sh "+TMP+"xa_in.gz "+TMP+"xa_in.tsv.gz"
	#---unzip-----
	cmd5 = "less "+TMP+"xa_in.tsv.gz > "+TMP+"samp"+c1+"_"+c2+"_CADD_anno.txt"
	######-----some clean up----
	cmd6 = "rm "+TMP+"*_in* ; rm "+TMP+"xa*"
	cmd7 = "rm "+TMP+"CADD_"+c2+"_*"
	######-----
	g1 =  subprocess.check_output(cmd1, shell=True)
	g2 =  subprocess.check_output(cmd2, shell=True)
	g3 =  subprocess.check_output(cmd3, shell=True)
	g4 =  subprocess.check_output(cmd4, shell=True)
	g5 =  subprocess.check_output(cmd5, shell=True)
	g6 =  subprocess.check_output(cmd6, shell=True)
	g7 =  subprocess.check_output(cmd7, shell=True)

def srchdb(RR):
	try:
		True
		cmdFls1 = "grep -P -w '"+str(RR[8]+"\t"+RR[9]+"\t([^\t]+)\t"+RR[10])+"' "+TMP+"samp"+c1+"_"+c2+"_CADD_anno.txt"	
		grepout =  subprocess.check_output(cmdFls1, shell=True)
	except subprocess.CalledProcessError:
		False
		grepout = "."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"\n"
	return grepout

owrite.write(str("annovar_chr	annovar_pos_1	annovar_ref	annovar_alt	allele_freq	nb_homoz	nb_heteroz	quality	position	ref	alt	rsID	DP	QD	ReadPosRankSum	SOR	")+Joinall_gvcfs1+str("	annovar_loc	annovar_gene_name	annovar_effect(ExonicFunc)	annovar_exonic_effects(AAChange)	Clinvar	ExAC_ALL	1000genome	ESP(exome sequencing project) 	SIFT(mean a variant with score<0.05 is predicted as deleterious)	 PolyPhen2:HVAR	 PolyPhen2:HDIV	GERP++_annotation(Higher the score, the more conserved the site)	CADD_Chrom	Pos	Ref	Anc	Alt	Type	Length	isTv	isDerived	AnnoType	Consequence	ConsScore	ConsDetail	GC	CpG	mapAbility20bp	mapAbility35bp	scoreSegDup	priPhCons	mamPhCons	verPhCons	priPhyloP	mamPhyloP	verPhyloP	GerpN	GerpS	GerpRS	GerpRSpval	bStatistic	mutIndex	dnaHelT	dnaMGW	dnaProT	dnaRoll	mirSVR-Score	mirSVR-E	mirSVR-Aln	targetScan	fitCons	cHmmTssA	cHmmTssAFlnk	cHmmTxFlnk	cHmmTx	cHmmTxWk	cHmmEnhG	cHmmEnh	cHmmZnfRpts	cHmmHet	cHmmTssBiv	cHmmBivFlnk	cHmmEnhBiv	cHmmReprPC	cHmmReprPCWk	cHmmQuies	EncExp	EncH3K27Ac	EncH3K4Me1	EncH3K4Me3	EncNucleo	EncOCC	EncOCCombPVal	EncOCDNasePVal	EncOCFairePVal	EncOCpolIIPVal	EncOCctcfPVal	EncOCmycPVal	EncOCDNaseSig	EncOCFaireSig	EncOCpolIISig	EncOCctcfSig	EncOCmycSig	Segway	tOverlapMotifs	motifDist	motifECount	motifEName	motifEHIPos	motifEScoreChng	TFBS	TFBSPeaks	TFBSPeaksMax	isKnownVariant	ESP_AF	ESP_AFR	ESP_EUR	TG_AF	TG_ASN	TG_AMR	TG_AFR	TG_EUR	minDistTSS	minDistTSE	GeneID	FeatureID	CCDS	GeneName	cDNApos	relcDNApos	CDSpos	relCDSpos	protPos	relProtPos	Domain	Dst2Splice	Dst2SplType	Exon	Intron	oAA	nAA	Grantham	PolyPhenCat	PolyPhenVal	SIFTcat	SIFTval	CADD_RawScore	PHRED")+"\n")

out = open(TMP+"xa", "w")
totlns= len(infile3)
count1 = 0; count2 = 0;count3 = 0; dit1={}

while (count1 < totlns):
	#-----file1------
	hom = 0; het = 0;
	#if infile1[count1][0]!="#":
	elm1 = infile1[count1].strip()
	####-----[9:14] is for 5 individuals , [9:15] is for 6 individuals,[9:21] is for 12 individuals. change here if num of individuals increases.
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
	lst2=[];lst22=[];lst3=[];lst4=[]
	for xx in range(len(bb)):
		####---Heterozygous-check----
		lst2.append(re.compile("[|/]").split(el1[xx])[0] != re.compile("[|/]").split(el1[xx])[1])
		####---Homozygous-check----
		lst22.append((re.compile("[|/]").split(el1[xx])[0] == "1") and (re.compile("[|/]").split(el1[xx])[1] == "1"))
		lst3.append(elm2[xx][0:3]+str("	")+elm2[xx][4:])
	#Joinall_gvcfs2 = '&&'.join(lst2)
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

	#----check hom & het in all-------------------
	#if (((ss3[0] == "1") and (ss3[1] == "1")) or ((ss5[0] == "1") and (ss5[1] == "1")) ):
	#if (((ss3[0] == "1") and (ss3[1] == "1"))):
	#if (((ss3[0] == "1") and (ss3[1] == "1")) or ((ss5[0] == "1") and (ss5[1] == "1")) ):
	#	False
	#elif ((ss1[0] != ss1[1]) and (ss2[0] != ss2[1]) and (ss3[0] == ss3[1]) and (ss4[0] != ss4[1]) and (ss5[0] == ss5[1])):
	#elif ((ss1[0] != ss1[1]) and (ss2[0] != ss2[1]) and (ss3[0] == ss3[1])):
	#if ((ss1[0] != ss1[1]) and (ss2[0] != ss2[1]) and (ss3[0] != ss3[1]) and (ss4[0] != ss4[1])):

	###---If all samples are HET use as below or else specify each sample whether true or false...
	#if sum(lst2)==len(bb):
	###--Filter for Dominant Homozygous cases (1/1) # order matters recheck!!!
	if (lst22[0]==True):
		False
	###--Filter for Dominant Heterozygous cases (0/1 or 1/0 == True ; 1/1 == False) # order matters recheck!!!
	elif (lst2[0]==False and lst2[1]==True and lst2[2]==True):
		True
		#----without rsIDs;MAF<1%;rm_syn;exonic----
		# and ((dit1['exac02'] == '.') or (dit1['exac02'] <= "0.01")) # and (dit1['Func.refGene'] == 'exonic') or (dit1['Func.refGene'] == 'ncRNA_exonic')
		#and ((dit1['ExonicFunc.refGene'] == '.') or (dit1['ExonicFunc.refGene'] != 'synonymous_SNV'))
		###########----FILTERS-----	
		#####------Input-1-----------------------------
		if (b3=="1"):
			True
			#####------SNPDB+MAF+RM_SYN+EXONIC-----------------------------
			if ((dit1['snp138'] == '.') and ((dit1['1000g2014oct_all'] == '.') or (dit1['1000g2014oct_all'] <= "0.01")) and (dit1['ExonicFunc.refGene'] != 'synonymous_SNV') and (dit1['Func.refGene'] == 'exonic')):
				True
				#----CADD output----------------
				out.write(str(elm4[0]+"\t"+elm4[1]+"\t"+elm4[2]+"\t"+elm4[3]+"\t"+elm4[4])+"\n")#+"\t"+elm4[5]
				owrite.write(str(glm2[0])+"\t"+str(glm2[2])+"\t"+str(glm2[3])+"\t"+str(glm2[4])+"\t"+str(glm2[5])+"\t"+str(hom)+"\t"+str(het)+"\t"+str(glm2[6])+"\t"+str(glm3)+"\t"+str(glm4)+"\t"+dit1['snp138']+"\t"+str(dit1['DP'])+"\t"+str(dit1['QD'])+"\t"+str(dit1['ReadPosRankSum'])+"\t"+str(dit1['SOR'])+"\t"+Joinall_gvcfs3.rstrip('\n')+"\t"+dit1['Func.refGene']+"\t"+dit1['Gene.refGene']+"\t"+dit1['ExonicFunc.refGene']+"\t"+dit1['AAChange.refGene']+"\t"+dit1['clinvar_20150629']+"\t"+dit1['ExAC_ALL']+"\t"+dit1['1000g2014oct_all']+"\t"+dit1['esp6500siv2_all']+"\t"+dit1['ljb23_sift']+"\t"+dit1['ljb23_pp2hvar']+"\t"+dit1['ljb23_pp2hdiv']+"\t"+dit1['ljb23_gerp++']+"\n") #+"\t"+str(srchcadd)) #+"\n")#
				count3 = count3 + 1
				
		#####------Input-2-----------------------------		
		elif (b3=="2"):
			True
			#####------MAF+RM_SYN+EXONIC----------------------------- 
			if (((dit1['1000g2014oct_all'] <= "0.01") or (dit1['1000g2014oct_all'] == '.')) and (dit1['ExonicFunc.refGene'] != 'synonymous_SNV') and (dit1['Func.refGene'] == 'exonic') ):
				True
				#----CADD output----------------
				out.write(str(elm4[0]+"\t"+elm4[1]+"\t"+elm4[2]+"\t"+elm4[3]+"\t"+elm4[4])+"\n")#+"\t"+elm4[5]
				owrite.write(str(glm2[0])+"\t"+str(glm2[2])+"\t"+str(glm2[3])+"\t"+str(glm2[4])+"\t"+str(glm2[5])+"\t"+str(hom)+"\t"+str(het)+"\t"+str(glm2[6])+"\t"+str(glm3)+"\t"+str(glm4)+"\t"+dit1['snp138']+"\t"+str(dit1['DP'])+"\t"+str(dit1['QD'])+"\t"+str(dit1['ReadPosRankSum'])+"\t"+str(dit1['SOR'])+"\t"+Joinall_gvcfs3.rstrip('\n')+"\t"+dit1['Func.refGene']+"\t"+dit1['Gene.refGene']+"\t"+dit1['ExonicFunc.refGene']+"\t"+dit1['AAChange.refGene']+"\t"+dit1['clinvar_20150629']+"\t"+dit1['ExAC_ALL']+"\t"+dit1['1000g2014oct_all']+"\t"+dit1['esp6500siv2_all']+"\t"+dit1['ljb23_sift']+"\t"+dit1['ljb23_pp2hvar']+"\t"+dit1['ljb23_pp2hdiv']+"\t"+dit1['ljb23_gerp++']+"\n") #+"\t"+str(srchcadd)) #+"\n")#
				count3 = count3 + 1

		#####------Input-3----------------------------- 
		elif (b3=="3"):
			True
			#####------SNPDB+MAF+NON_EXONIC----------------------------- 
			if ((dit1['snp138'] == '.') and ((dit1['1000g2014oct_all'] <= "0.01") or (dit1['1000g2014oct_all'] == '.')) and (dit1['Func.refGene'] != 'exonic')):
				True
				#----CADD output----------------
				out.write(str(elm4[0]+"\t"+elm4[1]+"\t"+elm4[2]+"\t"+elm4[3]+"\t"+elm4[4])+"\n")#+"\t"+elm4[5]
				owrite.write(str(glm2[0])+"\t"+str(glm2[2])+"\t"+str(glm2[3])+"\t"+str(glm2[4])+"\t"+str(glm2[5])+"\t"+str(hom)+"\t"+str(het)+"\t"+str(glm2[6])+"\t"+str(glm3)+"\t"+str(glm4)+"\t"+dit1['snp138']+"\t"+str(dit1['DP'])+"\t"+str(dit1['QD'])+"\t"+str(dit1['ReadPosRankSum'])+"\t"+str(dit1['SOR'])+"\t"+Joinall_gvcfs3.rstrip('\n')+"\t"+dit1['Func.refGene']+"\t"+dit1['Gene.refGene']+"\t"+dit1['ExonicFunc.refGene']+"\t"+dit1['AAChange.refGene']+"\t"+dit1['clinvar_20150629']+"\t"+dit1['ExAC_ALL']+"\t"+dit1['1000g2014oct_all']+"\t"+dit1['esp6500siv2_all']+"\t"+dit1['ljb23_sift']+"\t"+dit1['ljb23_pp2hvar']+"\t"+dit1['ljb23_pp2hdiv']+"\t"+dit1['ljb23_gerp++']+"\n") #+"\t"+str(srchcadd)) #+"\n")#
				count3 = count3 + 1
		count2 = count2 + 1
	count1 = count1 + 1

print "Total_Number of lines","\t",count1
print "After homozygous and heterozygous check","\t",count2
print "Varients without rsIDs","\t",count3

out.close()
owrite.close()

#-----Use dict1-------
eachCADD  = createCADD(count2)	#input count3 is not necessary.
#-------------Combine both ANNOVAR and CADD outs------------------------------------
#out1 = open(TMP+"samp"+c1+"_"+c2+"_Variants_GATK_ANNOVAR_CADD_homhet.txt","w")

if (b3=="1"):
	True
	#print "Hi this is working"
	out1 = open(AN_CADD+"samp"+c1+"_"+c2+"_Variants_GATK_ANNOVAR_CADD_het_MAF1%_nosnp138_nosyn_exonic.txt","w")
elif (b3=="2"):
	True
	out1 = open(AN_CADD+"samp"+c1+"_"+c2+"_Variants_GATK_ANNOVAR_CADD_het_MAF1%_nosyn_exonic.txt","w")
elif (b3=="3"):
	True
	out1 = open(AN_CADD+"samp"+c1+"_"+c2+"_Variants_GATK_ANNOVAR_CADD_het_MAF1%_nosnp138_non-exonic.txt","w")

#out2 = open(TMP+"samp"+c1+"_"+c2+"_Variants_GATK_head","w")
fl1 = open(TMP+"samp"+c1+"_"+c2+"_variants_Parser_ANNOVAR_OUT.txt",'r').readlines()

fl2 = fl1[1:]
fl3 = fl1[0]
#out2.write(str(fl3))
#out2.close()
out1.write(str(fl3))

for g in fl2:
	gg = g.strip().split("\t")
	#-----Use dict2-------
	CADDout = srchdb(gg)
	out1.write(str(g.rstrip())+"\t"+str(CADDout))#+"\n")
out1.close()

cmd11 = "rm "+TMP+"samp"+c1+"_"+c2+"_variants_Parser_ANNOVAR_OUT.txt"
cmd12 = "rm "+TMP+"samp"+c1+"_"+c2+"_CADD_anno.txt"

g11 =  subprocess.check_output(cmd11, shell=True)
g12 =  subprocess.check_output(cmd12, shell=True)

#-------------------------------------------------
#
#

