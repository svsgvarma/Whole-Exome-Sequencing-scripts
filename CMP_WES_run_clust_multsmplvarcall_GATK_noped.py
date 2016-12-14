#!/usr/bin/env python

import os
import sys
#import paramik
import subprocess
import threading
import subprocess32 #local install for check_output (not in python 2.X)
import fnmatch
import re
import tempfile
import argparse

### python script.py 1:(--sampl comname) 2:(--sampl ID) 3:(--chr num) 4:(--software dir) -p 5:(--cores) -C 6:(--runtime class) --EXEC
### python CMP_WES_run_clust_multsmplvarcall_GATK_noped.py '1,2,3,4,5' '17200,17671,18288,19208,19237' 1 -p 8 -C C --EXEC

DEBUG = 0
EXEC  = 1

def launch(cmdList,classe, jobname, nbProc, MODE):
    cmd_tmpfile = tempfile.NamedTemporaryFile()
    cmd_tmpfile.write( "# @ job_type = serial\n" +\
    	"# @ class = " + classe + " \n" +\
        "# @ resources = ConsumableCpus(1) ConsumableMemory(20Gb)\n" +\
        "# @ output=log/log_"  + jobname + "\n" + \
        "# @ error=log/error_" + jobname + "\n" + \
        "# @ job_name = " + jobname + "\n" +\
        "# @ notification = error\n\n")
    cmd_tmpfile.write("# @ queue\n")
    for cmd in cmdList:
    	cmd_tmpfile.write(cmd+"\n")
    cmd_tmpfile.flush()

    if MODE == DEBUG:
        os.system( "cat " + cmd_tmpfile.name )
    if MODE == EXEC:
        os.system( "llsubmit " + cmd_tmpfile.name )

parser = argparse.ArgumentParser()
parser.add_argument("a")
parser.add_argument("b")
parser.add_argument("c")
#parser.add_argument("d", type=int)
parser.add_argument("-p", dest="nbProc", type=int, required=False, default=1)
parser.add_argument("-C", dest="classe", choices=["A","B","C","D"], required=False, default="A")
parser.add_argument("--EXEC",action='store_true', default=False)
args = parser.parse_args()

a=sys.argv[1]
b=sys.argv[2]
c=sys.argv[3]
#d=sys.argv[4]

PICARD="/san2/home/varma/softwares"+c+"/picard-tools-1.119/"                            ########------ValidateSamFile.jar #SortSam.jar  #MarkDuplicates.jar
GATK="/san2/home/varma/softwares"+c+"/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar"
SAMTOOLS="/san2/home/varma/softwares"+c+"/samtools-1.1/samtools"
BCFTOOLS="/san2/home/varma/softwares"+c+"/bcftools-1.2/bcftools"
JAVA="/san2/home/varma/softwares"+c+"/jdk1.8.0_11/bin/java"
HG19="/san2/home/varma/softwares"+c+"/human_g1k_v37.fasta"
SOFT="/san2/home/varma/softwares"+c+"/"

INPUT_BAM="/san2/home/varma/proj_exome/CMP_WES_BAM/"
WORK_DIR="/san2/home/varma/proj_exome/work_WES/"
TMP="/san2/home/varma/proj_exome/work_WES/TMP/"
PEDS="/san2/home/varma/proj_exome/PEDS/"

print (args)
MODE= args.EXEC
print (MODE)
aa= list(map(int, a.split(',')))
bb= list(map(int, b.split(',')))

def magic(numList):
    s = ''.join(map(str, numList))
    return int(s)
a1 = magic(aa)
b1 = magic(bb)

i = a1; j = b1;
jobname = a+"_"+b+"_"+c
print "Jobname: ", jobname
cmdList = [];lst=[]


for xx in range(len(bb)):
        lst.append("--variant "+str(WORK_DIR)+"samp"+str(aa[xx])+"_withoutgel_"+str(bb[xx])+"_sort_mark_indelreln_recal_SNPnINDELcalls.g.vcf")
Joinall_gvcfs = ' '.join(lst)


########################################################------
#######-----Call variants using GATK GenotypeGVCFs------
#--variant "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls.g.vcf
#cmdFls15 = JAVA+" -Xmx10g -jar "+GATK+" -T GenotypeGVCFs -nt 10 -R "+HG19+" -ped "+PEDS+"samp"+str(i)+"_withoutgel_"+str(j)+"_pedigreeIDS.ped --standard_min_confidence_threshold_for_calling 30 --standard_min_confidence_threshold_for_emitting 30 --variant "+WORK_DIR+"fam"+str(i)+"_samp"+str(bb[0])+"_sort_mark_indelreln_recal_SNPnINDELcalls.g.vcf --variant "+WORK_DIR+"fam"+str(i)+"_samp"+str(bb[1])+"_sort_mark_indelreln_recal_SNPnINDELcalls.g.vcf --variant "+WORK_DIR+"fam"+str(i)+"_samp"+str(bb[2])+"_sort_mark_indelreln_recal_SNPnINDELcalls.g.vcf --variant "+WORK_DIR+"fam"+str(i)+"_samp"+str(bb[3])+"_sort_mark_indelreln_recal_SNPnINDELcalls.g.vcf -o "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls.vcf"

cmdFls15 = JAVA+" -Xmx10g -jar "+GATK+" -T GenotypeGVCFs -nt 10 -R "+HG19+" --standard_min_confidence_threshold_for_calling 30 --standard_min_confidence_threshold_for_emitting 30 "+Joinall_gvcfs+" -o "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls.vcf"

#######------Run VQSR----
#######------2.Build the SNP recalibration model-----
#removed for exome -an MQ

cmdFls16 = JAVA+" -Xmx10g -jar "+GATK+" -T VariantRecalibrator -R "+HG19+" -input "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 "+SOFT+"hapmap_3.3.b37.vcf -resource:omni,known=false,training=true,truth=true,prior=12.0 "+SOFT+"1000G_omni2.5.b37.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 "+SOFT+"1000G_phase1.snps.high_confidence.b37.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "+SOFT+"dbsnp_138.b37.vcf -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode SNP -recalFile "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls_SNP.recal -tranchesFile "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls_SNP.tranches -rscriptFile "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls_SNP_plots.R"

#######------3.Apply the desired level of recalibration to the SNPs in the call set------
cmdFls17 = JAVA+" -Xmx10g -jar "+GATK+" -T ApplyRecalibration -R "+HG19+" -input "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls.vcf -mode SNP --ts_filter_level 99.5 -recalFile "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls_SNP.recal -tranchesFile "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls_SNP.tranches -o "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls_snps_raw_indels.vcf"

#######------5.Build the Indel recalibration model--------
cmdFls18 = JAVA+" -Xmx10g -jar "+GATK+" -T VariantRecalibrator -R "+HG19+" -input "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls_snps_raw_indels.vcf -resource:mills,known=false,training=true,truth=true,prior=12.0 "+SOFT+"Mills_and_1000G_gold_standard.indels.b37.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "+SOFT+"dbsnp_138.b37.vcf -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode INDEL --maxGaussians 4 -recalFile "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls_INDEL.recal -tranchesFile "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls_INDEL.tranches -rscriptFile "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls_INDEL_plots.R"

######-------6.Apply the desired level of recalibration to the Indels in the call set---------
cmdFls19 = JAVA+" -Xmx10g -jar "+GATK+" -T ApplyRecalibration -R "+HG19+" -input "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls_snps_raw_indels.vcf -mode INDEL --ts_filter_level 99.0 -recalFile "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls_INDEL.recal -tranchesFile "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls_INDEL.tranches -o "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_variants_SNPnINDELcalls.vcf"

######-------7.PhasebyTransmission-----
#cmdFls20 = JAVA+" -Xmx10g -jar "+GATK+" -T PhaseByTransmission -R "+HG19+" -ped "+PEDS+"samp"+str(i)+"_withoutgel_"+str(j)+"_pedigreeIDS.ped -V "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_variants_SNPnINDELcalls.vcf -o "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_variants_GATK_SNPnINDELcalls.vcf"

cmdFls_clean1 = "rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls*"

#; rm "+WORK_DIR+"fam"+str(i)+"_samp_*_sort_mark_indelreln_recal_SNPnINDELcalls.g.vcf*";
#rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls.g.vcf* ;  rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls.g.vcf* ;  rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls.g.vcf* ;  rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls.g.vcf* ;  rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls.g.vcf*" #; rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_variants_SNPnINDELcalls.vcf*"

######--------8.SAmtools---------
#samtools mpileup -ugf <ref.fa> <sample1.bam> <sample2.bam> <sample3.bam> | bcftools call -vmO z -o <study.vcf.gz>
#cmdFls21 = SAMTOOLS+" mpileup -go "+WORK_DIR+"fam12345_"+str(j)+"_variants_samtools_SNPnINDELcalls.bcf -f "+HG19+" "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_indelreln_markdup_recal.bam "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_indelreln_markdup_recal.bam "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_indelreln_markdup_recal.bam "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_indelreln_markdup_recal.bam "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_indelreln_markdup_recal.bam"
#cmdFls22 = BCFTOOLS+" call -vmO z -o "+WORK_DIR+"fam12345_"+str(j)+"_variants_samtools_SNPnINDELcalls.vcf.gz "+WORK_DIR+"fam12345_"+str(j)+"_variants_samtools_SNPnINDELcalls.bcf"
#cmdFls23 = BCFTOOLS+" filter -O z -o "+WORK_DIR+"fam12345_"+str(j)+"_variants_samtools_SNPnINDELcalls_filtered.vcf.gz -s LOWQUAL -i'%QUAL>10' "+WORK_DIR+"fam12345_"+str(j)+"_variants_samtools_SNPnINDELcalls.vcf.gz"
#cmdFls_clean2 = "rm "+WORK_DIR+"fam12345_"+str(j)+"_variants_samtools_SNPnINDELcalls.bcf ; rm "+WORK_DIR+"fam12345_"+str(j)+"_variants_samtools_SNPnINDELcalls.vcf.gz"

######--------9.Varscan---------
######----Perform Cross-sample Variant Calling

##--------1: Run SAMtools mpileup on all BAM files in a single command:
#samtools mpileup -B -q 1 -f reference.fasta sample1.bam sample2.bam sample3.bam > cohort.mpileup
#cmdFls24 = SAMTOOLS+" mpileup -B -q 1 -f "+HG19+" "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_indelreln_markdup_recal.bam "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_indelreln_markdup_recal.bam "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_indelreln_markdup_recal.bam "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_indelreln_markdup_recal.bam "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_indelreln_markdup_recal.bam > "+WORK_DIR+"fam12345_"+str(j)+"_variants_Varscan.mpileup"

##--------2:Run VarScan mpileup2snp to call SNVs:
#java -jar VarScan.jar mpileup2snp cohort.mpileup --vcf-sample-list samples.txt --min-coverage 10 --min-var-freq 0.20 --p-value 0.05 --output-vcf 1 >cohort.varScan.snp.vcf
#cmdFls25 = JAVA+" -jar "+VARSCAN+" mpileup2snp "+WORK_DIR+"fam12345_"+str(j)+"_variants_Varscan.mpileup --vcf-sample-list "+WORK_DIR+"fam12345_"+str(j)+"_samples.txt --min-coverage 10 --min-var-freq 0.20 --p-value 0.05 --output-vcf 1 > "+WORK_DIR+"fam12345_"+str(j)+"_variants_Varscan_SNP.vcf"

##--------3:Run VarScan mpileup2indel to call indels:
#java -jar VarScan.jar mpileup2indel cohort.mpileup --vcf-sample-list samples.txt --min-coverage 10 --min-var-freq 0.10 --p-value 0.10 --output-vcf 1 > cohort.varScan.indel.vcf
#cmdFls26 = JAVA+" -jar "+VARSCAN+" mpileup2indel "+WORK_DIR+"fam12345_"+str(j)+"_variants_Varscan.mpileup --vcf-sample-list "+WORK_DIR+"fam12345_"+str(j)+"_samples.txt --min-coverage 10 --min-var-freq 0.10 --p-value 0.10 --output-vcf 1 > "+WORK_DIR+"fam12345_"+str(j)+"_variants_Varscan_INDEL.vcf"

#######-------10. Some clean up-----
#cmdFls_clean3 = "rm "+WORK_DIR+"fam12345_"+str(j)+"_variants_Varscan.mpileup ; rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_indelreln_markdup_recal.ba* ;  rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_indelreln_markdup_recal.ba* ;  rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_indelreln_markdup_recal.ba* ;  rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_indelreln_markdup_recal.ba* ;  rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_indelreln_markdup_recal.ba*"

#######------Call variants using GATK VariantRecalibrator and ApplyRecalibration (VQSR) -----
#cmdFls_test = "touch /san2/home/varma/WORK/testfls/new_emptyfl"+str(pers_no)+"_"+str(i)+".txt"
#cmdList.append(cmd1)
#cmdList.append(cmdFls)
#cmdLiINPUT=st.append(cmdFls3)

cmdList.append(cmdFls15)
cmdList.append(cmdFls16)
cmdList.append(cmdFls17)
cmdList.append(cmdFls18)
cmdList.append(cmdFls19)
#cmdList.append(cmdFls20)
cmdList.append(cmdFls_clean1)

#cmdList.append(cmdFls21)
#cmdList.append(cmdFls22)
#cmdList.append(cmdFls23)
#cmdList.append(cmdFls_clean2)

#cmdList.append(cmdFls24)
#cmdList.append(cmdFls25)
#cmdList.append(cmdFls26)
#cmdList.append(cmdFls_clean3)

#cmdList.append(cmdFls_test)
print ("job has been done...")
#cmdDirec = "touch /san2/home/varma/WORK/testfls/new_emptyfl1.txt"
#cmdList.append(cmdDirec)

launch(cmdList,args.classe, jobname, str(args.nbProc), MODE)

