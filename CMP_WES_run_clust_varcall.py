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

### python script.py 1:(--sampl family_no) 2:(--sampl ID_no) 3:(--software dir) -p 4:(--cores) -C 5:(--runtime class) --EXEC
### python CMP_WES_run_clust_varcall.py 2 17671 2 -p 8 -C C --EXEC

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
parser.add_argument("b", type=int)
#parser.add_argument("c")
parser.add_argument("d", type=int)
parser.add_argument("-p", dest="nbProc", type=int, required=False, default=1)
parser.add_argument("-C", dest="classe", choices=["A","B","C","D"], required=False, default="A")
parser.add_argument("--EXEC",action='store_true', default=False)
args = parser.parse_args()

a=sys.argv[1]
b=sys.argv[2]
#c=sys.argv[3]
d=sys.argv[3]

PICARD="/san2/home/varma/softwares"+d+"/picard-tools-1.119/" 				########------ValidateSamFile.jar #SortSam.jar  #MarkDuplicates.jar
GATK="/san2/home/varma/softwares"+d+"/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar"
SAMTOOLS="/san2/home/varma/softwares"+d+"/samtools-1.1/samtools"
BCFTOOLS="/san2/home/varma/softwares"+d+"/bcftools-1.2/bcftools"
JAVA="/san2/home/varma/softwares"+d+"/jdk1.8.0_11/bin/java"
HG19="/san2/home/varma/softwares"+d+"/human_g1k_v37.fasta"
SOFT="/san2/home/varma/softwares"+d+"/"

INPUT_BAM="/san2/home/varma/proj_exome/CMP_WES_BAM/"
WORK_DIR="/san2/home/varma/proj_exome/work_WES/"
TMP="/san2/home/varma/proj_exome/work_WES/TMP/"

print (args)
MODE= args.EXEC
print (MODE)

i = a; j = b;
jobname = a+"_"+b+"_"+d
print "Jobname: ", jobname
cmdList = []


######## Sorts reads according to coordinate, first 2 lanes (all lanes....)------------
#cmdFls3 = JAVA+" -Xmx8g -jar "+PICARD+"SortSam.jar INPUT="+INPUT_BAM_2+" OUTPUT="+WORK_DIR+"2_sort.bam TMP_DIR="+TMP+" SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=10000000"
#cmdFls3 = SAMTOOLS+" sort -m 4G "+INPUT_BAM+"V6_nochr_"+str(i)+".bam "+WORK_DIR+"V6_"+str(i)+"_sort.bam"

######## Index BAM-files
#$SAMTOOLS index ''$WORK_DIR''V6_"+str(i)+"_sort.bam
#or
#cmdFls4 = JAVA+" -Xmx4g -jar "+PICARD+"BuildBamIndex.jar INPUT="+WORK_DIR+"V6_"+str(i)+"_nochr_sort.bam"

######## additional FLAGSTAT
#SAMTOOLS flagstat ''$WORK_DIR''V60_sort.bam > ''$WORK_DIR''V60_sort_flagstat.stats
#or
#cmdFls5 = JAVA+" -Xmx4g -jar "+GATK+" -T FlagStat -R "+HG19+" -I "+WORK_DIR+"V6_"+str(i)+"_nochr_sort.bam -o "+WORK_DIR+"V6_"+str(i)+"_sort_flag.txt"

#####################------------
########----GATK-----
########----Use picard to mark duplicate reads (Mark Duplicates) (first 1 lane or all lanes....)----------------

cmdFls6 = JAVA+" -Xmx10g -XX:ParallelGCThreads=10 -jar "+PICARD+"MarkDuplicates.jar INPUT="+INPUT_BAM+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort.bam OUTPUT="+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark.bam METRICS_FILE="+TMP+"tmp_"+str(i)+"_"+str(j)+"_file REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=True AS=true" #MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=500"

########----Index BAM-files----
#cmdFls7 = JAVA+" -Xmx8g -jar "+PICARD+"BuildBamIndex.jar INPUT="+WORK_DIR+"V6_"+str(i)+"_"+str(j)+"_sort_indelreln_markdup.bam"

### Local realignment around indels
######## Indel intervals marking------------------------

cmdFls8 = JAVA+" -Xmx10g -jar "+GATK+" -T RealignerTargetCreator -dcov 1600 -nt 10 -log "+TMP+"tmp_"+str(i)+"_"+str(j)+"_file -R "+HG19+" -I "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark.bam -o "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_realn.intervals --known "+SOFT+"1000G_phase1.indels.b37.vcf --known "+SOFT+"Mills_and_1000G_gold_standard.indels.b37.vcf"

######## Indel intervals Realignment (Local realignment around indels)---------------
cmdFls9 = JAVA+" -Xmx10g -jar "+GATK+" -T IndelRealigner -dcov 1600 -log "+TMP+"tmp_"+str(i)+"_"+str(j)+"_file -R "+HG19+" -I "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark.bam -targetIntervals "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_realn.intervals -o "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln.bam" 

# -fixMisencodedQuals"

########----Base quality score recalibration
########----Perform quality recalibration with GATK (Base quality score recalibration:CountCovariates)----------------
# Removed -nct 8 and used -rf BadCigar

cmdFls10 = JAVA+" -Xmx10g -jar "+GATK+" -T BaseRecalibrator -dcov 1600 -nct 10 -rf BadCigar -l INFO -log "+TMP+"tmp_"+str(i)+"_"+str(j)+"_file -R "+HG19+" -I "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln.bam -o "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal.grp -knownSites "+SOFT+"dbsnp_138.b37.vcf -knownSites "+SOFT+"1000G_phase1.indels.b37.vcf -knownSites "+SOFT+"Mills_and_1000G_gold_standard.indels.b37.vcf -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate"

cmdFls11 = JAVA+" -Xmx10g -jar "+GATK+" -T PrintReads -dcov 1600 -nct 10 -l INFO -R "+HG19+" -I "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln.bam -o "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal.bam -BQSR "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal.grp"

#######-----Call variants using GATK either HaplotypeCaller or Unified Genotyper-----
#-T UnifiedGenotyper --read_filter BadCigar -glm BOTH

cmdFls12 = JAVA+" -Xmx16g -jar "+GATK+" -T HaplotypeCaller -dcov 1600 -nct 10 --read_filter BadCigar -stand_call_conf 30 -stand_emit_conf 10 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -log "+TMP+"tmp_"+str(i)+"_"+str(j)+"_file -R "+HG19+" -I "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal.bam --dbsnp "+SOFT+"dbsnp_138.b37.vcf -o "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls.g.vcf"

#######-----Some clean up-----
# ; rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal.ba*

cmdFls13 = "rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark.ba* ; rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln.ba*; rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal.grp; rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_realn.intervals ; rm "+TMP+"tmp_"+str(i)+"_"+str(j)+"_file ; rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal.ba*"


#; rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal.ba*
#; rm "+WORK_DIR+"samp"+str(i)+"_withoutgel_"+str(j)+"_sort_mark_indelreln_recal_SNPnINDELcalls*"


#######-----Some clean up-----
#cmdFls20 = "rm "+WORK_DIR+"V6_"+str(i)+"_"+str(j)+"_sort_indelreln_markdup_recal_SNPnINDELcalls*"

#######-----Call variants using GATK VariantRecalibrator and ApplyRecalibration (VQSR) -----
#cmdFls_test = "touch /san2/home/varma/WORK/testfls/new_emptyfl"+str(pers_no)+"_"+str(i)+".txt"
#cmdList.append(cmd1)
#cmdList.append(cmdFls)
#cmdList.append(cmdFls3)
#cmdList.append(cmdFls4)
#cmdList.append(cmdFls5)

cmdList.append(cmdFls6)
#cmdList.append(cmdFls7)
cmdList.append(cmdFls8)
cmdList.append(cmdFls9)
cmdList.append(cmdFls10)
cmdList.append(cmdFls11)
cmdList.append(cmdFls12)
cmdList.append(cmdFls13)

#cmdList.append(cmdFls_test)
print ("job has been done...")
#cmdDirec = "touch /san2/home/varma/WORK/testfls/new_emptyfl1.txt"
#cmdList.append(cmdDirec)

launch(cmdList,args.classe, jobname, str(args.nbProc), MODE)

