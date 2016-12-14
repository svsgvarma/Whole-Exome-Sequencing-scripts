#!/bin/bash

CON="/home/varma/softwares/annovar/convert2annovar.pl"
VAR="/home/varma/softwares/annovar/annotate_variation.pl"
HUMANDB="/home/varma/softwares/annovar/humandb/"
CADD="/home/varma/softwares/annovar/CADD/"
F1="/home/varma/proj_exome/annotation_parser_vcf/"

######---- split files for CADD scores------
WDIR="/home/varma/proj_exome/work_WES/"
TMP="/home/varma/proj_exome/work_WES/TMP/"
INPUT="/home/varma/proj_exome/DVT2_WES_VCF/"
OUTPUT="/home/varma/proj_exome/DVT2_WES_ANNOTATION_ANNOVAR/"

#####----INPUT------

#./runwfall_annotation_norm_ANNOVAR_GATK_forallsamples.bash 123456789101112 FSCD1FSCD2FSCD3FSCD4FSCD5FSCD6FSCD7FSCD8FSCD9FSCD10FSCD11FSCD12
#./runwfall_annotation_norm_ANNOVAR_GATK_forallsamples.bash 123456 FSCD1FSCD2FSCD3FSCD4FSCD5FSCD6
#./runwfall_annotation_norm_ANNOVAR_GATK_forallsamples.bash 789101112 FSCD7FSCD8FSCD9FSCD10FSCD11FSCD12

##############----examples------
#(4) Run either score.sh or score_anno.sh on a gzip or block-gezipped VCF file:
#bin/score.sh path-to-variants/myVariants.vcf.gz path-to-output/output.tsv.gz
#bin/score.sh path-to-variants/myVariants.vcf.gz path-to-output/output.tsv.gz
#EX:
#bin/score.sh test.vcf.gz test_output.tsv.gz
#bin/score_anno.sh  ex3.vcf.gz  ex3.tsv.gz
#/home/varma/softwares/CADD_v1.2/bin/score_anno.sh  /home/varma/softwares/CADD_v1.2/ex3.vcf.gz  /home/varma/softwares/CADD_v1.2/ex3.tsv.gz

###############----run work flow-------
#/home/varma/softwares/CADD_v1.2/bin/score_anno.sh /home/varma/proj/work_WGS/V6_12345_4_variants_norm_annovar.hg19_multianno.vcf.gz /home/varma/proj/work_WGS/V6_12345_4_variants_norm_annovar.hg19_multianno.tsv.gz

#cat -n V6_12345_10_variants_step2.vcf | less
#except header
#tail -n+145 V6_12345_10_variants_step2.vcf | less
#V6_12345_22_variants_SNPnINDELcalls_PhasebyT.vcf
#for ch in {1..22} #10 12 19 22 #{6..22}
#do
######-----preparing left-normalization vcf file-------

mkdir -p $TMP

/home/varma/softwares/bcftools-1.1/bcftools norm -m-both -o ''$TMP''samp''$1''_''$2''_variants_GATK_step1.vcf ''$INPUT''samp''$1''_''$2''_variants_SNPnINDELcalls.vcf 

/home/varma/softwares/bcftools-1.1/bcftools norm -f /home/varma/softwares/human_g1k_v37.fasta -o ''$OUTPUT''samp''$1''_''$2''_variants_GATK_step1_norm.vcf ''$TMP''samp''$1''_''$2''_variants_GATK_step1.vcf

######-----without CADD scores ---------------
/home/varma/softwares/annovar/table_annovar.pl ''$OUTPUT''samp''$1''_''$2''_variants_GATK_step1_norm.vcf /home/varma/softwares/annovar/humandb/ -buildver hg19 -out ''$OUTPUT''samp''$1''_''$2''_variants_GATK_step2_annovar -remove -protocol refGene,1000g2014oct_all,snp138,clinvar_20150629,exac03,esp6500siv2_all,ljb23_sift,ljb23_pp2hvar,ljb23_pp2hdiv,ljb23_gerp++ -operation g,f,f,f,f,f,f,f,f,f -nastring . -vcfinput

######-----CADD--database--annotation
#for header
#rm ''$TMP''xa*
#rm ''$TMP''CADD_''$ch''_*

rm ''$TMP''samp''$1''_''$2''_variants_GATK_step1.vcf

