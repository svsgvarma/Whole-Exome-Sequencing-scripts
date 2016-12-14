#!/bin/bash

##--Script-1
./runwfall_annotation_norm_ANNOVAR_GATK_forallsamples.bash 134 DBRIDMARDVIN


##--Script-2
./Pars_vcf_ANNOVAR.py 134 DBRIDMARDVIN DBRI,DMAR,DVIN

##--Script-3
#####----INPUT-----
##./Pars_vcf_withlocalCADD_GATK_forallsamples_addedDB.py 7 2153339 15,2,33,39 1

#####----INPUT1: #####------SNPDB+MAF+RM_SYN+EXONIC----------------------------------

./Pars_vcf_withlocalCADD_GATK_forallsamples_addedDB.py 134 DBRIDMARDVIN DBRI,DMAR,DVIN 1


#####----INPUT2: #####------MAF+RM_SYN+EXONIC-----------------------------

./Pars_vcf_withlocalCADD_GATK_forallsamples_addedDB.py 134 DBRIDMARDVIN DBRI,DMAR,DVIN 2


#####----INPUT3: #####------SNPDB+MAF+NON_EXONIC--------------------------

./Pars_vcf_withlocalCADD_GATK_forallsamples_addedDB.py 134 DBRIDMARDVIN DBRI,DMAR,DVIN 3



