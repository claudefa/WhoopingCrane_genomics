#!/bin/bash
java -Xmx8g -jar /home/qvw641/bin/snpEff/snpEff.jar GrusAme /projects/mjolnir1/people/qvw641/WhoopingCrane/VCF/Concat/WC_outgroups.singl.outmiss2.vcf.gz | bgzip > /projects/mjolnir1/people/qvw641/WhoopingCrane/Annotation/WC_outgroups.singl.outmiss.vcf.gz 
