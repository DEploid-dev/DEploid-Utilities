#!/usr/bin/env Rscript
rm(list=ls()); library("DEploid.utils")
# DESCRIPTION:
#
# USAGE:
#    ./dataExplore.r -vcf FILE -plaf FILE -o STRING
#    R --slave "--args -vcf FILE -plaf FILE -o STRING " < dataExplore.r > dataExplore.rout
#
# EXAMPLE:
#    utilities/dataExplore.r -vcf data/testData/PG0390-C.test.vcf.gz -plaf data/testData/labStrains.test.PLAF.txt -o PG0390-C
#    R --slave "--args -vcf data/testData/PG0390-C.test.vcf.gz -plaf data/testData/labStrains.test.PLAF.txt -o PG0390-C " < utilities/dataExplore.r
#    R --slave "--args -ref data/testData/PG0390-C.test.ref -alt data/testData/PG0390-C.test.alt -plaf data/testData/labStrains.test.PLAF.txt -o PG0390-C " < utilities/dataExplore.r


args = (commandArgs(TRUE))

myInput = fun_parse ( args )

if (myInput$helpBool){
    fun.print.help.explore()
}

myCoverageInfo = fun_extract_coverage ( myInput )

myPlafInfo = extractPLAF( myInput$plafFileName )

fun_dataExplore (myCoverageInfo, myPlafInfo, myInput$outPrefix, myInput$pdfBool, myInput$filter.threshold, myInput$filter.window)
