#!/usr/bin/env Rscript
library("DEploid.utils")
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

myInput = DEploid.utils:::fun_parse ( args )

if (myInput$helpBool){
    fun.print.help.explore()
}

myCoverageInfo = DEploid.utils:::fun_extract_coverage ( myInput )

myPlafInfo = extractPLAF( myInput$plafFileName )

DEploid.utils:::fun_dataExplore (myCoverageInfo, myPlafInfo, myInput$outPrefix, myInput$pdfBool, myInput$filter.threshold, myInput$filter.window)
