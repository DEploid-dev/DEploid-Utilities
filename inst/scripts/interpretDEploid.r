#!/usr/bin/env Rscript
library("DEploid.utils")
# DESCRIPTION:
#
# USAGE:
#    ./interpretDEploid.r -vcf FILE -plaf FILE -dEprefix STRING -o STRING
#    R --slave "--args -vcf FILE -plaf FILE -dEprefix STRING -o STRING " < utilities/interpretDEploid.r
#
# EXAMPLE:
#    ./interpretDEploid.r -vcf data/testData/PG0390-C.test.vcf -plaf data/testData/labStrains.test.PLAF.txt -dEprefix PG0390-CNopanel -o PG0390-CNopanel
#    R --slave "--args -vcf data/testData/PG0390-C.test.vcf -plaf data/testData/labStrains.test.PLAF.txt -dEprefix PG0390-CNopanel -o PG0390-CNopanel " < utilities/interpretDEploid.r
#    R --slave "--args -vcf data/testData/PG0390-C.test.vcf -plaf data/testData/labStrains.test.PLAF.txt -dEprefix PG0390-CNopanelExclude -o PG0390-CNopanelExclude -exclude data/testData/labStrains.test.exclude.txt " < utilities/interpretDEploid.r

if (!exists("dEploidRootDir")){
    print("dEploidRootDir undefined, try make dEploid again!")
}

args = (commandArgs(TRUE))
print(args)
myInput = DEploid.utils:::fun_parse ( args )

if (myInput$helpBool){
  DEploid.utils:::fun_print_help_interpret()
}

myCoverageInfo = DEploid.utils:::fun_extract_coverage ( myInput )

myPlafInfo = extractPLAF( myInput$plafFileName )

myExcludeInfo = DEploid.utils:::fun_extract_exclude (myInput$excludeFileName, myInput$excludeBool)

if (myInput$dEploid_v == "best") {
  DEploid.utils:::fun_interpretDEploid_best(myCoverageInfo, myPlafInfo, myInput$dEploidPrefix, myInput$outPrefix, myExcludeInfo, myInput$pdfBool)
  DEploid.utils:::fun_interpretDEploid_2 (myCoverageInfo, myInput$dEploidPrefix, myInput$outPrefix, myExcludeInfo, myInput$pdfBool, ringBool = FALSE, myInput$dEploid_v)
} else {
    if ( myInput$skip1Bool == FALSE ){
      DEploid.utils:::fun_interpretDEploid_1 (myCoverageInfo, myPlafInfo, myInput$dEploidPrefix, myInput$outPrefix, myExcludeInfo, myInput$pdfBool)
    }

  DEploid.utils:::fun_interpretDEploid_2 (myCoverageInfo, myInput$dEploidPrefix, myInput$outPrefix, myExcludeInfo, myInput$pdfBool)

  DEploid.utils:::fun_interpretDEploid_3 (myInput$dEploidPrefix, myInput$outPrefix, myInput$pdfBool, myInput$inbreedingBool)

  DEploid.utils:::fun_interpretDEploid_4 (myInput$dEploidPrefix, myInput$outPrefix, myInput$pdfBool)

    if (myInput$ringBool == TRUE){
      DEploid.utils:::fun_interpretDEploid_2 (myCoverageInfo, myInput$dEploidPrefix, myInput$outPrefix, myExcludeInfo, myInput$pdfBool, myInput$ringBool)
      DEploid.utils:::fun_interpretDEploid_3.ring (myInput$dEploidPrefix, myInput$outPrefix, myInput$pdfBool, myInput$inbreedingBool, myCoverageInfo, myExcludeInfo, myInput$ringDecreasingOrder, myInput$trackHeight, myInput$transformP)
    }

    if (myInput$ibdBool == TRUE){
        if ( myInput$skip1Bool == FALSE ){
          DEploid.utils:::fun_interpretDEploid_1 (myCoverageInfo, myPlafInfo, paste(myInput$dEploidPrefix, ".ibd", sep=""), paste(myInput$outPrefix, ".ibd", sep=""), myExcludeInfo, myInput$pdfBool)
        }

      DEploid.utils:::fun_interpretDEploid_2 (myCoverageInfo, paste(myInput$dEploidPrefix, ".ibd", sep=""), paste(myInput$outPrefix, ".ibd", sep=""), myExcludeInfo, myInput$pdfBool)

      DEploid.utils:::fun_interpretDEploid_3 (paste(myInput$dEploidPrefix, ".ibd", sep=""), paste(myInput$outPrefix, ".ibd", sep=""), myInput$pdfBool, myInput$inbreedingBool)

    }
}
