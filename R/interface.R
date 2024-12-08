fun_print_help <- function() {
  cat("    DEploid R utilities help\n")

  cat("
          Arguments:
               -help  --  Help. List the following content.
            -vcf STR  --  VCF file path.
            -ref STR  --  File path of reference allele count.
            -alt STR  --  File path of alternative allele count.
           -plaf STR  --  File path of population level allele frequencies.
        -exclude STR  --  File path of sites to be excluded.
              -o STR  --  Specify the file name prefix of the output.")
}

fun_print_help_interpret <- function() {
  fun_print_help()
  cat("
       -dEprefix STR  --  Specify DEploid output file prefix.\n")
  cat("
          Example:
          ./dEploid -vcf data/testData/PG0390-C.test.vcf \\
           -plaf data/testData/labStrains.test.PLAF.txt \\
           -o PG0390-CNopanel \\
           -noPanel
          ./utilities/dataExplore.r \\
          -vcf data/testData/PG0390-C.test.vcf.gz \\
          -plaf data/testData/labStrains.test.PLAF.txt \\
          -o PG0390-CNopanel \\
          -dEprefix PG0390-CNopanel\n\n")


  q(save = "no")
}

fun_print_help_explore <- function() {
  fun_print_help()
  cat("\n
          Example:
            ./utilities/dataExplore.r \\
             -vcf data/testData/PG0390-C.test.vcf.gz \\
             -plaf data/testData/labStrains.test.PLAF.txt \\
             -o PG0390-C\n\n")
  q(save = "no")
}

fun_parse <- function(args) {
  fun_local_checkAndIncreaseArgI <- function() {
    arg_i <- arg_i + 1
  }

  outPrefix <- "dataExplore"
  vcfFileName <- ""
  refFileName <- ""
  altFileName <- ""
  plafFileName <- ""
  excludeFileName <- ""
  dEploidPrefix <- ""
  excludeBool <- FALSE
  inbreedingBool <- FALSE
  ADFieldIndex <- 2
  pdfBool <- FALSE
  skip1Bool <- FALSE
  ibdBool <- FALSE
  helpBool <- FALSE
  ringBool <- FALSE
  ringDecreasingOrder <- TRUE
  transformP <- FALSE
  dEploid_v <- "classic"
  arg_i <- 1
  filter.window <- 10
  filter.threshold <- 0.995
  trackHeight <- 0.6
  while (arg_i <= length(args)) {
    argv <- args[arg_i]
    if (argv == "-vcf") {
      arg_i <- fun_local_checkAndIncreaseArgI()
      vcfFileName <- args[arg_i]
    } else if (argv == "-plaf") {
      arg_i <- fun_local_checkAndIncreaseArgI()
      plafFileName <- args[arg_i]
    } else if (argv == "-ref") {
      arg_i <- fun_local_checkAndIncreaseArgI()
      refFileName <- args[arg_i]
    } else if (argv == "-alt") {
      arg_i <- fun_local_checkAndIncreaseArgI()
      altFileName <- args[arg_i]
    } else if (argv == "-exclude") {
      arg_i <- fun_local_checkAndIncreaseArgI()
      excludeFileName <- args[arg_i]
      excludeBool <- TRUE
    } else if (argv == "-o") {
      arg_i <- fun_local_checkAndIncreaseArgI()
      outPrefix <- args[arg_i]
    } else if (argv == "-dEprefix") {
      arg_i <- fun_local_checkAndIncreaseArgI()
      dEploidPrefix <- args[arg_i]
    } else if (argv == "-inbreeding") {
      inbreedingBool <- TRUE
    } else if (argv == "-ADFieldIndex") {
      arg_i <- fun_local_checkAndIncreaseArgI()
      ADFieldIndex <- as.numeric(args[arg_i])
    } else if (argv == "-filter.threshold") {
      arg_i <- fun_local_checkAndIncreaseArgI()
      filter.threshold <- as.numeric(args[arg_i])
      # check range
      if (filter.threshold < 0 | filter.threshold > 1) {
        stop(paste("filter.threshold out of range [0, 1]:", filter.threshold))
      }
    } else if (argv == "-filter.window") {
      arg_i <- fun_local_checkAndIncreaseArgI()
      filter.window <- as.numeric(args[arg_i])
      # check range
      if (filter.window < 0) {
        stop(paste("filter.window cannot be negative:", filter.window))
      }
    } else if (argv == "-pdf") {
      pdfBool <- TRUE
    } else if (argv == "-skip1") {
      skip1Bool <- TRUE
    } else if (argv == "-ibd") {
      ibdBool <- TRUE
    } else if (argv == "-help") {
      helpBool <- TRUE
    } else if (argv == "-ring") {
      ringBool <- TRUE
    } else if (argv == "-reverseRing") {
      ringBool <- TRUE
      ringDecreasingOrder <- FALSE
    } else if (argv == "-best") {
      dEploid_v <- "best"
    } else if (argv == "-trackHeight") {
      arg_i <- fun_local_checkAndIncreaseArgI()
      trackHeight <- as.numeric(args[arg_i])
      # check range
      if (trackHeight < 0 | trackHeight > 1) {
        stop(paste("trackHeight out of range [0, 1]:", trackHeight))
      }
    } else if (argv == "-transformP") {
      transformP <- TRUE
    } else {
      stop(paste("Unknown flag:", argv))
    }

    arg_i <- arg_i + 1
  }

  if (length(args) == 0) {
    helpBool <- TRUE
  }

  #    if ( vcfFileName == "" || ( refFileName == "" && altFileName == "") ){
  #        stop ("Vcf File name not specified!")
  #    }
  #    cat ("vcfFileName: ", vcfFileName, "\n")
  #    cat ("plafFileName: ", plafFileName, "\n")

  return(list(
    vcfFileName = vcfFileName,
    refFileName = refFileName,
    altFileName = altFileName,
    plafFileName = plafFileName,
    outPrefix = outPrefix,
    dEploidPrefix = dEploidPrefix,
    excludeFileName = excludeFileName,
    excludeBool = excludeBool,
    ADFieldIndex = ADFieldIndex,
    pdfBool = pdfBool,
    inbreedingBool = inbreedingBool,
    skip1Bool = skip1Bool,
    ibdBool = ibdBool,
    helpBool = helpBool,
    filter.threshold = filter.threshold,
    filter.window = filter.window,
    ringBool = ringBool,
    ringDecreasingOrder = ringDecreasingOrder,
    trackHeight = trackHeight,
    transformP = transformP,
    dEploid_v = dEploid_v
  ))
}


fun_dEploidPrefix <- function(prefix, dEploid_v = "classic") {
  if (prefix == "") {
    stop("dEprefix ungiven!!!")
  }

  if (dEploid_v == "classic") {
    version_suffix <- "classic"
    return(list(
      propFileName = paste(prefix, ".", version_suffix, ".prop", sep = ""),
      hapFileName = paste(prefix, ".", version_suffix, ".hap", sep = ""),
      llkFileName = paste(prefix, ".", version_suffix, ".llk", sep = "")
      #                    dicLogFileName  = paste(prefix, "dic.log", sep = "")
    ))
  } else if (dEploid_v == "best") {
    return(list(
      propFileName.chooseK = paste(prefix, ".chooseK.prop", sep = ""),
      hapFileName.chooseK = paste(prefix, ".chooseK.hap", sep = ""),
      llkFileName.chooseK = paste(prefix, ".chooseK.llk", sep = ""),
      propFileName.ibd = paste(prefix, ".ibd.prop", sep = ""),
      hapFileName.ibd = paste(prefix, ".ibd.hap", sep = ""),
      llkFileName.ibd = paste(prefix, ".ibd.llk", sep = ""),
      hapFileName.final = paste(prefix, ".final.hap", sep = "")
      #                    dicLogFileName  = paste(prefix, "dic.log", sep = "")
    ))
  } else {
    return(NULL)
  }
}
