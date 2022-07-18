// RNA-seq Analysis Pipeline : custom functions 

def checkAlignmentPercent(prefix, logs) {
  def percentAligned = 0;
  logs.eachLine { line ->
    if ((matcher = line =~ /([\d\.]+) \+ ([\d\.]+) mapped \s*/)) {
      nbAligned = matcher[0][1]

    } else if ((matcher = line =~ /([\d\.]+) \+ ([\d\.]+) in total \s*/)) {
      nbTotal = matcher[0][1]
    }
  }
  percentAligned = nbAligned.toFloat() / nbTotal.toFloat() * 100
  if(percentAligned.toFloat() <= '2'.toFloat() ){
      log.info "#################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($prefix)    >> ${percentAligned}% <<"
      //skippedPoorAlignment << $prefix
      return false
  } else {
      log.info "          Passed alignment > ${prefix} >> ${percentAligned}% <<"
      return true
  }
}



def combineStrandness(input, strandness){

  def idx = input.toList().size()

  input
    .combine(strandness)
    .filter{it[0].id == it[2]}
    .map{ it ->
      meta = it[0]
      meta.strandness = it[3]
      return [meta, it[1]]
    }

}