// RNA-seq Analysis Pipeline : custom functions 

def checkAlignmentPercent(prefix, logs) {
  //log.info "### CHECK ${prefix} - ${logs}"
  def percentAligned = 0;
  logs.eachLine { line ->
    if ((matcher = line =~ /([\d\.]+) \+ ([\d\.]+) mapped \(\s*/)) {
      //log.info "## L7 ## ${matcher[0]}"
      nbAligned = matcher[0][1]
    } else if ((matcher = line =~ /([\d\.]+) \+ ([\d\.]+) in total \s*/)) {
      //log.info "## L9 ## ${matcher[0]}"
      nbTotal = matcher[0][1]
    }
  }
  percentAligned = nbAligned.toFloat() / nbTotal.toFloat() * 100
  percentAligned = percentAligned.round(3)
  if(percentAligned.toFloat() <= '2'.toFloat()){
      log.info "### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($prefix)    >> ${percentAligned}% <<"
      return false
  }else if (percentAligned.toFloat() > '100'.toFloat()){
      log.info "Error -  aligner reads in ${logs} > 100% >>  ${percentAligned}% (${nbAligned} ${nbTotal}) <<"
      System.exit(-1)
  }else {
      log.info "### Passed alignment > ${prefix} >> ${percentAligned}% <<"
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