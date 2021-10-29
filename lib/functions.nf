// RNA-seq Analysis Pipeline : custom functions 

def checkAlignmentPercent(prefix, logs) {
  def percent_aligned = 0;
  logs.eachLine { line ->
    if ((matcher = line =~ /([\d\.]+) \+ ([\d\.]+) mapped \s*/)) {
      nb_aligned = matcher[0][1]

    } else if ((matcher = line =~ /([\d\.]+) \+ ([\d\.]+) in total \s*/)) {
      nb_total = matcher[0][1]
    }
  }
  percent_aligned = nb_aligned.toFloat() / nb_total.toFloat() * 100
  if(percent_aligned.toFloat() <= '2'.toFloat() ){
      log.info "#################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($prefix)    >> ${percent_aligned}% <<"
      skipped_poor_alignment << $prefix
      return false
  } else {
      log.info "          Passed alignment > ${prefix} >> ${percent_aligned}% <<"
      return true
  }
}
