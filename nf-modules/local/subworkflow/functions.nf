// From nf-core
// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false
skippedPoorAlignment = []
def checkStarLog(logs) {
     def percentAligned = 0;
     logs.eachLine { line ->
       if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
         percentAligned = matcher[0][1]
       }else if ((matcher = line =~ /Uniquely mapped reads number\s*\|\s*([\d\.]+)/)) {
         numAligned = matcher[0][1]
       }
     }
     logname = logs.getBaseName() - 'Log.final'
     if(percentAligned.toFloat() <= '2'.toFloat() || numAligned.toInteger() <= 1000.toInteger() ){
         log.info "#################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percentAligned}% <<"
         skippedPoorAlignment << logname
         return false
     } else {
         log.info "          Passed alignment > star ($logname)   >> ${percentAligned}% <<"
         return true
     }
}
