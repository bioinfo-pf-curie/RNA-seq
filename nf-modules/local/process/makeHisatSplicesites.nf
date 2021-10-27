process makeHisatSplicesites {
     label 'hisat2'
     label 'minCpu'
     label 'lowMem'
     publishDir "${params.outDir}/mapping", mode: 'copy',
       saveAs: { filename ->
         if (params.saveAlignedIntermediates) filename
         else null
       }

     input:
     path gtf

     output:
     path "${gtf.baseName}.hisat2SpliceSites.txt", emit: alignmentSplicesites
     path ("versions.txt"), emit: versions

     script:
     """
     echo \$(hisat2 --version | awk 'NR==1{print "hisat2 "\$3}') > versions.txt
     hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2SpliceSites.txt
     """
  }