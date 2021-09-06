/* 
 * define the data analysis workflow 
 */

/* 
 * include requires tasks 
 */
include { markDuplicates } from '../process/markDuplicates'
include { dupradar } from '../process/dupradar'

workflow markdupFlow {
    // required inputs
    take:
      chBam
      chGtf
      chStrandedResults
    // workflow implementation
    main:
      // Duplicates
      markDuplicates(
        chBam
      )

      dupradar (
        markDuplicates.out.BamMd,
        chGtf.collect(),
        chStrandedResults
      )

    emit:
      chBamMd = chmarkDuplicates.out.BamMd
      chPicardResults = markDuplicates.out.PicardResults
      chPicardVersion = markDuplicates.out.PicardVersion
      chDupradarResults = dupradar.out.dupradarResults
}
