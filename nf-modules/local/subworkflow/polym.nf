/* 
 * define the data analysis workflow 
 */

/* 
 * include requires tasks 
 */
include { getPolym } from '../process/getPolym'
include { combinePolym } from '../process/combinePolym'

workflow polymFlow {
    // required inputs
    take:
      chFasta
      chPolymBed
      chBamMd
    // workflow implementation
    main:
      // Identito - polym
      getPolym(
        chFasta.collect(),
        chPolymBed.collect(),
        chBamMd
      )

      combinePolym(
        getPolym.out.clustPolym.collect()
      )

    emit:
      chBcftoolsVersion     = getPolym.out.version
      chClustPolymResults   = combinePolym.out.clustPolymResults
      chCombinePolymVersion = combinePolym.out.version
}
