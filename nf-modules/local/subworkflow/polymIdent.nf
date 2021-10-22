/* 
 * define the data analysis workflow 
 */

/* 
 * include requires tasks 
 */
include { getPolym } from '../process/getPolym'
include { combinePolym } from '../process/combinePolym'
include { identito } from '../process/identito'
include { combineIndentito } from '../process/combineIndentito'


workflow polymIdentitoFlow {
    // required inputs
    take:
      chFasta
      chFastaFai
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
      // Identito - monitoring
      identito(
        chFasta.collect(),
        chFastaFai.collect(),
        chPolymBed.collect(),
        chBamMd
      )
      combineIndentito(
        identito.out.clustPolym.collect()
      )

    emit:
      chBcftoolsVersion         = getPolym.out.version
      chClustPolymResults       = combinePolym.out.clustPolymResults
      chClustIdentitoResults    = combineIndentito.out.clustIdentitoResults
      chCombinePolymVersion     = combinePolym.out.version
      chbcftoolsIdentitoVersion = identito.out.version
}
