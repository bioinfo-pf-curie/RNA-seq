/* 
 * define the data analysis workflow 
 */

/* 
 * include requires tasks 
 */
include { identito } from '../process/identito'
include { combineIndentito } from '../process/combineIndentito'

workflow identitoFlow {
    // required inputs
    take:
      chFasta
      chFastaFai
      chPolymBed
      chBamMd
      chSplan
    // workflow implementation
    main:
      // Identito
      identito(
        chFasta,
        chFastaFai,
        chPolymBed,
        chBamMd
      )

      combineIndentito(
        identito.out.clustPolym.collect(),
        chSplan.collect()
      )

    emit:
      chbcftoolsIdentitoVersion = identito.out.version
      chClustPolymResults       = combineIndentito.out.clustPolymResults
}
