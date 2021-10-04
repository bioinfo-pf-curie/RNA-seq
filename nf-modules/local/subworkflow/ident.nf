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
      polymsCh
      mdBamPolymCh
      //ou 
      //chPolymBed
      //chBamMd
    // workflow implementation
    main:
      // Identito
      identito(
        chFasta,
        chFastaFai,
        chPolymBed.collect(),
        chBamMd
      )

      combineIndentito(
        identito.out.clustPolym.collect()
      )

    emit:
      chbcftoolsIdentitoVersion = identito.out.version
      chClustPolymResults       = combineIndentito.out.clustPolymResults
}
