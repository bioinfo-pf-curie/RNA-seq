/* 
 * define the data analysis workflow 
 */

/* 
 * include requires tasks 
 */
include { fastqc } from '../process/fastqc'

workflow qcFlow {
    // required inputs
    take:
      chRawReads 
    // workflow implementation
    main:
      fastqc(chRawReads)
    emit:
      chFastqcResults = fastqc.out.mqc       // channel: [ path *_fastqc.{zip,html} ]
      chFastqcVersion = fastqc.out.version   // channel: [ path v_fastqc.txt ]
}
