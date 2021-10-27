/* 
 * define the data analysis workflow 
 */

/* 
 * include requires tasks 
 */
include { identitoPolym } from '../process/identitoPolym'
include { identitoCombine } from '../process/identitoCombine'

workflow identitoFlow {
    // required inputs
    take:
    bam
    fasta
    fai
    polymBed

    // workflow implementation
    main:
    chVersions = Channel.empty()

    identitoPolym(
      bam,
      fasta.collect(),
      fai.collect(),
      polymBed.collect()
    )
    chVersions = chVersions.mix(identitoPolym.out.versions)

    identitoCombine(
      identitoPolym.out.polyms.collect()
    )
    chVersions = chVersions.mix(identitoCombine.out.versions)

    emit:
    results = identitoCombine.out.results
    versions  = chVersions
}
