/* 
 * Identito monitoring
 */

include { identitoPolym } from '../process/identitoPolym'
include { identitoCombine } from '../process/identitoCombine'

workflow identitoFlow {
    take:
    bam
    fasta
    fai
    polymBed

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
