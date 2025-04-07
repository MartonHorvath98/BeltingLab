include { GET } from '../modules/get/main.nf'
include { ARCHIVE_EXTRACT} from '../subworkflows/archive_extract/main.nf'

workflow {
    // STEP 1: Download the CTAT archive    
    def ctat_exists = file(params.ctat_lib).exists()
    if (!ctat_exists) {
        ctat_ch = GET(params.ctat_url)
        ctat_ch = ARCHIVE_EXTRACT(ctat_ch.out.archive)
    } else {
        log.info("CTAT archive found in the specified path: ${params.ctat_lib}.")
    }

    //STEP 2: Download the 3DID database
    did_exists = file(params.did_lib).exists()
    if (!did_exists) {
        did_ch = GET(params.did_url)
        did_ch = ARCHIVE_EXTRACT(did_ch.out.archive)
    } else {
        log.info("3DID archive found in the specified path: ${params.did_lib}.")
    }


}