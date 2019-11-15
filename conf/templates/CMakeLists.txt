

set(AP_NF_EXECUTOR "pbs" CACHE STRING "executor used by nextflow (e.g. pbs, slurm, etc.)")
set(AP_ANNOTATION_PATH "/to/replace" CACHE PATH "directory for the annotations used by the pipeline")
set(AP_SINGULARITY_IMAGE_PATH "/to/replace" CACHE PATH "directory for the singularity images used by the pipeline (when profile containers is used)")
set(AP_TOOLS_PATH "/to/replace" CACHE PATH "directory for the tools used by the pipeline (when profile toolsPath is used)")



if(${AP_USE_CONTAINERS})
    if(IS_DIRECTORY ${AP_SINGULARITY_IMAGE_PATH})
        message("OK - AP_SINGULARITY_IMAGE_PATH ${AP_SINGULARITY_IMAGE_PATH} exists")
        configure_file(containers.config.in containers.config @ONLY)
        install(FILES containers.config DESTINATION ${CMAKE_INSTALL_PREFIX}/conf)
    else()
        message(SEND_ERROR "AP_SINGULARITY_IMAGE_PATH ${AP_SINGULARITY_IMAGE_PATH} does not exist")
    endif()
endif()

if(${AP_USE_TOOLSPATH})
    if(IS_DIRECTORY ${AP_TOOLS_PATH})
        message("OK - AP_TOOLS_PATH ${AP_TOOLS_PATH} exists")
        configure_file(tools-path.config.in tools-path.config @ONLY)
        install(FILES tools-path.config DESTINATION ${CMAKE_INSTALL_PREFIX}/conf)
    else()
        message(SEND_ERROR "AP_TOOLS_PATH ${AP_TOOLS_PATH} does not exist")
    endif()
endif()

if(IS_DIRECTORY ${AP_ANNOTATION_PATH})
    message("OK - AP_ANNOTATION_PATH ${AP_ANNOTATION_PATH} exists")
    configure_file(igenomes.config.in igenomes.config @ONLY)
    install(FILES igenomes.config DESTINATION ${CMAKE_INSTALL_PREFIX}/conf)
else()
    message(SEND_ERROR "AP_ANNOTATION_PATH ${AP_ANNOTATION_PATH} does not exist")
endif()

configure_file(cluster.config.in cluster.config @ONLY)
install(FILES cluster.config DESTINATION ${CMAKE_INSTALL_PREFIX}/conf)
