#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/bactmap
========================================================================================
    Github : https://github.com/nf-core/bactmap
    Website: https://nf-co.re/bactmap
    Slack  : https://nfcore.slack.com/channels/bactmap
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

workflow NFCORE_BACTMAP {

    //
    // WORKFLOW: Run main nf-core/bactmap analysis pipeline
    //
    include { BACTMAP } from './workflows/bactmap'
    BACTMAP ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_BACTMAP ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
