import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test workflow
include { MINIMAP2_ALIGNMENT } from '/Users/ajv37/Documents/Pipelines/bactmap/subworkflows/local/minimap2_alignment/tests/../main.nf'

// define custom rules for JSON that will be generated.
def jsonOutput =
    new JsonGenerator.Options()
        .addConverter(Path) { value -> value.toAbsolutePath().toString() } // Custom converter for Path. Only filename
        .build()

def jsonWorkflowOutput = new JsonGenerator.Options().excludeNulls().build()

workflow {

    // run dependencies
    

    // workflow mapping
    def input = []
    
                input[0] = [
                    [ id:'test', single_end:false ], // meta map
                    file(params.modules_testdata_base_path + 'genomics/sarscov2/nanopore/fastq/test.fastq.gz', checkIfExists: true),
                ]
                input[1] = [
                    [ id:'genome' ],
                    file(params.modules_testdata_base_path + 'genomics/sarscov2/genome/genome.fasta', checkIfExists: true),
                ]
                
    //----

    //run workflow
    MINIMAP2_ALIGNMENT(*input)
    
    if (MINIMAP2_ALIGNMENT.output){

        // consumes all named output channels and stores items in a json file
        for (def name in MINIMAP2_ALIGNMENT.out.getNames()) {
            serializeChannel(name, MINIMAP2_ALIGNMENT.out.getProperty(name), jsonOutput)
        }	  
    
        // consumes all unnamed output channels and stores items in a json file
        def array = MINIMAP2_ALIGNMENT.out as Object[]
        for (def i = 0; i < array.length ; i++) {
            serializeChannel(i, array[i], jsonOutput)
        }    	

    }
}


def serializeChannel(name, channel, jsonOutput) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: {
            list.add(it)
        },
        onComplete: {
              def map = new HashMap()
              map[_name] = list
              def filename = "${params.nf_test_output}/output_${_name}.json"
              new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}


workflow.onComplete {

    def result = [
        success: workflow.success,
        exitStatus: workflow.exitStatus,
        errorMessage: workflow.errorMessage,
        errorReport: workflow.errorReport
    ]
    new File("${params.nf_test_output}/workflow.json").text = jsonWorkflowOutput.toJson(result)
    
}
