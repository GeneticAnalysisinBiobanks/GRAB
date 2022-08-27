version 1.0

workflow grab_ReadGeno {
	
	input {

    	File GenoFile
    	File GenoMarkerIndexFile
    	File GenoSampleFile
    	File IDsToIncludeFile
    	String OutputPrefix
    }

    call ReadGeno {
    	input : GenoFile = GenoFile, GenoMarkerIndexFile = GenoMarkerIndexFile, GenoSampleFile = GenoSampleFile, IDsToIncludeFile = IDsToIncludeFile, OutputPrefix = OutputPrefix
    }

    output {

        File GenoMatFile = ReadGeno.GenoMatFile
        File markerInfoFile = ReadGeno.markerInfoFile

    }	

}

task ReadGeno {
	input {

    	File GenoFile
    	File GenoMarkerIndexFile
    	File GenoSampleFile
    	File IDsToIncludeFile
    	String OutputPrefix
    }

    command <<<
    	set -euo pipefail
    	Rscript /grab/GRAB.ReadGeno.New.R  \
     		--GenoFile=~{GenoFile} \
     		--GenoMarkerIndexFile=~{GenoMarkerIndexFile} \
     		--GenoSampleFile=~{GenoSampleFile} \
    		--IDsToIncludeFile=~{IDsToIncludeFile} \
    		--OutputPrefix=~{OutputPrefix}  		
    >>>

    output {

        File GenoMatFile = OutputPrefix + ".GenoMat.txt"
        File markerInfoFile = OutputPrefix + ".markerInfo.txt"

    }

    runtime {
        docker: "dx://UK_Biobank_WES:/docker_images/grab_0.0.3.3.tar.gz"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}

