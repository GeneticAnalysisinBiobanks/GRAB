version 1.0

workflow grab_Region {
	
	input {

    	File objNullFile
    	File GenoFile
    	File GenoMarkerIndexFile
    	File GenoSampleFile
    	File GroupFile
    	File SparseGRMFile
    	File ControlFile
    	String MaxMAFVec
    	String annoVec
    	String OutputPrefix
    	
    }

    call Region {
    	input : objNullFile = objNullFile, GenoFile = GenoFile, GenoMarkerIndexFile = GenoMarkerIndexFile, GenoSampleFile = GenoSampleFile, GroupFile = GroupFile, SparseGRMFile = SparseGRMFile, ControlFile = ControlFile, MaxMAFVec = MaxMAFVec, annoVec = annoVec, OutputPrefix = OutputPrefix
    }

    output {

        File OutputFile = Region.OutputFile
        File indexFile = Region.indexFile
        File markerInfo = Region.markerInfo
        File otherMarkerInfo = Region.otherMarkerInfo

    }	

}

task Region {
	input {

    	File objNullFile
    	File GenoFile
    	File GenoMarkerIndexFile
    	File GenoSampleFile
    	File GroupFile
    	File SparseGRMFile
    	File ControlFile
    	String MaxMAFVec
    	String annoVec
    	String OutputPrefix

    }

    command <<<
    	set -euo pipefail
    	Rscript /grab/GRAB.Region.R  \
     		--objNullFile=~{objNullFile} \
     		--GenoFile=~{GenoFile} \
     		--GenoMarkerIndexFile=~{GenoMarkerIndexFile} \
     		--GenoSampleFile=~{GenoSampleFile} \
     		--GroupFile=~{GroupFile} \
    		--SparseGRMFile=~{SparseGRMFile} \
    		--ControlFile=~{ControlFile} \
                --MaxMAFVec=~{MaxMAFVec} \
                --annoVec=~{annoVec} \
    		--OutputPrefix=~{OutputPrefix}  		
    >>>

    output {

        File OutputFile = OutputPrefix
        File indexFile = OutputPrefix + ".index"
        File markerInfo = OutputPrefix + ".markerInfo"
        File otherMarkerInfo = OutputPrefix + ".otherMarkerInfo"

    }

    runtime {
        docker: "dx://UK_Biobank_WES:/docker_images/grab_0.0.3.3.tar.gz"
        dx_instance_type: "mem2_ssd1_v2_x2"
    }
}

