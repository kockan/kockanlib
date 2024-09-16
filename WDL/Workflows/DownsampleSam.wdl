version 1.0

import "../Tasks/Downsample.wdl" as Downsample

workflow DownsampleSam {
    input {
        File inputCram
        File refFasta
        File refFastaIndex
		Float probability
        String strategy
    }

    call Downsample.Downsample {
        input:
            inputCram = inputCram,
            refFasta = refFasta,
            refFastaIndex = refFastaIndex,
            probability = probability,
            strategy = strategy
    }

    output {
        File downsampledBam = Downsample.downsampledBam
    }
}
