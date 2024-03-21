version 1.0

import "../Tasks/Downsample.wdl" as Downsample

workflow DownsampleStepsofTen {
    input {
        File inputCram
        File refFasta
        File refFastaIndex
        String strategy
    }

    call Downsample.Downsample as Downsample_0_9{
        input:
            inputCram = inputCram,
            refFasta = refFasta,
            refFastaIndex = refFastaIndex,
            probability = 0.9,
            strategy = strategy
    }

    output {
        File downsampledBam_0_9 = Downsample_0_9.downsampledBam
    }
}