version 1.0

task Downsample {
    input {
        File inputCram
        File refFasta
        File refFastaIndex
        Float probability
        String strategy
    }

    String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.27.5"
    Int preemptible = 0
    Int cpu = 16
    Int memoryMB = 32000
    Int javaMemorySize = memoryMB - 16000
    Int maxHeap = memoryMB - 8000
    Int diskSize = ceil(3.5 * size(inputCram, "GiB") + 256)

    String outputBasename = sub(sub(basename(inputCram), "\\.bam$", ""), "\\.cram$", "")

    command <<<
        java -Xms~{javaMemorySize}m -Xmx~{maxHeap}m -jar /usr/picard/picard.jar \
            DownsampleSam \
            -I ~{inputCram} \
            -O ~{outputBasename}.downsampled.bam \
            -STRATEGY ~{strategy} \
            -P ~{probability} \
            -CREATE_INDEX false \
            -REFERENCE_SEQUENCE ~{refFasta}
    >>>

    runtime {
        docker: docker
        preemptible: preemptible
        memory: "~{memoryMB} MiB"
        cpu: cpu
        disks: "local-disk " + diskSize + " HDD"
    }

    output {
        File downsampledBam = outputBasename + ".downsampled.bam"
    }
}