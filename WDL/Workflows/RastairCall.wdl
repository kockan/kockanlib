version 1.0

task RastairCallTask {
    input {
        File bam
        File bai
        File reference
        File reference_index
        String nOT
        String nOB
        Int min_mapping_quality
        Int min_base_quality

        Int cpu = 2
        Int memory_gb = 32
        Int disk_size_gb = 256
    }

    String output_basename = sub(basename(bam), "\\.bam$", "")

    command <<<
        rastair call \
        --nOT ~{nOT} \
        --nOB ~{nOB} \
        --min-mapq ~{min_mapping_quality} \
        --min-baseq ~{min_base_quality} \
        -- --cpgs-only \
        --fasta-file ~{reference} \
        ~{bam} > ~{output_basename}.tsv
    >>>

    output {
        File rastair_methylation_calls = "~{output_basename}.tsv"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/rastair@sha256:09bb60205344cd7a841675fb0ad7736eb35fa6bc0eb9d46f58f1aab5f8242f75"
    }
}

workflow RastairCall {
    input {
        File bam
        File bai
        File reference
        File reference_index
        String nOT = "0,0,5,0"
        String nOB = "0,0,5,0"
        Int min_mapping_quality = 1
        Int min_base_quality = 10
    }

    call RastairCallTask {
        input:
            bam = bam,
            bai = bai,
            reference = reference,
            reference_index = reference_index,
            nOT = nOT,
            nOB = nOB,
            min_mapping_quality = min_mapping_quality,
            min_base_quality = min_base_quality
    }

    output {
        File rastair_methylation_calls = RastairCallTask.rastair_methylation_calls
    }
}