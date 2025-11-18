version 1.0

task SamtoolsBedcovTask {
    input {
        File bam
        File bai
        File dedup_bam
        File dedup_bai
        File bed
        String output_basename
        Int min_mapping_quality = 1

        Int? cpu = 2
        Int? memory_gb = 32
        Int? disk_size_gb = 512
    }

    command <<<
        samtools bedcov -Q ~{min_mapping_quality} -c ~{bed} ~{bam} > ~{output_basename}.raw.bedcov.tsv
        samtools bedcov -Q ~{min_mapping_quality} -c ~{bed} ~{dedup_bam} > ~{output_basename}.dedup.bedcov.tsv
    >>>

    output {
        File raw_bedcov = "~{output_basename}.raw.bedcov.tsv"
        File dedup_bedcov = "~{output_basename}.dedup.bedcov.tsv"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/tiwih@sha256:dbf6bc2b69c8b94d3b14451b60182fa4d1651f90781b7a60f9dd7cad2fe4c578"
    }
}

workflow SamtoolsBedcov {
    input {
        File bam
        File bai
        File dedup_bam
        File dedup_bai
        File bed
        String output_basename
    }

    call SamtoolsBedcovTask {
        input:
            bam = bam,
            bai = bai,
            dedup_bam = dedup_bam,
            dedup_bai = dedup_bai,
            bed = bed,
            output_basename = output_basename
    }

    output {
        File raw_bedcov = SamtoolsBedcovTask.raw_bedcov
        File dedup_bedcov = SamtoolsBedcovTask.dedup_bedcov
    }
}