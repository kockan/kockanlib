version 1.0

task SamtoolsView {
    input {
        File bam
        File bai
        String contig

        Int cpu = 1
        Int memory_gb = 16
        Int disk_size_gb = 128
    }

    command <<<
        samtools view -c ~{bam} ~{contig} > output.txt
    >>>

    output {
        Int num_reads = read_int("output.txt")
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/tiwih@sha256:dbf6bc2b69c8b94d3b14451b60182fa4d1651f90781b7a60f9dd7cad2fe4c578"
    }
}

workflow SamtoolsNumReadsOnContig {
    input {
        File bam
        File bai
        String contig
    }

    call SamtoolsView {
        input:
            bam = bam,
            bai = bai,
            contig = contig
    }

    output {
        Int num_reads = SamtoolsView.num_reads
    }
}