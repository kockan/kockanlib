version 1.0

task CustomConsensusFilter {
    input {
        File simplex_bam
        File duplex_bam
        File simplex_bam_index
        File duplex_bam_index

        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = 512
    }

    command <<<
        set -e
        python3 <<CODE

        import sys
        import pysam
        from collections import Counter

        ctr_simplex = Counter()
        ctr_duplex = Counter()

        infile_simplex = pysam.AlignmentFile("~{simplex_bam}", "rb")
        infile_duplex = pysam.AlignmentFile("~{duplex_bam}", "rb")

        for read in infile_simplex.fetch("HPV16_Ref"):
            ctr_simplex[read.get_tag("cD")] += 1

        for read in infile_duplex.fetch("HPV16_Ref"):
            ctr_duplex[read.get_tag("cD")] += 1

        infile_simplex.close()
        infile_duplex.close()

        if len(ctr_duplex) == 0:
            count_filter_a = 0
        else:
            count_filter_a = sum(ctr_simplex.values()) - ctr_simplex[1]

        if sum(ctr_simplex.values()) == ctr_simplex[1]:
            count_filter_b = 0
        else:
            count_filter_b = sum(ctr_simplex.values()) - ctr_simplex[1]

        with open("custom_consensus_filter_a.txt", 'w') as f:
            f.write(str(count_filter_a))

        with open("custom_consensus_filter_b.txt", 'w') as f:
            f.write(str(count_filter_b))

        CODE
    >>>

    output {
        Int custom_consensus_filter_a_count = read_int("custom_consensus_filter_a.txt")
        Int custom_consensus_filter_b_count = read_int("custom_consensus_filter_b.txt")
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/simple_pysam@sha256:a302f9efe0bf1d4f9998ee1e9dda406223454ccaea0b5619046742221c1d2a74"
    }
}

workflow CustomConsensusFilter {
    input {
        File simplex_bam
        File duplex_bam
        File simplex_bam_index
        File duplex_bam_index
    }

    call CustomConsensusFilter {
        input:
            simplex_bam = simplex_bam,
            duplex_bam = duplex_bam,
            simplex_bam_index = simplex_bam_index,
            duplex_bam_index = duplex_bam_index
    }

    output {
        Int custom_consensus_filter_a_count = CustomConsensusFilter.custom_consensus_filter_a_count
        Int custom_consensus_filter_b_count = CustomConsensusFilter.custom_consensus_filter_b_count
    }
}