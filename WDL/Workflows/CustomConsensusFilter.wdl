version 1.0

task CustomConsensusFilter {
    input {
        String sample_id
        File simplex_bam
        File duplex_bam
        File simplex_bam_index
        File duplex_bam_index
        File regions

        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = 512
    }

    command <<<
        set -e
        python3 <<CODE

        import pysam
        from collections import Counter

        infile_simplex = pysam.AlignmentFile("~{simplex_bam}", "rb")
        infile_duplex = pysam.AlignmentFile("~{duplex_bam}", "rb")

        region_list = []
        with open("~{regions}", 'r') as f:
            region_list = [line.strip() for line in f]

        outfile = open("~{sample_id}.consensus_read_filter.tsv", 'w')

        for region in region_list:
            ctr_simplex = Counter()
            ctr_duplex = Counter()

            for read in infile_simplex.fetch(region):
                ctr_simplex[read.get_tag("cD")] += 1

            for read in infile_duplex.fetch(region):
                ctr_duplex[read.get_tag("cD")] += 1

            if len(ctr_duplex) == 0:
                count_filter_a = 0
            else:
                count_filter_a = sum(ctr_simplex.values()) - ctr_simplex[1]

            if sum(ctr_simplex.values()) == ctr_simplex[1]:
                count_filter_b = 0
            else:
                count_filter_b = sum(ctr_simplex.values()) - ctr_simplex[1]

            if len(ctr_duplex) > 0:
                count_filter_c = sum(ctr_simplex.values())
                count_filter_d = sum(ctr_simplex.values())
                count_filter_e = sum(ctr_simplex.values())
                count_filter_f = sum(ctr_simplex.values())
                count_filter_g = sum(ctr_simplex.values())
            else:
                count_filter_c = sum(v for k, v in ctr_simplex.items() if k >= 3)
                count_filter_d = sum(v for k, v in ctr_simplex.items() if k >= 5)
                count_filter_e = sum(v for k, v in ctr_simplex.items() if k >= 10)
                count_filter_f = sum(v for k, v in ctr_simplex.items() if k >= 20)
                count_filter_g = sum(v for k, v in ctr_simplex.items() if k >= 50)

            outfile.write(region + "\t" + str(count_filter_a) + "\t" + str(count_filter_b) + "\t" + str(count_filter_c) + "\t")
            outfile.write(str(count_filter_d) + "\t" + str(count_filter_e) + "\t" + str(count_filter_f) + "\t" + str(count_filter_g) + "\n")

        outfile.close()

        infile_simplex.close()
        infile_duplex.close()

        CODE
    >>>

    output {
        File custom_consensus_filter = "~{sample_id}.consensus_read_filter.tsv"
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
        String sample_id
        File simplex_bam
        File duplex_bam
        File simplex_bam_index
        File duplex_bam_index
        File regions
    }

    call CustomConsensusFilter {
        input:
            sample_id = sample_id,
            simplex_bam = simplex_bam,
            duplex_bam = duplex_bam,
            simplex_bam_index = simplex_bam_index,
            duplex_bam_index = duplex_bam_index,
            regions = regions
    }

    output {
        File custom_consensus_filter = CustomConsensusFilter.custom_consensus_filter
    }
}