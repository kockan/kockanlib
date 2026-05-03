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
        outfile_fs_metrics = open("~{sample_id}.fs_metrics.tsv", 'w')

        for region in region_list:
            tokens = region.split('\t')
            chromosome = tokens[0]
            start = int(tokens[1])
            end = int(tokens[2])

            ctr_simplex = Counter()
            ctr_duplex = Counter()

            for read in infile_simplex.fetch(chromosome, start, end):
                if not read.is_unmapped and not read.is_duplicate and not read.is_secondary and not read.is_qcfail:
                    ctr_simplex[read.get_tag("cD")] += 1

            for read in infile_duplex.fetch(chromosome, start, end):
                if not read.is_unmapped and not read.is_duplicate and not read.is_secondary and not read.is_qcfail:
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

            num_simplex_reads = sum(ctr_simplex.values())
            num_duplex_reads = sum(ctr_duplex.values())
            num_simplex_reads_fs_g_1 = sum(v for k, v in ctr_simplex.items() if k > 1)
            num_simplex_reads_fs_geq_3 = sum(v for k, v in ctr_simplex.items() if k >= 3)
            num_simplex_reads_fs_geq_5 = sum(v for k, v in ctr_simplex.items() if k >= 5)
            num_simplex_reads_fs_geq_10 = sum(v for k, v in ctr_simplex.items() if k >= 10)

            mean_simplex_depth = 0
            total_depth = 0
            positions = 0

            for pileupcolumn in infile_simplex.pileup(chromosome, start, end, stepper = "all", truncate = False):
                total_depth += pileupcolumn.nsegments
                positions += 1

            if positions > 0:
                mean_simplex_depth = total_depth / positions
            else:
                mean_simplex_depth = 0

            outfile.write(chromosome + ":" + str(start) + "-" + str(end) + "\t" + str(count_filter_a) + "\t" + str(count_filter_b) + "\t" + str(count_filter_c) + "\t")
            outfile.write(str(count_filter_d) + "\t" + str(count_filter_e) + "\t" + str(count_filter_f) + "\t" + str(count_filter_g) + "\n")

            outfile_fs_metrics.write(chromosome + ":" + str(start) + "-" + str(end) + "\t" + str(num_simplex_reads) + "\t" + str(num_duplex_reads) + "\t")
            outfile_fs_metrics.write(str(num_simplex_reads_fs_g_1) + "\t" + str(num_simplex_reads_fs_geq_3) + "\t" + str(num_simplex_reads_fs_geq_5) + "\t" + str(num_simplex_reads_fs_geq_10) + "\t")
            outfile_fs_metrics.write(str(mean_simplex_depth) + "\n")

        outfile.close()
        outfile_fs_metrics.close()

        infile_simplex.close()
        infile_duplex.close()

        CODE
    >>>

    output {
        File custom_consensus_filter = "~{sample_id}.consensus_read_filter.tsv"
        File fs_metrics = "~{sample_id}.fs_metrics.tsv"
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
        File fs_metrics = CustomConsensusFilter.fs_metrics
    }
}