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
        outfile_hpv_mean_depths_filtered = open("~{sample_id}.hpv_mean_depths_filtered.tsv", 'w')

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

            if chromosome.startswith("HPV"):
                total_depth_fs_g_1 = 0
                total_depth_fs_geq_3 = 0
                total_depth_fs_geq_5 = 0
                total_depth_fs_geq_10 = 0
                num_positions = 0

                for pileupcolumn in infile_simplex.pileup(chromosome, start, end, stepper = "all", truncate = False, max_depth = 1000000, ignore_overlaps = True):
                    count_fs_g_1 = 0
                    count_fs_geq_3 = 0
                    count_fs_geq_5 = 0
                    count_fs_geq_10 = 0

                    for pileupread in pileupcolumn.pileups:
                        if pileupread.alignment.get_tag("cD") > 1:
                            count_fs_g_1 += 1
                        if pileupread.alignment.get_tag("cD") >= 3:
                            count_fs_geq_3 += 1
                        if pileupread.alignment.get_tag("cD") >= 5:
                            count_fs_geq_5 += 1
                        if pileupread.alignment.get_tag("cD") >= 10:
                            count_fs_geq_10 += 1

                    total_depth_fs_g_1 += count_fs_g_1
                    total_depth_fs_geq_3 += count_fs_geq_3
                    total_depth_fs_geq_5 += count_fs_geq_5
                    total_depth_fs_geq_10 += count_fs_geq_10
                    num_positions += 1

                mean_simplex_depth_fs_g_1 = 0
                mean_simplex_depth_fs_geq_3 = 0
                mean_simplex_depth_fs_geq_5 = 0
                mean_simplex_depth_fs_geq_10 = 0
                if num_positions > 0:
                    mean_simplex_depth_fs_g_1 = total_depth_fs_g_1 / num_positions
                    mean_simplex_depth_fs_geq_3 = total_depth_fs_geq_3 / num_positions
                    mean_simplex_depth_fs_geq_5 = total_depth_fs_geq_5 / num_positions
                    mean_simplex_depth_fs_geq_10 = total_depth_fs_geq_10 / num_positions
                outfile_hpv_mean_depths_filtered.write(chromosome + "\t" + str(mean_simplex_depth_fs_g_1) + "\t" + str(mean_simplex_depth_fs_geq_3) + "\t" + str(mean_simplex_depth_fs_geq_5) + "\t" + str(mean_simplex_depth_fs_geq_10) + "\n")

            outfile.write(chromosome + ":" + str(start) + "-" + str(end) + "\t" + str(count_filter_a) + "\t" + str(count_filter_b) + "\t" + str(count_filter_c) + "\t")
            outfile.write(str(count_filter_d) + "\t" + str(count_filter_e) + "\t" + str(count_filter_f) + "\t" + str(count_filter_g) + "\n")

            outfile_fs_metrics.write(chromosome + ":" + str(start) + "-" + str(end) + "\t" + str(num_simplex_reads) + "\t" + str(num_duplex_reads) + "\t")
            outfile_fs_metrics.write(str(num_simplex_reads_fs_g_1) + "\t" + str(num_simplex_reads_fs_geq_3) + "\t" + str(num_simplex_reads_fs_geq_5) + "\t" + str(num_simplex_reads_fs_geq_10) + "\t")
            outfile_fs_metrics.write(str(mean_simplex_depth) + "\n")

        outfile.close()
        outfile_fs_metrics.close()
        outfile_hpv_mean_depths_filtered.close()

        infile_simplex.close()
        infile_duplex.close()

        CODE
    >>>

    output {
        File custom_consensus_filter = "~{sample_id}.consensus_read_filter.tsv"
        File fs_metrics = "~{sample_id}.fs_metrics.tsv"
        File hpv_mean_depths_filtered = "~{sample_id}.hpv_mean_depths_filtered.tsv"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/simple_pysam@sha256:a302f9efe0bf1d4f9998ee1e9dda406223454ccaea0b5619046742221c1d2a74"
    }
}

task SummarizeStats {
    input {
        String sample_id
        String top_hpv_contig
        File fs_stats
        File gapdh_regions
        File fp_regions
        File hpv_mean_depths_filtered

        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = 512
    }

    command <<<
        set -e
        python3 <<CODE

        gapdh_region_list = []
        with open("~{gapdh_regions}", 'r') as f:
            gapdh_region_list = [line.strip() for line in f]

        fp_region_list = []
        with open("~{fp_regions}", 'r') as f:
            fp_region_list = [line.strip() for line in f]

        hpv_simplex_reads = 0
        hpv_duplex_reads = 0
        hpv_simplex_reads_fs_g_1 = 0
        hpv_simplex_reads_fs_geq_3 = 0
        hpv_simplex_reads_fs_geq_5 = 0
        hpv_simplex_reads_fs_geq_10 = 0
        gapdh_simplex_reads = 0
        gapdh_duplex_reads = 0
        hg38_simplex_reads = 0
        hg38_duplex_reads = 0
        hg38_simplex_reads_fp_only = 0
        hg38_duplex_reads_fp_only = 0
        mean_simplex_depth_hg38_fp_only = 0.0
        mean_simplex_depth_hpv_fs_g_1 = 0.0
        mean_simplex_depth_hpv_fs_geq_3 = 0.0
        mean_simplex_depth_hpv_fs_geq_5 = 0.0
        mean_simplex_depth_hpv_fs_geq_10 = 0.0

        with open("~{hpv_mean_depths_filtered}", 'r') as f:
            for line in f:
                line = line.rstrip()
                columns = line.split('\t')

                region = columns[0]
                if region.startswith("~{top_hpv_contig}"):
                    mean_simplex_depth_hpv_fs_g_1 = float(columns[1])
                    mean_simplex_depth_hpv_fs_geq_3 = float(columns[2])
                    mean_simplex_depth_hpv_fs_geq_5 = float(columns[3])
                    mean_simplex_depth_hpv_fs_geq_10 = float(columns[4])

        with open("~{fs_stats}", 'r') as f:
            for line in f:
                line = line.rstrip()
                columns = line.split('\t')

                region = columns[0]
                if region.startswith("~{top_hpv_contig}"):
                    hpv_simplex_reads = int(columns[1])
                    hpv_duplex_reads = int(columns[2])
                    hpv_simplex_reads_fs_g_1 = int(columns[3])
                    hpv_simplex_reads_fs_geq_3 = int(columns[4])
                    hpv_simplex_reads_fs_geq_5 = int(columns[5])
                    hpv_simplex_reads_fs_geq_10 = int(columns[6])

                if region in gapdh_region_list:
                    gapdh_simplex_reads = gapdh_simplex_reads + int(columns[1])
                    gapdh_duplex_reads = gapdh_duplex_reads + int(columns[2])

                if region in fp_region_list:
                    hg38_simplex_reads_fp_only = hg38_simplex_reads_fp_only + int(columns[1])
                    hg38_duplex_reads_fp_only = hg38_duplex_reads_fp_only + int(columns[2])
                    mean_simplex_depth_hg38_fp_only = mean_simplex_depth_hg38_fp_only + float(columns[7])

                if not region.startswith("HPV"):
                    hg38_simplex_reads = hg38_simplex_reads + int(columns[1])
                    hg38_duplex_reads = hg38_duplex_reads + int(columns[2])

        mean_simplex_depth_hg38_fp_only = mean_simplex_depth_hg38_fp_only / len(fp_region_list)

        outfile = open("~{sample_id}.fs_stats_summary.tsv", 'w')

        outfile.write("sample_id" + "\t" + "hpv_simplex_reads" + "\t" + "hpv_duplex_reads" + "\t" + "hpv_simplex_reads_fs_g_1" + "\t")
        outfile.write("hpv_simplex_reads_fs_geq_3" + "\t" + "hpv_simplex_reads_fs_geq_5" + "\t" + "hpv_simplex_reads_fs_geq_10" + "\t")
        outfile.write("hg38_simplex_reads" + "\t" + "hg38_duplex_reads" + "\t" + "hg38_simplex_reads_fp_only" + "\t" + "hg38_duplex_reads_fp_only" + "\t")
        outfile.write("mean_simplex_depth_hg38_fp_only" + "\t" + "gapdh_simplex_reads" + "\t" + "gapdh_duplex_reads" + "\t")
        outfile.write("mean_simplex_depth_hpv_fs_g_1" + "\t" + "mean_simplex_depth_hpv_fs_geq_3" + "\t" + "mean_simplex_depth_hpv_fs_geq_5" + "\t" + "mean_simplex_depth_hpv_fs_geq_10" + "\n")

        outfile.write("~{sample_id}" + "\t" + str(hpv_simplex_reads) + "\t" + str(hpv_duplex_reads) + "\t" + str(hpv_simplex_reads_fs_g_1) + "\t")
        outfile.write(str(hpv_simplex_reads_fs_geq_3) + "\t" + str(hpv_simplex_reads_fs_geq_5) + "\t" + str(hpv_simplex_reads_fs_geq_10) + "\t")
        outfile.write(str(hg38_simplex_reads) + "\t" + str(hg38_duplex_reads) + "\t")
        outfile.write(str(hg38_simplex_reads_fp_only) + "\t" + str(hg38_duplex_reads_fp_only) + "\t" + str(mean_simplex_depth_hg38_fp_only) + "\t")
        outfile.write(str(gapdh_simplex_reads) + "\t" + str(gapdh_duplex_reads) + "\t")
        outfile.write(str(mean_simplex_depth_hpv_fs_g_1) + "\t" + str(mean_simplex_depth_hpv_fs_geq_3) + "\t" + str(mean_simplex_depth_hpv_fs_geq_5) + "\t" + str(mean_simplex_depth_hpv_fs_geq_10) + "\n")

        outfile.close()

        CODE
    >>>

    output {
        File fs_stats_summary = "~{sample_id}.fs_stats_summary.tsv"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/simple_pysam@sha256:a302f9efe0bf1d4f9998ee1e9dda406223454ccaea0b5619046742221c1d2a74"
    }
}

task GetCoverageStatsFiltered {
    input {
        String sample_id
        String top_hpv_contig
        File simplex_bam
        File simplex_bam_index

        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = 512
    }

    command <<<
        samtools view -h -e '[cD] > 1' ~{simplex_bam} ~{top_hpv_contig} | samtools coverage - > ~{sample_id}.fs_g_1.coverage.txt
        samtools view -h -e '[cD] >= 3' ~{simplex_bam} ~{top_hpv_contig} | samtools coverage - > ~{sample_id}.fs_geq_3.coverage.txt
        samtools view -h -e '[cD] >= 5' ~{simplex_bam} ~{top_hpv_contig} | samtools coverage - > ~{sample_id}.fs_geq_5.coverage.txt
        samtools view -h -e '[cD] >= 10' ~{simplex_bam} ~{top_hpv_contig} | samtools coverage - > ~{sample_id}.fs_geq_10.coverage.txt
    >>>

    output {
        File fs_g_1_coverage = "~{sample_id}.fs_g_1.coverage.txt"
        File fs_geq_3_coverage = "~{sample_id}.fs_geq_3.coverage.txt"
        File fs_geq_5_coverage = "~{sample_id}.fs_geq_5.coverage.txt"
        File fs_geq_10_coverage = "~{sample_id}.fs_geq_10.coverage.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/samtools@sha256:d92861a76f16d9f4ee2d89ee52af46d7beb27764675682a410b3cd847f3d5a17"
    }
}
workflow CustomConsensusFilter {
    input {
        String sample_id
        String top_hpv_contig
        File simplex_bam
        File duplex_bam
        File simplex_bam_index
        File duplex_bam_index
        File regions
        File gapdh_regions
        File fp_regions
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

    call SummarizeStats {
        input:
            sample_id = sample_id,
            top_hpv_contig = top_hpv_contig,
            fs_stats = CustomConsensusFilter.fs_metrics,
            gapdh_regions = gapdh_regions,
            fp_regions = fp_regions,
            hpv_mean_depths_filtered = CustomConsensusFilter.hpv_mean_depths_filtered
    }

    call GetCoverageStatsFiltered {
        input:
            sample_id = sample_id,
            top_hpv_contig = top_hpv_contig,
            simplex_bam = simplex_bam,
            simplex_bam_index = simplex_bam_index
    }

    output {
        File custom_consensus_filter = CustomConsensusFilter.custom_consensus_filter
        File fs_metrics = CustomConsensusFilter.fs_metrics
        File fs_stats_summary = SummarizeStats.fs_stats_summary
        File fs_g_1_coverage = GetCoverageStatsFiltered.fs_g_1_coverage
        File fs_geq_3_coverage = GetCoverageStatsFiltered.fs_geq_3_coverage
        File fs_geq_5_coverage = GetCoverageStatsFiltered.fs_geq_5_coverage
        File fs_geq_10_coverage = GetCoverageStatsFiltered.fs_geq_10_coverage
    }
}