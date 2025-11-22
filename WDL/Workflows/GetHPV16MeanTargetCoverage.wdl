version 1.0

task GetHPV16MeanTargetCoverageTask {
    input {
        File raw_per_target_hs_metrics
        File dedup_per_target_hs_metrics

        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = 512
    }

    command <<<
        set -e
        python3 <<CODE

        import sys

        outfile_raw = open("raw_hpv16_mtc.txt", 'w')
        outfile_dedup = open("dedup_hpv16_mtc.txt", 'w')

        with open("~{raw_per_target_hs_metrics}", 'r') as infile:
            header = infile.readline()

            for line in infile:
                line = line.rstrip()
                columns = line.split('\t')
                if columns[0] == "HPV16_Ref":
                    # columns[6] is mean_coverage for that target
                    outfile_raw.write(columns[6])

        with open("~{dedup_per_target_hs_metrics}", 'r') as infile:
            header = infile.readline()

            for line in infile:
                line = line.rstrip()
                columns = line.split('\t')
                if columns[0] == "HPV16_Ref":
                    # columns[6] is mean_coverage for that target
                    outfile_dedup.write(columns[6])

        outfile_raw.close()
        outfile_dedup.close()
        CODE
    >>>

    output {
        Float raw_hpv16_mtc = read_float("raw_hpv16_mtc.txt")
        Float dedup_hpv16_mtc = read_float("dedup_hpv16_mtc.txt")
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/simple_pysam@sha256:a302f9efe0bf1d4f9998ee1e9dda406223454ccaea0b5619046742221c1d2a74"
    }
}

workflow GetHPV16MeanTargetCoverage {
    input {
        File raw_per_target_hs_metrics
        File dedup_per_target_hs_metrics
    }

    call GetHPV16MeanTargetCoverageTask {
        input:
            raw_per_target_hs_metrics = raw_per_target_hs_metrics,
            dedup_per_target_hs_metrics = dedup_per_target_hs_metrics
    }

    output {
        Float raw_hpv16_mtc = GetHPV16MeanTargetCoverageTask.raw_hpv16_mtc
        Float dedup_hpv16_mtc = GetHPV16MeanTargetCoverageTask.dedup_hpv16_mtc
    }
}