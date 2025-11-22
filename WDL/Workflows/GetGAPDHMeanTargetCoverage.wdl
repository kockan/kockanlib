version 1.0

task GetGAPDHMeanTargetCoverageTask {
    input {
        File raw_per_target_hs_metrics

        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = 512
    }

    command <<<
        set -e
        python3 <<CODE

        import sys

        outfile_raw = open("raw_gapdh_mtc.txt", 'w')

        gapdh_targets_length_total = 0
        gapdh_coverage_total = 0.0

        with open("~{raw_per_target_hs_metrics}", 'r') as infile:
            header = infile.readline()

            for line in infile:
                line = line.rstrip()
                columns = line.split('\t')
                if columns[4] == "GAPDH":
                    # columns[3]: length (int)
                    # columns[4]: name (string)
                    # columns[6]: mean_coverage (float)
                    length = int(columns[3])
                    mean_coverage = float(columns[6])

                    gapdh_targets_length_total = gapdh_targets_length_total + length
                    gapdh_coverage_total = gapdh_coverage_total + (length * mean_coverage)

        gapdh_mtc = gapdh_coverage_total / gapdh_targets_length_total
        outfile_raw.write(str(gapdh_mtc))
        outfile_raw.close()
        CODE
    >>>

    output {
        Float raw_gapdh_mtc = read_float("raw_gapdh_mtc.txt")
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/simple_pysam@sha256:a302f9efe0bf1d4f9998ee1e9dda406223454ccaea0b5619046742221c1d2a74"
    }
}

workflow GetGAPDHMeanTargetCoverage {
    input {
        File raw_per_target_hs_metrics
    }

    call GetGAPDHMeanTargetCoverageTask {
        input:
            raw_per_target_hs_metrics = raw_per_target_hs_metrics
    }

    output {
        Float raw_gapdh_mtc = GetGAPDHMeanTargetCoverageTask.raw_gapdh_mtc
    }
}