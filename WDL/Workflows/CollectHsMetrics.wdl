version 1.0

task CollectHsMetricsTask {
    input {
        File bam
        File bai
        File reference
        File reference_fai
        File reference_dict
        File bait_interval_list
        File target_interval_list
        String bait_set_name
        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((3 * size(bam, "GiB")) + 50)
    }

    String prefix = basename(bam, ".sorted.bam")

    command <<<
        gatk CollectHsMetrics \
        --BAIT_SET_NAME ~{bait_set_name} \
        --BAIT_INTERVALS ~{bait_interval_list} \
        --TARGET_INTERVALS ~{target_interval_list} \
        --INPUT ~{bam} \
        --OUTPUT ~{prefix}.hs_metrics.txt \
        --METRIC_ACCUMULATION_LEVEL ALL_READS \
        --REFERENCE_SEQUENCE ~{reference} \
        --COVERAGE_CAP 100000 \
		--PER_TARGET_COVERAGE ~{prefix}.per_target_coverage.txt \
        --PER_BASE_COVERAGE ~{prefix}.per_base_coverage.txt \
        --VALIDATION_STRINGENCY LENIENT
    >>>

    output {
        File hs_metrics = "~{prefix}.hs_metrics.txt"
		File per_target_coverage = "~{prefix}.per_target_coverage.txt"
        File per_base_coverage = "~{prefix}.per_base_coverage.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

workflow CollectHsMetrics {
    input {
        File bam
        File bai
        File reference
        File reference_fai
        File reference_dict
        File bait_interval_list
        File target_interval_list
        String bait_set_name
    }

    call CollectHsMetricsTask {
        input:
			bam = bam,
			bai = bai,
			reference = reference,
			reference_fai = reference_fai,
			reference_dict = reference_dict,
			bait_interval_list = bait_interval_list,
			target_interval_list = target_interval_list,
			bait_set_name = bait_set_name
    }

    output {
        File hs_metrics = CollectHsMetricsTask.hs_metrics
		File per_target_coverage = CollectHsMetricsTask.per_target_coverage
        File per_base_coverage = CollectHsMetricsTask.per_base_coverage
    }
}
