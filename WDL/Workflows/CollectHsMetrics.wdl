version 1.0

task CollectHsMetricsTask {
    input {
        File bam
        File bai
        File dedup_bam
        File dedup_bai
        File reference
        File reference_fai
        File reference_dict
        File hpv_bait_interval_list
        File hpv_target_interval_list
        File hg38_bait_interval_list
        File hg38_target_interval_list
        String bait_set_name
        String output_basename

        Int? cpu = 2
        Int? memory_gb = 32
        Int? disk_size_gb = 512
    }

    command <<<
        gatk CollectHsMetrics \
        --BAIT_SET_NAME ~{bait_set_name} \
        --BAIT_INTERVALS ~{hpv_bait_interval_list} \
        --TARGET_INTERVALS ~{hpv_target_interval_list} \
        --INPUT ~{bam} \
        --OUTPUT ~{output_basename}.raw.hpv.hs_metrics.txt \
        --METRIC_ACCUMULATION_LEVEL ALL_READS \
        --REFERENCE_SEQUENCE ~{reference} \
        --COVERAGE_CAP 10000000 \
        --PER_TARGET_COVERAGE ~{output_basename}.raw.hpv.per_target_coverage.txt \
        --VALIDATION_STRINGENCY LENIENT

        gatk CollectHsMetrics \
        --BAIT_SET_NAME ~{bait_set_name} \
        --BAIT_INTERVALS ~{hpv_bait_interval_list} \
        --TARGET_INTERVALS ~{hpv_target_interval_list} \
        --INPUT ~{dedup_bam} \
        --OUTPUT ~{output_basename}.dedup.hpv.hs_metrics.txt \
        --METRIC_ACCUMULATION_LEVEL ALL_READS \
        --REFERENCE_SEQUENCE ~{reference} \
        --COVERAGE_CAP 10000000 \
        --PER_TARGET_COVERAGE ~{output_basename}.dedup.hpv.per_target_coverage.txt \
        --VALIDATION_STRINGENCY LENIENT

        gatk CollectHsMetrics \
        --BAIT_SET_NAME ~{bait_set_name} \
        --BAIT_INTERVALS ~{hg38_bait_interval_list} \
        --TARGET_INTERVALS ~{hg38_target_interval_list} \
        --INPUT ~{bam} \
        --OUTPUT ~{output_basename}.raw.hg38.hs_metrics.txt \
        --METRIC_ACCUMULATION_LEVEL ALL_READS \
        --REFERENCE_SEQUENCE ~{reference} \
        --COVERAGE_CAP 10000000 \
        --PER_TARGET_COVERAGE ~{output_basename}.raw.hg38.per_target_coverage.txt \
        --VALIDATION_STRINGENCY LENIENT

        gatk CollectHsMetrics \
        --BAIT_SET_NAME ~{bait_set_name} \
        --BAIT_INTERVALS ~{hg38_bait_interval_list} \
        --TARGET_INTERVALS ~{hg38_target_interval_list} \
        --INPUT ~{dedup_bam} \
        --OUTPUT ~{output_basename}.dedup.hg38.hs_metrics.txt \
        --METRIC_ACCUMULATION_LEVEL ALL_READS \
        --REFERENCE_SEQUENCE ~{reference} \
        --COVERAGE_CAP 10000000 \
        --PER_TARGET_COVERAGE ~{output_basename}.dedup.hg38.per_target_coverage.txt \
        --VALIDATION_STRINGENCY LENIENT
    >>>

    output {
        File raw_hpv_hs_metrics = "~{output_basename}.raw.hpv.hs_metrics.txt"
        File raw_hpv_per_target_coverage = "~{output_basename}.raw.hpv.per_target_coverage.txt"
        File raw_hg38_hs_metrics = "~{output_basename}.raw.hg38.hs_metrics.txt"
        File raw_hg38_per_target_coverage = "~{output_basename}.raw.hg38.per_target_coverage.txt"

        File dedup_hpv_hs_metrics = "~{output_basename}.dedup.hpv.hs_metrics.txt"
        File dedup_hpv_per_target_coverage = "~{output_basename}.dedup.hpv.per_target_coverage.txt"
        File dedup_hg38_hs_metrics = "~{output_basename}.dedup.hg38.hs_metrics.txt"
        File dedup_hg38_per_target_coverage = "~{output_basename}.dedup.hg38.per_target_coverage.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/hds@sha256:56f964695f08ddb74e3a29c63c3bc902334c1ddd735735cc98ba6d6a4212285c"
    }
}

workflow CollectHsMetrics {
    input {
        File bam
        File bai
        File dedup_bam
        File dedup_bai
        File reference
        File reference_fai
        File reference_dict
        File hpv_bait_interval_list
        File hpv_target_interval_list
        File hg38_bait_interval_list
        File hg38_target_interval_list
        String bait_set_name
        String output_basename
    }

    call CollectHsMetricsTask {
        input:
            bam = bam,
            bai = bai,
            dedup_bam = dedup_bam,
            dedup_bai = dedup_bai,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            hpv_bait_interval_list = hpv_bait_interval_list,
            hpv_target_interval_list = hpv_target_interval_list,
            hg38_bait_interval_list = hg38_bait_interval_list,
            hg38_target_interval_list = hg38_target_interval_list,
            bait_set_name = bait_set_name,
            output_basename = output_basename
    }

    output {
        File raw_hpv_hs_metrics = CollectHsMetricsTask.raw_hpv_hs_metrics
        File raw_hpv_per_target_coverage = CollectHsMetricsTask.raw_hpv_per_target_coverage
        File raw_hg38_hs_metrics = CollectHsMetricsTask.raw_hg38_hs_metrics
        File raw_hg38_per_target_coverage = CollectHsMetricsTask.raw_hg38_per_target_coverage

        File dedup_hpv_hs_metrics = CollectHsMetricsTask.dedup_hpv_hs_metrics
        File dedup_hpv_per_target_coverage = CollectHsMetricsTask.dedup_hpv_per_target_coverage
        File dedup_hg38_hs_metrics = CollectHsMetricsTask.dedup_hg38_hs_metrics
        File dedup_hg38_per_target_coverage = CollectHsMetricsTask.dedup_hg38_per_target_coverage
    }
}