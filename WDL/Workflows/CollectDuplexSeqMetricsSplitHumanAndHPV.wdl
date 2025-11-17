version 1.0

task CollectDuplexSeqMetricsSplitHumanAndHPVTask {
    input {
        File bam
        File hpv_intervals
        File hg38_intervals

        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = ceil((3 * size(bam, "GiB")) + 50)
    }

    String prefix = basename(bam, ".umi_grouped.bam")

    command <<<
        fgbio CollectDuplexSeqMetrics \
        --input ~{bam} \
        --output ~{prefix}.hpv \
        --intervals ~{hpv_intervals}

        fgbio CollectDuplexSeqMetrics \
        --input ~{bam} \
        --output ~{prefix}.hg38 \
        --intervals ~{hg38_intervals}
    >>>

    output {
        File hpv_family_sizes = "~{prefix}.hpv.family_sizes.txt"
        File hpv_duplex_family_sizes = "~{prefix}.hpv.duplex_family_sizes.txt"
        File hpv_duplex_yield_metrics = "~{prefix}.hpv.duplex_yield_metrics.txt"
        File hpv_umi_counts = "~{prefix}.hpv.umi_counts.txt"
        File hpv_duplex_qc = "~{prefix}.hpv.duplex_qc.pdf"

        File hg38_family_sizes = "~{prefix}.hg38.family_sizes.txt"
        File hg38_duplex_family_sizes = "~{prefix}.hg38.duplex_family_sizes.txt"
        File hg38_duplex_yield_metrics = "~{prefix}.hg38.duplex_yield_metrics.txt"
        File hg38_umi_counts = "~{prefix}.hg38.umi_counts.txt"
        File hg38_duplex_qc = "~{prefix}.hg38.duplex_qc.pdf"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/fgbio@sha256:b6869a0ae243d9f1b183e4a986fbe0853df2a56a1c6d7c0fec2965b6d8a7af1d"
    }
}

workflow CollectDuplexSeqMetricsSplitHumanAndHPV {
    input {
        File bam
        File hpv_intervals
        File hg38_intervals
    }

    call CollectDuplexSeqMetricsSplitHumanAndHPVTask {
        input:
            bam = bam,
            hpv_intervals = hpv_intervals,
            hg38_intervals = hg38_intervals
    }

    output {
        File hpv_family_sizes = CollectDuplexSeqMetricsSplitHumanAndHPVTask.hpv_family_sizes
        File hpv_duplex_family_sizes = CollectDuplexSeqMetricsSplitHumanAndHPVTask.hpv_duplex_family_sizes
        File hpv_duplex_yield_metrics = CollectDuplexSeqMetricsSplitHumanAndHPVTask.hpv_duplex_yield_metrics
        File hpv_umi_counts = CollectDuplexSeqMetricsSplitHumanAndHPVTask.hpv_umi_counts
        File hpv_duplex_qc = CollectDuplexSeqMetricsSplitHumanAndHPVTask.hpv_duplex_qc

        File hg38_family_sizes = CollectDuplexSeqMetricsSplitHumanAndHPVTask.hg38_family_sizes
        File hg38_duplex_family_sizes = CollectDuplexSeqMetricsSplitHumanAndHPVTask.hg38_duplex_family_sizes
        File hg38_duplex_yield_metrics = CollectDuplexSeqMetricsSplitHumanAndHPVTask.hg38_duplex_yield_metrics
        File hg38_umi_counts = CollectDuplexSeqMetricsSplitHumanAndHPVTask.hg38_umi_counts
        File hg38_duplex_qc = CollectDuplexSeqMetricsSplitHumanAndHPVTask.hg38_duplex_qc
    }
}