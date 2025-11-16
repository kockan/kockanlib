version 1.0

task DragenCytosineReportToBismarkTask {
    input {
        File cytosine_report

        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = 512
    }

    String prefix = basename(cytosine_report, ".CX_report.txt.gz")

    command <<<
        set -e
        python3 <<CODE

        import sys
        import gzip

        # The Bismark output format has the following tab-delimited columns (but no actual header)
        # 0: Chromosome
        # 1: Start position
        # 2: End position
        # 3: Methylation Percentage
        # 4: Count C's Methylated
        # 5: Count C's Unmethylated
        outfile = open("~{prefix}.5base.bismark_converted.cov", 'w')

        primary_assembly = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]

        with gzip.open("~{cytosine_report}", 'rt') as infile:
            for line in infile:
                line = line.rstrip()

                # Illumina DRAGEN cytosine reports are tab delimited with no header
                #
                # File format:
                #
                # columns[0]: chromosome
                # columns[1]: position
                # columns[2]: strand
                # columns[3]: number of reads covering site, supporting methylated C's
                # columns[4]: number of reads covering site, supporting unmethylated C's
                # columns[5]: methylation context
                # columns[6]: trinucleotide content
                columns = line.split('\t')

                # Skip if not in primary assembly
                if columns[0] not in primary_assembly:
                    continue

                # We only care about CpG context
                if columns[5] != "CG":
                    continue

                chrom = columns[0]
                pos = columns[1]
                mcs = int(columns[3])
                umcs = int(columns[4])

                percent_methylated = 0.00

                coverage = mcs + umcs
                if coverage > 0:
                    percent_methylated = mcs * 100 / coverage

                outfile.write(chrom + "\t" + pos + "\t" + pos + "\t" + str(percent_methylated) + "\t" + str(mcs) + "\t" + str(umcs) + "\n")
        outfile.close()
        CODE
    >>>

    output {
        File bismark_cov = "~{prefix}.5base.bismark_converted.cov"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/methylation_misc@sha256:e33b932cd0adb0d7dacc4fa8c8134378e3e871ec2704e0574e462dacdc464d60"
    }
}

task BismarkCovToBedgraph {
    input {
        File bismark_cov

        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = 512
    }

    String prefix = basename(bismark_cov, ".5base.bismark_converted.cov")

    command <<<
        set -e
        python3 <<CODE

        import sys

        # The Bismark output format has the following tab-delimited columns (but no actual header)
        # 0: Chromosome
        # 1: Start position
        # 2: End position
        # 3: Methylation Percentage
        # 4: Count C's Methylated
        # 5: Count C's Unmethylated
        outfile = open("~{prefix}.5base.bismark_converted.bedgraph", 'w')

        with open("~{bismark_cov}", 'r') as infile:
            for line in infile:
                line = line.rstrip()

                columns = line.split('\t')

                chrom = columns[0]
                pos = int(columns[1])
                percent_methylated = float(columns[3])

                outfile.write(chrom + "\t" + str(pos - 1) + "\t" + str(pos) + "\t" + str(percent_methylated) + "\n")
        outfile.close()
        CODE
    >>>

    output {
        File bedgraph = "~{prefix}.5base.bismark_converted.bedgraph"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/methylation_misc@sha256:e33b932cd0adb0d7dacc4fa8c8134378e3e871ec2704e0574e462dacdc464d60"
    }
}

task BedgraphToBigwig {
    input {
        File bedgraph
        File chrom_sizes

        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = 512
    }

    String prefix = basename(bedgraph, ".5base.bismark_converted.bedgraph")

    command <<<
        bedGraphToBigWig ~{bedgraph} ~{chrom_sizes} ~{prefix}.5base.bismark_converted.bw
    >>>

    output {
        File bigwig = "~{prefix}.5base.bismark_converted.bw"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/methylation_misc@sha256:e33b932cd0adb0d7dacc4fa8c8134378e3e871ec2704e0574e462dacdc464d60"
    }
}

workflow DragenCytosineReportToBismark {
    input {
        File cytosine_report
        File chrom_sizes
    }

    call DragenCytosineReportToBismarkTask {
        input:
            cytosine_report = cytosine_report
    }

    call BismarkCovToBedgraph {
        input:
            bismark_cov = DragenCytosineReportToBismarkTask.bismark_cov
    }

    call BedgraphToBigwig {
        input:
            bedgraph = BismarkCovToBedgraph.bedgraph,
            chrom_sizes = chrom_sizes
    }

    output {
        File bismark_cov = DragenCytosineReportToBismarkTask.bismark_cov
        File bedgraph = BismarkCovToBedgraph.bedgraph
        File bigwig = BedgraphToBigwig.bigwig
    }
}