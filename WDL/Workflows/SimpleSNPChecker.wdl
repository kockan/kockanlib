version 1.0

task SimpleSNPCheckerTask {
    input {
        File bam
        File bai
        File reference
        File reference_fai
        File snp_bed

        Int? cpu = 2
        Int? memory_gb = 16
        Int? disk_size_gb = 512
    }

    String prefix = basename(bam, ".deduped.sorted.bam")

    command <<<
        set -e
        python3 <<CODE

        import sys
        import pysam

        MAX_PILEUP_DEPTH = 500000
        MIN_BASE_QUALITY = 20

        # ====================================
        # Needs three inputs:
        # - Aligned BAM (position sorted and indexed)
        # - Reference FASTA
        # - Target SNP positions in BED format
        # ====================================

        bam_file = pysam.AlignmentFile("~{bam}", "rb")
        fasta_file = pysam.FastaFile("~{reference}")

        fingerprinting_snps = []

        with open("~{snp_bed}", 'r') as infile:
            for line in infile:
                line = line.rstrip()
                columns = line.split('\t')
                fingerprinting_snps.append((columns[0], int(columns[1]), int(columns[2])))

        outfile = open("~{prefix}.fp_snp_vafs.tsv", 'w')

        for fingerprinting_snp in fingerprinting_snps:
            chromosome = fingerprinting_snp[0]
            start = fingerprinting_snp[1]
            end = fingerprinting_snp[2]

            reference_base = fasta_file.fetch(chromosome, start, end).upper()
            num_ref_bases = 0
            num_alt_bases = 0

            # pysam uses 0-based, half-open intervals
            for pileup_column in bam_file.pileup(chromosome, start, end, max_depth = MAX_PILEUP_DEPTH, min_base_quality = MIN_BASE_QUALITY):
                for pileup_read in pileup_column.pileups:
                    if not pileup_read.is_del and not pileup_read.is_refskip and pileup_column.pos == start:
                        base = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
                        if base == reference_base:
                            num_ref_bases = num_ref_bases + 1
                        else:
                            num_alt_bases = num_alt_bases + 1

                        coverage = num_ref_bases + num_alt_bases
                        if coverage > 100:
                            vaf = num_alt_bases / coverage
                            outfile.write(chromosome + ":" + str(start) + "\t" + "{:.3f}".format(vaf) + "\n")
        outfile.close()
        CODE
    >>>

    output {
        File fp_snp_vafs = "~{prefix}.fp_snp_vafs.tsv"
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size_gb} SSD"
        docker: "us-central1-docker.pkg.dev/broad-gp-hydrogen/hydrogen-dockers/kockan/simple_pysam@sha256:a302f9efe0bf1d4f9998ee1e9dda406223454ccaea0b5619046742221c1d2a74"
    }
}

workflow SimpleSNPChecker {
    input {
        File bam
        File bai
        File reference
        File reference_fai
        File snp_bed
    }

    call SimpleSNPCheckerTask {
        input:
            bam = bam,
            bai = bai,
            reference = reference,
            reference_fai = reference_fai,
            snp_bed = snp_bed
    }

    output {
        File fp_snp_vafs = SimpleSNPCheckerTask.fp_snp_vafs
    }
}