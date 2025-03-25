process fastq_scan {
    input:
    path fastq_file

    output:
    path "results/*"

    script:
    """
    fastq-scan --input $fastq_file --output results/
    """
}

