#!/usr/bin/env nextflow

nextflow.enable.dsl=2

COUNTS = params.counts
RNASEQC = params.rnaseqc
SAMPLE_LIST_GLOB = params.sample_list_glob

process pca {

    publishDir "${params.results}/pca"
    memory '75 GB'
    container 'docker.io/porchard/general:20230411143108'
    time '10h'
    cpus 1
    tag "${name}"

    input:
    tuple val(name), path(sample_list), path(gene_counts)

    output:
    tuple val(name), path("${name}.PC-scores.txt"), path("${name}.PC-variance-explained.txt")

    """
    pca.py --sample-list $sample_list $gene_counts
    mv PC-scores.txt ${name}.PC-scores.txt
    mv PC-variance-explained.txt ${name}.PC-variance-explained.txt
    """

}


process get_qc_outliers {

    publishDir "${params.results}/pca-outliers"
    memory '30 GB'
    container 'library://porchard/default/general:20220107'
    tag "${name}"

    input:
    tuple val(name), path(pc_scores)

    output:
    path("${name}.*")
    path("${name}.outlier-status.txt"), emit: outlier_status

    """
    pc-scores-to-absolute-deviation.py --normalize-by-mad $pc_scores > ${name}.ad.txt
    cut -f1-11 ${name}.ad.txt > ad.top.txt
    pc-scores-to-mahalanobis-distance.py --use-pcs 5 $pc_scores > ${name}.md.txt
    score-pca-outliers.py --normalized-absolute-deviation ad.top.txt --mahalanobis-distance ${name}.md.txt --mahalanobis-pcs 5 --max-normalized-absolute-deviation 5 --mahalanobis-pvalue-threshold 0.001 --prefix ${name}.
    """

}


process plot_pca {

    publishDir "${params.results}/plot-pca"
    memory '20 GB'
    container 'library://porchard/default/general:20220107'
    tag "${name}"

    input:
    tuple val(name), path(pc_scores), path(var_explained)

    output:
    path("${name}.pca.png")

    """
    plot-pca.py --panel-height 3 --panel-width 3 --plot-top 30 --id tor $pc_scores $var_explained ${name}.pca.png
    """

}

process rnaseqc_vs_pcs_heatmap {

    publishDir "${params.results}/pca-vs-qc-metrics"
    memory '20 GB'
    container 'library://porchard/default/general:20220107'
    tag "${name}"

    input:
    tuple val(name), path(pc_scores), path(var_explained), path(rnaseqc)

    output:
    path("${name}.heatmap.png")

    """
    plot-pcs-vs-covariates.py --pcs $pc_scores --variance-explained $var_explained --covariates $rnaseqc --index tor --number-pcs 30 --prefix ${name}.
    """

}

process rnaseqc_boxplots {

    publishDir "${params.results}/qc-metric-boxplots"
    memory '20 GB'
    container 'library://porchard/default/general:20220107'

    input:
    path(rnaseqc)
    path(outlier_status)

    output:
    path("rnaseqc.png")

    """
    rnaseqc-metric-boxplots.py rnaseqc.png $rnaseqc ${outlier_status.join(' ')}
    """

}


workflow {
    rnaseqc = Channel.fromPath(RNASEQC)
    counts = Channel.fromPath(COUNTS)
    groups = Channel.fromPath(SAMPLE_LIST_GLOB).map({it -> [it.getName().replaceAll('.txt', ''), it]}) // group name, sample list
    
    pca_output = groups.combine(counts) | pca
    plot_pca(pca_output)
    os = pca_output.map({it -> it[0..1]}) | get_qc_outliers
    rnaseqc_vs_pcs_heatmap(pca_output.combine(rnaseqc))
    rnaseqc_boxplots(rnaseqc, os.outlier_status.toSortedList())
}
