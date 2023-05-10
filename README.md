# Pre-xQTL scan cohort-level QC of RNA-seq data

This is for per-cohort, per-tissue level QC of RNA-seq data. The TOPMed IRC receives BAM files and RNA-SeQC metrics from the sequencing cores. For the sake of e/sQTL scans, we flag RNA-seq libraries that are gene expression outliers based on top gene expression PCs. Outlier libraries are not considered when later selecting samples to be used for e/sQTL scans.

## Input

This code requires the following:
* A matrix of gene counts (read counts from RNA-SeQC). This is to be in parquet format, with rows = gene IDs and columns = TOR IDs.
* A matrix of TPM-normalized gene expression (TPM from RNA-SeQC).  This is to be in parquet format, with rows = gene IDs and columns = TOR IDs.
* A matrix of RNA-SeQC metrics. This is to be in TSV format, with columns = metrics and rows = TOR IDs.
* List(s) of TOR IDs to analyze. Potentially multiple lists, in the case that multiple groups of samples should be analyzed independently (e.g., each list could correpond to a cohort/tissue combination). The file names (stripping off any .txt suffix) will be used to identify the groups, so a naming scheme such as `{cohort}_{tissue}.txt` should be used.

## Output

* Gene expression PCs
* A file noting whether each sample is/isn't an outlier in PC space
* Plots of QC metrics, top PCs, and the relationships between them

## Running

You must only have NextFlow (>= v. 21.04) and Singularity (v. 3) installed. NextFlow should be configured as appropriate for your computing platform. Once this has been done, you can run the following command:

```bash
nextflow run -resume --results /path/to/pipeline_results \
                        --counts /path/to/gene-counts.parquet \
                        --rnaseqc /path/to/rnaseqc-metrics.txt \
                        --sample_list_glob '/path/to/sample-lists/*' \
                        /path/to/qc.nf
```

Where `/path/to/sample-lists/*` is a shell glob capturing the file(s) that list TOR IDs to analyze. This should be in quotes as the glob itself should be passed to the pipeline.