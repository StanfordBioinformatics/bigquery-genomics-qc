# bigquery-genomics-qc
Genomics QC pipeline

These scripts run queries on genomics data stored in BigQuery that identify samples and positions that fail various quality control metrics.  Once failing samples and positions they are removed from the variantset set in config.py.

Any new queries added should fall under sample or variant level qc.  Sample level qc queries should return sample names that fail with the column header 'sample_id'.  Variant level qc queries should return genomics positions that fail qc with columns for chromosome, start, and end, with the column headers 'reference_name', 'start', and 'end', respectively.  Add query names to GenomicsQueries.py to include them in the analysis and add a file with the query to the sql directory.
