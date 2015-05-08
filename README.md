# bigquery-genomics-qc
Genomics QC pipeline on Google's BigQuery
Gregory McInnes
gmcinnes@stanford.edu

## Description
These scripts run queries on genomics data stored in BigQuery that identify samples and positions that fail various quality control metrics.  Once failing samples and positions they are removed from the variantset set in config.py.

Any new queries added should fall under sample or variant level qc.  Sample level qc queries should return sample names that fail with the column header 'sample_id'.  Variant level qc queries should return genomics positions that fail qc with columns for chromosome, start, and end, with the column headers 'reference_name', 'start', and 'end', respectively.  Add query names to GenomicsQueries.py to include them in the analysis and add a file with the query to the sql directory.

## Scripts
#### run_bq_qc.py
Run sample level and/or variant level qc on genomes in BigQuery.  

#### remove_variants.py
Remove samples and/or variants provided in a file from a Google Genomics dataset.

#### Config.py
Settings for access to Google Cloud. Both credentials and dataset information are stored here.

#### GenomicsQueries.py
Query configuration file.  Preset cutoff information is also stored in this file. 

#### BigQueryGenomicsQC.py
Class storing the majority of qc functions.

#### BigQueryClient.py
Client to access BigQuery API and execute basic functions.

#### GoogleGenomicsClient.py
Client to access Google Genomics API and execute basic functions.

 
