<!-- R Markdown Documentation, DO NOT EDIT THE PLAIN MARKDOWN VERSION OF THIS FILE -->

<!-- Copyright 2015 Google Inc. All rights reserved. -->

<!-- Licensed under the Apache License, Version 2.0 (the "License"); -->
<!-- you may not use this file except in compliance with the License. -->
<!-- You may obtain a copy of the License at -->

<!--     http://www.apache.org/licenses/LICENSE-2.0 -->

<!-- Unless required by applicable law or agreed to in writing, software -->
<!-- distributed under the License is distributed on an "AS IS" BASIS, -->
<!-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. -->
<!-- See the License for the specific language governing permissions and -->
<!-- limitations under the License. -->

# Variant-Level QC





* [Missingness Rate](#missingness-rate)
* [Blacklisted Variants](#blacklisted-variants)
* [Heterozygous Haplotype](#heterozygous-haplotype)
* [Ti/Tv by Genomic Window](#titv-by-genomic-window)
* [Ti/Tv by Depth](#titv-by-depth)
* [Ti/Tv by Alternate Allele Counts](#titv-by-alternate-allele-counts)
* [Hardy Weinberg Equilibrium](#hardy-weinberg-equilibrium)


```r
queryReplacements <- list("_THE_TABLE_"="va_aaa_pilot_data.all_genomes_gvcfs_20150514",
                          "_THE_EXPANDED_TABLE_"="va_aaa_pilot_data.all_genomes_expanded_vcfs_java3",
                          "_PATIENT_INFO_"="va_aaa_pilot_data.patient_info",
                          "_BLACKLISTED_TABLE_"="resources.blacklisted_positions")
```

## Missingness Rate

Identify all variants with a missingness rate greater than a specified cutoff.


```r
cutoff = list("_CUTOFF_"="0.1")
sortAndLimit <- list("#_ORDER_BY_" = "LIMIT 1000")
outputTable = 'qc_tables.variant_missingness'
result <- DisplayAndDispatchQuery("../sql/variant-level-missingness-fail.sql",
                                  project=project,
                                  replacements=c(cutoff, queryReplacements),
                                  outputTable=outputTable)
```

```
SELECT 
variant_id,
"variant_missingness" AS failure_reason,
missingness_rate,
FROM (
  SELECT
  variant_id,
  reference_name,
  start,
  end,
  reference_bases,
  alternate_bases,
  called_allele_count,
  1 - (called_allele_count)/sample_count AS missingness_rate,
  sample_count
  FROM (
    SELECT
    variant_id,
    reference_name,
    start,
    end,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
    SUM(call.genotype >= 0) WITHIN RECORD AS called_allele_count,
    FROM
    [va_aaa_pilot_data.all_genomes_expanded_vcfs_java3]
    # Optionally add clause here to limit the query to a particular
    # region of the genome.
    #_WHERE_
  ) AS g
  CROSS JOIN (
    SELECT
    COUNT(call.call_set_name) AS sample_count
    FROM (
      SELECT 
      call.call_set_name
      FROM
      [va_aaa_pilot_data.all_genomes_gvcfs]
      GROUP BY 
      call.call_set_name)) AS count )
WHERE
missingness_rate > 0.1
# Optionally add a clause here to sort and limit the results.
#_ORDER_BY_

Running query:   RUNNING  2.4s
Running query:   RUNNING  3.1s
Running query:   RUNNING  3.7s
Running query:   RUNNING  4.4s
Running query:   RUNNING  5.0s
Running query:   RUNNING  5.6s
Running query:   RUNNING  6.3s
Running query:   RUNNING  6.9s
Running query:   RUNNING  7.6s
Running query:   RUNNING  8.2s
Running query:   RUNNING  8.9s
Running query:   RUNNING  9.5s
Running query:   RUNNING 10.1s
Running query:   RUNNING 10.8s
Running query:   RUNNING 11.4s
Running query:   RUNNING 12.1s
Running query:   RUNNING 12.7s
Running query:   RUNNING 13.4s
Running query:   RUNNING 14.0s
Running query:   RUNNING 14.6s
Running query:   RUNNING 15.3s
Running query:   RUNNING 15.9s
Running query:   RUNNING 16.5s
Running query:   RUNNING 17.2s
Running query:   RUNNING 17.8s
Running query:   RUNNING 18.5s
Running query:   RUNNING 19.1s
Running query:   RUNNING 19.7s
Running query:   RUNNING 20.4s
Running query:   RUNNING 21.0s

Retrieving data:  3.1s
Retrieving data:  4.3s
Retrieving data:  5.9s
Retrieving data:  7.1s
Retrieving data:  8.3s
Retrieving data:  9.6s
Retrieving data: 11.7s
Retrieving data: 13.0s
```
Number of rows returned by this query: **100000**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri May 29 01:38:39 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> failure_reason </th> <th> missingness_rate </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTUYkMncLyCMs-GNmMze2Uw </td> <td> variant_missingness </td> <td align="right"> 0.74 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNhjaiJkiILrb8aW85dSamAE </td> <td> variant_missingness </td> <td align="right"> 0.99 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyMxj4m6oNILLekc-Goaq_Hw </td> <td> variant_missingness </td> <td align="right"> 0.90 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTYYyNLDGiCozICax4WkqgQ </td> <td> variant_missingness </td> <td align="right"> 0.97 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTYY3tXDGiDWkImGkeeJgHM </td> <td> variant_missingness </td> <td align="right"> 0.95 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTAYwdXhOyC7nOCPrp6hyKMB </td> <td> variant_missingness </td> <td align="right"> 0.38 </td> </tr>
   </table>

## Blacklisted Variants


```r
query <- "../sql/blacklisted-variants.sql"
sortAndLimit <- list("#_ORDER_BY_" = "LIMIT 1000")
outputTable = 'qc_tables.blacklisted'
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(queryReplacements),
                                  outputTable=outputTable)
```

```
SELECT
  seq.variant_id AS variant_id,
  "blacklisted" AS failure_reason,
  bl.Artifact_Type AS Artifact_Type 
FROM (
  SELECT
    variant_id,
    reference_name,
    start,
    end,
  FROM
    [va_aaa_pilot_data.all_genomes_expanded_vcfs_java3]) as seq
JOIN (
  SELECT 
    reference_name,
    start - 1 AS start,
    end,
    Artifact_Type
  FROM 
    [resources.blacklisted_positions]) AS bl
ON
  seq.reference_name = bl.reference_name
WHERE 
  seq.start >= bl.start AND
  seq.end <= bl.end
#_ORDER_BY_


Running query:   RUNNING  2.6s
Running query:   RUNNING  3.2s
Running query:   RUNNING  3.8s
Running query:   RUNNING  4.5s
Running query:   RUNNING  5.1s
Running query:   RUNNING  5.8s
Running query:   RUNNING  6.4s
Running query:   RUNNING  7.0s
Running query:   RUNNING  7.7s
Running query:   RUNNING  8.3s
Running query:   RUNNING  9.0s
Running query:   RUNNING  9.6s
Running query:   RUNNING 10.3s
Running query:   RUNNING 10.9s
Running query:   RUNNING 11.5s
Running query:   RUNNING 12.2s
Running query:   RUNNING 12.8s
Running query:   RUNNING 13.5s
Running query:   RUNNING 14.1s
Running query:   RUNNING 14.8s
Running query:   RUNNING 15.4s
Running query:   RUNNING 16.0s
Running query:   RUNNING 16.7s
Running query:   RUNNING 17.4s
Running query:   RUNNING 18.0s
Running query:   RUNNING 18.6s
Running query:   RUNNING 19.3s
Running query:   RUNNING 19.9s
Running query:   RUNNING 20.6s
Running query:   RUNNING 21.2s
Running query:   RUNNING 21.8s
Running query:   RUNNING 22.5s
Running query:   RUNNING 23.1s
Running query:   RUNNING 23.7s
Running query:   RUNNING 24.6s
Running query:   RUNNING 25.2s
Running query:   RUNNING 25.9s
Running query:   RUNNING 26.5s
Running query:   RUNNING 27.1s
Running query:   RUNNING 27.8s
Running query:   RUNNING 28.4s
Running query:   RUNNING 29.1s
Running query:   RUNNING 29.7s
Running query:   RUNNING 30.4s
Running query:   RUNNING 31.0s
Running query:   RUNNING 31.7s
Running query:   RUNNING 32.3s
Running query:   RUNNING 32.9s
Running query:   RUNNING 33.6s
Running query:   RUNNING 34.2s
Running query:   RUNNING 34.9s
Running query:   RUNNING 35.5s
Running query:   RUNNING 36.2s
Running query:   RUNNING 36.8s
Running query:   RUNNING 37.4s
Running query:   RUNNING 38.1s
Running query:   RUNNING 38.7s
Running query:   RUNNING 39.3s
Running query:   RUNNING 40.0s
Running query:   RUNNING 40.6s

Retrieving data:  3.2s
Retrieving data:  4.5s
Retrieving data:  6.2s
Retrieving data:  7.4s
Retrieving data:  8.5s
Retrieving data:  9.7s
Retrieving data: 10.9s
Retrieving data: 12.4s
```

Number of rows returned by this query: **100000**.

First few results
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri May 29 01:39:34 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> failure_reason </th> <th> Artifact_Type </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNxjK9MkbIJmz2tqq5aPZIA </td> <td> blacklisted </td> <td> centromeric_repeat </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNxjY9ckbIJvEw8zv-aa_4AE </td> <td> blacklisted </td> <td> centromeric_repeat </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNxja9skbILuo0aD8gqWVyQE </td> <td> blacklisted </td> <td> centromeric_repeat </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNxj-9skbIMi8oIXuz9blcA </td> <td> blacklisted </td> <td> centromeric_repeat </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNxiq98kbILq5taqRjNjPCQ </td> <td> blacklisted </td> <td> centromeric_repeat </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNxik-MkbIJ-Yk4Xp-4CY6QE </td> <td> blacklisted </td> <td> centromeric_repeat </td> </tr>
   </table>

## Heterozygous Haplotype
For each variant within the X and Y chromosome, identify heterozygous variants in male genomes.


```r
sortAndLimit <- list("_LIMIT_" = "LIMIT 1000")
outputTable = 'qc_tables.heterozygous_haplotype'
result <- DisplayAndDispatchQuery("../sql/sex-chromosome-heterozygous-haplotypes.sql",
                                  project=project,
                                  replacements=c(queryReplacements),
                                  outputTable=outputTable)
```

```
SELECT
  variant_id,
  sample_id,
  "heterozygous_haplotype" AS failure_reason,
FROM (
SELECT
  variant_id,
  reference_name,
  start,
  end,
  reference_bases,
  call.call_set_name AS sample_id,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype,
FROM(FLATTEN((
  [va_aaa_pilot_data.all_genomes_expanded_vcfs_java3]), call.call_set_name))
WHERE
  reference_name IN ('chrX', 'chrY')
OMIT
  call IF (2 > COUNT(call.genotype))
  OR EVERY(call.genotype <= 0)
  OR EVERY(call.genotype = 1)
  # Pseudoautosomal Region 1
  OR (reference_name = 'chrX'
    AND start > 60001
    AND end < 2699520)
  OR (reference_name = 'chrY'
    AND start > 10001
    AND end < 2649520)
  # Pseudoautosomal Region 2
  OR (reference_name = 'chrX'
    AND start > 155260560
    AND end < 155270560)
  OR (reference_name = 'chrY' 
    AND start > 59363566
    AND end < 59373566)) AS seq
JOIN (
  SELECT
    IlluminaID,
    SEX
  FROM
    [va_aaa_pilot_data.patient_info] ) AS info
ON
  seq.sample_id = info.IlluminaID
WHERE
  SEX = 'M'
# Optionally add a clause here to sort and limit the results.
#_ORDER_BY_

Running query:   RUNNING  2.6s
Running query:   RUNNING  3.2s
Running query:   RUNNING  3.9s
Running query:   RUNNING  4.5s
Running query:   RUNNING  5.1s
Running query:   RUNNING  5.8s
Running query:   RUNNING  7.2s
Running query:   RUNNING  7.8s
Running query:   RUNNING  8.5s
Running query:   RUNNING  9.1s
Running query:   RUNNING  9.8s
Running query:   RUNNING 10.5s
Running query:   RUNNING 11.1s
Running query:   RUNNING 11.7s
Running query:   RUNNING 12.4s
Running query:   RUNNING 13.0s
Running query:   RUNNING 13.7s
Running query:   RUNNING 14.3s
Running query:   RUNNING 14.9s
Running query:   RUNNING 15.6s
Running query:   RUNNING 16.2s
Running query:   RUNNING 16.8s
Running query:   RUNNING 17.5s
Running query:   RUNNING 18.1s
Running query:   RUNNING 18.8s
Running query:   RUNNING 19.4s
Running query:   RUNNING 20.0s
Running query:   RUNNING 20.7s
Running query:   RUNNING 21.3s
Running query:   RUNNING 22.0s
Running query:   RUNNING 22.6s
Running query:   RUNNING 23.2s
Running query:   RUNNING 23.9s
Running query:   RUNNING 24.5s
Running query:   RUNNING 25.1s
Running query:   RUNNING 25.8s
Running query:   RUNNING 26.4s
Running query:   RUNNING 27.1s
Running query:   RUNNING 27.7s
Running query:   RUNNING 28.3s
Running query:   RUNNING 29.0s
Running query:   RUNNING 29.6s
Running query:   RUNNING 30.2s
Running query:   RUNNING 30.9s
Running query:   RUNNING 31.5s
Running query:   RUNNING 32.1s
Running query:   RUNNING 32.8s
Running query:   RUNNING 33.4s
Running query:   RUNNING 34.1s
Running query:   RUNNING 34.7s
Running query:   RUNNING 35.3s
Running query:   RUNNING 36.0s
Running query:   RUNNING 36.6s
Running query:   RUNNING 37.2s
Running query:   RUNNING 37.9s
Running query:   RUNNING 38.6s
Running query:   RUNNING 39.2s
Running query:   RUNNING 39.8s
Running query:   RUNNING 40.5s
Running query:   RUNNING 41.1s
Running query:   RUNNING 41.8s
Running query:   RUNNING 42.4s
Running query:   RUNNING 43.0s
Running query:   RUNNING 43.7s
Running query:   RUNNING 44.3s
Running query:   RUNNING 44.9s
Running query:   RUNNING 45.6s
Running query:   RUNNING 46.2s
Running query:   RUNNING 46.9s
Running query:   RUNNING 47.5s
Running query:   RUNNING 48.1s
Running query:   RUNNING 48.8s
Running query:   RUNNING 49.4s
Running query:   RUNNING 50.0s
Running query:   RUNNING 50.7s
Running query:   RUNNING 51.3s
Running query:   RUNNING 52.0s
Running query:   RUNNING 52.6s
Running query:   RUNNING 53.2s
Running query:   RUNNING 53.9s
Running query:   RUNNING 54.5s
Running query:   RUNNING 55.3s
Running query:   RUNNING 55.9s
Running query:   RUNNING 56.5s
Running query:   RUNNING 57.2s
Running query:   RUNNING 57.9s
Running query:   RUNNING 58.5s
Running query:   RUNNING 59.1s
Running query:   RUNNING 59.8s
Running query:   RUNNING 60.4s
Running query:   RUNNING 61.0s
Running query:   RUNNING 61.7s
Running query:   RUNNING 62.3s
Running query:   RUNNING 62.9s
Running query:   RUNNING 63.6s
Running query:   RUNNING 64.2s
Running query:   RUNNING 64.8s
Running query:   RUNNING 65.5s
Running query:   RUNNING 66.5s

Retrieving data:  2.1s
Retrieving data:  4.1s
Retrieving data:  5.4s
Retrieving data:  7.1s
Retrieving data:  8.5s
Retrieving data: 10.2s
Retrieving data: 11.5s
Retrieving data: 13.3s
Retrieving data: 15.4s
```
Number of rows returned by this query: **100000**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri May 29 01:41:00 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> sample_id </th> <th> failure_reason </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBiwjtc6IMSTg8mc_ti7LA </td> <td> LP6005038-DNA_B02 </td> <td> heterozygous_haplotype </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBiwjtc6IMSTg8mc_ti7LA </td> <td> LP6005243-DNA_E06 </td> <td> heterozygous_haplotype </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBiMkdc6IMuXtLKq5IO5oAE </td> <td> LP6005144-DNA_D11 </td> <td> heterozygous_haplotype </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBjalNc6IOCh9P-Q_ebUFA </td> <td> LP6005692-DNA_G01 </td> <td> heterozygous_haplotype </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBiyu-kVILDA8-LAg-jr0wE </td> <td> LP6005038-DNA_B06 </td> <td> heterozygous_haplotype </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBi-u-kVIKazkNSwteO2Jg </td> <td> LP6005692-DNA_B02 </td> <td> heterozygous_haplotype </td> </tr>
   </table>

## Ti/Tv By Genomic Window

```r
query <- "../sql/titv-by-genomic-window-fail.sql"
sortAndLimit <- list("#_ORDER_BY_" = "LIMIT 1000")
max <- 3.0
min <- 1.5
cutoffs <- list("_MAX_" = max, "_MIN_" = min)
outputTable = 'qc_tables.titv_genomic_window'
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(queryReplacements, cutoffs),
                                  outputTable=outputTable)
```

```
SELECT
var.variant_id AS variant_id,
titv,
"titv_by_genomic_window" AS failure_reason,
FROM (
  SELECT
  variant_id,
  reference_name,
  start,
  end,
  INTEGER(FLOOR(start / 100000)) AS window,
  FROM
  [va_aaa_pilot_data.all_genomes_expanded_vcfs_java3]
  OMIT call IF EVERY(call.genotype <= 0)
  # Optionally add clause here to limit the query to a particular
  # region of the genome.
  #_WHERE_ 
) AS var
JOIN (
  SELECT
  reference_name,
  window,
  window_start,
  transitions,
  transversions,
  titv,
  num_variants_in_window,
  FROM (
    SELECT
    reference_name,
    window,
    window * 100000 AS window_start,
    transitions,
    transversions,
    transitions/transversions AS titv,
    num_variants_in_window,
    FROM (
      SELECT
      reference_name,
      window,
      SUM(mutation IN ('A->G', 'G->A', 'C->T', 'T->C')) AS transitions,
      SUM(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                       'A->T', 'T->A', 'C->G', 'G->C')) AS transversions,
      COUNT(mutation) AS num_variants_in_window
      FROM (
        SELECT
        reference_name,
        INTEGER(FLOOR(start / 100000)) AS window,
        CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
        COUNT(alternate_bases) WITHIN RECORD AS num_alts,
        FROM
        [va_aaa_pilot_data.all_genomes_expanded_vcfs_java3]
        # Optionally add clause here to limit the query to a particular
        # region of the genome.
        WHERE reference_name = 'chr22'
        HAVING
        # Skip 1/2 genotypes _and non-SNP variants
        num_alts = 1
        AND reference_bases IN ('A','C','G','T')
        AND alternate_bases IN ('A','C','G','T'))
      GROUP BY
      reference_name,
      window))
  WHERE
  titv > 3 OR
  titv < 1.5) as win
ON 
var.window = win.window
GROUP BY 
variant_id,
var.reference_name,
win.window_start,
titv,
#_ORDER_BY_


Running query:   RUNNING  2.5s
Running query:   RUNNING  3.2s
Running query:   RUNNING  3.8s
Running query:   RUNNING  4.5s
Running query:   RUNNING  5.1s
Running query:   RUNNING  5.8s
Running query:   RUNNING  6.4s
Running query:   RUNNING  7.0s
Running query:   RUNNING  7.7s
Running query:   RUNNING  8.3s
Running query:   RUNNING  8.9s
Running query:   RUNNING  9.6s
Running query:   RUNNING 10.3s
Running query:   RUNNING 10.9s
Running query:   RUNNING 11.5s
Running query:   RUNNING 12.2s
Running query:   RUNNING 12.8s
Running query:   RUNNING 13.5s
Running query:   RUNNING 14.1s
Running query:   RUNNING 14.7s
Running query:   RUNNING 15.4s
Running query:   RUNNING 16.0s
Running query:   RUNNING 16.6s
Running query:   RUNNING 17.3s
Running query:   RUNNING 17.9s
Running query:   RUNNING 18.5s
Running query:   RUNNING 19.2s
Running query:   RUNNING 19.8s
Running query:   RUNNING 20.5s
Running query:   RUNNING 21.1s
Running query:   RUNNING 21.8s
Running query:   RUNNING 22.4s
Running query:   RUNNING 23.1s
Running query:   RUNNING 23.7s
Running query:   RUNNING 24.3s

Retrieving data:  2.2s
Retrieving data:  3.4s
Retrieving data:  5.1s
Retrieving data:  6.2s
Retrieving data:  7.3s
Retrieving data:  8.5s
Retrieving data:  9.9s
Retrieving data: 11.3s
Retrieving data: 12.9s
```

Number of rows returned by this query: **100000**.

First few results
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri May 29 01:41:42 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> titv </th> <th> failure_reason </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTMY9vbsCSD1-rb0j4f9_6cB </td> <td align="right"> 1.09 </td> <td> titv_by_genomic_window </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTMYovfsCSD7tvuynrCh4Fc </td> <td align="right"> 1.09 </td> <td> titv_by_genomic_window </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTMY6vjsCSCL-d3AroH_5IEB </td> <td align="right"> 1.09 </td> <td> titv_by_genomic_window </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTMYrvrsCSDA19uqquDAsVQ </td> <td align="right"> 1.09 </td> <td> titv_by_genomic_window </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTMY4PrsCSDHie_9wOq1pJsB </td> <td align="right"> 1.09 </td> <td> titv_by_genomic_window </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTMYgPvsCSCS68bysbzw1-0B </td> <td align="right"> 1.09 </td> <td> titv_by_genomic_window </td> </tr>
   </table>

## Ti/Tv By Depth

We want to identify all the regions of the genome where the Ti/Tv ratio is outside of the expected range.  Another method we can use to do this is calculating the transition-transversion ratio by depth of coverage.  

```r
query <- "../sql/titv-by-depth-fail.sql"
max <- 3
min <- 1.5
cutoffs <- list("_MAX_" = max, "_MIN_" = min)
sortAndLimit <- list("#_ORDER_BY_" = "LIMIT 1000")
outputTable = 'qc_tables.titv_depth'
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(queryReplacements, cutoffs),
                                  outputTable)
```

```
SELECT 
  variant_id,
  titv.titv_ratio,
  "titv_by_depth" AS failure_reason,
FROM(
SELECT
  call.call_set_name AS sample_id,
  variant_id,
  call.DP AS depth,
FROM(FLATTEN((
  [va_aaa_pilot_data.all_genomes_expanded_vcfs_java3]), call.call_set_name))) AS var
JOIN (
SELECT
  sample_id,
  titv_ratio,
  depth,
  FROM (
    SELECT
    call.call_set_name AS sample_id,
    (transitions/transversions) AS titv_ratio,
    call.DP AS depth,
    FROM (
      SELECT
      call.call_set_name,
      SUM(mutation IN ('A->G', 'G->A', 'C->T', 'T->C')) AS transitions,
      SUM(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                       'A->T', 'T->A', 'C->G', 'G->C')) AS transversions,
      call.DP,
      FROM (
        SELECT
        call.call_set_name,
        CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
        COUNT(alternate_bases) WITHIN RECORD AS num_alts,
        call.DP
        FROM (
          SELECT
          call.call_set_name,
          reference_bases,
          GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
          call.genotype,
          call.DP,
          FROM
          [va_aaa_pilot_data.all_genomes_expanded_vcfs_java3]
          # Optionally add clause here to limit the query to a particular
          # region of the genome.
          #_WHERE_  
        )
        WHERE
        call.DP is not null
        HAVING
        # Skip 1/2 genotypes _and non-SNP variants
        num_alts = 1
        AND reference_bases IN ('A','C','G','T')
        AND alternate_bases IN ('A','C','G','T'))
      GROUP BY 
      call.call_set_name,
      call.DP,)
    WHERE
    transversions > 0
    GROUP BY
    sample_id,
    titv_ratio,
    depth,)
WHERE
titv_ratio > 3
OR titv_ratio < 1.5) AS titv
ON
  var.sample_id = titv.sample_id
  AND var.depth = titv.depth
#_ORDER_BY_
Running query:   RUNNING  2.0s
Running query:   RUNNING  2.7s
Running query:   RUNNING  3.3s
Running query:   RUNNING  3.9s
Running query:   RUNNING  4.6s
Running query:   RUNNING  5.2s
Running query:   RUNNING  5.8s
Running query:   RUNNING  6.5s
Running query:   RUNNING  7.1s
Running query:   RUNNING  7.8s
Running query:   RUNNING  8.4s
Running query:   RUNNING  9.0s
Running query:   RUNNING  9.7s
Running query:   RUNNING 10.3s
Running query:   RUNNING 10.9s
Running query:   RUNNING 11.6s
Running query:   RUNNING 12.3s
Running query:   RUNNING 12.9s
Running query:   RUNNING 13.6s
Running query:   RUNNING 14.2s
Running query:   RUNNING 14.9s
Running query:   RUNNING 15.5s
Running query:   RUNNING 16.2s
Running query:   RUNNING 16.8s
Running query:   RUNNING 17.4s
Running query:   RUNNING 18.1s
Running query:   RUNNING 18.7s
Running query:   RUNNING 19.3s
Running query:   RUNNING 20.0s
Running query:   RUNNING 20.6s
Running query:   RUNNING 21.3s
Running query:   RUNNING 21.9s
Running query:   RUNNING 22.5s
Running query:   RUNNING 23.2s
Running query:   RUNNING 23.8s
Running query:   RUNNING 24.4s
Running query:   RUNNING 25.1s
Running query:   RUNNING 25.7s
Running query:   RUNNING 26.3s
Running query:   RUNNING 27.0s
Running query:   RUNNING 27.6s
Running query:   RUNNING 28.2s
Running query:   RUNNING 28.9s
Running query:   RUNNING 29.5s
Running query:   RUNNING 30.1s
Running query:   RUNNING 30.8s
Running query:   RUNNING 31.4s
Running query:   RUNNING 32.1s
Running query:   RUNNING 32.7s
Running query:   RUNNING 33.3s
Running query:   RUNNING 34.0s
Running query:   RUNNING 34.6s
Running query:   RUNNING 35.2s
Running query:   RUNNING 35.9s
Running query:   RUNNING 36.5s
Running query:   RUNNING 37.2s
Running query:   RUNNING 37.8s
Running query:   RUNNING 38.4s
Running query:   RUNNING 39.1s
Running query:   RUNNING 39.7s
Running query:   RUNNING 40.4s
Running query:   RUNNING 41.0s
Running query:   RUNNING 41.6s
Running query:   RUNNING 42.3s
Running query:   RUNNING 42.9s
Running query:   RUNNING 43.5s
Running query:   RUNNING 44.2s
Running query:   RUNNING 44.8s
Running query:   RUNNING 45.4s
Running query:   RUNNING 46.1s
Running query:   RUNNING 46.7s
Running query:   RUNNING 47.4s
Running query:   RUNNING 48.0s
Running query:   RUNNING 48.6s
Running query:   RUNNING 49.3s
Running query:   RUNNING 49.9s
Running query:   RUNNING 50.5s
Running query:   RUNNING 51.2s
Running query:   RUNNING 51.8s
Running query:   RUNNING 53.9s
Running query:   RUNNING 54.5s
Running query:   RUNNING 55.2s
Running query:   RUNNING 55.8s
Running query:   RUNNING 56.4s
Running query:   RUNNING 57.1s
Running query:   RUNNING 57.7s
Running query:   RUNNING 58.4s
Running query:   RUNNING 59.0s
Running query:   RUNNING 59.7s
Running query:   RUNNING 60.3s
Running query:   RUNNING 60.9s
Running query:   RUNNING 61.6s
Running query:   RUNNING 62.2s
Running query:   RUNNING 62.9s
Running query:   RUNNING 63.5s
Running query:   RUNNING 64.2s
Running query:   RUNNING 64.8s
Running query:   RUNNING 65.4s
Running query:   RUNNING 66.1s
Running query:   RUNNING 66.8s
Running query:   RUNNING 67.4s
Running query:   RUNNING 68.0s
Running query:   RUNNING 68.7s
Running query:   RUNNING 69.3s

Retrieving data:  3.4s
Retrieving data:  5.4s
Retrieving data:  7.2s
Retrieving data:  8.7s
Retrieving data: 10.2s
Retrieving data: 11.9s
Retrieving data: 13.3s
Retrieving data: 15.0s
```

Number of rows returned by this query: **100000**.

First few results
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri May 29 01:43:11 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> titv_titv_ratio </th> <th> failure_reason </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNxjp8adLILDflpnJurzi8gE </td> <td align="right"> 1.22 </td> <td> titv_by_depth </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNxjt8adLIPn8rrXMwsiJKQ </td> <td align="right"> 1.25 </td> <td> titv_by_depth </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNxjt8adLIPn8rrXMwsiJKQ </td> <td align="right"> 1.49 </td> <td> titv_by_depth </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyMxjur6YRIN-qjqbMzpW1iwE </td> <td align="right"> 1.41 </td> <td> titv_by_depth </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyMhjGgbISIObS0bvW6vHAuwE </td> <td align="right"> 1.44 </td> <td> titv_by_depth </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyMhjGgbISIObS0bvW6vHAuwE </td> <td align="right"> 1.43 </td> <td> titv_by_depth </td> </tr>
   </table>


## Ti/Tv By Alternate Allele Counts
Collect all the alternate allele counts that are outside our desired range.

```r
query <- "../sql/titv-by-allternate-allele-fail.sql"
max <- 3
min <- 1.5
cutoffs <- list("_MAX_" = max, "_MIN_" = min)
sortAndLimit <- list("#_ORDER_BY_" = "LIMIT 1000")
#result <- DisplayAndDispatchQuery(query,
#                                  project=project,
#                                  replacements=c(queryReplacements, cutoffs, sortAndLimit))
```

Get the variant ids for all the failed groups


I'll write this query later.  There are no variants that fail this qc step for our current dataset.

## Hardy Weinberg Equilibrium

Here we want to identify the variants that are out of Hardy Weinberg Equilibrium.  We want to remove the top 0.05 quantile of variants, so first we have to define what the cutoff for the chi squared value should be.

```r
quantile <- list("_QUANTILE_" = 1999) # <- Define quantile by number. 
                                  # The 1999th quantile selects the value that partitions the top 0.05% of values, 
                                  # assuming there are 2000 quantiles.
result <- DisplayAndDispatchQuery("../sql/hwe-quantile.sql",
                                  project=project,
                                  replacements=c(queryReplacements, quantile))
```

```
# Get chisq quantile cutoff
SELECT 
  quantile,
  row_num,
FROM(
  SELECT
    quantile,
    ROW_NUMBER() OVER (ORDER BY quantile ASC) row_num,
  FROM (
    SELECT
      QUANTILES(chisq, 2000) AS quantile
    FROM js(
      (SELECT
        reference_name,
        start,
        reference_bases,
        alternate_bases,
        hom_ref AS obs_hom1,
        het AS obs_het,
        hom_alt AS obs_hom2,
        hom_ref + het + hom_alt AS sample_count,
       FROM (
         SELECT
          reference_name,
          start,
          END,
          reference_bases,
          GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
          COUNT(alternate_bases) WITHIN RECORD AS num_alts,
          SUM(EVERY(0 = call.genotype)) WITHIN call AS hom_ref,
          SUM(EVERY(1 = call.genotype)) WITHIN call AS hom_alt,
          SUM(SOME(0 = call.genotype)
             AND SOME(1 = call.genotype)) WITHIN call AS het,
         FROM
          [va_aaa_pilot_data.all_genomes_expanded_vcfs_java3]
         # Optionally add a clause here to limit the query to a particular
         # region of the genome.
         #_WHERE_
         HAVING
         # Skip 1/2 genotypes
          num_alts = 1
       )),
      // Start javascript function
      // Input Columns
      reference_name, start, reference_bases, alternate_bases, obs_hom1, obs_hom2, obs_het, sample_count,
      // Output Schema
      "[{name: 'reference_name', type: 'string'},
      {name: 'start', type: 'integer'},
      {name: 'reference_bases', type: 'string'},
      {name: 'alternate_bases', type: 'string'},
      {name: 'obs_hom1', type: 'integer'},
      {name: 'obs_het', type: 'integer'},
      {name: 'obs_hom2', type: 'integer'},
      {name: 'e_hom1', type: 'integer'},
      {name: 'e_het', type: 'integer'},
      {name: 'e_hom2', type: 'integer'},        
      {name: 'chisq', type: 'float'}]",
      // Function
      "function(r, emit) {
        var e_hom1 = Math.pow((r.obs_hom1 + (r.obs_het/2)) / r.sample_count, 2) * r.sample_count;
        var e_het = 2 * ((r.obs_hom1 + (r.obs_het/2)) / r.sample_count) * ((r.obs_hom2 + (r.obs_het/2)) / r.sample_count) * r.sample_count;
        var e_hom2 = Math.pow((r.obs_hom2 + (r.obs_het/2)) / r.sample_count, 2) * r.sample_count;
        var chisq = (Math.pow(r.obs_hom1 - e_hom1, 2) / e_hom1) + (Math.pow(r.obs_het - e_het, 2) / e_het) + (Math.pow(r.obs_hom2 - e_hom2, 2) / e_hom2);
        emit({
          reference_name: r.reference_name,
          start: r.start,
          reference_bases: r.reference_bases,
          alternate_bases: r.alternate_bases,
          obs_hom1: r.obs_hom1,
          obs_hom2: r.obs_hom2,
          obs_het: r.obs_het,
          e_hom1: e_hom1,
          e_hom2: e_hom2,
          e_het: e_het,
          chisq: chisq
        })
      }"
      )))
WHERE 
  row_num = 1999

Running query:   RUNNING  2.6s
Running query:   RUNNING  3.2s
Running query:   RUNNING  4.1s
Running query:   RUNNING  4.8s
Running query:   RUNNING  5.4s
Running query:   RUNNING  6.1s
Running query:   RUNNING  6.7s
Running query:   RUNNING  7.3s
Running query:   RUNNING  8.0s
Running query:   RUNNING  8.6s
Running query:   RUNNING  9.3s
Running query:   RUNNING  9.9s
Running query:   RUNNING 10.6s
Running query:   RUNNING 11.2s
Running query:   RUNNING 11.9s
Running query:   RUNNING 12.5s
Running query:   RUNNING 13.1s
Running query:   RUNNING 13.8s
Running query:   RUNNING 14.4s
Running query:   RUNNING 15.1s
Running query:   RUNNING 15.7s
Running query:   RUNNING 16.3s
Running query:   RUNNING 17.0s
Running query:   RUNNING 17.6s
Running query:   RUNNING 18.2s
Running query:   RUNNING 18.9s
Running query:   RUNNING 19.5s
Running query:   RUNNING 20.1s
Running query:   RUNNING 20.8s
Running query:   RUNNING 21.4s
Running query:   RUNNING 22.1s
Running query:   RUNNING 22.7s
Running query:   RUNNING 23.3s
Running query:   RUNNING 24.0s
Running query:   RUNNING 24.6s
Running query:   RUNNING 25.3s
Running query:   RUNNING 25.9s
Running query:   RUNNING 26.5s
Running query:   RUNNING 27.2s
Running query:   RUNNING 27.8s
Running query:   RUNNING 28.5s
Running query:   RUNNING 29.2s
Running query:   RUNNING 29.8s
Running query:   RUNNING 30.5s
Running query:   RUNNING 31.1s
Running query:   RUNNING 31.7s
Running query:   RUNNING 32.4s
Running query:   RUNNING 33.0s
Running query:   RUNNING 33.6s
Running query:   RUNNING 34.3s
Running query:   RUNNING 34.9s
Running query:   RUNNING 35.6s
Running query:   RUNNING 36.2s
Running query:   RUNNING 36.8s
Running query:   RUNNING 37.5s
Running query:   RUNNING 38.1s
Running query:   RUNNING 38.7s
Running query:   RUNNING 39.4s
Running query:   RUNNING 40.0s
Running query:   RUNNING 40.7s
Running query:   RUNNING 41.3s
Running query:   RUNNING 41.9s
Running query:   RUNNING 42.6s
Running query:   RUNNING 43.2s
Running query:   RUNNING 43.9s
Running query:   RUNNING 44.5s
Running query:   RUNNING 45.1s
Running query:   RUNNING 45.8s
Running query:   RUNNING 46.4s
Running query:   RUNNING 47.1s
Running query:   RUNNING 47.7s
Running query:   RUNNING 48.4s
Running query:   RUNNING 49.0s
Running query:   RUNNING 49.7s
Running query:   RUNNING 50.3s
Running query:   RUNNING 51.0s
Running query:   RUNNING 51.6s
Running query:   RUNNING 52.3s
Running query:   RUNNING 52.9s
Running query:   RUNNING 53.5s
Running query:   RUNNING 54.2s
```

Displaying the results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri May 29 01:44:08 2015 -->
<table border=1>
<tr> <th> quantile </th> <th> row_num </th>  </tr>
  <tr> <td align="right"> 663.30 </td> <td align="right"> 1999 </td> </tr>
   </table>

Determine the cutoffs:

```r
maxChiSq = result$quantile
```
Cutoff: 663.3034678

Determine which genomes are outside our desired range

```r
values = list("_CUTOFF_" = maxChiSq)
sortAndLimit <- list("#_ORDER_BY_" = "LIMIT 1000")
outputTable = 'qc_tables.hardy_weinberg'
result <- DisplayAndDispatchQuery("../sql/hwe-fail.sql",
                                  project=project,
                                  replacements=c(queryReplacements, values),
                                  outputTable=outputTable)
```

```
# Get all variants that have a chi squared value above a definited limit
SELECT
  variant_id,
  chisq,
  "hardy_weinberg" AS failure_reason,
FROM js(
    (SELECT
      variant_id,
      reference_name,
      start,
      reference_bases,
      alternate_bases,
      hom_ref AS obs_hom1,
      het AS obs_het,
      hom_alt AS obs_hom2,
      hom_ref + het + hom_alt AS sample_count,
    FROM (
      SELECT
        variant_id,
        reference_name,
        start,
        END,
        reference_bases,
        GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
        COUNT(alternate_bases) WITHIN RECORD AS num_alts,
        SUM(EVERY(0 = call.genotype)) WITHIN call AS hom_ref,
        SUM(EVERY(1 = call.genotype)) WITHIN call AS hom_alt,
        SUM(SOME(0 = call.genotype)
          AND SOME(1 = call.genotype)) WITHIN call AS het,
      FROM
        [va_aaa_pilot_data.all_genomes_expanded_vcfs_java3]
      # Optionally add a clause here to limit the query to a particular
      # region of the genome.
      #_WHERE_
      HAVING
        # Skip 1/2 genotypes
        num_alts = 1
        )),
        // Start javascript function
        // Input Columns
        variant_id, reference_name, start, reference_bases, alternate_bases, obs_hom1, obs_hom2, obs_het, sample_count,
        // Output Schema
        "[{name: 'variant_id', type: 'string'},
        {name: 'reference_name', type: 'string'},
        {name: 'start', type: 'integer'},
        {name: 'reference_bases', type: 'string'},
        {name: 'alternate_bases', type: 'string'},
        {name: 'obs_hom1', type: 'integer'},
        {name: 'obs_het', type: 'integer'},
        {name: 'obs_hom2', type: 'integer'},
        {name: 'e_hom1', type: 'integer'},
        {name: 'e_het', type: 'integer'},
        {name: 'e_hom2', type: 'integer'},        
        {name: 'chisq', type: 'float'}]",
        // Function
        "function(r, emit) {
          var e_hom1 = Math.pow((r.obs_hom1 + (r.obs_het/2)) / r.sample_count, 2) * r.sample_count;
          var e_het = 2 * ((r.obs_hom1 + (r.obs_het/2)) / r.sample_count) * ((r.obs_hom2 + (r.obs_het/2)) / r.sample_count) * r.sample_count;
          var e_hom2 = Math.pow((r.obs_hom2 + (r.obs_het/2)) / r.sample_count, 2) * r.sample_count;
          var chisq = (Math.pow(r.obs_hom1 - e_hom1, 2) / e_hom1) + (Math.pow(r.obs_het - e_het, 2) / e_het) + (Math.pow(r.obs_hom2 - e_hom2, 2) / e_hom2);
          emit({
            variant_id: r.variant_id,
            reference_name: r.reference_name,
            start: r.start,
            reference_bases: r.reference_bases,
            alternate_bases: r.alternate_bases,
            obs_hom1: r.obs_hom1,
            obs_hom2: r.obs_hom2,
            obs_het: r.obs_het,
            e_hom1: e_hom1,
            e_hom2: e_hom2,
            e_het: e_het,
            chisq: chisq
          })
        }"
      )
WHERE
  chisq > 663.303467797427
#_ORDER_BY_


Running query:   RUNNING  2.6s
Running query:   RUNNING  3.2s
Running query:   RUNNING  3.9s
Running query:   RUNNING  4.5s
Running query:   RUNNING  5.2s
Running query:   RUNNING  5.8s
Running query:   RUNNING  6.4s
Running query:   RUNNING  7.1s
Running query:   RUNNING  7.7s
Running query:   RUNNING  8.4s
Running query:   RUNNING  9.0s
Running query:   RUNNING  9.6s
Running query:   RUNNING 10.3s
Running query:   RUNNING 10.9s
Running query:   RUNNING 11.5s
Running query:   RUNNING 12.2s
Running query:   RUNNING 12.8s
Running query:   RUNNING 13.5s
Running query:   RUNNING 14.1s
Running query:   RUNNING 14.7s
Running query:   RUNNING 15.4s
Running query:   RUNNING 16.0s
Running query:   RUNNING 16.7s
Running query:   RUNNING 17.3s
Running query:   RUNNING 17.9s
Running query:   RUNNING 18.9s
Running query:   RUNNING 19.6s
Running query:   RUNNING 20.2s
Running query:   RUNNING 20.8s
Running query:   RUNNING 21.5s
Running query:   RUNNING 22.1s
Running query:   RUNNING 22.8s
Running query:   RUNNING 23.5s
Running query:   RUNNING 24.1s
Running query:   RUNNING 24.8s
Running query:   RUNNING 25.4s
Running query:   RUNNING 26.0s
Running query:   RUNNING 26.7s
Running query:   RUNNING 27.3s
Running query:   RUNNING 27.9s
Running query:   RUNNING 28.6s
Running query:   RUNNING 29.2s
Running query:   RUNNING 29.9s
Running query:   RUNNING 30.5s
Running query:   RUNNING 31.1s
Running query:   RUNNING 31.8s
Running query:   RUNNING 32.4s
Running query:   RUNNING 33.0s
Running query:   RUNNING 33.7s
Running query:   RUNNING 34.3s
Running query:   RUNNING 35.0s
Running query:   RUNNING 35.6s
Running query:   RUNNING 36.2s
Running query:   RUNNING 36.9s
Running query:   RUNNING 37.5s
Running query:   RUNNING 38.2s
Running query:   RUNNING 38.8s
Running query:   RUNNING 39.4s
Running query:   RUNNING 40.1s
Running query:   RUNNING 40.7s
Running query:   RUNNING 41.4s
Running query:   RUNNING 42.0s
Running query:   RUNNING 42.7s
Running query:   RUNNING 43.3s
Running query:   RUNNING 43.9s
Running query:   RUNNING 44.6s
Running query:   RUNNING 45.2s
Running query:   RUNNING 45.8s
Running query:   RUNNING 46.5s
Running query:   RUNNING 47.1s
Running query:   RUNNING 47.8s
Running query:   RUNNING 48.4s
Running query:   RUNNING 49.0s
Running query:   RUNNING 49.7s
Running query:   RUNNING 50.3s
Running query:   RUNNING 51.0s
Running query:   RUNNING 51.6s
Running query:   RUNNING 52.2s
Running query:   RUNNING 52.9s
Running query:   RUNNING 53.5s
Running query:   RUNNING 54.2s
Running query:   RUNNING 54.8s
Running query:   RUNNING 55.4s
Running query:   RUNNING 56.1s
Running query:   RUNNING 56.7s
Running query:   RUNNING 57.4s
Running query:   RUNNING 58.0s
Running query:   RUNNING 58.7s
Running query:   RUNNING 59.3s
Running query:   RUNNING 59.9s
Running query:   RUNNING 60.6s
Running query:   RUNNING 61.3s
Running query:   RUNNING 61.9s
Running query:   RUNNING 62.6s
Running query:   RUNNING 63.2s
Running query:   RUNNING 63.8s
Running query:   RUNNING 64.5s
Running query:   RUNNING 65.1s
Running query:   RUNNING 65.8s
Running query:   RUNNING 66.4s
Running query:   RUNNING 67.0s
Running query:   RUNNING 67.7s
Running query:   RUNNING 68.3s
Running query:   RUNNING 68.9s
Running query:   RUNNING 69.6s
Running query:   RUNNING 70.2s
Running query:   RUNNING 70.9s
Running query:   RUNNING 71.5s
Running query:   RUNNING 72.1s
Running query:   RUNNING 72.8s
Running query:   RUNNING 73.4s
Running query:   RUNNING 74.1s
Running query:   RUNNING 74.7s
Running query:   RUNNING 75.4s
Running query:   RUNNING 76.0s
Running query:   RUNNING 76.7s
Running query:   RUNNING 77.3s
Running query:   RUNNING 77.9s
Running query:   RUNNING 78.6s
Running query:   RUNNING 79.3s
Running query:   RUNNING 79.9s
Running query:   RUNNING 80.5s
Running query:   RUNNING 81.2s
Running query:   RUNNING 81.8s
Running query:   RUNNING 82.5s
Running query:   RUNNING 83.1s
Running query:   RUNNING 83.8s
Running query:   RUNNING 84.4s
Running query:   RUNNING 85.1s
Running query:   RUNNING 85.7s
Running query:   RUNNING 86.3s
Running query:   RUNNING 87.0s
Running query:   RUNNING 87.6s
Running query:   RUNNING 88.3s
Running query:   RUNNING 88.9s
Running query:   RUNNING 89.5s
Running query:   RUNNING 90.2s
Running query:   RUNNING 90.8s
Running query:   RUNNING 91.4s
Running query:   RUNNING 92.1s
Running query:   RUNNING 92.7s
Running query:   RUNNING 93.4s
Running query:   RUNNING 94.0s
Running query:   RUNNING 94.6s
Running query:   RUNNING 95.3s
Running query:   RUNNING 95.9s
Running query:   RUNNING 96.5s
Running query:   RUNNING 97.2s
Running query:   RUNNING 97.8s
Running query:   RUNNING 98.5s
Running query:   RUNNING 99.1s
Running query:   RUNNING 99.8s
Running query:   RUNNING 100.4s
Running query:   RUNNING 101.0s
Running query:   RUNNING 101.7s
Running query:   RUNNING 102.3s
Running query:   RUNNING 102.9s
Running query:   RUNNING 103.6s
Running query:   RUNNING 104.2s
Running query:   RUNNING 104.8s
Running query:   RUNNING 105.5s
Running query:   RUNNING 106.1s
Running query:   RUNNING 106.8s
Running query:   RUNNING 107.4s
Running query:   RUNNING 108.0s
Running query:   RUNNING 108.7s
Running query:   RUNNING 109.4s
Running query:   RUNNING 110.0s
Running query:   RUNNING 110.6s
Running query:   RUNNING 111.3s
Running query:   RUNNING 111.9s
Running query:   RUNNING 112.5s
Running query:   RUNNING 113.2s
Running query:   RUNNING 113.8s
Running query:   RUNNING 114.5s
Running query:   RUNNING 115.1s
Running query:   RUNNING 115.8s
Running query:   RUNNING 116.4s
Running query:   RUNNING 117.1s
Running query:   RUNNING 117.7s
Running query:   RUNNING 118.3s
Running query:   RUNNING 119.0s
Running query:   RUNNING 119.6s
Running query:   RUNNING 120.2s
Running query:   RUNNING 120.9s
Running query:   RUNNING 121.5s
Running query:   RUNNING 122.2s
Running query:   RUNNING 122.8s
Running query:   RUNNING 123.5s
Running query:   RUNNING 124.2s
Running query:   RUNNING 124.9s
Running query:   RUNNING 125.5s
Running query:   RUNNING 126.2s
Running query:   RUNNING 126.8s
Running query:   RUNNING 127.5s
Running query:   RUNNING 128.1s
Running query:   RUNNING 128.8s
Running query:   RUNNING 129.5s
Running query:   RUNNING 130.1s
Running query:   RUNNING 130.7s
Running query:   RUNNING 131.4s
Running query:   RUNNING 132.0s
Running query:   RUNNING 132.6s
Running query:   RUNNING 133.3s
Running query:   RUNNING 133.9s
Running query:   RUNNING 134.6s
Running query:   RUNNING 135.2s
Running query:   RUNNING 135.9s
Running query:   RUNNING 136.5s
Running query:   RUNNING 137.2s
Running query:   RUNNING 137.8s
Running query:   RUNNING 138.5s
Running query:   RUNNING 139.1s
Running query:   RUNNING 139.7s
Running query:   RUNNING 140.4s
Running query:   RUNNING 141.0s
Running query:   RUNNING 141.7s
Running query:   RUNNING 142.3s
Running query:   RUNNING 143.0s
Running query:   RUNNING 143.6s
Running query:   RUNNING 144.2s
Running query:   RUNNING 144.9s
Running query:   RUNNING 145.5s
Running query:   RUNNING 146.2s
Running query:   RUNNING 146.8s
Running query:   RUNNING 147.4s
Running query:   RUNNING 148.1s
Running query:   RUNNING 148.7s
Running query:   RUNNING 149.4s
Running query:   RUNNING 150.0s
Running query:   RUNNING 150.6s
Running query:   RUNNING 151.3s
Running query:   RUNNING 151.9s
Running query:   RUNNING 152.6s
Running query:   RUNNING 153.2s
Running query:   RUNNING 153.9s
Running query:   RUNNING 154.5s
Running query:   RUNNING 155.1s
Running query:   RUNNING 155.8s
Running query:   RUNNING 156.4s
Running query:   RUNNING 157.1s
Running query:   RUNNING 157.7s
Running query:   RUNNING 158.4s
Running query:   RUNNING 159.0s
Running query:   RUNNING 159.6s
Running query:   RUNNING 160.3s
Running query:   RUNNING 160.9s
Running query:   RUNNING 161.6s
Running query:   RUNNING 162.2s
Running query:   RUNNING 162.8s
Running query:   RUNNING 163.5s
Running query:   RUNNING 164.1s
Running query:   RUNNING 164.7s
Running query:   RUNNING 165.4s
Running query:   RUNNING 166.0s
Running query:   RUNNING 166.7s
Running query:   RUNNING 167.3s
Running query:   RUNNING 167.9s
Running query:   RUNNING 168.6s
Running query:   RUNNING 169.2s
Running query:   RUNNING 169.9s
Running query:   RUNNING 170.5s
Running query:   RUNNING 171.2s
Running query:   RUNNING 171.8s
Running query:   RUNNING 172.4s
Running query:   RUNNING 173.1s
Running query:   RUNNING 173.7s
Running query:   RUNNING 174.3s
Running query:   RUNNING 175.0s
Running query:   RUNNING 175.6s
Running query:   RUNNING 176.3s
Running query:   RUNNING 176.9s
Running query:   RUNNING 177.5s
Running query:   RUNNING 178.2s
Running query:   RUNNING 178.8s
Running query:   RUNNING 179.4s
Running query:   RUNNING 180.1s
Running query:   RUNNING 180.7s
Running query:   RUNNING 181.4s
Running query:   RUNNING 182.0s
Running query:   RUNNING 182.6s
Running query:   RUNNING 183.3s
Running query:   RUNNING 183.9s
Running query:   RUNNING 184.6s
Running query:   RUNNING 185.2s
Running query:   RUNNING 185.8s
Running query:   RUNNING 186.5s
Running query:   RUNNING 187.1s
Running query:   RUNNING 187.8s
Running query:   RUNNING 188.4s
Running query:   RUNNING 189.0s
Running query:   RUNNING 189.7s
Running query:   RUNNING 190.3s

Retrieving data:  2.0s
```

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri May 29 01:47:24 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> chisq </th> <th> failure_reason </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBjiwqwJINj4hrDW9uC9qQE </td> <td align="right"> 687.30 </td> <td> hardy_weinberg </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBj1w6wJIIjmgdqrkczy2wE </td> <td align="right"> 880.00 </td> <td> hardy_weinberg </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNBjivrwCIKHL-ND5vLiFUw </td> <td align="right"> 881.00 </td> <td> hardy_weinberg </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBjhpI4_IPbu9N6a2YTzxwE </td> <td align="right"> 881.00 </td> <td> hardy_weinberg </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBjDhrYWIP3N0OOsj52_Dw </td> <td align="right"> 881.00 </td> <td> hardy_weinberg </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBjlookNIOrX7K-EzOHlCw </td> <td align="right"> 752.18 </td> <td> hardy_weinberg </td> </tr>
   </table>







