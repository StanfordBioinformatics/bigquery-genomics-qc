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
sampleData <- read.csv("./data/patient_info.csv")
sampleInfo <- select(sampleData, call_call_set_name=Catalog.ID, gender=Gender)
```

## Missingness Rate

Identify all variants with a missingness rate greater than a specified cutoff.


```r
cutoff = list("_CUTOFF_"="0.1")
sortAndLimit <- list("#_ORDER_BY_" = "LIMIT 1000")
result <- DisplayAndDispatchQuery("../sql/variant-level-missingness-fail.sql",
                                  project=project,
                                  replacements=c(cutoff,
                                                 queryReplacements,
                                                 sortAndLimit))
```

```
SELECT 
variant_id,
reference_name,
start,
end,
called_allele_count,
sample_count,
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
LIMIT 1000
```
Number of rows returned by this query: **1000**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Thu May 28 22:34:29 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> reference_name </th> <th> start </th> <th> end </th> <th> called_allele_count </th> <th> sample_count </th> <th> missingness_rate </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNBjGoqcBILa5z5Xskdv5GA </td> <td> chr4 </td> <td align="right"> 2740550 </td> <td align="right"> 2740553 </td> <td align="right">  10 </td> <td align="right"> 480 </td> <td align="right"> 0.98 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNBj1oqcBIKv21NDO_NqVbA </td> <td> chr4 </td> <td align="right"> 2740597 </td> <td align="right"> 2740602 </td> <td align="right">   2 </td> <td align="right"> 480 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTQYnpb7GCDKq-ipvN7F3js </td> <td> chr14 </td> <td align="right"> 52349726 </td> <td align="right"> 52349727 </td> <td align="right">   4 </td> <td align="right"> 480 </td> <td align="right"> 0.99 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNhjV5vEbIJ6ApI3vwPyEngE </td> <td> chr6 </td> <td align="right"> 58487637 </td> <td align="right"> 58487638 </td> <td align="right"> 276 </td> <td align="right"> 480 </td> <td align="right"> 0.43 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNhi56PEbIK3T6Lrv562apwE </td> <td> chr6 </td> <td align="right"> 58487865 </td> <td align="right"> 58487866 </td> <td align="right"> 194 </td> <td align="right"> 480 </td> <td align="right"> 0.60 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNhjc6PEbIObnqPiO4MSaFQ </td> <td> chr6 </td> <td align="right"> 58487900 </td> <td align="right"> 58487901 </td> <td align="right">  74 </td> <td align="right"> 480 </td> <td align="right"> 0.85 </td> </tr>
   </table>

## Blacklisted Variants


```r
query <- "../sql/blacklisted-variants.sql"
sortAndLimit <- list("#_ORDER_BY_" = "LIMIT 1000")
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(queryReplacements, sortAndLimit))
```

```
SELECT
  seq.variant_id AS variant_id,
  seq.reference_name AS reference_name,
  seq.start AS start,
  seq.end AS end,
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
LIMIT 1000
```

Number of rows returned by this query: **1000**.

First few results
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Thu May 28 22:34:31 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> reference_name </th> <th> start </th> <th> end </th> <th> Artifact_Type </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTEYi72EGCDEwbzL1Obo43U </td> <td> chr11 </td> <td align="right"> 50405003 </td> <td align="right"> 50405004 </td> <td> centromeric_repeat </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTEY2r2EGCChw63qz72392Y </td> <td> chr11 </td> <td align="right"> 50405082 </td> <td align="right"> 50405083 </td> <td> centromeric_repeat </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTEYir6EGCCJ8eya08SY2JwB </td> <td> chr11 </td> <td align="right"> 50405130 </td> <td align="right"> 50405131 </td> <td> centromeric_repeat </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTEY8r6EGCDJwKyomc7o6tAB </td> <td> chr11 </td> <td align="right"> 50405234 </td> <td align="right"> 50405235 </td> <td> centromeric_repeat </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTEY-76EGCCHy4Do4LOrmXw </td> <td> chr11 </td> <td align="right"> 50405243 </td> <td align="right"> 50405244 </td> <td> centromeric_repeat </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTEYtL-EGCCfy5q83uuftskB </td> <td> chr11 </td> <td align="right"> 50405300 </td> <td align="right"> 50405301 </td> <td> centromeric_repeat </td> </tr>
   </table>

## Heterozygous Haplotype
For each variant within the X and Y chromosome, identify heterozygous variants in male genomes.


```r
sortAndLimit <- "LIMIT 1000"
result <- DisplayAndDispatchQuery("../sql/sex-chromosome-heterozygous-haplotypes.sql",
                                  project=project,
                                  replacements=c("#_ORDER_BY_"=sortAndLimit,
                                                 queryReplacements))
```

```
SELECT
  variant_id,
  reference_name,
  start,
  end,
  reference_bases,
  sample_id,
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
LIMIT 1000
```
Number of rows returned by this query: **1000**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Thu May 28 22:34:34 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> reference_name </th> <th> start </th> <th> end </th> <th> reference_bases </th> <th> sample_id </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWRis4rgGINPtubrox5zX6wE </td> <td> chrY </td> <td align="right"> 13513004 </td> <td align="right"> 13513005 </td> <td> A </td> <td> LP6005243-DNA_E03 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWRis4rgGINPtubrox5zX6wE </td> <td> chrY </td> <td align="right"> 13513004 </td> <td align="right"> 13513005 </td> <td> A </td> <td> LP6005243-DNA_E06 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWRis4rgGINPtubrox5zX6wE </td> <td> chrY </td> <td align="right"> 13513004 </td> <td align="right"> 13513005 </td> <td> A </td> <td> LP6005051-DNA_D10 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWRis4rgGINPtubrox5zX6wE </td> <td> chrY </td> <td align="right"> 13513004 </td> <td align="right"> 13513005 </td> <td> A </td> <td> LP6005051-DNA_F10 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWRis4rgGINPtubrox5zX6wE </td> <td> chrY </td> <td align="right"> 13513004 </td> <td align="right"> 13513005 </td> <td> A </td> <td> LP6005051-DNA_C12 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWRis4rgGINPtubrox5zX6wE </td> <td> chrY </td> <td align="right"> 13513004 </td> <td align="right"> 13513005 </td> <td> A </td> <td> LP6005692-DNA_B08 </td> </tr>
   </table>

## Ti/Tv By Genomic Window

```r
query <- "../sql/titv-by-genomic-window-fail.sql"
sortAndLimit <- list("#_ORDER_BY_" = "LIMIT 1000")
max <- 3.0
min <- 1.5
cutoffs <- list("_MAX_" = max, "_MIN_" = min)
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(queryReplacements, cutoffs, sortAndLimit))
```

```
SELECT
var.variant_id AS variant_id,
var.reference_name AS reference_name,
win.window_start AS window_start,
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
reference_name,
window_start,
LIMIT 1000
```

Number of rows returned by this query: **1000**.

First few results
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Thu May 28 22:34:36 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> reference_name </th> <th> window_start </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTcY-OCdGCD6p5PL_72AvNkB </td> <td> chr17 </td> <td align="right"> 50800000 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTcYq-KdGCCh9KfbhPm_4yE </td> <td> chr17 </td> <td align="right"> 50800000 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTcYv-OdGCDk-b_0xOmlqYUB </td> <td> chr17 </td> <td align="right"> 50800000 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTcY6-OdGCDppvaNjJi93F4 </td> <td> chr17 </td> <td align="right"> 50800000 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTcYreSdGCCPoOTNtZTTz3I </td> <td> chr17 </td> <td align="right"> 50800000 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTcYxuSdGCDZlpGl5MuN0VE </td> <td> chr17 </td> <td align="right"> 50800000 </td> </tr>
   </table>

## Ti/Tv By Depth

```r
query <- "../sql/titv-by-depth-fail.sql"
max <- 3
min <- 1.5
cutoffs <- list("_MAX_" = max, "_MIN_" = min)
sortAndLimit <- list("#_ORDER_BY_" = "LIMIT 1000")
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(queryReplacements, cutoffs, sortAndLimit))
```

```
SELECT
*
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
OR titv_ratio < 1.5
LIMIT 1000
```

Number of rows returned by this query: **1000**.

First few results
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Thu May 28 22:34:39 2015 -->
<table border=1>
<tr> <th> sample_id </th> <th> titv_ratio </th> <th> depth </th>  </tr>
  <tr> <td> LP6005051-DNA_A11 </td> <td align="right"> 1.49 </td> <td align="right">   1 </td> </tr>
  <tr> <td> LP6005038-DNA_C06 </td> <td align="right"> 1.47 </td> <td align="right">   6 </td> </tr>
  <tr> <td> LP6005038-DNA_H04 </td> <td align="right"> 1.49 </td> <td align="right">   7 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 1.42 </td> <td align="right">   4 </td> </tr>
  <tr> <td> LP6005038-DNA_E05 </td> <td align="right"> 1.49 </td> <td align="right">   5 </td> </tr>
  <tr> <td> LP6005038-DNA_E01 </td> <td align="right"> 1.41 </td> <td align="right">   3 </td> </tr>
   </table>

Generate query replacement to retreive variant ids for each sample/depth pair that is outside our cutoffs.

```r
# Format
#  (call.call_set_name = 'LP6005692-DNA_B08' AND call.DP = 112)
result$query = paste("(call.call_set_name = '", result$sample_id, "' AND call.DP = ", result$depth, ")", sep='')
sampleReplacement = list("_SAMPLE_DEPTH_" = paste(result$query, collapse = ' OR '))
```

Run query to retreive all variant ids for variants failing ti/tv by depth

```r
query <- "../sql/titv-by-depth-fail-variants.sql"
sortAndLimit <- list("#_ORDER_BY_" = "LIMIT 1000")
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(queryReplacements, sampleReplacement, sortAndLimit))
```

```
SELECT
  variant_id,
  call.call_set_name,
FROM
  [va_aaa_pilot_data.all_genomes_expanded_vcfs_java3]
WHERE
  # A list of samples and the depths of interest
  # Example:
  #   (call.call_set_name = 'LP6005692-DNA_B08'
  #      AND call.DP = 112)
  #    OR (call.call_set_name = 'LP6005692-DNA_C10'
  #      AND call.DP = 125)
  (call.call_set_name = 'LP6005051-DNA_A11' AND call.DP = 1) OR (call.call_set_name = 'LP6005038-DNA_C06' AND call.DP = 6) OR (call.call_set_name = 'LP6005038-DNA_H04' AND call.DP = 7) OR (call.call_set_name = 'LP6005038-DNA_A03' AND call.DP = 4) OR (call.call_set_name = 'LP6005038-DNA_E05' AND call.DP = 5) OR (call.call_set_name = 'LP6005038-DNA_E01' AND call.DP = 3) OR (call.call_set_name = 'LP6005038-DNA_H06' AND call.DP = 5) OR (call.call_set_name = 'LP6005144-DNA_C12' AND call.DP = 5) OR (call.call_set_name = 'LP6005693-DNA_E01' AND call.DP = 94) OR (call.call_set_name = 'LP6005038-DNA_C05' AND call.DP = 6) OR (call.call_set_name = 'LP6005144-DNA_A06' AND call.DP = 4) OR (call.call_set_name = 'LP6005038-DNA_E06' AND call.DP = 4) OR (call.call_set_name = 'LP6005051-DNA_C12' AND call.DP = 2) OR (call.call_set_name = 'LP6005144-DNA_C11' AND call.DP = 4) OR (call.call_set_name = 'LP6005051-DNA_F11' AND call.DP = 2) OR (call.call_set_name = 'LP6005051-DNA_E11' AND call.DP = 2) OR (call.call_set_name = 'LP6005051-DNA_F04' AND call.DP = 2) OR (call.call_set_name = 'LP6005038-DNA_D05' AND call.DP = 4) OR (call.call_set_name = 'LP6005051-DNA_H04' AND call.DP = 3) OR (call.call_set_name = 'LP6005144-DNA_D05' AND call.DP = 2) OR (call.call_set_name = 'LP6005144-DNA_A06' AND call.DP = 3) OR (call.call_set_name = 'LP6005051-DNA_F04' AND call.DP = 75) OR (call.call_set_name = 'LP6005692-DNA_A06' AND call.DP = 93) OR (call.call_set_name = 'LP6005144-DNA_A08' AND call.DP = 68) OR (call.call_set_name = 'LP6005144-DNA_C09' AND call.DP = 77) OR (call.call_set_name = 'LP6005692-DNA_B12' AND call.DP = 100) OR (call.call_set_name = 'LP6005051-DNA_E10' AND call.DP = 86) OR (call.call_set_name = 'LP6005693-DNA_B01' AND call.DP = 115) OR (call.call_set_name = 'LP6005692-DNA_D06' AND call.DP = 98) OR (call.call_set_name = 'LP6005144-DNA_E11' AND call.DP = 95) OR (call.call_set_name = 'LP6005051-DNA_A03' AND call.DP = 104) OR (call.call_set_name = 'LP6005144-DNA_H03' AND call.DP = 94) OR (call.call_set_name = 'LP6005051-DNA_B04' AND call.DP = 90) OR (call.call_set_name = 'LP6005144-DNA_H04' AND call.DP = 81) OR (call.call_set_name = 'LP6005144-DNA_H05' AND call.DP = 102) OR (call.call_set_name = 'LP6005051-DNA_G06' AND call.DP = 84) OR (call.call_set_name = 'LP6005692-DNA_G02' AND call.DP = 93) OR (call.call_set_name = 'LP6005243-DNA_C08' AND call.DP = 95) OR (call.call_set_name = 'LP6005038-DNA_D06' AND call.DP = 85) OR (call.call_set_name = 'LP6005243-DNA_A02' AND call.DP = 86) OR (call.call_set_name = 'LP6005038-DNA_D04' AND call.DP = 105) OR (call.call_set_name = 'LP6005038-DNA_B08' AND call.DP = 126) OR (call.call_set_name = 'LP6005692-DNA_F03' AND call.DP = 78) OR (call.call_set_name = 'LP6005243-DNA_D02' AND call.DP = 85) OR (call.call_set_name = 'LP6005038-DNA_B05' AND call.DP = 83) OR (call.call_set_name = 'LP6005144-DNA_C08' AND call.DP = 63) OR (call.call_set_name = 'LP6005243-DNA_D01' AND call.DP = 93) OR (call.call_set_name = 'LP6005144-DNA_E05' AND call.DP = 96) OR (call.call_set_name = 'LP6005144-DNA_F08' AND call.DP = 67) OR (call.call_set_name = 'LP6005051-DNA_C06' AND call.DP = 77) OR (call.call_set_name = 'LP6005243-DNA_H07' AND call.DP = 87) OR (call.call_set_name = 'LP6005692-DNA_E01' AND call.DP = 82) OR (call.call_set_name = 'LP6005144-DNA_H09' AND call.DP = 81) OR (call.call_set_name = 'LP6005038-DNA_H08' AND call.DP = 85) OR (call.call_set_name = 'LP6005051-DNA_A12' AND call.DP = 84) OR (call.call_set_name = 'LP6005243-DNA_C10' AND call.DP = 83) OR (call.call_set_name = 'LP6005051-DNA_A07' AND call.DP = 96) OR (call.call_set_name = 'LP6005038-DNA_E07' AND call.DP = 84) OR (call.call_set_name = 'LP6005692-DNA_E07' AND call.DP = 81) OR (call.call_set_name = 'LP6005692-DNA_A05' AND call.DP = 73) OR (call.call_set_name = 'LP6005051-DNA_A03' AND call.DP = 76) OR (call.call_set_name = 'LP6005144-DNA_C11' AND call.DP = 7) OR (call.call_set_name = 'LP6005038-DNA_C03' AND call.DP = 2) OR (call.call_set_name = 'LP6005038-DNA_D06' AND call.DP = 4) OR (call.call_set_name = 'LP6005144-DNA_C09' AND call.DP = 5) OR (call.call_set_name = 'LP6005051-DNA_H11' AND call.DP = 3) OR (call.call_set_name = 'LP6005051-DNA_B11' AND call.DP = 2) OR (call.call_set_name = 'LP6005051-DNA_A06' AND call.DP = 1) OR (call.call_set_name = 'LP6005051-DNA_B04' AND call.DP = 4) OR (call.call_set_name = 'LP6005038-DNA_E01' AND call.DP = 5) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 250) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 244) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 242) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 248) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 245) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 249) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 245) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 244) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 247) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 249) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 219) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 170) OR (call.call_set_name = 'LP6005051-DNA_F03' AND call.DP = 102) OR (call.call_set_name = 'LP6005144-DNA_D10' AND call.DP = 105) OR (call.call_set_name = 'LP6005243-DNA_C02' AND call.DP = 92) OR (call.call_set_name = 'LP6005051-DNA_G03' AND call.DP = 86) OR (call.call_set_name = 'LP6005243-DNA_G04' AND call.DP = 106) OR (call.call_set_name = 'LP6005051-DNA_F05' AND call.DP = 87) OR (call.call_set_name = 'LP6005144-DNA_H11' AND call.DP = 97) OR (call.call_set_name = 'LP6005692-DNA_G05' AND call.DP = 101) OR (call.call_set_name = 'LP6005692-DNA_G06' AND call.DP = 101) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 161) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 207) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 229) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 250) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 193) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 204) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 212) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 224) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 187) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 209) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 105) OR (call.call_set_name = 'LP6005693-DNA_B01' AND call.DP = 121) OR (call.call_set_name = 'LP6005144-DNA_D01' AND call.DP = 105) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 155) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 159) OR (call.call_set_name = 'LP6005144-DNA_F12' AND call.DP = 110) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 97) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 168) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 187) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 111) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 223) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 139) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 242) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 153) OR (call.call_set_name = 'LP6005243-DNA_B05' AND call.DP = 97) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 231) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 125) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 228) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 126) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 182) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 169) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 126) OR (call.call_set_name = 'LP6005051-DNA_E08' AND call.DP = 74) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 127) OR (call.call_set_name = 'LP6005038-DNA_D07' AND call.DP = 109) OR (call.call_set_name = 'LP6005038-DNA_E02' AND call.DP = 128) OR (call.call_set_name = 'LP6005038-DNA_E03' AND call.DP = 176) OR (call.call_set_name = 'LP6005144-DNA_D10' AND call.DP = 85) OR (call.call_set_name = 'LP6005038-DNA_E01' AND call.DP = 77) OR (call.call_set_name = 'LP6005692-DNA_G10' AND call.DP = 99) OR (call.call_set_name = 'LP6005051-DNA_F01' AND call.DP = 77) OR (call.call_set_name = 'LP6005144-DNA_B01' AND call.DP = 71) OR (call.call_set_name = 'LP6005144-DNA_C12' AND call.DP = 1) OR (call.call_set_name = 'LP6005038-DNA_G04' AND call.DP = 104) OR (call.call_set_name = 'LP6005051-DNA_E04' AND call.DP = 77) OR (call.call_set_name = 'LP6005038-DNA_E03' AND call.DP = 101) OR (call.call_set_name = 'LP6005051-DNA_G03' AND call.DP = 94) OR (call.call_set_name = 'LP6005144-DNA_G10' AND call.DP = 108) OR (call.call_set_name = 'LP6005038-DNA_E03' AND call.DP = 5) OR (call.call_set_name = 'LP6005038-DNA_G06' AND call.DP = 5) OR (call.call_set_name = 'LP6005038-DNA_A03' AND call.DP = 5) OR (call.call_set_name = 'LP6005144-DNA_C11' AND call.DP = 5) OR (call.call_set_name = 'LP6005038-DNA_B03' AND call.DP = 5) OR (call.call_set_name = 'LP6005038-DNA_F06' AND call.DP = 4) OR (call.call_set_name = 'LP6005051-DNA_B04' AND call.DP = 3) OR (call.call_set_name = 'LP6005144-DNA_H05' AND call.DP = 5) OR (call.call_set_name = 'LP6005144-DNA_H04' AND call.DP = 4) OR (call.call_set_name = 'LP6005051-DNA_A04' AND call.DP = 2) OR (call.call_set_name = 'LP6005038-DNA_H03' AND call.DP = 1) OR (call.call_set_name = 'LP6005051-DNA_D11' AND call.DP = 1) OR (call.call_set_name = 'LP6005243-DNA_F08' AND call.DP = 1) OR (call.call_set_name = 'LP6005051-DNA_E05' AND call.DP = 77) OR (call.call_set_name = 'LP6005051-DNA_D03' AND call.DP = 74) OR (call.call_set_name = 'LP6005243-DNA_D02' AND call.DP = 82) OR (call.call_set_name = 'LP6005051-DNA_G12' AND call.DP = 75) OR (call.call_set_name = 'LP6005051-DNA_B12' AND call.DP = 73) OR (call.call_set_name = 'LP6005051-DNA_E06' AND call.DP = 103) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 104) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 103) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 80) OR (call.call_set_name = 'LP6005243-DNA_C02' AND call.DP = 83) OR (call.call_set_name = 'LP6005038-DNA_F05' AND call.DP = 66) OR (call.call_set_name = 'LP6005051-DNA_E11' AND call.DP = 67) OR (call.call_set_name = 'LP6005038-DNA_A03' AND call.DP = 74) OR (call.call_set_name = 'LP6005038-DNA_F11' AND call.DP = 71) OR (call.call_set_name = 'LP6005038-DNA_F01' AND call.DP = 84) OR (call.call_set_name = 'LP6005038-DNA_D08' AND call.DP = 92) OR (call.call_set_name = 'LP6005038-DNA_D09' AND call.DP = 70) OR (call.call_set_name = 'LP6005243-DNA_A06' AND call.DP = 92) OR (call.call_set_name = 'LP6005692-DNA_C12' AND call.DP = 88) OR (call.call_set_name = 'LP6005144-DNA_A09' AND call.DP = 73) OR (call.call_set_name = 'LP6005243-DNA_G08' AND call.DP = 1) OR (call.call_set_name = 'LP6005243-DNA_C10' AND call.DP = 1) OR (call.call_set_name = 'LP6005692-DNA_H06' AND call.DP = 1) OR (call.call_set_name = 'LP6005038-DNA_G01' AND call.DP = 5) OR (call.call_set_name = 'LP6005144-DNA_B09' AND call.DP = 4) OR (call.call_set_name = 'LP6005051-DNA_B04' AND call.DP = 2) OR (call.call_set_name = 'LP6005051-DNA_D12' AND call.DP = 1) OR (call.call_set_name = 'LP6005038-DNA_B03' AND call.DP = 4) OR (call.call_set_name = 'LP6005038-DNA_G04' AND call.DP = 6) OR (call.call_set_name = 'LP6005038-DNA_F03' AND call.DP = 4) OR (call.call_set_name = 'LP6005144-DNA_H01' AND call.DP = 72) OR (call.call_set_name = 'LP6005038-DNA_G03' AND call.DP = 5) OR (call.call_set_name = 'LP6005144-DNA_B12' AND call.DP = 5) OR (call.call_set_name = 'LP6005038-DNA_B03' AND call.DP = 3) OR (call.call_set_name = 'LP6005144-DNA_E07' AND call.DP = 5) OR (call.call_set_name = 'LP6005144-DNA_B12' AND call.DP = 2) OR (call.call_set_name = 'LP6005051-DNA_A10' AND call.DP = 1) OR (call.call_set_name = 'LP6005144-DNA_B04' AND call.DP = 2) OR (call.call_set_name = 'LP6005038-DNA_B03' AND call.DP = 1) OR (call.call_set_name = 'LP6005243-DNA_A04' AND call.DP = 1) OR (call.call_set_name = 'LP6005051-DNA_C04' AND call.DP = 1) OR (call.call_set_name = 'LP6005051-DNA_B11' AND call.DP = 1) OR (call.call_set_name = 'LP6005051-DNA_C06' AND call.DP = 2) OR (call.call_set_name = 'LP6005051-DNA_H11' AND call.DP = 1) OR (call.call_set_name = 'LP6005243-DNA_G03' AND call.DP = 1) OR (call.call_set_name = 'LP6005038-DNA_G03' AND call.DP = 2) OR (call.call_set_name = 'LP6005051-DNA_D10' AND call.DP = 1) OR (call.call_set_name = 'LP6005038-DNA_C03' AND call.DP = 1) OR (call.call_set_name = 'LP6005051-DNA_B06' AND call.DP = 1) OR (call.call_set_name = 'LP6005051-DNA_H10' AND call.DP = 1) OR (call.call_set_name = 'LP6005038-DNA_G06' AND call.DP = 4) OR (call.call_set_name = 'LP6005051-DNA_A12' AND call.DP = 2) OR (call.call_set_name = 'LP6005144-DNA_C11' AND call.DP = 3) OR (call.call_set_name = 'LP6005038-DNA_H06' AND call.DP = 4) OR (call.call_set_name = 'LP6005051-DNA_B04' AND call.DP = 1) OR (call.call_set_name = 'LP6005051-DNA_C11' AND call.DP = 2) OR (call.call_set_name = 'LP6005051-DNA_D06' AND call.DP = 3) OR (call.call_set_name = 'LP6005144-DNA_B12' AND call.DP = 3) OR (call.call_set_name = 'LP6005144-DNA_B09' AND call.DP = 2) OR (call.call_set_name = 'LP6005038-DNA_H04' AND call.DP = 5) OR (call.call_set_name = 'LP6005038-DNA_H03' AND call.DP = 2) OR (call.call_set_name = 'LP6005038-DNA_E03' AND call.DP = 2) OR (call.call_set_name = 'LP6005051-DNA_C04' AND call.DP = 3) OR (call.call_set_name = 'LP6005144-DNA_H08' AND call.DP = 4) OR (call.call_set_name = 'LP6005243-DNA_F09' AND call.DP = 107) OR (call.call_set_name = 'LP6005038-DNA_B10' AND call.DP = 93) OR (call.call_set_name = 'LP6005692-DNA_G02' AND call.DP = 85) OR (call.call_set_name = 'LP6005038-DNA_H12' AND call.DP = 116) OR (call.call_set_name = 'LP6005144-DNA_H09' AND call.DP = 119) OR (call.call_set_name = 'LP6005051-DNA_G05' AND call.DP = 74) OR (call.call_set_name = 'LP6005243-DNA_D11' AND call.DP = 92) OR (call.call_set_name = 'LP6005692-DNA_A05' AND call.DP = 91) OR (call.call_set_name = 'LP6005144-DNA_D08' AND call.DP = 115) OR (call.call_set_name = 'LP6005051-DNA_G04' AND call.DP = 71) OR (call.call_set_name = 'LP6005051-DNA_D05' AND call.DP = 123) OR (call.call_set_name = 'LP6005038-DNA_A10' AND call.DP = 79) OR (call.call_set_name = 'LP6005243-DNA_E04' AND call.DP = 133) OR (call.call_set_name = 'LP6005243-DNA_F10' AND call.DP = 107) OR (call.call_set_name = 'LP6005051-DNA_B10' AND call.DP = 97) OR (call.call_set_name = 'LP6005038-DNA_H09' AND call.DP = 107) OR (call.call_set_name = 'LP6005243-DNA_F03' AND call.DP = 117) OR (call.call_set_name = 'LP6005051-DNA_B09' AND call.DP = 106) OR (call.call_set_name = 'LP6005038-DNA_E05' AND call.DP = 114) OR (call.call_set_name = 'LP6005693-DNA_B01' AND call.DP = 154) OR (call.call_set_name = 'LP6005144-DNA_E11' AND call.DP = 118) OR (call.call_set_name = 'LP6005243-DNA_C08' AND call.DP = 119) OR (call.call_set_name = 'LP6005051-DNA_E01' AND call.DP = 108) OR (call.call_set_name = 'LP6005144-DNA_D12' AND call.DP = 108) OR (call.call_set_name = 'LP6005144-DNA_G11' AND call.DP = 82) OR (call.call_set_name = 'LP6005693-DNA_A01' AND call.DP = 155) OR (call.call_set_name = 'LP6005144-DNA_G12' AND call.DP = 122) OR (call.call_set_name = 'LP6005051-DNA_F10' AND call.DP = 116) OR (call.call_set_name = 'LP6005144-DNA_H06' AND call.DP = 96) OR (call.call_set_name = 'LP6005038-DNA_A12' AND call.DP = 110) OR (call.call_set_name = 'LP6005051-DNA_C03' AND call.DP = 117) OR (call.call_set_name = 'LP6005051-DNA_B12' AND call.DP = 72) OR (call.call_set_name = 'LP6005038-DNA_C07' AND call.DP = 84) OR (call.call_set_name = 'LP6005692-DNA_A06' AND call.DP = 151) OR (call.call_set_name = 'LP6005243-DNA_G01' AND call.DP = 148) OR (call.call_set_name = 'LP6005243-DNA_D04' AND call.DP = 111) OR (call.call_set_name = 'LP6005243-DNA_H07' AND call.DP = 112) OR (call.call_set_name = 'LP6005144-DNA_D11' AND call.DP = 145) OR (call.call_set_name = 'LP6005243-DNA_E08' AND call.DP = 91) OR (call.call_set_name = 'LP6005051-DNA_G09' AND call.DP = 88) OR (call.call_set_name = 'LP6005243-DNA_F07' AND call.DP = 102) OR (call.call_set_name = 'LP6005692-DNA_D10' AND call.DP = 143) OR (call.call_set_name = 'LP6005038-DNA_G12' AND call.DP = 152) OR (call.call_set_name = 'LP6005038-DNA_D09' AND call.DP = 103) OR (call.call_set_name = 'LP6005144-DNA_C09' AND call.DP = 107) OR (call.call_set_name = 'LP6005038-DNA_H12' AND call.DP = 107) OR (call.call_set_name = 'LP6005051-DNA_F12' AND call.DP = 138) OR (call.call_set_name = 'LP6005038-DNA_F07' AND call.DP = 115) OR (call.call_set_name = 'LP6005692-DNA_F03' AND call.DP = 117) OR (call.call_set_name = 'LP6005692-DNA_E06' AND call.DP = 155) OR (call.call_set_name = 'LP6005144-DNA_E05' AND call.DP = 104) OR (call.call_set_name = 'LP6005144-DNA_H03' AND call.DP = 134) OR (call.call_set_name = 'LP6005038-DNA_A05' AND call.DP = 108) OR (call.call_set_name = 'LP6005692-DNA_D01' AND call.DP = 151) OR (call.call_set_name = 'LP6005144-DNA_E11' AND call.DP = 141) OR (call.call_set_name = 'LP6005038-DNA_E12' AND call.DP = 139) OR (call.call_set_name = 'LP6005144-DNA_B04' AND call.DP = 76) OR (call.call_set_name = 'LP6005038-DNA_G10' AND call.DP = 113) OR (call.call_set_name = 'LP6005038-DNA_H04' AND call.DP = 116) OR (call.call_set_name = 'LP6005144-DNA_F07' AND call.DP = 162) OR (call.call_set_name = 'LP6005038-DNA_C08' AND call.DP = 157) OR (call.call_set_name = 'LP6005144-DNA_C07' AND call.DP = 120) OR (call.call_set_name = 'LP6005038-DNA_B06' AND call.DP = 149) OR (call.call_set_name = 'LP6005692-DNA_F01' AND call.DP = 178) OR (call.call_set_name = 'LP6005051-DNA_A06' AND call.DP = 141) OR (call.call_set_name = 'LP6005692-DNA_F06' AND call.DP = 118) OR (call.call_set_name = 'LP6005144-DNA_B09' AND call.DP = 104) OR (call.call_set_name = 'LP6005051-DNA_E07' AND call.DP = 108) OR (call.call_set_name = 'LP6005038-DNA_D07' AND call.DP = 86) OR (call.call_set_name = 'LP6005692-DNA_D07' AND call.DP = 133) OR (call.call_set_name = 'LP6005144-DNA_D05' AND call.DP = 165) OR (call.call_set_name = 'LP6005144-DNA_D07' AND call.DP = 135) OR (call.call_set_name = 'LP6005144-DNA_G11' AND call.DP = 80) OR (call.call_set_name = 'LP6005144-DNA_C02' AND call.DP = 99) OR (call.call_set_name = 'LP6005693-DNA_C01' AND call.DP = 137) OR (call.call_set_name = 'LP6005243-DNA_E06' AND call.DP = 99) OR (call.call_set_name = 'LP6005051-DNA_F08' AND call.DP = 97) OR (call.call_set_name = 'LP6005051-DNA_C07' AND call.DP = 105) OR (call.call_set_name = 'LP6005038-DNA_B04' AND call.DP = 179) OR (call.call_set_name = 'LP6005038-DNA_A10' AND call.DP = 77) OR (call.call_set_name = 'LP6005051-DNA_C06' AND call.DP = 159) OR (call.call_set_name = 'LP6005038-DNA_G10' AND call.DP = 124) OR (call.call_set_name = 'LP6005051-DNA_E09' AND call.DP = 104) OR (call.call_set_name = 'LP6005038-DNA_A06' AND call.DP = 175) OR (call.call_set_name = 'LP6005038-DNA_H03' AND call.DP = 129) OR (call.call_set_name = 'LP6005051-DNA_A11' AND call.DP = 124) OR (call.call_set_name = 'LP6005243-DNA_H05' AND call.DP = 118) OR (call.call_set_name = 'LP6005692-DNA_F04' AND call.DP = 127) OR (call.call_set_name = 'LP6005038-DNA_A11' AND call.DP = 115) OR (call.call_set_name = 'LP6005243-DNA_E05' AND call.DP = 151) OR (call.call_set_name = 'LP6005051-DNA_B11' AND call.DP = 139) OR (call.call_set_name = 'LP6005038-DNA_F03' AND call.DP = 82) OR (call.call_set_name = 'LP6005243-DNA_C10' AND call.DP = 137) OR (call.call_set_name = 'LP6005243-DNA_F08' AND call.DP = 114) OR (call.call_set_name = 'LP6005693-DNA_A02' AND call.DP = 103) OR (call.call_set_name = 'LP6005144-DNA_C10' AND call.DP = 110) OR (call.call_set_name = 'LP6005692-DNA_E08' AND call.DP = 125) OR (call.call_set_name = 'LP6005692-DNA_C07' AND call.DP = 99) OR (call.call_set_name = 'LP6005692-DNA_E07' AND call.DP = 99) OR (call.call_set_name = 'LP6005144-DNA_B10' AND call.DP = 137) OR (call.call_set_name = 'LP6005051-DNA_D10' AND call.DP = 72) OR (call.call_set_name = 'LP6005144-DNA_B12' AND call.DP = 157) OR (call.call_set_name = 'LP6005243-DNA_A05' AND call.DP = 98) OR (call.call_set_name = 'LP6005051-DNA_H10' AND call.DP = 103) OR (call.call_set_name = 'LP6005038-DNA_G01' AND call.DP = 90) OR (call.call_set_name = 'LP6005692-DNA_B09' AND call.DP = 90) OR (call.call_set_name = 'LP6005144-DNA_C12' AND call.DP = 103) OR (call.call_set_name = 'LP6005051-DNA_B11' AND call.DP = 134) OR (call.call_set_name = 'LP6005243-DNA_F06' AND call.DP = 151) OR (call.call_set_name = 'LP6005692-DNA_C12' AND call.DP = 141) OR (call.call_set_name = 'LP6005051-DNA_H05' AND call.DP = 91) OR (call.call_set_name = 'LP6005051-DNA_E08' AND call.DP = 142) OR (call.call_set_name = 'LP6005692-DNA_E06' AND call.DP = 142) OR (call.call_set_name = 'LP6005038-DNA_F03' AND call.DP = 79) OR (call.call_set_name = 'LP6005243-DNA_G03' AND call.DP = 153) OR (call.call_set_name = 'LP6005243-DNA_H02' AND call.DP = 124) OR (call.call_set_name = 'LP6005243-DNA_A07' AND call.DP = 108) OR (call.call_set_name = 'LP6005051-DNA_E07' AND call.DP = 109) OR (call.call_set_name = 'LP6005051-DNA_B01' AND call.DP = 97) OR (call.call_set_name = 'LP6005243-DNA_H10' AND call.DP = 100) OR (call.call_set_name = 'LP6005243-DNA_C11' AND call.DP = 81) OR (call.call_set_name = 'LP6005038-DNA_E07' AND call.DP = 100) OR (call.call_set_name = 'LP6005051-DNA_A06' AND call.DP = 137) OR (call.call_set_name = 'LP6005051-DNA_A08' AND call.DP = 105) OR (call.call_set_name = 'LP6005243-DNA_C06' AND call.DP = 161) OR (call.call_set_name = 'LP6005144-DNA_H08' AND call.DP = 121) OR (call.call_set_name = 'LP6005038-DNA_E11' AND call.DP = 110) OR (call.call_set_name = 'LP6005243-DNA_A11' AND call.DP = 143) OR (call.call_set_name = 'LP6005144-DNA_H05' AND call.DP = 152) OR (call.call_set_name = 'LP6005051-DNA_E03' AND call.DP = 164) OR (call.call_set_name = 'LP6005144-DNA_D01' AND call.DP = 168) OR (call.call_set_name = 'LP6005144-DNA_G05' AND call.DP = 166) OR (call.call_set_name = 'LP6005038-DNA_A09' AND call.DP = 146) OR (call.call_set_name = 'LP6005051-DNA_B09' AND call.DP = 109) OR (call.call_set_name = 'LP6005692-DNA_C09' AND call.DP = 103) OR (call.call_set_name = 'LP6005243-DNA_A04' AND call.DP = 125) OR (call.call_set_name = 'LP6005051-DNA_H03' AND call.DP = 190) OR (call.call_set_name = 'LP6005051-DNA_G05' AND call.DP = 79) OR (call.call_set_name = 'LP6005038-DNA_C07' AND call.DP = 99) OR (call.call_set_name = 'LP6005692-DNA_H07' AND call.DP = 121) OR (call.call_set_name = 'LP6005144-DNA_A03' AND call.DP = 130) OR (call.call_set_name = 'LP6005038-DNA_D07' AND call.DP = 93) OR (call.call_set_name = 'LP6005692-DNA_C09' AND call.DP = 121) OR (call.call_set_name = 'LP6005243-DNA_D10' AND call.DP = 131) OR (call.call_set_name = 'LP6005038-DNA_H12' AND call.DP = 132) OR (call.call_set_name = 'LP6005144-DNA_B10' AND call.DP = 146) OR (call.call_set_name = 'LP6005038-DNA_E04' AND call.DP = 155) OR (call.call_set_name = 'LP6005038-DNA_C09' AND call.DP = 170) OR (call.call_set_name = 'LP6005692-DNA_E07' AND call.DP = 135) OR (call.call_set_name = 'LP6005051-DNA_F04' AND call.DP = 119) OR (call.call_set_name = 'LP6005038-DNA_C09' AND call.DP = 176) OR (call.call_set_name = 'LP6005692-DNA_E09' AND call.DP = 87) OR (call.call_set_name = 'LP6005692-DNA_E08' AND call.DP = 121) OR (call.call_set_name = 'LP6005144-DNA_B04' AND call.DP = 77) OR (call.call_set_name = 'LP6005243-DNA_B07' AND call.DP = 75) OR (call.call_set_name = 'LP6005692-DNA_H09' AND call.DP = 109) OR (call.call_set_name = 'LP6005243-DNA_E01' AND call.DP = 127) OR (call.call_set_name = 'LP6005144-DNA_B04' AND call.DP = 70) OR (call.call_set_name = 'LP6005051-DNA_D06' AND call.DP = 65) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 106) OR (call.call_set_name = 'LP6005051-DNA_C04' AND call.DP = 110) OR (call.call_set_name = 'LP6005051-DNA_B08' AND call.DP = 124) OR (call.call_set_name = 'LP6005038-DNA_H09' AND call.DP = 111) OR (call.call_set_name = 'LP6005243-DNA_C05' AND call.DP = 94) OR (call.call_set_name = 'LP6005243-DNA_F10' AND call.DP = 121) OR (call.call_set_name = 'LP6005051-DNA_A11' AND call.DP = 134) OR (call.call_set_name = 'LP6005693-DNA_A03' AND call.DP = 85) OR (call.call_set_name = 'LP6005038-DNA_H06' AND call.DP = 125) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 144) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 115) OR (call.call_set_name = 'LP6005038-DNA_F08' AND call.DP = 164) OR (call.call_set_name = 'LP6005144-DNA_B10' AND call.DP = 171) OR (call.call_set_name = 'LP6005692-DNA_F03' AND call.DP = 155) OR (call.call_set_name = 'LP6005692-DNA_F04' AND call.DP = 174) OR (call.call_set_name = 'LP6005051-DNA_G09' AND call.DP = 176) OR (call.call_set_name = 'LP6005038-DNA_H06' AND call.DP = 148) OR (call.call_set_name = 'LP6005692-DNA_F09' AND call.DP = 145) OR (call.call_set_name = 'LP6005243-DNA_F02' AND call.DP = 182) OR (call.call_set_name = 'LP6005144-DNA_E09' AND call.DP = 100) OR (call.call_set_name = 'LP6005243-DNA_F07' AND call.DP = 158) OR (call.call_set_name = 'LP6005038-DNA_G03' AND call.DP = 120) OR (call.call_set_name = 'LP6005038-DNA_C09' AND call.DP = 179) OR (call.call_set_name = 'LP6005051-DNA_B05' AND call.DP = 148) OR (call.call_set_name = 'LP6005243-DNA_B04' AND call.DP = 139) OR (call.call_set_name = 'LP6005051-DNA_C02' AND call.DP = 103) OR (call.call_set_name = 'LP6005051-DNA_D05' AND call.DP = 177) OR (call.call_set_name = 'LP6005038-DNA_F07' AND call.DP = 143) OR (call.call_set_name = 'LP6005243-DNA_C02' AND call.DP = 122) OR (call.call_set_name = 'LP6005051-DNA_B10' AND call.DP = 152) OR (call.call_set_name = 'LP6005144-DNA_H08' AND call.DP = 159) OR (call.call_set_name = 'LP6005051-DNA_A10' AND call.DP = 154) OR (call.call_set_name = 'LP6005051-DNA_E08' AND call.DP = 220) OR (call.call_set_name = 'LP6005243-DNA_A10' AND call.DP = 108) OR (call.call_set_name = 'LP6005243-DNA_G02' AND call.DP = 133) OR (call.call_set_name = 'LP6005051-DNA_G05' AND call.DP = 117) OR (call.call_set_name = 'LP6005051-DNA_H12' AND call.DP = 126) OR (call.call_set_name = 'LP6005051-DNA_E03' AND call.DP = 225) OR (call.call_set_name = 'LP6005051-DNA_C04' AND call.DP = 146) OR (call.call_set_name = 'LP6005038-DNA_A07' AND call.DP = 104) OR (call.call_set_name = 'LP6005144-DNA_H06' AND call.DP = 144) OR (call.call_set_name = 'LP6005038-DNA_E01' AND call.DP = 151) OR (call.call_set_name = 'LP6005051-DNA_B06' AND call.DP = 182) OR (call.call_set_name = 'LP6005038-DNA_E09' AND call.DP = 105) OR (call.call_set_name = 'LP6005243-DNA_E06' AND call.DP = 121) OR (call.call_set_name = 'LP6005243-DNA_D09' AND call.DP = 115) OR (call.call_set_name = 'LP6005051-DNA_F09' AND call.DP = 93) OR (call.call_set_name = 'LP6005051-DNA_B11' AND call.DP = 179) OR (call.call_set_name = 'LP6005038-DNA_B06' AND call.DP = 159) OR (call.call_set_name = 'LP6005243-DNA_D05' AND call.DP = 112) OR (call.call_set_name = 'LP6005038-DNA_H04' AND call.DP = 142) OR (call.call_set_name = 'LP6005038-DNA_F02' AND call.DP = 136) OR (call.call_set_name = 'LP6005038-DNA_E04' AND call.DP = 175) OR (call.call_set_name = 'LP6005051-DNA_B01' AND call.DP = 122) OR (call.call_set_name = 'LP6005692-DNA_E11' AND call.DP = 151) OR (call.call_set_name = 'LP6005038-DNA_B08' AND call.DP = 129) OR (call.call_set_name = 'LP6005038-DNA_B05' AND call.DP = 151) OR (call.call_set_name = 'LP6005692-DNA_F11' AND call.DP = 168) OR (call.call_set_name = 'LP6005051-DNA_G06' AND call.DP = 103) OR (call.call_set_name = 'LP6005144-DNA_G04' AND call.DP = 152) OR (call.call_set_name = 'LP6005144-DNA_D07' AND call.DP = 189) OR (call.call_set_name = 'LP6005692-DNA_F01' AND call.DP = 207) OR (call.call_set_name = 'LP6005144-DNA_G10' AND call.DP = 99) OR (call.call_set_name = 'LP6005243-DNA_H08' AND call.DP = 154) OR (call.call_set_name = 'LP6005144-DNA_B10' AND call.DP = 152) OR (call.call_set_name = 'LP6005243-DNA_A05' AND call.DP = 120) OR (call.call_set_name = 'LP6005038-DNA_E12' AND call.DP = 149) OR (call.call_set_name = 'LP6005051-DNA_G09' AND call.DP = 159) OR (call.call_set_name = 'LP6005144-DNA_H12' AND call.DP = 115) OR (call.call_set_name = 'LP6005038-DNA_B11' AND call.DP = 144) OR (call.call_set_name = 'LP6005051-DNA_B01' AND call.DP = 109) OR (call.call_set_name = 'LP6005038-DNA_B05' AND call.DP = 142) OR (call.call_set_name = 'LP6005038-DNA_H04' AND call.DP = 129) OR (call.call_set_name = 'LP6005051-DNA_B05' AND call.DP = 114) OR (call.call_set_name = 'LP6005693-DNA_A01' AND call.DP = 196) OR (call.call_set_name = 'LP6005144-DNA_H07' AND call.DP = 167) OR (call.call_set_name = 'LP6005692-DNA_F01' AND call.DP = 195) OR (call.call_set_name = 'LP6005692-DNA_F04' AND call.DP = 154) OR (call.call_set_name = 'LP6005692-DNA_E12' AND call.DP = 167) OR (call.call_set_name = 'LP6005243-DNA_B09' AND call.DP = 167) OR (call.call_set_name = 'LP6005038-DNA_H08' AND call.DP = 170) OR (call.call_set_name = 'LP6005243-DNA_C10' AND call.DP = 153) OR (call.call_set_name = 'LP6005692-DNA_E11' AND call.DP = 153) OR (call.call_set_name = 'LP6005692-DNA_B08' AND call.DP = 159) OR (call.call_set_name = 'LP6005038-DNA_F10' AND call.DP = 109) OR (call.call_set_name = 'LP6005038-DNA_G03' AND call.DP = 106) OR (call.call_set_name = 'LP6005144-DNA_D11' AND call.DP = 172) OR (call.call_set_name = 'LP6005051-DNA_E08' AND call.DP = 207) OR (call.call_set_name = 'LP6005144-DNA_H09' AND call.DP = 169) OR (call.call_set_name = 'LP6005243-DNA_E07' AND call.DP = 112) OR (call.call_set_name = 'LP6005038-DNA_E04' AND call.DP = 181) OR (call.call_set_name = 'LP6005243-DNA_G02' AND call.DP = 117) OR (call.call_set_name = 'LP6005038-DNA_D02' AND call.DP = 171) OR (call.call_set_name = 'LP6005243-DNA_F03' AND call.DP = 133) OR (call.call_set_name = 'LP6005051-DNA_A06' AND call.DP = 159) OR (call.call_set_name = 'LP6005038-DNA_D09' AND call.DP = 108) OR (call.call_set_name = 'LP6005051-DNA_D11' AND call.DP = 113) OR (call.call_set_name = 'LP6005243-DNA_A09' AND call.DP = 162) OR (call.call_set_name = 'LP6005692-DNA_F02' AND call.DP = 136) OR (call.call_set_name = 'LP6005144-DNA_B11' AND call.DP = 123) OR (call.call_set_name = 'LP6005038-DNA_G10' AND call.DP = 167) OR (call.call_set_name = 'LP6005144-DNA_H11' AND call.DP = 126) OR (call.call_set_name = 'LP6005051-DNA_B09' AND call.DP = 148) OR (call.call_set_name = 'LP6005144-DNA_D06' AND call.DP = 125) OR (call.call_set_name = 'LP6005144-DNA_C12' AND call.DP = 101) OR (call.call_set_name = 'LP6005144-DNA_C09' AND call.DP = 125) OR (call.call_set_name = 'LP6005144-DNA_G02' AND call.DP = 132) OR (call.call_set_name = 'LP6005038-DNA_E06' AND call.DP = 141) OR (call.call_set_name = 'LP6005692-DNA_G04' AND call.DP = 126) OR (call.call_set_name = 'LP6005038-DNA_E02' AND call.DP = 201) OR (call.call_set_name = 'LP6005243-DNA_F10' AND call.DP = 133) OR (call.call_set_name = 'LP6005243-DNA_A06' AND call.DP = 115) OR (call.call_set_name = 'LP6005243-DNA_F08' AND call.DP = 129) OR (call.call_set_name = 'LP6005051-DNA_G08' AND call.DP = 96) OR (call.call_set_name = 'LP6005692-DNA_A03' AND call.DP = 136) OR (call.call_set_name = 'LP6005051-DNA_F09' AND call.DP = 80) OR (call.call_set_name = 'LP6005692-DNA_F09' AND call.DP = 115) OR (call.call_set_name = 'LP6005243-DNA_A10' AND call.DP = 78) OR (call.call_set_name = 'LP6005038-DNA_C03' AND call.DP = 88) OR (call.call_set_name = 'LP6005051-DNA_G02' AND call.DP = 133) OR (call.call_set_name = 'LP6005051-DNA_G06' AND call.DP = 98) OR (call.call_set_name = 'LP6005051-DNA_F07' AND call.DP = 74) OR (call.call_set_name = 'LP6005038-DNA_G02' AND call.DP = 115) OR (call.call_set_name = 'LP6005051-DNA_B03' AND call.DP = 84) OR (call.call_set_name = 'LP6005243-DNA_E08' AND call.DP = 116) OR (call.call_set_name = 'LP6005038-DNA_E05' AND call.DP = 173) OR (call.call_set_name = 'LP6005243-DNA_A05' AND call.DP = 118) OR (call.call_set_name = 'LP6005051-DNA_H07' AND call.DP = 132) OR (call.call_set_name = 'LP6005038-DNA_A11' AND call.DP = 156) OR (call.call_set_name = 'LP6005243-DNA_H09' AND call.DP = 150) OR (call.call_set_name = 'LP6005038-DNA_E11' AND call.DP = 120) OR (call.call_set_name = 'LP6005144-DNA_B07' AND call.DP = 101) OR (call.call_set_name = 'LP6005144-DNA_F11' AND call.DP = 167) OR (call.call_set_name = 'LP6005051-DNA_H12' AND call.DP = 106) OR (call.call_set_name = 'LP6005051-DNA_G07' AND call.DP = 138) OR (call.call_set_name = 'LP6005038-DNA_A03' AND call.DP = 161) OR (call.call_set_name = 'LP6005051-DNA_E12' AND call.DP = 109) OR (call.call_set_name = 'LP6005692-DNA_C12' AND call.DP = 173) OR (call.call_set_name = 'LP6005051-DNA_C02' AND call.DP = 77) OR (call.call_set_name = 'LP6005038-DNA_H12' AND call.DP = 123) OR (call.call_set_name = 'LP6005144-DNA_H01' AND call.DP = 97) OR (call.call_set_name = 'LP6005038-DNA_B04' AND call.DP = 209) OR (call.call_set_name = 'LP6005144-DNA_C09' AND call.DP = 112) OR (call.call_set_name = 'LP6005243-DNA_A01' AND call.DP = 133) OR (call.call_set_name = 'LP6005051-DNA_F02' AND call.DP = 136) OR (call.call_set_name = 'LP6005243-DNA_E05' AND call.DP = 175) OR (call.call_set_name = 'LP6005051-DNA_B08' AND call.DP = 113) OR (call.call_set_name = 'LP6005038-DNA_C11' AND call.DP = 132) OR (call.call_set_name = 'LP6005144-DNA_B05' AND call.DP = 95) OR (call.call_set_name = 'LP6005144-DNA_F12' AND call.DP = 146) OR (call.call_set_name = 'LP6005243-DNA_F03' AND call.DP = 137) OR (call.call_set_name = 'LP6005243-DNA_F08' AND call.DP = 126) OR (call.call_set_name = 'LP6005693-DNA_A02' AND call.DP = 100) OR (call.call_set_name = 'LP6005692-DNA_D11' AND call.DP = 117) OR (call.call_set_name = 'LP6005038-DNA_F07' AND call.DP = 138) OR (call.call_set_name = 'LP6005243-DNA_G08' AND call.DP = 102) OR (call.call_set_name = 'LP6005038-DNA_F10' AND call.DP = 99) OR (call.call_set_name = 'LP6005144-DNA_G03' AND call.DP = 107) OR (call.call_set_name = 'LP6005051-DNA_F10' AND call.DP = 165) OR (call.call_set_name = 'LP6005144-DNA_G05' AND call.DP = 167) OR (call.call_set_name = 'LP6005692-DNA_F03' AND call.DP = 105) OR (call.call_set_name = 'LP6005243-DNA_F04' AND call.DP = 189) OR (call.call_set_name = 'LP6005144-DNA_C12' AND call.DP = 108) OR (call.call_set_name = 'LP6005692-DNA_E12' AND call.DP = 168) OR (call.call_set_name = 'LP6005692-DNA_D10' AND call.DP = 193) OR (call.call_set_name = 'LP6005038-DNA_H07' AND call.DP = 121) OR (call.call_set_name = 'LP6005051-DNA_D01' AND call.DP = 139) OR (call.call_set_name = 'LP6005144-DNA_E03' AND call.DP = 142) OR (call.call_set_name = 'LP6005692-DNA_E06' AND call.DP = 197) OR (call.call_set_name = 'LP6005051-DNA_D11' AND call.DP = 104) OR (call.call_set_name = 'LP6005692-DNA_A02' AND call.DP = 178) OR (call.call_set_name = 'LP6005051-DNA_G12' AND call.DP = 137) OR (call.call_set_name = 'LP6005692-DNA_B11' AND call.DP = 227) OR (call.call_set_name = 'LP6005144-DNA_E10' AND call.DP = 160) OR (call.call_set_name = 'LP6005243-DNA_A09' AND call.DP = 152) OR (call.call_set_name = 'LP6005051-DNA_B11' AND call.DP = 174) OR (call.call_set_name = 'LP6005038-DNA_B09' AND call.DP = 119) OR (call.call_set_name = 'LP6005038-DNA_E01' AND call.DP = 108) OR (call.call_set_name = 'LP6005051-DNA_A01' AND call.DP = 132) OR (call.call_set_name = 'LP6005692-DNA_A06' AND call.DP = 173) OR (call.call_set_name = 'LP6005038-DNA_D02' AND call.DP = 149) OR (call.call_set_name = 'LP6005692-DNA_F02' AND call.DP = 127) OR (call.call_set_name = 'LP6005692-DNA_D06' AND call.DP = 123) OR (call.call_set_name = 'LP6005144-DNA_D08' AND call.DP = 147) OR (call.call_set_name = 'LP6005038-DNA_F02' AND call.DP = 132) OR (call.call_set_name = 'LP6005038-DNA_B05' AND call.DP = 141) OR (call.call_set_name = 'LP6005038-DNA_H04' AND call.DP = 130) OR (call.call_set_name = 'LP6005144-DNA_D01' AND call.DP = 172) OR (call.call_set_name = 'LP6005144-DNA_F06' AND call.DP = 111) OR (call.call_set_name = 'LP6005144-DNA_F07' AND call.DP = 200) OR (call.call_set_name = 'LP6005038-DNA_B07' AND call.DP = 106) OR (call.call_set_name = 'LP6005144-DNA_F08' AND call.DP = 101) OR (call.call_set_name = 'LP6005038-DNA_G10' AND call.DP = 162) OR (call.call_set_name = 'LP6005243-DNA_B01' AND call.DP = 208) OR (call.call_set_name = 'LP6005038-DNA_G11' AND call.DP = 93) OR (call.call_set_name = 'LP6005038-DNA_B06' AND call.DP = 180) OR (call.call_set_name = 'LP6005692-DNA_B12' AND call.DP = 221) OR (call.call_set_name = 'LP6005243-DNA_H04' AND call.DP = 249) OR (call.call_set_name = 'LP6005051-DNA_G01' AND call.DP = 132) OR (call.call_set_name = 'LP6005051-DNA_H05' AND call.DP = 111) OR (call.call_set_name = 'LP6005038-DNA_C01' AND call.DP = 113) OR (call.call_set_name = 'LP6005144-DNA_D05' AND call.DP = 172) OR (call.call_set_name = 'LP6005038-DNA_C09' AND call.DP = 171) OR (call.call_set_name = 'LP6005051-DNA_C09' AND call.DP = 157) OR (call.call_set_name = 'LP6005243-DNA_B09' AND call.DP = 150) OR (call.call_set_name = 'LP6005243-DNA_D10' AND call.DP = 113) OR (call.call_set_name = 'LP6005038-DNA_C08' AND call.DP = 141) OR (call.call_set_name = 'LP6005051-DNA_F08' AND call.DP = 114) OR (call.call_set_name = 'LP6005243-DNA_C02' AND call.DP = 102) OR (call.call_set_name = 'LP6005038-DNA_A05' AND call.DP = 110) OR (call.call_set_name = 'LP6005144-DNA_G01' AND call.DP = 145) OR (call.call_set_name = 'LP6005051-DNA_G11' AND call.DP = 250) OR (call.call_set_name = 'LP6005051-DNA_A06' AND call.DP = 150) OR (call.call_set_name = 'LP6005243-DNA_B03' AND call.DP = 167) OR (call.call_set_name = 'LP6005144-DNA_A04' AND call.DP = 153) OR (call.call_set_name = 'LP6005692-DNA_E11' AND call.DP = 150) OR (call.call_set_name = 'LP6005051-DNA_C07' AND call.DP = 129) OR (call.call_set_name = 'LP6005051-DNA_C04' AND call.DP = 118) OR (call.call_set_name = 'LP6005692-DNA_A05' AND call.DP = 103) OR (call.call_set_name = 'LP6005038-DNA_G03' AND call.DP = 110) OR (call.call_set_name = 'LP6005038-DNA_E07' AND call.DP = 126) OR (call.call_set_name = 'LP6005144-DNA_B12' AND call.DP = 193) OR (call.call_set_name = 'LP6005051-DNA_C06' AND call.DP = 171) OR (call.call_set_name = 'LP6005144-DNA_G06' AND call.DP = 115) OR (call.call_set_name = 'LP6005144-DNA_A06' AND call.DP = 121) OR (call.call_set_name = 'LP6005243-DNA_F06' AND call.DP = 180) OR (call.call_set_name = 'LP6005038-DNA_A09' AND call.DP = 173) OR (call.call_set_name = 'LP6005144-DNA_E05' AND call.DP = 115) OR (call.call_set_name = 'LP6005693-DNA_B01' AND call.DP = 155) OR (call.call_set_name = 'LP6005051-DNA_E01' AND call.DP = 126) OR (call.call_set_name = 'LP6005144-DNA_G11' AND call.DP = 110) OR (call.call_set_name = 'LP6005243-DNA_D01' AND call.DP = 170) OR (call.call_set_name = 'LP6005038-DNA_G12' AND call.DP = 180) OR (call.call_set_name = 'LP6005038-DNA_C04' AND call.DP = 123) OR (call.call_set_name = 'LP6005243-DNA_H05' AND call.DP = 150) OR (call.call_set_name = 'LP6005038-DNA_H02' AND call.DP = 155) OR (call.call_set_name = 'LP6005243-DNA_B06' AND call.DP = 146) OR (call.call_set_name = 'LP6005692-DNA_B01' AND call.DP = 143) OR (call.call_set_name = 'LP6005051-DNA_B05' AND call.DP = 143) OR (call.call_set_name = 'LP6005243-DNA_B05' AND call.DP = 111) OR (call.call_set_name = 'LP6005692-DNA_C05' AND call.DP = 99) OR (call.call_set_name = 'LP6005144-DNA_G02' AND call.DP = 110) OR (call.call_set_name = 'LP6005051-DNA_H07' AND call.DP = 143) OR (call.call_set_name = 'LP6005243-DNA_E06' AND call.DP = 105) OR (call.call_set_name = 'LP6005243-DNA_D02' AND call.DP = 109) OR (call.call_set_name = 'LP6005144-DNA_C08' AND call.DP = 109) OR (call.call_set_name = 'LP6005144-DNA_B08' AND call.DP = 82) OR (call.call_set_name = 'LP6005144-DNA_F03' AND call.DP = 105) OR (call.call_set_name = 'LP6005051-DNA_A08' AND call.DP = 121) OR (call.call_set_name = 'LP6005243-DNA_F03' AND call.DP = 138) OR (call.call_set_name = 'LP6005692-DNA_H05' AND call.DP = 108) OR (call.call_set_name = 'LP6005038-DNA_G06' AND call.DP = 90) OR (call.call_set_name = 'LP6005692-DNA_D11' AND call.DP = 109) OR (call.call_set_name = 'LP6005243-DNA_F05' AND call.DP = 132) OR (call.call_set_name = 'LP6005038-DNA_F07' AND call.DP = 130) OR (call.call_set_name = 'LP6005051-DNA_F01' AND call.DP = 84) OR (call.call_set_name = 'LP6005243-DNA_G08' AND call.DP = 95) OR (call.call_set_name = 'LP6005144-DNA_H01' AND call.DP = 87) OR (call.call_set_name = 'LP6005038-DNA_E03' AND call.DP = 106) OR (call.call_set_name = 'LP6005051-DNA_H06' AND call.DP = 110) OR (call.call_set_name = 'LP6005038-DNA_A02' AND call.DP = 106) OR (call.call_set_name = 'LP6005243-DNA_E01' AND call.DP = 156) OR (call.call_set_name = 'LP6005243-DNA_F10' AND call.DP = 126) OR (call.call_set_name = 'LP6005038-DNA_F10' AND call.DP = 93) OR (call.call_set_name = 'LP6005243-DNA_B04' AND call.DP = 135) OR (call.call_set_name = 'LP6005051-DNA_G02' AND call.DP = 158) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 105) OR (call.call_set_name = 'LP6005144-DNA_H04' AND call.DP = 100) OR (call.call_set_name = 'LP6005243-DNA_A01' AND call.DP = 138) OR (call.call_set_name = 'LP6005692-DNA_B02' AND call.DP = 106) OR (call.call_set_name = 'LP6005243-DNA_A02' AND call.DP = 103) OR (call.call_set_name = 'LP6005144-DNA_G02' AND call.DP = 114) OR (call.call_set_name = 'LP6005051-DNA_C05' AND call.DP = 106) OR (call.call_set_name = 'LP6005038-DNA_B10' AND call.DP = 115) OR (call.call_set_name = 'LP6005051-DNA_D07' AND call.DP = 132) OR (call.call_set_name = 'LP6005144-DNA_D02' AND call.DP = 95) OR (call.call_set_name = 'LP6005038-DNA_A10' AND call.DP = 99) OR (call.call_set_name = 'LP6005243-DNA_F03' AND call.DP = 143) OR (call.call_set_name = 'LP6005038-DNA_G06' AND call.DP = 96) OR (call.call_set_name = 'LP6005243-DNA_G08' AND call.DP = 101) OR (call.call_set_name = 'LP6005051-DNA_A12' AND call.DP = 89) OR (call.call_set_name = 'LP6005038-DNA_F03' AND call.DP = 94) OR (call.call_set_name = 'LP6005051-DNA_B12' AND call.DP = 85) OR (call.call_set_name = 'LP6005051-DNA_E11' AND call.DP = 108) OR (call.call_set_name = 'LP6005243-DNA_F10' AND call.DP = 128) OR (call.call_set_name = 'LP6005038-DNA_F10' AND call.DP = 101) OR (call.call_set_name = 'LP6005051-DNA_H07' AND call.DP = 133) OR (call.call_set_name = 'LP6005051-DNA_A04' AND call.DP = 80) OR (call.call_set_name = 'LP6005243-DNA_F04' AND call.DP = 193) OR (call.call_set_name = 'LP6005051-DNA_E07' AND call.DP = 115) OR (call.call_set_name = 'LP6005243-DNA_A05' AND call.DP = 126) OR (call.call_set_name = 'LP6005692-DNA_H11' AND call.DP = 100) OR (call.call_set_name = 'LP6005051-DNA_F06' AND call.DP = 112) OR (call.call_set_name = 'LP6005692-DNA_B11' AND call.DP = 225) OR (call.call_set_name = 'LP6005038-DNA_E06' AND call.DP = 145) OR (call.call_set_name = 'LP6005144-DNA_E03' AND call.DP = 156) OR (call.call_set_name = 'LP6005692-DNA_E06' AND call.DP = 180) OR (call.call_set_name = 'LP6005051-DNA_F04' AND call.DP = 116) OR (call.call_set_name = 'LP6005038-DNA_F11' AND call.DP = 73) OR (call.call_set_name = 'LP6005692-DNA_G02' AND call.DP = 92) OR (call.call_set_name = 'LP6005038-DNA_B11' AND call.DP = 169) OR (call.call_set_name = 'LP6005144-DNA_D12' AND call.DP = 105) OR (call.call_set_name = 'LP6005243-DNA_A09' AND call.DP = 109) OR (call.call_set_name = 'LP6005038-DNA_B02' AND call.DP = 100) OR (call.call_set_name = 'LP6005038-DNA_B09' AND call.DP = 110) OR (call.call_set_name = 'LP6005038-DNA_E01' AND call.DP = 109) OR (call.call_set_name = 'LP6005051-DNA_G04' AND call.DP = 123) OR (call.call_set_name = 'LP6005051-DNA_H03' AND call.DP = 229) OR (call.call_set_name = 'LP6005692-DNA_D12' AND call.DP = 217) OR (call.call_set_name = 'LP6005038-DNA_D10' AND call.DP = 213) OR (call.call_set_name = 'LP6005051-DNA_B01' AND call.DP = 118) OR (call.call_set_name = 'LP6005243-DNA_D02' AND call.DP = 120) OR (call.call_set_name = 'LP6005693-DNA_D01' AND call.DP = 194) OR (call.call_set_name = 'LP6005144-DNA_G12' AND call.DP = 151) OR (call.call_set_name = 'LP6005144-DNA_D01' AND call.DP = 199) OR (call.call_set_name = 'LP6005051-DNA_G10' AND call.DP = 166) OR (call.call_set_name = 'LP6005038-DNA_G08' AND call.DP = 196) OR (call.call_set_name = 'LP6005038-DNA_C11' AND call.DP = 114) OR (call.call_set_name = 'LP6005051-DNA_B06' AND call.DP = 153) OR (call.call_set_name = 'LP6005038-DNA_B01' AND call.DP = 84) OR (call.call_set_name = 'LP6005038-DNA_B10' AND call.DP = 122) OR (call.call_set_name = 'LP6005144-DNA_H07' AND call.DP = 141) OR (call.call_set_name = 'LP6005038-DNA_D07' AND call.DP = 79) OR (call.call_set_name = 'LP6005144-DNA_B05' AND call.DP = 79) OR (call.call_set_name = 'LP6005692-DNA_F04' AND call.DP = 150) OR (call.call_set_name = 'LP6005038-DNA_G11' AND call.DP = 97) OR (call.call_set_name = 'LP6005051-DNA_C06' AND call.DP = 174) OR (call.call_set_name = 'LP6005692-DNA_F07' AND call.DP = 163) OR (call.call_set_name = 'LP6005692-DNA_H05' AND call.DP = 118) OR (call.call_set_name = 'LP6005243-DNA_H02' AND call.DP = 142) OR (call.call_set_name = 'LP6005051-DNA_D02' AND call.DP = 133) OR (call.call_set_name = 'LP6005144-DNA_H02' AND call.DP = 127) OR (call.call_set_name = 'LP6005692-DNA_H09' AND call.DP = 122) OR (call.call_set_name = 'LP6005144-DNA_B10' AND call.DP = 163) OR (call.call_set_name = 'LP6005692-DNA_D03' AND call.DP = 184) OR (call.call_set_name = 'LP6005051-DNA_A03' AND call.DP = 182) OR (call.call_set_name = 'LP6005144-DNA_D05' AND call.DP = 171) OR (call.call_set_name = 'LP6005038-DNA_C09' AND call.DP = 154) OR (call.call_set_name = 'LP6005144-DNA_C01' AND call.DP = 115) OR (call.call_set_name = 'LP6005051-DNA_E06' AND call.DP = 107) OR (call.call_set_name = 'LP6005243-DNA_F07' AND call.DP = 120) OR (call.call_set_name = 'LP6005051-DNA_F01' AND call.DP = 98) OR (call.call_set_name = 'LP6005038-DNA_H08' AND call.DP = 135) OR (call.call_set_name = 'LP6005144-DNA_B06' AND call.DP = 85) OR (call.call_set_name = 'LP6005051-DNA_B07' AND call.DP = 165) OR (call.call_set_name = 'LP6005051-DNA_C12' AND call.DP = 131) OR (call.call_set_name = 'LP6005243-DNA_H06' AND call.DP = 125) OR (call.call_set_name = 'LP6005243-DNA_C01' AND call.DP = 141) OR (call.call_set_name = 'LP6005144-DNA_A01' AND call.DP = 140) OR (call.call_set_name = 'LP6005051-DNA_E09' AND call.DP = 154) OR (call.call_set_name = 'LP6005051-DNA_H01' AND call.DP = 115) OR (call.call_set_name = 'LP6005243-DNA_A11' AND call.DP = 146) OR (call.call_set_name = 'LP6005144-DNA_B09' AND call.DP = 125) OR (call.call_set_name = 'LP6005243-DNA_E09' AND call.DP = 108) OR (call.call_set_name = 'LP6005144-DNA_A04' AND call.DP = 159) OR (call.call_set_name = 'LP6005692-DNA_F11' AND call.DP = 164) OR (call.call_set_name = 'LP6005692-DNA_E11' AND call.DP = 148) OR (call.call_set_name = 'LP6005051-DNA_E10' AND call.DP = 78) OR (call.call_set_name = 'LP6005243-DNA_D06' AND call.DP = 98) OR (call.call_set_name = 'LP6005243-DNA_D03' AND call.DP = 170) OR (call.call_set_name = 'LP6005051-DNA_C04' AND call.DP = 117) OR (call.call_set_name = 'LP6005243-DNA_F10' AND call.DP = 115) OR (call.call_set_name = 'LP6005038-DNA_E07' AND call.DP = 103) OR (call.call_set_name = 'LP6005038-DNA_F05' AND call.DP = 117) OR (call.call_set_name = 'LP6005243-DNA_F02' AND call.DP = 170) OR (call.call_set_name = 'LP6005243-DNA_B04' AND call.DP = 127) OR (call.call_set_name = 'LP6005144-DNA_B12' AND call.DP = 176) OR (call.call_set_name = 'LP6005051-DNA_F12' AND call.DP = 153) OR (call.call_set_name = 'LP6005144-DNA_H08' AND call.DP = 139) OR (call.call_set_name = 'LP6005038-DNA_C05' AND call.DP = 164) OR (call.call_set_name = 'LP6005243-DNA_F01' AND call.DP = 208) OR (call.call_set_name = 'LP6005243-DNA_C05' AND call.DP = 115) OR (call.call_set_name = 'LP6005144-DNA_H05' AND call.DP = 184) OR (call.call_set_name = 'LP6005144-DNA_H09' AND call.DP = 159) OR (call.call_set_name = 'LP6005692-DNA_B02' AND call.DP = 113) OR (call.call_set_name = 'LP6005243-DNA_C04' AND call.DP = 95) OR (call.call_set_name = 'LP6005051-DNA_B05' AND call.DP = 161) OR (call.call_set_name = 'LP6005051-DNA_F10' AND call.DP = 141) OR (call.call_set_name = 'LP6005144-DNA_H03' AND call.DP = 144) OR (call.call_set_name = 'LP6005243-DNA_H09' AND call.DP = 141) OR (call.call_set_name = 'LP6005243-DNA_H05' AND call.DP = 145) OR (call.call_set_name = 'LP6005243-DNA_G01' AND call.DP = 206) OR (call.call_set_name = 'LP6005692-DNA_F06' AND call.DP = 152) OR (call.call_set_name = 'LP6005243-DNA_D05' AND call.DP = 109) OR (call.call_set_name = 'LP6005243-DNA_E04' AND call.DP = 168) OR (call.call_set_name = 'LP6005051-DNA_C03' AND call.DP = 172) OR (call.call_set_name = 'LP6005243-DNA_G03' AND call.DP = 150) OR (call.call_set_name = 'LP6005051-DNA_E08' AND call.DP = 193) OR (call.call_set_name = 'LP6005692-DNA_A08' AND call.DP = 106) OR (call.call_set_name = 'LP6005692-DNA_A06' AND call.DP = 152) OR (call.call_set_name = 'LP6005038-DNA_D05' AND call.DP = 114) OR (call.call_set_name = 'LP6005038-DNA_H12' AND call.DP = 115) OR (call.call_set_name = 'LP6005051-DNA_E01' AND call.DP = 151) OR (call.call_set_name = 'LP6005051-DNA_D11' AND call.DP = 109) OR (call.call_set_name = 'LP6005144-DNA_A08' AND call.DP = 102) OR (call.call_set_name = 'LP6005692-DNA_D06' AND call.DP = 115) OR (call.call_set_name = 'LP6005243-DNA_B07' AND call.DP = 95) OR (call.call_set_name = 'LP6005051-DNA_E05' AND call.DP = 114) OR (call.call_set_name = 'LP6005038-DNA_E02' AND call.DP = 162) OR (call.call_set_name = 'LP6005243-DNA_E01' AND call.DP = 139) OR (call.call_set_name = 'LP6005692-DNA_H07' AND call.DP = 124) OR (call.call_set_name = 'LP6005051-DNA_E11' AND call.DP = 96) OR (call.call_set_name = 'LP6005038-DNA_D02' AND call.DP = 144) OR (call.call_set_name = 'LP6005243-DNA_B06' AND call.DP = 132) OR (call.call_set_name = 'LP6005692-DNA_H06' AND call.DP = 123) OR (call.call_set_name = 'LP6005051-DNA_E03' AND call.DP = 186) OR (call.call_set_name = 'LP6005038-DNA_B11' AND call.DP = 167) OR (call.call_set_name = 'LP6005038-DNA_H02' AND call.DP = 134) OR (call.call_set_name = 'LP6005051-DNA_D11' AND call.DP = 92) OR (call.call_set_name = 'LP6005051-DNA_D01' AND call.DP = 122) OR (call.call_set_name = 'LP6005051-DNA_B08' AND call.DP = 107) OR (call.call_set_name = 'LP6005243-DNA_H05' AND call.DP = 120) OR (call.call_set_name = 'LP6005038-DNA_D10' AND call.DP = 203) OR (call.call_set_name = 'LP6005144-DNA_H11' AND call.DP = 103) OR (call.call_set_name = 'LP6005051-DNA_H04' AND call.DP = 124) OR (call.call_set_name = 'LP6005692-DNA_F02' AND call.DP = 133) OR (call.call_set_name = 'LP6005051-DNA_A11' AND call.DP = 151) OR (call.call_set_name = 'LP6005144-DNA_A01' AND call.DP = 120) OR (call.call_set_name = 'LP6005692-DNA_E12' AND call.DP = 128) OR (call.call_set_name = 'LP6005144-DNA_D05' AND call.DP = 174) OR (call.call_set_name = 'LP6005243-DNA_D03' AND call.DP = 146) OR (call.call_set_name = 'LP6005243-DNA_D10' AND call.DP = 122) OR (call.call_set_name = 'LP6005051-DNA_B05' AND call.DP = 139) OR (call.call_set_name = 'LP6005144-DNA_A07' AND call.DP = 105) OR (call.call_set_name = 'LP6005038-DNA_B02' AND call.DP = 94) OR (call.call_set_name = 'LP6005051-DNA_H08' AND call.DP = 112) OR (call.call_set_name = 'LP6005038-DNA_A06' AND call.DP = 166) OR (call.call_set_name = 'LP6005692-DNA_B08' AND call.DP = 152) OR (call.call_set_name = 'LP6005051-DNA_B06' AND call.DP = 145) OR (call.call_set_name = 'LP6005051-DNA_D10' AND call.DP = 81) OR (call.call_set_name = 'LP6005051-DNA_D02' AND call.DP = 115) OR (call.call_set_name = 'LP6005692-DNA_F06' AND call.DP = 149) OR (call.call_set_name = 'LP6005692-DNA_H10' AND call.DP = 116) OR (call.call_set_name = 'LP6005243-DNA_E06' AND call.DP = 111) OR (call.call_set_name = 'LP6005693-DNA_D01' AND call.DP = 197) OR (call.call_set_name = 'LP6005038-DNA_G09' AND call.DP = 75) OR (call.call_set_name = 'LP6005692-DNA_D10' AND call.DP = 184) OR (call.call_set_name = 'LP6005038-DNA_G12' AND call.DP = 158) OR (call.call_set_name = 'LP6005243-DNA_H10' AND call.DP = 106) OR (call.call_set_name = 'LP6005243-DNA_E09' AND call.DP = 113) OR (call.call_set_name = 'LP6005051-DNA_H06' AND call.DP = 98) OR (call.call_set_name = 'LP6005692-DNA_E06' AND call.DP = 170) OR (call.call_set_name = 'LP6005038-DNA_D02' AND call.DP = 134) OR (call.call_set_name = 'LP6005051-DNA_C09' AND call.DP = 131) OR (call.call_set_name = 'LP6005051-DNA_G01' AND call.DP = 130) OR (call.call_set_name = 'LP6005243-DNA_C01' AND call.DP = 99) OR (call.call_set_name = 'LP6005692-DNA_A06' AND call.DP = 141) OR (call.call_set_name = 'LP6005038-DNA_E01' AND call.DP = 89) OR (call.call_set_name = 'LP6005144-DNA_F07' AND call.DP = 192) OR (call.call_set_name = 'LP6005038-DNA_A10' AND call.DP = 82) OR (call.call_set_name = 'LP6005243-DNA_E04' AND call.DP = 165) OR (call.call_set_name = 'LP6005038-DNA_D05' AND call.DP = 108) OR (call.call_set_name = 'LP6005693-DNA_A01' AND call.DP = 161) OR (call.call_set_name = 'LP6005051-DNA_F10' AND call.DP = 143) OR (call.call_set_name = 'LP6005051-DNA_F08' AND call.DP = 96) OR (call.call_set_name = 'LP6005038-DNA_C05' AND call.DP = 154) OR (call.call_set_name = 'LP6005692-DNA_B03' AND call.DP = 104) OR (call.call_set_name = 'LP6005038-DNA_H01' AND call.DP = 138) OR (call.call_set_name = 'LP6005144-DNA_F08' AND call.DP = 107) OR (call.call_set_name = 'LP6005051-DNA_H05' AND call.DP = 92) OR (call.call_set_name = 'LP6005243-DNA_G03' AND call.DP = 156) OR (call.call_set_name = 'LP6005051-DNA_B01' AND call.DP = 117) OR (call.call_set_name = 'LP6005692-DNA_H06' AND call.DP = 125) OR (call.call_set_name = 'LP6005144-DNA_F06' AND call.DP = 105) OR (call.call_set_name = 'LP6005243-DNA_H01' AND call.DP = 110) OR (call.call_set_name = 'LP6005243-DNA_H06' AND call.DP = 105) OR (call.call_set_name = 'LP6005144-DNA_C12' AND call.DP = 89) OR (call.call_set_name = 'LP6005038-DNA_E09' AND call.DP = 83) OR (call.call_set_name = 'LP6005144-DNA_A04' AND call.DP = 142) OR (call.call_set_name = 'LP6005144-DNA_A06' AND call.DP = 115) OR (call.call_set_name = 'LP6005038-DNA_F03' AND call.DP = 93) OR (call.call_set_name = 'LP6005038-DNA_C08' AND call.DP = 137) OR (call.call_set_name = 'LP6005051-DNA_F02' AND call.DP = 138) OR (call.call_set_name = 'LP6005051-DNA_G04' AND call.DP = 106) OR (call.call_set_name = 'LP6005243-DNA_H02' AND call.DP = 125) OR (call.call_set_name = 'LP6005038-DNA_C01' AND call.DP = 128) OR (call.call_set_name = 'LP6005693-DNA_A03' AND call.DP = 67) OR (call.call_set_name = 'LP6005243-DNA_B03' AND call.DP = 166) OR (call.call_set_name = 'LP6005051-DNA_E10' AND call.DP = 81) OR (call.call_set_name = 'LP6005038-DNA_F05' AND call.DP = 116) OR (call.call_set_name = 'LP6005051-DNA_B09' AND call.DP = 126) OR (call.call_set_name = 'LP6005144-DNA_D12' AND call.DP = 97) OR (call.call_set_name = 'LP6005692-DNA_F03' AND call.DP = 127) OR (call.call_set_name = 'LP6005051-DNA_C04' AND call.DP = 126) OR (call.call_set_name = 'LP6005051-DNA_G12' AND call.DP = 150) OR (call.call_set_name = 'LP6005692-DNA_H07' AND call.DP = 119) OR (call.call_set_name = 'LP6005692-DNA_F01' AND call.DP = 183) OR (call.call_set_name = 'LP6005243-DNA_G08' AND call.DP = 112) OR (call.call_set_name = 'LP6005051-DNA_B07' AND call.DP = 141) OR (call.call_set_name = 'LP6005051-DNA_D03' AND call.DP = 95) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 92) OR (call.call_set_name = 'LP6005693-DNA_B01' AND call.DP = 146) OR (call.call_set_name = 'LP6005144-DNA_B04' AND call.DP = 89) OR (call.call_set_name = 'LP6005692-DNA_G01' AND call.DP = 125) OR (call.call_set_name = 'LP6005051-DNA_C06' AND call.DP = 178) OR (call.call_set_name = 'LP6005051-DNA_D06' AND call.DP = 122) OR (call.call_set_name = 'LP6005243-DNA_G03' AND call.DP = 199) OR (call.call_set_name = 'LP6005243-DNA_D04' AND call.DP = 147) OR (call.call_set_name = 'LP6005051-DNA_H05' AND call.DP = 133) OR (call.call_set_name = 'LP6005051-DNA_B05' AND call.DP = 177) OR (call.call_set_name = 'LP6005051-DNA_C03' AND call.DP = 170) OR (call.call_set_name = 'LP6005038-DNA_D10' AND call.DP = 250) OR (call.call_set_name = 'LP6005243-DNA_G10' AND call.DP = 131) OR (call.call_set_name = 'LP6005038-DNA_E01' AND call.DP = 158) OR (call.call_set_name = 'LP6005038-DNA_H01' AND call.DP = 187) OR (call.call_set_name = 'LP6005692-DNA_B03' AND call.DP = 160) OR (call.call_set_name = 'LP6005038-DNA_G07' AND call.DP = 147) OR (call.call_set_name = 'LP6005051-DNA_H12' AND call.DP = 168) OR (call.call_set_name = 'LP6005038-DNA_G05' AND call.DP = 159) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 151) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 183) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 150) OR (call.call_set_name = 'LP6005051-DNA_C03' AND call.DP = 208) OR (call.call_set_name = 'LP6005243-DNA_F02' AND call.DP = 225) OR (call.call_set_name = 'LP6005693-DNA_A03' AND call.DP = 138) OR (call.call_set_name = 'LP6005243-DNA_B03' AND call.DP = 237) OR (call.call_set_name = 'LP6005038-DNA_H01' AND call.DP = 217) OR (call.call_set_name = 'LP6005051-DNA_D06' AND call.DP = 167) OR (call.call_set_name = 'LP6005038-DNA_H10' AND call.DP = 166) OR (call.call_set_name = 'LP6005144-DNA_F07' AND call.DP = 250) OR (call.call_set_name = 'LP6005051-DNA_A07' AND call.DP = 192) OR (call.call_set_name = 'LP6005243-DNA_C10' AND call.DP = 234) OR (call.call_set_name = 'LP6005051-DNA_F04' AND call.DP = 187) OR (call.call_set_name = 'LP6005038-DNA_C10' AND call.DP = 189) OR (call.call_set_name = 'LP6005051-DNA_E02' AND call.DP = 156) OR (call.call_set_name = 'LP6005051-DNA_F02' AND call.DP = 215) OR (call.call_set_name = 'LP6005692-DNA_E01' AND call.DP = 182) OR (call.call_set_name = 'LP6005038-DNA_E01' AND call.DP = 168) OR (call.call_set_name = 'LP6005051-DNA_B03' AND call.DP = 157) OR (call.call_set_name = 'LP6005051-DNA_F05' AND call.DP = 171) OR (call.call_set_name = 'LP6005051-DNA_H02' AND call.DP = 154) OR (call.call_set_name = 'LP6005243-DNA_F06' AND call.DP = 233) OR (call.call_set_name = 'LP6005038-DNA_F06' AND call.DP = 142) OR (call.call_set_name = 'LP6005243-DNA_A07' AND call.DP = 220) OR (call.call_set_name = 'LP6005051-DNA_H12' AND call.DP = 175) OR (call.call_set_name = 'LP6005051-DNA_A06' AND call.DP = 228) OR (call.call_set_name = 'LP6005038-DNA_F05' AND call.DP = 182) OR (call.call_set_name = 'LP6005051-DNA_H11' AND call.DP = 118) OR (call.call_set_name = 'LP6005692-DNA_F11' AND call.DP = 244) OR (call.call_set_name = 'LP6005051-DNA_D01' AND call.DP = 229) OR (call.call_set_name = 'LP6005051-DNA_F03' AND call.DP = 223) OR (call.call_set_name = 'LP6005051-DNA_C04' AND call.DP = 188) OR (call.call_set_name = 'LP6005038-DNA_E06' AND call.DP = 193) OR (call.call_set_name = 'LP6005038-DNA_E07' AND call.DP = 196) OR (call.call_set_name = 'LP6005038-DNA_A05' AND call.DP = 191) OR (call.call_set_name = 'LP6005051-DNA_B01' AND call.DP = 202) OR (call.call_set_name = 'LP6005692-DNA_G05' AND call.DP = 158) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 155) OR (call.call_set_name = 'LP6005144-DNA_C11' AND call.DP = 139) OR (call.call_set_name = 'LP6005038-DNA_H03' AND call.DP = 233) OR (call.call_set_name = 'LP6005693-DNA_A01' AND call.DP = 250) OR (call.call_set_name = 'LP6005038-DNA_F07' AND call.DP = 213) OR (call.call_set_name = 'LP6005243-DNA_C04' AND call.DP = 164) OR (call.call_set_name = 'LP6005144-DNA_G03' AND call.DP = 217) OR (call.call_set_name = 'LP6005038-DNA_G06' AND call.DP = 196) OR (call.call_set_name = 'LP6005051-DNA_F01' AND call.DP = 156) OR (call.call_set_name = 'LP6005243-DNA_H04' AND call.DP = 250) OR (call.call_set_name = 'LP6005243-DNA_A06' AND call.DP = 201) OR (call.call_set_name = 'LP6005038-DNA_G10' AND call.DP = 250) OR (call.call_set_name = 'LP6005243-DNA_D04' AND call.DP = 203) OR (call.call_set_name = 'LP6005051-DNA_D08' AND call.DP = 179) OR (call.call_set_name = 'LP6005038-DNA_C04' AND call.DP = 239) OR (call.call_set_name = 'LP6005243-DNA_B01' AND call.DP = 250) OR (call.call_set_name = 'LP6005144-DNA_D08' AND call.DP = 229) OR (call.call_set_name = 'LP6005051-DNA_C03' AND call.DP = 229) OR (call.call_set_name = 'LP6005692-DNA_G07' AND call.DP = 237) OR (call.call_set_name = 'LP6005144-DNA_D04' AND call.DP = 220) OR (call.call_set_name = 'LP6005692-DNA_B07' AND call.DP = 164) OR (call.call_set_name = 'LP6005038-DNA_D09' AND call.DP = 225) OR (call.call_set_name = 'LP6005038-DNA_E01' AND call.DP = 210) OR (call.call_set_name = 'LP6005692-DNA_D02' AND call.DP = 150) OR (call.call_set_name = 'LP6005692-DNA_B03' AND call.DP = 216) OR (call.call_set_name = 'LP6005051-DNA_C05' AND call.DP = 214) OR (call.call_set_name = 'LP6005038-DNA_H04' AND call.DP = 231) OR (call.call_set_name = 'LP6005144-DNA_G09' AND call.DP = 212) OR (call.call_set_name = 'LP6005692-DNA_F04' AND call.DP = 248) OR (call.call_set_name = 'LP6005692-DNA_E01' AND call.DP = 200) OR (call.call_set_name = 'LP6005038-DNA_F07' AND call.DP = 208) OR (call.call_set_name = 'LP6005051-DNA_G09' AND call.DP = 250) OR (call.call_set_name = 'LP6005144-DNA_A12' AND call.DP = 211) OR (call.call_set_name = 'LP6005051-DNA_D06' AND call.DP = 207) OR (call.call_set_name = 'LP6005051-DNA_D02' AND call.DP = 223) OR (call.call_set_name = 'LP6005038-DNA_D01' AND call.DP = 232) OR (call.call_set_name = 'LP6005243-DNA_C08' AND call.DP = 217) OR (call.call_set_name = 'LP6005144-DNA_B11' AND call.DP = 221) OR (call.call_set_name = 'LP6005038-DNA_G04' AND call.DP = 250) OR (call.call_set_name = 'LP6005693-DNA_B01' AND call.DP = 249) OR (call.call_set_name = 'LP6005051-DNA_B03' AND call.DP = 196) OR (call.call_set_name = 'LP6005144-DNA_B03' AND call.DP = 219) OR (call.call_set_name = 'LP6005038-DNA_C03' AND call.DP = 190) OR (call.call_set_name = 'LP6005038-DNA_C01' AND call.DP = 205) OR (call.call_set_name = 'LP6005051-DNA_F01' AND call.DP = 189) OR (call.call_set_name = 'LP6005051-DNA_H02' AND call.DP = 210) OR (call.call_set_name = 'LP6005144-DNA_A11' AND call.DP = 164) OR (call.call_set_name = 'LP6005692-DNA_G04' AND call.DP = 196) OR (call.call_set_name = 'LP6005692-DNA_C02' AND call.DP = 168) OR (call.call_set_name = 'LP6005038-DNA_H10' AND call.DP = 176) OR (call.call_set_name = 'LP6005243-DNA_H04' AND call.DP = 247) OR (call.call_set_name = 'LP6005693-DNA_B01' AND call.DP = 248) OR (call.call_set_name = 'LP6005038-DNA_F09' AND call.DP = 243) OR (call.call_set_name = 'LP6005692-DNA_B03' AND call.DP = 207) OR (call.call_set_name = 'LP6005051-DNA_C05' AND call.DP = 207) OR (call.call_set_name = 'LP6005144-DNA_G06' AND call.DP = 206) OR (call.call_set_name = 'LP6005693-DNA_A01' AND call.DP = 249) OR (call.call_set_name = 'LP6005692-DNA_D03' AND call.DP = 249) OR (call.call_set_name = 'LP6005144-DNA_H08' AND call.DP = 227) OR (call.call_set_name = 'LP6005243-DNA_A11' AND call.DP = 250) OR (call.call_set_name = 'LP6005038-DNA_D01' AND call.DP = 223) OR (call.call_set_name = 'LP6005144-DNA_A11' AND call.DP = 158) OR (call.call_set_name = 'LP6005692-DNA_F11' AND call.DP = 248) OR (call.call_set_name = 'LP6005243-DNA_F10' AND call.DP = 219) OR (call.call_set_name = 'LP6005038-DNA_F10' AND call.DP = 205) OR (call.call_set_name = 'LP6005038-DNA_G10' AND call.DP = 248) OR (call.call_set_name = 'LP6005692-DNA_H02' AND call.DP = 167) OR (call.call_set_name = 'LP6005038-DNA_G04' AND call.DP = 239) OR (call.call_set_name = 'LP6005051-DNA_F09' AND call.DP = 153) OR (call.call_set_name = 'LP6005144-DNA_G09' AND call.DP = 190) OR (call.call_set_name = 'LP6005243-DNA_D07' AND call.DP = 185) OR (call.call_set_name = 'LP6005144-DNA_B08' AND call.DP = 175) OR (call.call_set_name = 'LP6005051-DNA_E07' AND call.DP = 232) OR (call.call_set_name = 'LP6005051-DNA_A06' AND call.DP = 243) OR (call.call_set_name = 'LP6005144-DNA_A12' AND call.DP = 200) OR (call.call_set_name = 'LP6005038-DNA_F07' AND call.DP = 175) OR (call.call_set_name = 'LP6005051-DNA_B03' AND call.DP = 173) OR (call.call_set_name = 'LP6005051-DNA_H02' AND call.DP = 191)
LIMIT 1000
```

First few results
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Thu May 28 22:34:43 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> call_call_set_name </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTYYxa_YHyChkfm3vL_J8_UB </td> <td> LP6005051-DNA_D12 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTYY0q_YHyCoh4O5gKnYlVM </td> <td> LP6005051-DNA_D12 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTYY1a_YHyDA7_fkkdKghv0B </td> <td> LP6005051-DNA_D12 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTYY3K_YHyDQ7eCY4OWi_Xk </td> <td> LP6005051-DNA_D12 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTYY3a_YHyDI_frmo-aI0UU </td> <td> LP6005051-DNA_D12 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTYY76_YHyC2n_b3-oCHh7oB </td> <td> LP6005051-DNA_D12 </td> </tr>
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

Running query:   RUNNING  2.5s
Running query:   RUNNING  3.2s
Running query:   RUNNING  3.8s
Running query:   RUNNING  4.4s
Running query:   RUNNING  5.0s
Running query:   RUNNING  5.7s
Running query:   RUNNING  6.3s
Running query:   RUNNING  6.9s
Running query:   RUNNING  7.6s
Running query:   RUNNING  8.2s
Running query:   RUNNING  8.8s
Running query:   RUNNING  9.5s
Running query:   RUNNING 10.1s
Running query:   RUNNING 10.7s
Running query:   RUNNING 11.4s
Running query:   RUNNING 12.0s
Running query:   RUNNING 12.6s
Running query:   RUNNING 13.2s
Running query:   RUNNING 13.9s
Running query:   RUNNING 14.5s
Running query:   RUNNING 15.1s
Running query:   RUNNING 15.7s
Running query:   RUNNING 16.4s
Running query:   RUNNING 17.0s
Running query:   RUNNING 17.6s
Running query:   RUNNING 18.3s
Running query:   RUNNING 18.9s
Running query:   RUNNING 19.5s
Running query:   RUNNING 20.2s
Running query:   RUNNING 20.8s
Running query:   RUNNING 21.4s
Running query:   RUNNING 22.0s
Running query:   RUNNING 22.7s
Running query:   RUNNING 23.3s
Running query:   RUNNING 23.9s
Running query:   RUNNING 24.5s
Running query:   RUNNING 25.2s
Running query:   RUNNING 25.8s
Running query:   RUNNING 26.4s
Running query:   RUNNING 27.1s
Running query:   RUNNING 27.7s
Running query:   RUNNING 28.3s
Running query:   RUNNING 29.0s
Running query:   RUNNING 29.7s
Running query:   RUNNING 30.3s
Running query:   RUNNING 30.9s
Running query:   RUNNING 31.5s
Running query:   RUNNING 32.2s
Running query:   RUNNING 33.0s
Running query:   RUNNING 33.7s
Running query:   RUNNING 34.3s
Running query:   RUNNING 34.9s
Running query:   RUNNING 35.5s
Running query:   RUNNING 36.2s
Running query:   RUNNING 36.8s
Running query:   RUNNING 37.4s
```

Displaying the results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Thu May 28 22:35:24 2015 -->
<table border=1>
<tr> <th> quantile </th> <th> row_num </th>  </tr>
  <tr> <td align="right"> 657.00 </td> <td align="right"> 1999 </td> </tr>
   </table>

Determine the cutoffs:

```r
maxChiSq = result$quantile
```
Cutoff: 657

Determine which genomes are outside our desired range

```r
values = list("_CUTOFF_" = maxChiSq)
sortAndLimit <- list("#_ORDER_BY_" = "LIMIT 1000")
result <- DisplayAndDispatchQuery("../sql/hwe-fail.sql",
                                  project=project,
                                  replacements=c(queryReplacements, values, sortAndLimit))
```

```
# Get all variants that have a chi squared value above a definited limit
SELECT
  variant_id,
  chisq,
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
  chisq > 657
LIMIT 1000


Running query:   RUNNING  2.5s
Running query:   RUNNING  3.2s
Running query:   RUNNING  3.8s
Running query:   RUNNING  4.4s
Running query:   RUNNING  5.1s
Running query:   RUNNING  5.7s
Running query:   RUNNING  6.3s
Running query:   RUNNING  7.0s
Running query:   RUNNING  7.6s
Running query:   RUNNING  8.2s
Running query:   RUNNING  8.9s
Running query:   RUNNING  9.5s
Running query:   RUNNING 10.1s
Running query:   RUNNING 10.7s
Running query:   RUNNING 11.4s
Running query:   RUNNING 12.0s
Running query:   RUNNING 12.7s
Running query:   RUNNING 13.4s
Running query:   RUNNING 14.0s
Running query:   RUNNING 14.6s
Running query:   RUNNING 15.3s
Running query:   RUNNING 15.9s
Running query:   RUNNING 16.6s
Running query:   RUNNING 17.2s
Running query:   RUNNING 17.8s
Running query:   RUNNING 18.5s
Running query:   RUNNING 19.1s
Running query:   RUNNING 19.7s
Running query:   RUNNING 20.4s
Running query:   RUNNING 21.1s
Running query:   RUNNING 21.7s
Running query:   RUNNING 22.4s
Running query:   RUNNING 23.0s
Running query:   RUNNING 23.6s
Running query:   RUNNING 24.2s
Running query:   RUNNING 24.9s
Running query:   RUNNING 25.5s
Running query:   RUNNING 26.1s
Running query:   RUNNING 26.8s
Running query:   RUNNING 27.4s
Running query:   RUNNING 28.0s
Running query:   RUNNING 28.6s
Running query:   RUNNING 29.3s
Running query:   RUNNING 29.9s
Running query:   RUNNING 30.5s
Running query:   RUNNING 31.2s
Running query:   RUNNING 31.8s
Running query:   RUNNING 32.5s
Running query:   RUNNING 33.1s
Running query:   RUNNING 33.7s
Running query:   RUNNING 34.3s
Running query:   RUNNING 35.0s
Running query:   RUNNING 35.6s
Running query:   RUNNING 36.2s
Running query:   RUNNING 36.9s
Running query:   RUNNING 37.5s
Running query:   RUNNING 38.1s
Running query:   RUNNING 38.8s
Running query:   RUNNING 39.4s
Running query:   RUNNING 40.1s
Running query:   RUNNING 40.8s
Running query:   RUNNING 41.4s
Running query:   RUNNING 42.0s
Running query:   RUNNING 42.7s
Running query:   RUNNING 43.3s
Running query:   RUNNING 43.9s
Running query:   RUNNING 44.6s
Running query:   RUNNING 45.2s
Running query:   RUNNING 45.8s
Running query:   RUNNING 46.4s
Running query:   RUNNING 47.1s
Running query:   RUNNING 47.7s
Running query:   RUNNING 48.3s
Running query:   RUNNING 49.0s
Running query:   RUNNING 49.6s
Running query:   RUNNING 50.2s
Running query:   RUNNING 50.8s
Running query:   RUNNING 51.5s
Running query:   RUNNING 52.2s
```

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Thu May 28 22:36:20 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> chisq </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBjo5rNCIIPIr7SCv8_QHA </td> <td align="right"> 877.00 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBiZ57NCIPm9gsCurJ-0jQE </td> <td align="right"> 878.00 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBjl57NCILHDy5yRh8rbUg </td> <td align="right"> 877.00 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBio6LNCIKG_l9SU9ZXYrgE </td> <td align="right"> 878.00 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyMhjdiqJzIP_HpKrtzev4Dg </td> <td align="right"> 851.00 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBisip8RINDgkvqZot_Qfg </td> <td align="right"> 881.00 </td> </tr>
   </table>







