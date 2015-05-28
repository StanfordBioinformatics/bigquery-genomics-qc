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
* [Ti/Tv by Depth](#titv-by-depth)
* [Blacklisted Variants](#blacklisted-variants)
* [Heterozygous Haplotype](#heterozygous-haplotype)
* [Ti/Tv by Genomic Window](#titv-by-genomic-window)
* [Ti/Tv by Dpeth](#titv-by-depth)
* [Hardy Weinberg Equilibrium](#hardy-weinberg-equilibrium)


```r
queryReplacements <- list("_THE_TABLE_"="va_aaa_pilot_data.all_genomes_gvcfs_20150514",
                          "_THE_EXPANDED_TABLE_"="va_aaa_pilot_data.all_genomes_expanded_vcfs_java3",
                          "_BLACKLISTED_TABLE_"="resources.blacklisted_positions")
sampleData <- read.csv("./data/patient_info.csv")
sampleInfo <- select(sampleData, call_call_set_name=Catalog.ID, gender=Gender)

# To run this against other public data, source in one of the dataset helpers.  For example:
# source("./rHelpers/pgpCGIOnlyDataset.R")
```

## Missingness Rate

For each variant, compute the missingness rate.  This query can be used to identify variants with a poor call rate.


```r
cutoff = list("_CUTOFF_"="0.9")
result <- DisplayAndDispatchQuery("../sql/variant-level-missingness-fail.sql",
                                  project=project,
                                  replacements=c(cutoff,
                                                 queryReplacements))
```

```
SELECT 
variant_id,
reference_name,
start,
end,
missingness_rate,
FROM (
  SELECT
  variant_id,
  reference_name,
  start,
  end,
  reference_bases,
  alternate_bases,
  no_calls,
  all_calls,
  (no_calls/all_calls) AS no_call_rate,
  1 - (all_calls-no_calls)/sample_count AS missingness_rate,
  sample_count
  FROM (
    SELECT
    variant_id,
    reference_name,
    start,
    end,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
    SUM(call.genotype == -1) WITHIN RECORD AS no_calls,
    COUNT(call.genotype) WITHIN RECORD AS all_calls,
    FROM
    [va_aaa_pilot_data.all_genomes_expanded_vcfs_java3]
    # Optionally add clause here to limit the query to a particular
    # region of the genome.
    #_WHERE_
  ) as.g
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
missingness_rate > 0.9
# Optionally add a clause here to sort and limit the results.
#_ORDER_BY_

Running query:   RUNNING  2.5s
Running query:   RUNNING  3.1s
Running query:   RUNNING  3.8s
Running query:   RUNNING  4.4s
Running query:   RUNNING  5.0s
```

```
Error: Response too large to return. Consider setting allowLargeResults to true in your job configuration. For more details, see https://cloud.google.com/bigquery/querying-data#largequeryresults

 responseTooLarge. Response too large to return. Consider setting allowLargeResults to true in your job configuration. For more details, see https://cloud.google.com/bigquery/querying-data#largequeryresults
```
Number of rows returned by this query: **100000**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Thu May 28 13:12:33 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> reference_name </th> <th> start </th> <th> reference_bases </th> <th> alternate_bases </th> <th> obs_hom1 </th> <th> obs_het </th> <th> obs_hom2 </th> <th> e_hom1 </th> <th> e_het </th> <th> e_hom2 </th> <th> chisq </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTUY1pHyLCCbx8GAh67pswU </td> <td> chr15 </td> <td align="right"> 94144726 </td> <td> C </td> <td> T </td> <td align="right"> 442 </td> <td align="right">  78 </td> <td align="right"> 359 </td> <td align="right"> 263 </td> <td align="right"> 435 </td> <td align="right"> 180 </td> <td align="right"> 592.38 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTUY7ZLyLCDct--Ys4zl4JIB </td> <td> chr15 </td> <td align="right"> 94144877 </td> <td> C </td> <td> T </td> <td align="right"> 442 </td> <td align="right">  79 </td> <td align="right"> 358 </td> <td align="right"> 263 </td> <td align="right"> 435 </td> <td align="right"> 179 </td> <td align="right"> 589.01 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTUYiZPyLCDCpebr0s3KmY4B </td> <td> chr15 </td> <td align="right"> 94144905 </td> <td> A </td> <td> G </td> <td align="right"> 506 </td> <td align="right"> 218 </td> <td align="right"> 155 </td> <td align="right"> 430 </td> <td align="right"> 369 </td> <td align="right">  79 </td> <td align="right"> 147.68 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBiyk-NBILWysJCYwbnFXw </td> <td> chrX </td> <td align="right"> 137939378 </td> <td> G </td> <td> A </td> <td align="right"> 418 </td> <td align="right">   4 </td> <td align="right">  19 </td> <td align="right"> 399 </td> <td align="right">  39 </td> <td align="right">   0 </td> <td align="right"> 357.21 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBiJlONBIPmqxc30_emK2gE </td> <td> chrX </td> <td align="right"> 137939465 </td> <td> C </td> <td> T </td> <td align="right"> 438 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right"> 437 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 439.00 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTAYrqK6KCDU5K-Dk9rN1vsB </td> <td> chr10 </td> <td align="right"> 84840750 </td> <td> T </td> <td> TA </td> <td align="right">   0 </td> <td align="right"> 221 </td> <td align="right"> 117 </td> <td align="right">  36 </td> <td align="right"> 148 </td> <td align="right"> 153 </td> <td align="right"> 79.74 </td> </tr>
   </table>

## Blacklisted Variants


```r
query <- "../sql/blacklisted-variants.sql"
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(queryReplacements))
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
Running query:   RUNNING  2.5s
Running query:   RUNNING  3.2s
Running query:   RUNNING  3.8s
Running query:   RUNNING  4.4s
Running query:   RUNNING  5.0s

Retrieving data:  3.8s
Retrieving data:  6.2s
Retrieving data:  8.0s
Retrieving data:  9.7s
Retrieving data: 11.0s
Retrieving data: 12.2s
Retrieving data: 13.3s
Retrieving data: 14.5s
Retrieving data: 15.8s
```

Number of rows returned by this query: **100000**.

First few results
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Thu May 28 13:12:59 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> reference_name </th> <th> start </th> <th> end </th> <th> Artifact_Type </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTIY7c7FECDHzY7mhr7Kzi8 </td> <td> chr12 </td> <td align="right"> 34695021 </td> <td align="right"> 34695022 </td> <td> centromeric_repeat </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTIYgc_FECCV283G08OrxR0 </td> <td> chr12 </td> <td align="right"> 34695041 </td> <td align="right"> 34695042 </td> <td> centromeric_repeat </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTIYk9DFECCr4qvI-qXomMcB </td> <td> chr12 </td> <td align="right"> 34695187 </td> <td align="right"> 34695188 </td> <td> centromeric_repeat </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTIYpdDFECDQ4Nntl_zC5-MB </td> <td> chr12 </td> <td align="right"> 34695205 </td> <td align="right"> 34695206 </td> <td> centromeric_repeat </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTIYotHFECDiwZjrlsCs_Z8B </td> <td> chr12 </td> <td align="right"> 34695330 </td> <td align="right"> 34695331 </td> <td> centromeric_repeat </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTIYhtLFECDgiPLk8qzd0Rg </td> <td> chr12 </td> <td align="right"> 34695430 </td> <td align="right"> 34695431 </td> <td> centromeric_repeat </td> </tr>
   </table>

## Heterozygous Haplotype
For each variant within the X and Y chromosome, identify heterozygous variants in male genomes.

First we use our sample information to determine which genomes are male.  

```r
maleSampleIds <- paste("'", filter(sampleInfo, gender == "M")$call_call_set_name, "'", sep="", collapse=",")
```


```r
sortAndLimit <- "ORDER BY reference_name, start, alternate_bases, call.call_set_name LIMIT 1000"
result <- DisplayAndDispatchQuery("../sql/sex-chromosome-heterozygous-haplotypes.sql",
                                  project=project,
                                  replacements=c("_MALE_SAMPLE_IDS_"=maleSampleIds,
                                                 "#_ORDER_BY_"=sortAndLimit,
                                                 queryReplacements))
```

```
# Retrieve heterozygous haplotype calls on chromosomes X and Y.
SELECT
  variant_id,
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  call.call_set_name,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype,
FROM
  [va_aaa_pilot_data.all_genomes_expanded_vcfs_java3]
WHERE
  reference_name IN ('chrX', 'chrY')
OMIT
  call if (2 > COUNT(call.genotype))
  OR EVERY(call.genotype <= 0)
  OR EVERY(call.genotype = 1)
HAVING call.call_set_name IN ('LP6005038-DNA_A01','LP6005038-DNA_A02','LP6005038-DNA_A04','LP6005038-DNA_A06','LP6005038-DNA_A07','LP6005038-DNA_A08','LP6005038-DNA_A09','LP6005038-DNA_A10','LP6005038-DNA_A11','LP6005038-DNA_A12','LP6005038-DNA_B01','LP6005038-DNA_B02','LP6005038-DNA_B03','LP6005038-DNA_B04','LP6005038-DNA_B05','LP6005038-DNA_B06','LP6005038-DNA_B07','LP6005038-DNA_B08','LP6005038-DNA_B09','LP6005038-DNA_B10','LP6005038-DNA_B11','LP6005038-DNA_B12','LP6005038-DNA_C01','LP6005038-DNA_C03','LP6005038-DNA_C04','LP6005038-DNA_C05','LP6005038-DNA_C07','LP6005038-DNA_C08','LP6005038-DNA_C09','LP6005038-DNA_C10','LP6005038-DNA_C11','LP6005038-DNA_D01','LP6005038-DNA_D02','LP6005038-DNA_D04','LP6005038-DNA_D05','LP6005038-DNA_D06','LP6005038-DNA_D07','LP6005038-DNA_D08','LP6005038-DNA_D09','LP6005038-DNA_D10','LP6005038-DNA_D11','LP6005038-DNA_E01','LP6005038-DNA_E02','LP6005038-DNA_E03','LP6005038-DNA_E04','LP6005038-DNA_E05','LP6005038-DNA_E06','LP6005038-DNA_E07','LP6005038-DNA_E08','LP6005038-DNA_E09','LP6005038-DNA_E10','LP6005038-DNA_E11','LP6005038-DNA_E12','LP6005038-DNA_F01','LP6005038-DNA_F02','LP6005038-DNA_F03','LP6005038-DNA_F04','LP6005038-DNA_F05','LP6005038-DNA_F06','LP6005038-DNA_F07','LP6005038-DNA_F08','LP6005038-DNA_F09','LP6005038-DNA_F10','LP6005038-DNA_F11','LP6005038-DNA_F12','LP6005038-DNA_G01','LP6005038-DNA_G02','LP6005038-DNA_G03','LP6005038-DNA_G04','LP6005038-DNA_G05','LP6005038-DNA_G06','LP6005038-DNA_G07','LP6005038-DNA_G08','LP6005038-DNA_G09','LP6005038-DNA_G10','LP6005038-DNA_G11','LP6005038-DNA_G12','LP6005038-DNA_H01','LP6005038-DNA_H02','LP6005038-DNA_H03','LP6005038-DNA_H05','LP6005038-DNA_H06','LP6005038-DNA_H07','LP6005038-DNA_H08','LP6005038-DNA_H09','LP6005038-DNA_H10','LP6005038-DNA_H11','LP6005038-DNA_H12','LP6005051-DNA_A01','LP6005051-DNA_A02','LP6005051-DNA_A03','LP6005051-DNA_A04','LP6005051-DNA_A05','LP6005051-DNA_A06','LP6005051-DNA_A07','LP6005051-DNA_A09','LP6005051-DNA_A10','LP6005051-DNA_A11','LP6005051-DNA_B01','LP6005051-DNA_B02','LP6005051-DNA_B03','LP6005051-DNA_B04','LP6005051-DNA_B05','LP6005051-DNA_B06','LP6005051-DNA_B07','LP6005051-DNA_B08','LP6005051-DNA_B09','LP6005051-DNA_B10','LP6005051-DNA_B12','LP6005051-DNA_C01','LP6005051-DNA_C02','LP6005051-DNA_C03','LP6005051-DNA_C04','LP6005051-DNA_C05','LP6005051-DNA_C06','LP6005051-DNA_C07','LP6005051-DNA_C08','LP6005051-DNA_C09','LP6005051-DNA_C10','LP6005051-DNA_C12','LP6005051-DNA_D01','LP6005051-DNA_D02','LP6005051-DNA_D03','LP6005051-DNA_D04','LP6005051-DNA_D05','LP6005051-DNA_D06','LP6005051-DNA_D08','LP6005051-DNA_D09','LP6005051-DNA_D10','LP6005051-DNA_D11','LP6005051-DNA_D12','LP6005051-DNA_E01','LP6005051-DNA_E03','LP6005051-DNA_E04','LP6005051-DNA_E05','LP6005051-DNA_E06','LP6005051-DNA_E07','LP6005051-DNA_E08','LP6005051-DNA_E09','LP6005051-DNA_E10','LP6005051-DNA_E11','LP6005051-DNA_E12','LP6005051-DNA_F02','LP6005051-DNA_F03','LP6005051-DNA_F04','LP6005051-DNA_F05','LP6005051-DNA_F06','LP6005051-DNA_F07','LP6005051-DNA_F08','LP6005051-DNA_F09','LP6005051-DNA_F10','LP6005051-DNA_F11','LP6005051-DNA_F12','LP6005051-DNA_G01','LP6005051-DNA_G02','LP6005051-DNA_G03','LP6005051-DNA_G04','LP6005051-DNA_G05','LP6005051-DNA_G06','LP6005051-DNA_G08','LP6005051-DNA_G09','LP6005051-DNA_G11','LP6005051-DNA_G12','LP6005051-DNA_H01','LP6005051-DNA_H02','LP6005051-DNA_H03','LP6005051-DNA_H04','LP6005051-DNA_H05','LP6005051-DNA_H06','LP6005051-DNA_H07','LP6005051-DNA_H08','LP6005051-DNA_H09','LP6005051-DNA_H11','LP6005051-DNA_H12','LP6005144-DNA_A03','LP6005144-DNA_A04','LP6005144-DNA_A05','LP6005144-DNA_A07','LP6005144-DNA_A08','LP6005144-DNA_A09','LP6005144-DNA_A10','LP6005144-DNA_A11','LP6005144-DNA_A12','LP6005144-DNA_B01','LP6005144-DNA_B02','LP6005144-DNA_B03','LP6005144-DNA_B04','LP6005144-DNA_B05','LP6005144-DNA_B07','LP6005144-DNA_B08','LP6005144-DNA_B09','LP6005144-DNA_B10','LP6005144-DNA_B11','LP6005144-DNA_B12','LP6005144-DNA_C01','LP6005144-DNA_C02','LP6005144-DNA_C03','LP6005144-DNA_C04','LP6005144-DNA_C05','LP6005144-DNA_C06','LP6005144-DNA_C07','LP6005144-DNA_C08','LP6005144-DNA_C09','LP6005144-DNA_C10','LP6005144-DNA_D01','LP6005144-DNA_D02','LP6005144-DNA_D03','LP6005144-DNA_D04','LP6005144-DNA_D05','LP6005144-DNA_D06','LP6005144-DNA_D07','LP6005144-DNA_D08','LP6005144-DNA_D09','LP6005144-DNA_D10','LP6005144-DNA_D11','LP6005144-DNA_D12','LP6005144-DNA_E01','LP6005144-DNA_E02','LP6005144-DNA_E03','LP6005144-DNA_E05','LP6005144-DNA_E06','LP6005144-DNA_E07','LP6005144-DNA_E08','LP6005144-DNA_E09','LP6005144-DNA_E10','LP6005144-DNA_E11','LP6005144-DNA_E12','LP6005144-DNA_F01','LP6005144-DNA_F03','LP6005144-DNA_F04','LP6005144-DNA_F05','LP6005144-DNA_F07','LP6005144-DNA_F08','LP6005144-DNA_F09','LP6005144-DNA_F10','LP6005144-DNA_F12','LP6005144-DNA_G02','LP6005144-DNA_G03','LP6005144-DNA_G04','LP6005144-DNA_G05','LP6005144-DNA_G06','LP6005144-DNA_G07','LP6005144-DNA_G08','LP6005144-DNA_G09','LP6005144-DNA_G10','LP6005144-DNA_G11','LP6005144-DNA_G12','LP6005144-DNA_H01','LP6005144-DNA_H02','LP6005144-DNA_H03','LP6005144-DNA_H04','LP6005144-DNA_H05','LP6005144-DNA_H06','LP6005144-DNA_H07','LP6005144-DNA_H08','LP6005144-DNA_H09','LP6005144-DNA_H10','LP6005144-DNA_H11','LP6005144-DNA_H12','LP6005243-DNA_A01','LP6005243-DNA_A03','LP6005243-DNA_A04','LP6005243-DNA_A05','LP6005243-DNA_A06','LP6005243-DNA_A07','LP6005243-DNA_A08','LP6005243-DNA_A10','LP6005243-DNA_A11','LP6005243-DNA_B01','LP6005243-DNA_B03','LP6005243-DNA_B04','LP6005243-DNA_B05','LP6005243-DNA_B06','LP6005243-DNA_B07','LP6005243-DNA_B08','LP6005243-DNA_B09','LP6005243-DNA_B10','LP6005243-DNA_B11','LP6005243-DNA_C03','LP6005243-DNA_C04','LP6005243-DNA_C05','LP6005243-DNA_C06','LP6005243-DNA_C07','LP6005243-DNA_C08','LP6005243-DNA_C09','LP6005243-DNA_C11','LP6005243-DNA_D03','LP6005243-DNA_D04','LP6005243-DNA_D05','LP6005243-DNA_D06','LP6005243-DNA_D07','LP6005243-DNA_D09','LP6005243-DNA_D10','LP6005243-DNA_D11','LP6005243-DNA_E01','LP6005243-DNA_E02','LP6005243-DNA_E03','LP6005243-DNA_E04','LP6005243-DNA_E05','LP6005243-DNA_E06','LP6005243-DNA_E07','LP6005243-DNA_E08','LP6005243-DNA_E09','LP6005243-DNA_E10','LP6005243-DNA_E11','LP6005243-DNA_F01','LP6005243-DNA_F03','LP6005243-DNA_F04','LP6005243-DNA_F05','LP6005243-DNA_F06','LP6005243-DNA_F07','LP6005243-DNA_F08','LP6005243-DNA_F09','LP6005243-DNA_F10','LP6005243-DNA_G01','LP6005243-DNA_G02','LP6005243-DNA_G03','LP6005243-DNA_G04','LP6005243-DNA_G05','LP6005243-DNA_G06','LP6005243-DNA_G07','LP6005243-DNA_G09','LP6005243-DNA_G10','LP6005243-DNA_H02','LP6005243-DNA_H03','LP6005243-DNA_H04','LP6005243-DNA_H05','LP6005243-DNA_H06','LP6005243-DNA_H07','LP6005243-DNA_H08','LP6005243-DNA_H09','LP6005243-DNA_H10','LP6005692-DNA_A01','LP6005692-DNA_A02','LP6005692-DNA_A03','LP6005692-DNA_A05','LP6005692-DNA_A06','LP6005692-DNA_A08','LP6005692-DNA_A09','LP6005692-DNA_A10','LP6005692-DNA_B01','LP6005692-DNA_B02','LP6005692-DNA_B03','LP6005692-DNA_B05','LP6005692-DNA_B07','LP6005692-DNA_B08','LP6005692-DNA_B09','LP6005692-DNA_B11','LP6005692-DNA_B12','LP6005692-DNA_C01','LP6005692-DNA_C02','LP6005692-DNA_C03','LP6005692-DNA_C05','LP6005692-DNA_C06','LP6005692-DNA_C07','LP6005692-DNA_C08','LP6005692-DNA_C09','LP6005692-DNA_C10','LP6005692-DNA_C12','LP6005692-DNA_D01','LP6005692-DNA_D02','LP6005692-DNA_D03','LP6005692-DNA_D05','LP6005692-DNA_D06','LP6005692-DNA_D07','LP6005692-DNA_D08','LP6005692-DNA_D10','LP6005692-DNA_D11','LP6005692-DNA_D12','LP6005692-DNA_E01','LP6005692-DNA_E02','LP6005692-DNA_E05','LP6005692-DNA_E06','LP6005692-DNA_E07','LP6005692-DNA_E08','LP6005692-DNA_E09','LP6005692-DNA_E10','LP6005692-DNA_E11','LP6005692-DNA_E12','LP6005692-DNA_F01','LP6005692-DNA_F02','LP6005692-DNA_F03','LP6005692-DNA_F04','LP6005692-DNA_F06','LP6005692-DNA_F07','LP6005692-DNA_F08','LP6005692-DNA_F09','LP6005692-DNA_F11','LP6005692-DNA_G01','LP6005692-DNA_G02','LP6005692-DNA_G03','LP6005692-DNA_G04','LP6005692-DNA_G05','LP6005692-DNA_G06','LP6005692-DNA_G07','LP6005692-DNA_G08','LP6005692-DNA_G09','LP6005692-DNA_G10','LP6005692-DNA_G11','LP6005692-DNA_H01','LP6005692-DNA_H02','LP6005692-DNA_H03','LP6005692-DNA_H04','LP6005692-DNA_H05','LP6005692-DNA_H06','LP6005692-DNA_H07','LP6005692-DNA_H09','LP6005692-DNA_H10','LP6005692-DNA_H11','LP6005693-DNA_A01','LP6005693-DNA_A02','LP6005693-DNA_A03','LP6005693-DNA_B01','LP6005693-DNA_C01','LP6005693-DNA_D01','LP6005693-DNA_E01','LP6005693-DNA_F01')
# Optionally add a clause here to sort and limit the results.
ORDER BY reference_name, start, alternate_bases, call.call_set_name LIMIT 1000

Running query:   RUNNING  2.6s
Running query:   RUNNING  3.2s
Running query:   RUNNING  3.9s
Running query:   RUNNING  4.5s
Running query:   RUNNING  5.2s
Running query:   RUNNING  5.9s
Running query:   RUNNING  6.5s
Running query:   RUNNING  7.1s
Running query:   RUNNING  7.8s
Running query:   RUNNING  8.4s
Running query:   RUNNING  9.0s
Running query:   RUNNING  9.6s
Running query:   RUNNING 10.2s
Running query:   RUNNING 10.9s
Running query:   RUNNING 11.5s
Running query:   RUNNING 12.1s
Running query:   RUNNING 12.8s
Running query:   RUNNING 13.4s
Running query:   RUNNING 14.0s
Running query:   RUNNING 14.6s
```
Number of rows returned by this query: **1000**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Thu May 28 13:13:16 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> reference_name </th> <th> start </th> <th> end </th> <th> reference_bases </th> <th> alternate_bases </th> <th> call_call_set_name </th> <th> genotype </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBjt36QBILSAtP6PhLTJLA </td> <td> chrX </td> <td align="right"> 2699245 </td> <td align="right"> 2699246 </td> <td> C </td> <td> A </td> <td> LP6005243-DNA_E07 </td> <td> 0,1 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBiK4KQBIJ-ggPnz7P7nhgE </td> <td> chrX </td> <td align="right"> 2699274 </td> <td align="right"> 2699275 </td> <td> T </td> <td> G </td> <td> LP6005243-DNA_A06 </td> <td> 0,1 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBjV4KQBIOW6iOn9mK6GHw </td> <td> chrX </td> <td align="right"> 2699349 </td> <td align="right"> 2699350 </td> <td> A </td> <td> T </td> <td> LP6005243-DNA_A06 </td> <td> 0,1 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBi54aQBILX1q4fYtpGd9wE </td> <td> chrX </td> <td align="right"> 2699449 </td> <td align="right"> 2699450 </td> <td> A </td> <td> C </td> <td> LP6005243-DNA_A06 </td> <td> 0,1 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBi54aQBILX1q4fYtpGd9wE </td> <td> chrX </td> <td align="right"> 2699449 </td> <td align="right"> 2699450 </td> <td> A </td> <td> C </td> <td> LP6005692-DNA_F04 </td> <td> 0,1 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBip5KQBIOzBmYaX-qX2qAE </td> <td> chrX </td> <td align="right"> 2699817 </td> <td align="right"> 2699818 </td> <td> A </td> <td> G </td> <td> LP6005692-DNA_F04 </td> <td> 0,1 </td> </tr>
   </table>

## Ti/Tv By Genomic Window

```r
query <- "../sql/titv-by-genomic-window-fail.sql"
max <- 3.0
min <- 1.5
cutoffs <- list("_MAX_" = max, "_MIN_" = min)
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(queryReplacements, cutoffs))
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
#_ORDER_
Retrieving data:  2.6s
Retrieving data:  3.4s
Retrieving data:  4.3s
Retrieving data:  5.2s
Retrieving data:  6.1s
Retrieving data:  6.9s
Retrieving data:  7.7s
```

Number of rows returned by this query: **100000**.

First few results
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Thu May 28 13:13:27 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> reference_name </th> <th> window_start </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTcY-OCdGCD6p5PL_72AvNkB </td> <td> chr17 </td> <td align="right"> 50800000 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTcYq-KdGCCh9KfbhPm_4yE </td> <td> chr17 </td> <td align="right"> 50800000 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTcYv-OdGCDk-b_0xOmlqYUB </td> <td> chr17 </td> <td align="right"> 50800000 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTcY6-OdGCDppvaNjJi93F4 </td> <td> chr17 </td> <td align="right"> 50800000 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTcYreSdGCCPoOTNtZTTz3I </td> <td> chr17 </td> <td align="right"> 50800000 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTcYxuSdGCDZlpGl5MuN0VE </td> <td> chr17 </td> <td align="right"> 50800000 </td> </tr>
   </table>

## Identify variants with average dpeth outside of defined range

```r
query <- "../sql/variant-depth-fail.sql"
max <- 150
min <- 10
cutoffs <- list("_MAX_" = max, "_MIN_" = min)
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(queryReplacements, cutoffs))
```

```
SELECT
variant_id,
reference_name,
start,
depth
FROM(
  SELECT
  variant_id,
  reference_name,
  start,
  ROUND(AVG(call.DP)) AS depth
  FROM
  [va_aaa_pilot_data.all_genomes_expanded_vcfs_java3]
  GROUP EACH BY 
  variant_id,
  reference_name,
  start)
WHERE
depth < 10 OR
depth > 150

Retrieving data:  2.4s
Retrieving data:  3.5s
Retrieving data:  4.6s
Retrieving data:  6.0s
Retrieving data:  7.2s
Retrieving data:  8.8s
Retrieving data:  9.7s
Retrieving data: 10.8s
```

Number of rows returned by this query: **100000**.

First few results
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Thu May 28 13:13:40 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> reference_name </th> <th> start </th> <th> depth </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWBiBkqdCIJ2M0qKRmcmZNw </td> <td> chrX </td> <td align="right"> 139053313 </td> <td align="right"> 9.00 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNxi56MNIII7xqp2CuJup-AE </td> <td> chr7 </td> <td align="right"> 152106041 </td> <td align="right"> 194.00 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIOY2hyVW5fZ2wwMDAyMTQYk6gBIKSGksLZkbyP_gE </td> <td> chrUn_gl000214 </td> <td align="right"> 21523 </td> <td align="right"> 190.00 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTAY0sqqFCCYof3bsKDd40M </td> <td> chr10 </td> <td align="right"> 42640722 </td> <td align="right"> 185.00 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTAYjO-bFCCSgfyCu-WxrawB </td> <td> chr10 </td> <td align="right"> 42399628 </td> <td align="right"> 243.00 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyWRj338MGIO789LHa76-JCg </td> <td> chrY </td> <td align="right"> 13692919 </td> <td align="right"> 4.00 </td> </tr>
   </table>

## Hardy Weinberg Equilibrium

Here we want to identify the variants that are out of Hardy Weinberg Equilibrium.  We want to remove the top 0.05 quantile of variants, so first we have to define what the cutoff for the chi squared value should be.

```r
quantile <- list("_QUANTILE_" = 19) # <- Define quantile by number. 
                                  # The 19th quantile selects the value that partitions the top 5% of values, 
                                  # assuming there are 20 quantiles.
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
    QUANTILES(chisq, 20) AS quantile
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
  row_num = 19

Running query:   RUNNING  2.5s
Running query:   RUNNING  3.4s
Running query:   RUNNING  4.0s
Running query:   RUNNING  4.6s
Running query:   RUNNING  5.2s
Running query:   RUNNING  5.9s
Running query:   RUNNING  6.5s
Running query:   RUNNING  7.1s
Running query:   RUNNING  7.8s
Running query:   RUNNING  8.4s
Running query:   RUNNING  9.0s
Running query:   RUNNING  9.6s
Running query:   RUNNING 10.2s
Running query:   RUNNING 10.9s
Running query:   RUNNING 11.5s
Running query:   RUNNING 12.1s
Running query:   RUNNING 12.7s
Running query:   RUNNING 13.3s
Running query:   RUNNING 14.0s
Running query:   RUNNING 14.6s
Running query:   RUNNING 15.2s
Running query:   RUNNING 15.8s
Running query:   RUNNING 16.4s
Running query:   RUNNING 17.1s
Running query:   RUNNING 17.7s
Running query:   RUNNING 18.3s
Running query:   RUNNING 18.9s
Running query:   RUNNING 19.5s
Running query:   RUNNING 20.2s
Running query:   RUNNING 20.8s
Running query:   RUNNING 21.5s
Running query:   RUNNING 22.1s
Running query:   RUNNING 22.7s
Running query:   RUNNING 23.3s
Running query:   RUNNING 24.0s
Running query:   RUNNING 24.6s
Running query:   RUNNING 25.2s
Running query:   RUNNING 25.9s
Running query:   RUNNING 26.5s
Running query:   RUNNING 27.1s
Running query:   RUNNING 27.7s
Running query:   RUNNING 28.3s
Running query:   RUNNING 29.0s
Running query:   RUNNING 29.6s
Running query:   RUNNING 30.2s
Running query:   RUNNING 30.8s
Running query:   RUNNING 31.5s
Running query:   RUNNING 32.1s
Running query:   RUNNING 32.7s
Running query:   RUNNING 33.3s
Running query:   RUNNING 33.9s
Running query:   RUNNING 34.6s
Running query:   RUNNING 35.2s
Running query:   RUNNING 35.8s
Running query:   RUNNING 36.5s
Running query:   RUNNING 37.1s
Running query:   RUNNING 37.7s
Running query:   RUNNING 38.3s
Running query:   RUNNING 39.0s
Running query:   RUNNING 39.6s
Running query:   RUNNING 40.2s
Running query:   RUNNING 40.8s
Running query:   RUNNING 41.5s
Running query:   RUNNING 42.1s
Running query:   RUNNING 42.7s
Running query:   RUNNING 43.3s
Running query:   RUNNING 44.0s
Running query:   RUNNING 44.6s
Running query:   RUNNING 45.2s
Running query:   RUNNING 45.8s
Running query:   RUNNING 46.4s
Running query:   RUNNING 47.1s
Running query:   RUNNING 47.7s
Running query:   RUNNING 48.3s
Running query:   RUNNING 48.9s
Running query:   RUNNING 49.5s
Running query:   RUNNING 50.2s
Running query:   RUNNING 50.8s
Running query:   RUNNING 51.4s
Running query:   RUNNING 52.0s
Running query:   RUNNING 52.7s
Running query:   RUNNING 53.3s
Running query:   RUNNING 53.9s
Running query:   RUNNING 54.6s
Running query:   RUNNING 55.2s
Running query:   RUNNING 55.8s
Running query:   RUNNING 56.4s
Running query:   RUNNING 57.0s
Running query:   RUNNING 57.7s
Running query:   RUNNING 58.3s
Running query:   RUNNING 58.9s
Running query:   RUNNING 59.5s
Running query:   RUNNING 60.2s
Running query:   RUNNING 60.8s
Running query:   RUNNING 61.4s
Running query:   RUNNING 62.0s
Running query:   RUNNING 62.6s
Running query:   RUNNING 63.2s
Running query:   RUNNING 63.9s
Running query:   RUNNING 64.5s
Running query:   RUNNING 65.1s
Running query:   RUNNING 65.7s
Running query:   RUNNING 66.3s
Running query:   RUNNING 66.9s
Running query:   RUNNING 67.6s
Running query:   RUNNING 68.2s
Running query:   RUNNING 68.8s
Running query:   RUNNING 69.4s
Running query:   RUNNING 70.1s
Running query:   RUNNING 70.7s
Running query:   RUNNING 71.3s
Running query:   RUNNING 72.0s
Running query:   RUNNING 72.6s
Running query:   RUNNING 73.2s
Running query:   RUNNING 73.9s
Running query:   RUNNING 74.5s
Running query:   RUNNING 75.1s
Running query:   RUNNING 75.7s
Running query:   RUNNING 76.3s
```

Displaying the results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Thu May 28 13:14:59 2015 -->
<table border=1>
<tr> <th> quantile </th> <th> row_num </th>  </tr>
  <tr> <td align="right"> 36.68 </td> <td align="right">  19 </td> </tr>
   </table>

Determine the cutoffs:

```r
maxChiSq = result$quantile
```
Cutoff: 36.6805293

Determine which genomes are outside our desired range

```r
values = list("_CUTOFF_" = maxChiSq)
result <- DisplayAndDispatchQuery("../sql/hwe-fail.sql",
                                  project=project,
                                  replacements=c(queryReplacements, values))
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
  chisq > 36.6805293005671
GROUP BY
  variant_id,
  reference_name,
  start,
  reference_bases,
  alternate_bases,
  obs_hom1,
  obs_het,
  obs_hom2,
  e_hom1,
  e_het,
  e_hom2,
  chisq

Running query:   RUNNING  2.5s
Running query:   RUNNING  3.1s
Running query:   RUNNING  3.7s
Running query:   RUNNING  4.4s
Running query:   RUNNING  5.0s
Running query:   RUNNING  5.6s
Running query:   RUNNING  6.3s
Running query:   RUNNING  6.9s
Running query:   RUNNING  7.5s
Running query:   RUNNING  8.5s
Running query:   RUNNING  9.1s
Running query:   RUNNING  9.7s
Running query:   RUNNING 10.3s
Running query:   RUNNING 11.0s
Running query:   RUNNING 11.6s
Running query:   RUNNING 12.2s
Running query:   RUNNING 12.8s
Running query:   RUNNING 13.5s
Running query:   RUNNING 14.1s
Running query:   RUNNING 14.7s
Running query:   RUNNING 15.4s
Running query:   RUNNING 16.0s
Running query:   RUNNING 16.7s
Running query:   RUNNING 17.4s
Running query:   RUNNING 19.6s
Running query:   RUNNING 20.2s
Running query:   RUNNING 20.8s
Running query:   RUNNING 21.4s
Running query:   RUNNING 22.1s
Running query:   RUNNING 22.7s
Running query:   RUNNING 23.4s
Running query:   RUNNING 24.6s
Running query:   RUNNING 25.2s
Running query:   RUNNING 25.8s
Running query:   RUNNING 26.4s
Running query:   RUNNING 27.1s
Running query:   RUNNING 27.7s
Running query:   RUNNING 28.3s
Running query:   RUNNING 28.9s
Running query:   RUNNING 29.5s
Running query:   RUNNING 30.2s
Running query:   RUNNING 30.8s
Running query:   RUNNING 31.4s
Running query:   RUNNING 32.0s
Running query:   RUNNING 32.6s
Running query:   RUNNING 33.3s
Running query:   RUNNING 33.9s
Running query:   RUNNING 34.5s
Running query:   RUNNING 35.1s
Running query:   RUNNING 35.8s
Running query:   RUNNING 36.6s
Running query:   RUNNING 37.3s
Running query:   RUNNING 37.9s
Running query:   RUNNING 38.5s

Retrieving data:  2.0s
Retrieving data:  2.8s
Retrieving data:  3.5s
Retrieving data:  4.2s
Retrieving data:  4.9s
Retrieving data:  5.5s
Retrieving data:  6.2s
```

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Thu May 28 13:15:47 2015 -->
<table border=1>
<tr> <th> variant_id </th> <th> chisq </th>  </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTYYjsX8DSCKgrHKx86pwTM </td> <td align="right"> 47.59 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTYYvMX8DSCtivOzqa2ugW8 </td> <td align="right"> 343.18 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTYY58X8DSCpnfWt_b3LqNgB </td> <td align="right"> 54.83 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTYY-MX8DSDih96Okf6ioxI </td> <td align="right"> 295.78 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIFY2hyMTIY3P_TLCDjnaGo8NGZkPoB </td> <td align="right"> 153.03 </td> </tr>
  <tr> <td> CJ_JqPj1p5_yIRIEY2hyNRiu-r1BIKbDqZj1gZf1qQE </td> <td align="right"> 58.27 </td> </tr>
   </table>







