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
          [_THE_EXPANDED_TABLE_]
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
titv_ratio > _MAX_ 
OR titv_ratio < _MIN_
#_ORDER_BY_

