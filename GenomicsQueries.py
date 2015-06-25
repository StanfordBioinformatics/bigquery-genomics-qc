from config import Config

class Queries(object):
    # Sample level qc queries
    GENDER_CHECK = "gender-check-fail.sql"
    GENOTYPING_CONCORDANCE = "genotyping-conconcordance-fail.sql"
    HETEROZYGOSITY = "heterozygous-call-counts-fail.sql"
    HETEROZYGOSITY_METRICS = "heterozygous-call-counts-metrics.sql"
    HETEROZYGOSITY_MAIN = "heterozygous-calls-count.sql"
    INBREEDING_COEFFICIENT = "inbreeding-coefficient-fail.sql"
    INBREEDING_COEFFICIENT_METRICS = "inbreeding-coefficient-metrics.sql"
    INBREEDING_COEFFICIENT_MAIN = "homozygous-variants.sql"
    MISSINGNESS_SAMPLE_LEVEL = "missingness-sample-level-fail.sql"
    PRIVATE_VARIANTS = "private-variants-fail.sql"
    PRIVATE_VARIANT_METRICS = "private-variants-metrics.sql"
    PRIVATE_VARIANT_MAIN = "private-variants.sql"

    # Variant level qc queries
    BLACKLISTED = "blacklisted-variants.sql"
    HARDY_WEINBERG = "hwe-fail.sql"
    HARDY_WEINBERG_METRICS = "hwe-quantile.sql"
    HETERZYGOUS_HAPLOTYPE = "sex-chromosome-heterozygous-haplotypes.sql"
    MISSINGNESS_VARIANT_LEVEL = "variant-level-missingness-fail.sql"
    TITV_DEPTH = "titv-by-depth-fail.sql"
    TITV_GENOMIC_WINDOW = "titv-by-genomic-window-fail.sql"

    SAMPLE_LEVEL_QC_QUERIES = [GENDER_CHECK,
                               GENOTYPING_CONCORDANCE,
                               INBREEDING_COEFFICIENT,
                               MISSINGNESS_SAMPLE_LEVEL,
                               PRIVATE_VARIANTS,
                               HETEROZYGOSITY]

    VARIANT_LEVEL_QC_QUERIES = [BLACKLISTED,
                                HARDY_WEINBERG,
                                HETERZYGOUS_HAPLOTYPE,
                                MISSINGNESS_VARIANT_LEVEL,
                                TITV_DEPTH,
                                TITV_GENOMIC_WINDOW]

    VARIANT_LEVEL_REMOVAL = {
        BLACKLISTED: 'position',
        HARDY_WEINBERG: 'position',
        HETERZYGOUS_HAPLOTYPE: 'sample_call',
        MISSINGNESS_VARIANT_LEVEL: 'position',
        TITV_DEPTH: 'sample_call',
        TITV_GENOMIC_WINDOW: 'position'}

    AVERAGE_STDDEV = {
        PRIVATE_VARIANTS: PRIVATE_VARIANT_METRICS,
        INBREEDING_COEFFICIENT: INBREEDING_COEFFICIENT_METRICS,
        HETEROZYGOSITY: HETEROZYGOSITY_METRICS}

    CUTOFF = {
        HARDY_WEINBERG: HARDY_WEINBERG_METRICS
    }

    MAIN_QUERY = {
        PRIVATE_VARIANTS: PRIVATE_VARIANT_MAIN,
        PRIVATE_VARIANT_METRICS: PRIVATE_VARIANT_MAIN,
        INBREEDING_COEFFICIENT: INBREEDING_COEFFICIENT_MAIN,
        INBREEDING_COEFFICIENT_METRICS: INBREEDING_COEFFICIENT_MAIN,
        HETEROZYGOSITY: HETEROZYGOSITY_MAIN,
        HETEROZYGOSITY_METRICS: HETEROZYGOSITY_MAIN}

    PRESET_CUTOFFS = {
        GENOTYPING_CONCORDANCE: {"_CUTOFF_": "0.95"},
        GENDER_CHECK: {"_MALE_CUTOFF_": "0.2",
                       "_FEMALE_CUTOFF_": "0.5"},
        MISSINGNESS_SAMPLE_LEVEL: {"_CUTOFF_": "0.9"},
        MISSINGNESS_VARIANT_LEVEL: {"_CUTOFF_": "0.9"},
        TITV_DEPTH: {"_MAX_": "3.0",
                              "_MIN_": "1.5"},
        TITV_GENOMIC_WINDOW: {"_MAX_": "3.0",
                              "_MIN_": "1.5"},
        BLACKLISTED: {"_BLACKLISTED_TABLE_": Config.BLACKLISTED_TABLE},
        HARDY_WEINBERG_METRICS: {"_QUANTILE_": "1999"}

    }
