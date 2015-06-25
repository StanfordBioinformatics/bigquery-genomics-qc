
class Queries(object):
    GENDER_CHECK = "gender-check-fail.sql"
    GENOTYPING_CONCORDANCE = "genotyping-conconcordance-fail.sql"
    INBREEDING_COEFFICIENT = "inbreeding-coefficient-fail.sql"
    INBREEDING_COEFFICIENT_METRICS = "inbreeding-coefficient-metrics.sql"
    INBREEDING_COEFFICIENT_MAIN = "homozygous-variants.sql"
    MISSINGNESS_SAMPLE_LEVEL = "missingness-sample-level-fail.sql"
    PRIVATE_VARIANTS = "private-variants-fail.sql"
    PRIVATE_VARIANT_METRICS = "private-variants-metrics.sql"
    PRIVATE_VARIANT_MAIN = "private-variants.sql"
    VARIANT_DEPTH = "variant-depth-fail.sql"
    MISSINGNESS_VARIANT_LEVEL = "variant-level-missingness-fail.sql"
    HETEROZYGOSITY = "heterozygous-call-counts-fail.sql"
    HETEROZYGOSITY_METRICS = "heterozygous-call-counts-metrics.sql"
    HETEROZYGOSITY_MAIN = "heterozygous-calls-count.sql"

    SAMPLE_LEVEL_QC_QUERIES = [GENDER_CHECK,
                               #GENOTYPING_CONCORDANCE,
                               INBREEDING_COEFFICIENT,
                               MISSINGNESS_SAMPLE_LEVEL,
                               PRIVATE_VARIANTS,
                               HETEROZYGOSITY
    ]

    VARIANT_LEVEL_QC_QUERIES = [VARIANT_DEPTH,
                                MISSINGNESS_VARIANT_LEVEL]

    VARIANT_LEVEL_REMOVAL = {
        VARIANT_DEPTH: 'sample_call',
        MISSINGNESS_VARIANT_LEVEL: 'position'
    }

    AVERAGE_STDDEV = {
        PRIVATE_VARIANTS: PRIVATE_VARIANT_METRICS,
        INBREEDING_COEFFICIENT: INBREEDING_COEFFICIENT_METRICS,
        HETEROZYGOSITY: HETEROZYGOSITY_METRICS
    }

    MAIN_QUERY = {
        PRIVATE_VARIANTS: PRIVATE_VARIANT_MAIN,
        PRIVATE_VARIANT_METRICS: PRIVATE_VARIANT_MAIN,
        INBREEDING_COEFFICIENT: INBREEDING_COEFFICIENT_MAIN,
        INBREEDING_COEFFICIENT_METRICS: INBREEDING_COEFFICIENT_MAIN,
        HETEROZYGOSITY: HETEROZYGOSITY_MAIN,
        HETEROZYGOSITY_METRICS: HETEROZYGOSITY_MAIN
    }

    PRESET_CUTOFFS = {
        GENOTYPING_CONCORDANCE: {"_CUTOFF_": "0.95"},
        GENDER_CHECK: {"_MALE_CUTOFF_": "0.2",
                       "_FEMALE_CUTOFF_": "0.5"},
        MISSINGNESS_SAMPLE_LEVEL: {"_CUTOFF_": "0.9"},
        MISSINGNESS_VARIANT_LEVEL: {"_CUTOFF_": "0.9"},
        VARIANT_DEPTH: {"_MAX_VALUE_": "70",
                        "_MIN_VALUE_": "8"}
    }
