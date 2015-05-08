import argparse
from BigQueryGenomicsQC import GenomicsQC

def RunQC(sample_level=False, variant_level=False, verbose=False):
    qc = GenomicsQC(verbose=verbose)
    if sample_level is True:
        qc.sample_qc()
    if variant_level is True:
        qc.variant_qc()

def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'This script runs qc queries on genomics data stored in BigQuery.  Failing samples and variants'
                      'will be removed from the variantset.  Sample level qc and variant level qc can be run at the '
                      'same time.  See config.py to set dataset and table names.')

    parser.add_argument("--sample_qc", default=False,
                                help="Run sample level qc.")
    parser.add_argument("--variant_qc", default=False,
                                help="Run variant level qc.")
    parser.add_argument("--verbose", default=False,
                                help="Logs will be very detailed.")

    options = parser.parse_args()
    if options.sample_qc is False and options.variant_qc is False:
        print "Exiting, no qc specified.\nSpecify sample qc, variant qc, or both.\n--sample_qc and/or --variant_qc"
        exit(0)
    return options

if __name__ == "__main__":
    options = parse_command_line()
    RunQC(sample_level=options.sample_qc, variant_level=options.variant_qc, verbose=options.verbose)