#!/usr/bin/python
import argparse
import json
from BigQueryGenomicsQC import GenomicsQC
from BigQueryGenomicsQC import QCSteps

def RunQC(verbose= False, sample_level=False, remove_samples=False, variant_level=False, client_secrets=None, project_number=None,
          dataset=None, variant_table=None, expanded_table=None, poll=False, qc_step=None):


    qc = GenomicsQC(verbose=verbose, client_secrets=client_secrets, project_number=project_number, dataset=dataset,
                    variant_table=variant_table, expanded_table=expanded_table)

    if sample_level is True:
        qc.sample_qc(remove=remove_samples)
    if variant_level is True:
        qc.variant_qc(poll)

    if qc_step is not None:
        qc.custom_list(qc_step)

def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'This script runs qc queries on genomics data stored in BigQuery.  Failing samples and variants'
                      'will be removed from the variantset.  Sample level qc and variant level qc can be run at the '
                      'same time.  See config.py to set dataset and table names.')

    parser.add_argument("--sample_qc", action='store_true', default=False,
                                help="Run sample level qc.")
    parser.add_argument("--variant_qc", action='store_true', default=False,
                                help="Run variant level qc.")
    parser.add_argument("--qc", nargs='+', default=None,
                                help="Pass a list of qc steps to run.  Valid qc steps include: gender_check, "
                                     "genotyping_concordance, heterozygosity_rate, inbreeding_coefficient, "
                                     "missingness_rate, singletons, blacklisted, hardy_weinberg, heterozygous_haplotype,"
                                     "titv_by_depth, titv_by_genomic_window.")
    parser.add_argument("--remove-samples", action='store_true', default=False,
                                help="Remove samples from dataset.")
    parser.add_argument("--variant_table", default=None,
                                help="OPTIONAL. Variant table to query. Defaults to value in config.py")
    parser.add_argument("--expanded_table", default=None,
                                help="OPTIONAL. Expanded variant table to query. Defaults to value in config.py")
    parser.add_argument("--project_number", default=None,
                                help="OPTIONAL. Google Cloud project number. Defaults to value in config.py")
    parser.add_argument("--dataset", default=None,
                                help="OPTIONAL. Variant store dataset. Defaults to value in config.py")
    parser.add_argument("--client_secrets", default=None,
                                help="OPTIONAL. client_secrets.json. Defaults to value in config.py")
    parser.add_argument("--verbose", action='store_true', default=False,
                                help="OPTIONAL. Logs will be very detailed. Kind of noise regardless thanks to Google"
                                     "API Client.")
    parser.add_argument("--poll", action='store_true', default=False,
                                help="OPTIONAL. Wait for submitted jobs to finish.  Only applies to variant qc.")

    options = parser.parse_args()

    if options.sample_qc is False and options.variant_qc is False and options.qc is None:
        print "Exiting, no qc specified.\nSpecify sample qc, variant qc, or both.\n--sample_qc and/or --variant_qc." \
              "Individual QC steps can be run using --qc"
        exit(0)
    return options

if __name__ == "__main__":
    options = parse_command_line()
    RunQC(sample_level=options.sample_qc, variant_level=options.variant_qc, verbose=options.verbose,
          client_secrets=options.client_secrets, project_number=options.project_number, dataset=options.dataset,
          variant_table=options.variant_table, expanded_table=options.expanded_table,
          remove_samples=options.remove_samples, poll=options.poll, qc_step=options.qc)