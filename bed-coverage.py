#!/usr/bin/python
import argparse
from config import Config
from GenomicsFunctions import GenomicsUtils

def GenomeCoverage(bed_file, sample_id=None, project_number=None, dataset=None, client_secrets=None, verbose=False):
    if project_number is None:
        project_number = Config.PROJECT_NUMBER
    if dataset is None:
        dataset = Config.DATASET
    if client_secrets is None:
        client_secrets = Config.CLIENT_SECRETS
    gf = GenomicsUtils(verbose, client_secrets, project_number, dataset)
    gf.bed_coverage(bed_file)



def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'This script removes ')

    parser.add_argument("--bed_file", default=None,
                                help="Bed file containing genomics regions for which to calculate coverage")
    parser.add_argument("--sample_id", default=None,
                                help="OPTIONAL. Sample id to compute coverages for.  Defaults to calculate coverage for all samples"
                                     "in the dataset.")
    parser.add_argument("--project_number", default=None,
                                help="OPTIONAL. Google Cloud project number. Defaults to value in config.py")
    parser.add_argument("--dataset", default=None,
                                help="OPTIONAL. Read store dataset to compute coverages. Defaults to value in config.py")
    parser.add_argument("--client_secrets", default=None,
                                help="OPTIONAL. client_secrets.json. Defaults to value in config.py")
    parser.add_argument("--verbose", action='store_true', default=False,
                                help="OPTIONAL. Logs will be very detailed. Kind of noise regardless thanks to Google"
                                     "API Client.")

    options = parser.parse_args()
    if options.bed_file is None:
        print "Exiting, no bed file provided.  Use --bed_file to define a bed file.\n"
        exit(0)

    return options

if __name__ == "__main__":
    options = parse_command_line()
    GenomeCoverage(verbose=options.verbose, sample_id=options.sample_id,
          client_secrets=options.client_secrets, bed_file=options.bed_file,
          project_number=options.project_number, dataset=options.dataset)