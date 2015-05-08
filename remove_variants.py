import argparse
from BigQueryGenomicsQC import GenomicsQC

def RunQC(verbose= False, sample_id=None, sample_file=None, variant_file=None, client_secrets=None,
          project_number=None, dataset=None, variant_table=None, expanded_table=None):

    qc = GenomicsQC(verbose=verbose, client_secrets=client_secrets, project_number=project_number, dataset=dataset,
                    variant_table=variant_table, expanded_table=expanded_table)

    if sample_id is not None:
        qc.remove_sample(sample_id)

    if sample_file is not None:
        qc.remove_samples_from_file(sample_file)

    if variant_file is not None:
        qc.remove_positions_from_file(variant_file)

def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'This script removes ')

    parser.add_argument("--sample_id", default=None,
                                help="Sample name to remove from dataset. To remove multiple samples use --sample_file")
    parser.add_argument("--sample_file", default=None,
                                help="File containing list of sample names to remove from dataset. One per line")
    parser.add_argument("--variant_file", default=None,
                                help="File containing list of positions to remove from dataset. BED format required.")
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
                                help="OPTIONAL. Logs will be very detailed.")

    options = parser.parse_args()
    if options.sample_id is None and options.sample_file is None and options.variant_file is None:
        print "Exiting, nothing to remove.\nSpecify samples, variants, or both."
        exit(0)
    return options

if __name__ == "__main__":
    options = parse_command_line()
    RunQC(verbose=options.verbose, sample_id=options.sample_id, sample_file=options.sample_file,
          variant_file=options.variant_file, client_secrets=options.client_secrets,
          project_number=options.project_number, dataset=options.dataset,
          variant_table=options.variant_table, expanded_table=options.expanded_table)