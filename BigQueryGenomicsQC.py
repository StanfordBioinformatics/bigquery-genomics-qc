import os
import sys
from GenomicsQueries import Queries
from BigQueryClient import BigQuery
from config import Config
from GoogleGenomicsClient import GoogleGenomicsClient
import logging

class GenomicsQC(object):
    def __init__(self, verbose=False, client_secrets=None, project_number=None, dataset=None, variant_table=None,
                 expanded_table=None):

        self.query_repo = Config.QUERY_REPO
        if variant_table is None:
            variant_table = Config.VARIANT_TABLE
        self.variant_table = variant_table

        if expanded_table is None:
            expanded_table = Config.EXPANDED_TABLE
        self.expanded_table = expanded_table

        if client_secrets is None:
            client_secrets = Config.CLIENT_SECRETS
        self.client_secrets_path = client_secrets

        if project_number is None:
            project_number = Config.PROJECT_NUMBER
        self.project_number = project_number

        if dataset is None:
            dataset = Config.DATASET
        self.dataset = dataset

        # Set up logging
        self.setup_log(verbose)

        # Set up API clients
        self.bq = BigQuery(project_number=self.project_number,
                           client_secrets=self.client_secrets_path)

        self.gg = GoogleGenomicsClient(client_secrets=self.client_secrets_path,
                                       project_number=self.project_number,
                                       dataset=self.dataset)

        # Empty dictionaries for failures
        self.failed_samples = {}
        self.failed_positions = {}
        self.failed_sample_calls = {}

    #### Specific types of QC functions ####
    # Sample level QC
    def sample_qc(self):
        logging.info("Running Sample Level QC")
        queries = Queries.SAMPLE_LEVEL_QC_QUERIES
        for q in queries:
            logging.debug(q)
            result = self.run_analysis(query_file=q, level='sample')
            self.collect_failed_samples(result, q)
        self.remove_failed_samples(self.failed_samples)

    # Variant level QC
    def variant_qc(self):
        logging.info("Running Variant Level QC")
        queries = Queries.VARIANT_LEVEL_QC_QUERIES
        for q in queries:
            logging.debug(q)
            failed = self.run_analysis(query_file=q, level='variant')
            # Determine if entire position needs to be removed or remove call for a specific sample
            if Queries.VARIANT_LEVEL_REMOVAL[q] == 'sample_call':
                self.collect_failed_sample_calls(result, query)
            elif Queries.VARIANT_LEVEL_REMOVAL[q] == 'position':
                self.collect_failed_positions(result, query)
        self.remove_failed_sample_calls(self.failed_sample_calls)
        self.remove_failed_positions(self.failed_positions)

    # Run all required analysis for each query
    def run_analysis(self, query_file, level):
        cutoffs = None
        # Check if this query requires cutoffs to be defined by average values
        if query_file in Queries.AVERAGE_STDDEV:
            prequery = Queries.AVERAGE_STDDEV[query_file]
            cutoffs = self.average_stddev_cutoffs(prequery)
        query = self.prepare_query(query_file, preset_cutoffs=cutoffs)
        result = self.bq.run(query)
        return result

    #### Query set up ####
    # Set up the query, read it in, apply substitutions
    def prepare_query(self, query_file, preset_cutoffs=None):
        raw_query = self.get_query(query_file)
        if preset_cutoffs is None:
            preset_cutoffs = self.get_preset_cutoffs(query_file)
        query = self.query_substitutions(raw_query, other=preset_cutoffs)
        return query

    # Read raw query in from file
    def get_query(self, file):
        path = os.path.join(self.query_repo, file)
        query = ''
        with open (path, "r") as f:
            query = f.read()
        return query

    # Apply any substitutions. Substitutions set on other must be in a dictionary
    def query_substitutions(self, query, other=None):
        replacements = {
            "_THE_TABLE_": Config.VARIANT_TABLE,
            "_THE_EXPANDED_TABLE_": Config.EXPANDED_TABLE,
            "_PATIENT_INFO_": Config.PATIENT_INFO,
            "_GENOTYPING_TABLE_": Config.GENOTYPING_TABLE
        }
        for r in replacements:
            query = query.replace(r, replacements[r])
        if other is not None:
            for r in other:
                query = query.replace(r, other[r])
        return query

    # Get preset cutoffs from query file
    def get_preset_cutoffs(self, query):
        cutoffs = None
        if query in Queries.PRESET_CUTOFFS:
            cutoffs = Queries.PRESET_CUTOFFS[query]
        return cutoffs

    # Run metrics query to define cutoffs based on average and standard deviation values
    def average_stddev_cutoffs(self, query_file):
        query = self.prepare_query(query_file)
        result = self.bq.run(query)
        average, stddev = self.get_average_stddev(result)
        max, min = self.calculate_max_min(average, stddev)
        substitutions = self.create_max_min_substitutions(max, min)
        return substitutions

    # Get average and standard deviation values from a parsed BigQuery result
    def get_average_stddev(self, result):
        for r in result:
            average = r['average']
            stddev = r['stddev']
            return average, stddev

    # Calculate maximum and minimum values based on average and standard deviation.
    # Cutoffs are defined as more than three standard deviations away from the average.
    def calculate_max_min(self, average, stddev):
        max = float(average) + (3 * float(stddev))
        min = float(average) - (3 * float(stddev))
        return max, min

    # Create maximum and minimum cutoff dictionary
    def create_max_min_substitutions(self, max, min):
        dict = {
            "_MAX_VALUE_": "%s" % max,
            "_MIN_VALUE_": "%s" % min,
        }
        return dict

    #### Functions to accumulate failures ####
    # Get failed sample ids
    def get_failed_samples(self, result):
        if result:
            failed_ids = []
            for r in result:
                failed_ids.append(r['sample_id'])
            return failed_ids
        return None

    # Get failed positions
    def get_failed_positions(self, result):
        if result:
            failed_positions = []
            for r in result:
                position = {
                    "chr": r['reference_name'],
                    "start": r['start'],
                    "end": r['end']
                }
                failed_positions.append(position)
            return failed_positions
        return None

    # Add failed samples to total
    def collect_failed_samples(self, result, query):
        if result is None:
            return
        for r in result:
            sample_id = r['sample_id']
            if sample_id in self.failed_samples:
                self.failed_samples[sample_id].append(query)
            else:
                self.failed_samples[sample_id] = [query]

    # Add failed positions to total
    def collect_failed_positions(self, result, query):
        if result is None:
            return
        for r in result:
            position = "/".join([r['reference_name'], r['start'], r['end']])
            if position in self.failed_positions:
                self.failed_positions[position].append(query)
            else:
                self.failed_positions[position] = [query]

    # Add failed sample calls to total
    def collect_failed_sample_calls(self, result, query):
        if result is None:
            return
        for r in result:
            position = "/".join([r['sample_id'], r['reference_name'], r['start'], r['end']])
            if position in self.failed_sample_calls:
                self.failed_sample_calls[position].append(query)
            else:
                self.failed_sample_calls[position] = [query]

    #### Functions for removing based on provided files ####
    def remove_samples_from_file(self, file):
        samples = self.sample_file_to_dict(file)
        self.remove_failed_samples(samples)

    def remove_positions_from_file(self, file):
        positions = self.position_file_to_dict(file)
        self.remove_failed_positions(positions)

    # Accumulate samples from file for removal
    def sample_file_to_dict(self, file):
        samples = {}
        for line in open(file) :
            samples[line.rstrip('\n')] = ''
        return samples

    # Accumulate positions from file for removal
    def position_file_to_dict(self, file):
        positions = {}
        for line in open(file) :
            p = re.sub(r'\t', "/", line.rstrip('\n'))
            positions[p] = ''
        return positions

    #### Functions to remove each type of failure ####
    # Remove samples from variant set given a dictionary of samples.  sample_id is dictionary key
    def remove_failed_samples(self, failed):
        for s in failed:
            self.remove_sample(s)

    # Remove positions from variant set given a dictionary of positions.  position is dictionary key
    def remove_failed_positions(self, failed):
        for p in failed:
            reference_name, start, end = p.split("/")
            self.remove_variant(reference_name, start, end)

    # Remove sample calls from variant set given a dictionary of sample positions.  sample-call is dictionary key
    def remove_failed_sample_calls(self, failed):
        for p in failed:
            sample_id, reference_name, start, end = p.split("/")
            self.remove_sample_call(sample_id, reference_name, start, end)

    #### Google Genomics functions ####
    # Remove a single sample from a variant set
    def remove_sample(self, sample_id):
        call_set_id = self.gg.get_call_set_id(sample_id)
        if call_set_id is None:
            logging.warn("Failed to retrieve call set id for %s." % sample_id)
            return None
        self.gg.delete_call_set(call_set_id)

    # Remove an entire position from a variant set
    def remove_variant(self, reference_name, start, end):
        variant_id = self.gg.get_variant_id(reference_name=reference_name, start=start, end=end)
        if variant_id is None:
            logging.warn("Failed to retrieve variant id for %s %s %s" % (reference_name, start, end))
        self.gg.delete_variant(variant_id)

    # Remove a variant of a specific call from a variant set
    def remove_sample_call(self, sample_id, reference_name, start, end):
        # todo
        return None

    #### Miscellaneous functions ####
    # Check if a path exists
    def check_path(self, file):
        if not os.path.exists(file):
            raise Exception("%s not found!" % file)

    def setup_log(self, verbose):
        log_level = ''
        if verbose is True:
            log_level = logging.DEBUG
        else:
            log_level = logging.INFO
        # todo write to file
        logging.basicConfig(format='%(asctime)s %(message)s', level=log_level)
