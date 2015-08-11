import os
import sys
from GenomicsQueries import Queries
from BigQueryClient import BigQuery
from config import Config
from GoogleGenomicsClient import GoogleGenomicsClient
import logging
import time
import json
from pyflow import WorkflowRunner

class GenomicsQC(object):
    def __init__(self, verbose=False, client_secrets=None, project_number=None, dataset=None, variant_table=None,
                 expanded_table=None):

        # Set global variables
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

        qc_dataset = Config.QC_DATASET
        project_name = Config.PROJECT_NAME

        # Set up logging
        self.date = time.strftime("%Y%m%d-%H%M%S")
        self.setup_log(verbose)

        # Set up API clients
        self.bq = BigQuery(project_number=self.project_number,
                           client_secrets=self.client_secrets_path,
                           project_name=project_name,
                           qc_dataset=qc_dataset)

        self.gg = GoogleGenomicsClient(client_secrets=self.client_secrets_path,
                                       project_number=self.project_number,
                                       dataset=self.dataset)

        self.queries = QCSteps(verbose=verbose, client_secrets=self.client_secrets_path,
                               project_number=self.project_number, dataset=self.dataset,
                               variant_table=self.variant_table, expanded_table=self.expanded_table)

    #### Specific types of QC functions ####
    # Sample level QC
    def sample_qc(self, remove=False):
        logging.info("Running Sample Level QC")
        failed_samples = self.queries.sample_level_qc()
        if remove is True:
            self.remove_failed_samples(self.failed_samples)
        self.print_removed(failed_samples, 'samples')

    # Variant level QC
    def variant_qc(self, poll=False):
        logging.info("Running Variant Level QC")
        self.queries.variant_level_qc()
        if poll is True:
            logging.debug("Waiting for query completion")
            self.poll_jobs(job_ids)

    # Execute a custom list of queries
    def custom_list(self, qc_list):
        for s in qc_list:
            try:
                method = getattr(self.queries, s)
                result = method()
                print json.dumps(result, sort_keys=True, indent=2)
            except Exception as e:
                print "%s not a valid qc method!" % s

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
        logging.debug("Removing failed samples.")
        for s in failed:
            self.remove_sample(s)
        self.print_removed(failed, "samples")

    # Remove positions from variant set given a dictionary of positions.  position is dictionary key
    def remove_failed_positions(self, failed):
        logging.debug("Removing failed positions.")
        for p in failed:
            reference_name, start, end = p.split("/")
            self.remove_variant(reference_name, start, end)
        self.print_removed(failed, "positions")

    # Remove sample calls from variant set given a dictionary of sample positions.  sample-call is dictionary key
    def remove_failed_sample_calls(self, failed):
        logging.debug("Removing failed calls.")
        for p in failed:
            sample_id, reference_name, start, end = p.split("/")
            self.remove_sample_call(sample_id, reference_name, start, end)
        self.print_removed(failed, "calls")

    #### Google Genomics functions ####
    # Remove a single sample from a variant set
    def remove_sample(self, sample_id):
        logging.debug("Removing sample: %s" % sample_id)
        call_set_id = self.gg.get_call_set_id(sample_id)

        if call_set_id is None:
            logging.error("Failed to retrieve call set id for %s." % sample_id)
            return None
        logging.debug("callset id: %s" % call_set_id)
        self.gg.delete_call_set(call_set_id)

    # Remove an entire position from a variant set
    def remove_variant(self, reference_name, start, end):
        logging.debug("Removing position: %s %s %s" % (reference_name, start, end))
        variant_id = self.gg.get_variant_id(reference_name=reference_name, start=start, end=end)
        if variant_id is None:
            logging.error("Failed to retrieve variant id for %s %s %s" % (reference_name, start, end))
        self.gg.delete_variant(variant_id)

    # Remove a variant of a specific call from a variant set
    def remove_sample_call(self, sample_id, reference_name, start, end):
        logging.debug("Removing call: %s %s %s %s" % (sample_id, reference_name, start, end))
        # todo figure this out
        return None

    #### Printing ####
    def print_removed(self, failed, type):
        file = "failed_%s.%s.tsv" % (type, self.date)
        f = open(file,'a')
        for r in failed:
            f.write("%s\t%s\n" % (r, ",".join(failed[r])))
        f.close()

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
        logging.basicConfig(filename='genomics-qc.%s.log' % self.date,
                            format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', level=log_level)

    # Check for job status given a set of BigQuery ids
    def poll_jobs(self, ids):
        print "Waiting for all queries to complete"
        for query in ids:
            id = ids[query]
            logging.debug("Entering polling for %s" % query)
            result = self.bq.poll_job(id)
            logging.info("Query complete: %s %s" % (query, result))
            if result != 'DONE':
                print "Query failed: %s %s" % (query, result)
        logging.debug("Polling complete")

'''
QCSteps Class

Self contained functions for executing QC queries on BigQuery.  Individual queries can be executed as well as all
sample level or variant level qcs.
'''
class QCSteps(object):
    def __init__(self, verbose=False, client_secrets=None, project_number=None, dataset=None, variant_table=None,
                 expanded_table=None):

        # Set global variables
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

        qc_dataset = Config.QC_DATASET
        project_name = Config.PROJECT_NAME

        # Set up API clients
        self.bq = BigQuery(project_number=self.project_number,
                           client_secrets=self.client_secrets_path,
                           project_name=project_name,
                           qc_dataset=qc_dataset)

        self.gg = GoogleGenomicsClient(client_secrets=self.client_secrets_path,
                                       project_number=self.project_number,
                                       dataset=self.dataset)

        self.failed_samples = {}

    '''
    Execute all sample level qc steps
    A dictionary containing the ids of samples failing any qc steps and the qcs they failed will be returned
    '''
    def sample_level_qc(self):
        self.__collect_failed_samples(self.gender_check(), "gender_check")
        #self.__collect_failed_samples(self.genotyping_concordance(), "genotype_concordance")
        self.__collect_failed_samples(self.heterozygosity_rate(), "heterozygosity_rate")
        self.__collect_failed_samples(self.inbreeding_coefficient(), "inbreeding_coefficient")
        self.__collect_failed_samples(self.missingness_rate(), "missingness_rate")
        self.__collect_failed_samples(self.singletons(), "private_variants")
        return self.failed_samples

    '''
    Execute all variant level qc steps
    As some of these queries can return very large results they are output to a table by default.
    A list of job ids is returned that can be used to check the query status.
    '''
    def variant_level_qc(self, save_to_table=True):
        ids = {}
        ids["blacklisted"] = self.blacklisted(save_to_table)
        ids["hardy_weinberg"] = self.hardy_weinberg(save_to_table)
        ids["heterozygous_haplotype"] = self.heterozygous_haplotype(save_to_table)
        ids["titv_by_alts"] = self.titv_by_alternate_allele_counts(save_to_table)
        ids["titv_by_depth"] = self.titv_by_depth(save_to_table)
        ids["titv_by_genomic_window"] = self.titv_by_genomic_window(save_to_table)
        return ids

    '''
    SAMPLE LEVEL QC Queries
    '''

    '''
    Gender Check

    Gender is inferred for each genome by calculating the heterozygosity rate on the X chromosome. Genomes who's
    inferred sex is different from that of the reported sex are removed from the cohort. Although it is possible for
    people to be genotypically male and phenotypically female, it is more likely that samples and phenotypic records
    were mislabeled.
    '''
    def gender_check(self):
        print "Running Gender Check"
        query_file = Queries.GENDER_CHECK
        query = self.__prepare_query(query_file)
        return self.bq.run(query, self.__query_name(query_file), False)

    '''
    Genotyping Concordance

    We next want to look at the concordance between SNPs called from the sequencing data and those called through the
    use genotyping. This allows us to identify samples that may have been mixed up in the laboratory.  Any genomes
    with a concordance less than 0.95 are returned
    '''
    def genotyping_concordance(self):
        print "Running Genotype Concordance"
        query_file = Queries.GENOTYPING_CONCORDANCE
        query = self.__prepare_query(query_file)
        return self.bq.run(query, self.__query_name(query_file), False)

    '''
    Heterozygosity Rate

    Heterozygosity rate is defined as the number of heterozygous calls in a genome. Genomes with a heterozygosity rate
    more than 3 standard deviations away from the mean are returned.
    '''
    def heterozygosity_rate(self):
        print "Running Heterozygosity Rate"
        prequery_file = Queries.HETEROZYGOSITY_METRICS
        cutoffs = self.__three_sigma_cutoffs(prequery_file)
        query_file = Queries.HETEROZYGOSITY
        query = self.__prepare_query(query_file, cutoffs)
        return self.bq.run(query, self.__query_name(query_file), False)

    '''
    Inbreeding Coefficient

    The inbreeding coefficient (F) is a measure of expected homozygosity rates vs observed homozygosity rates for
    individual genomes. Here, we calculate the inbreeding coefficient using the method-of-moments estimator. Genomes
    with an inbreeding coefficient more than 3 standard deviations away from the mean are removed from the cohort.
    '''
    def inbreeding_coefficient(self):
        print "Running Inbreeding Coefficient"
        prequery_file = Queries.INBREEDING_COEFFICIENT_METRICS
        cutoffs = self.__three_sigma_cutoffs(prequery_file)
        query_file = Queries.INBREEDING_COEFFICIENT
        query = self.__prepare_query(query_file, cutoffs)
        return self.bq.run(query, self.__query_name(query_file), False)

    '''
    Missingness Rate

    Missingess is defined as the proportion of sites found in the reference genome that are not called in a given
    genome. We calculate the missingness rate of each genome in our cohort in order to identify samples that are
    potentially low quality. If a sample has a high missingness rate it may be indicative of issues with sample
    preparation or sequencing. Genomes with a missingness rate greater than 0.1 are returned.
    '''
    def missingness_rate(self):
        print "Running Sample Level Missingness Rate"
        query_file = Queries.MISSINGNESS_SAMPLE_LEVEL
        query = self.__prepare_query(query_file)
        return self.bq.run(query, self.__query_name(query_file), False)

    '''
    Singleton Rate

    Singleton rate is defined as the number of variants that are unique to a genome. If a variant is found in only one
    genome in the cohort it is considered a singleton. Genomes with singleton rates more than 3 standard deviations away
    from the mean are returned.
    '''
    def singletons(self):
        print "Running Private Variants"
        prequery_file = Queries.PRIVATE_VARIANT_METRICS
        cutoffs = self.__three_sigma_cutoffs(prequery_file)
        query_file = Queries.PRIVATE_VARIANTS
        query = self.__prepare_query(query_file, cutoffs)
        return self.bq.run(query, self.__query_name(query_file), False)


    '''
    VARIANT LEVEL QC

    By default all output is directed to a BigQuery table as the results may be very large.
    '''
    '''
    Blacklisted Variants

    Identify all variants within our cohort that have been blacklisted. For more information on what variants are
    blacklisted and why see here: https://sites.google.com/site/anshulkundaje/projects/blacklists

    By default all output is directed to a BigQuery table as the results may be very large.  Output is directed to
    qc_tables.blacklisted_variants
    '''
    def blacklisted(self, save_to_table=False):
        print "Running Blacklisted"
        query_file = Queries.BLACKLISTED
        query = self.__prepare_query(query_file)
        return self.bq.run(query, self.__query_name(query_file), save_to_table)

    '''
    Hardy-Weinberg Equilibrium

    For each variant, compute the expected versus observed relationship between allele frequencies and genotype
    frequencies per the Hardy-Weinberg Equilibrium.

    By default all output is directed to a BigQuery table as the results may be very large.  Output is directed to
    qc_tables.hwe_fail
    '''
    def hardy_weinberg(self, save_to_table=False):
        print "Running Hardy-Weinberg Equilibrium"
        prequery_file = Queries.HARDY_WEINBERG_METRICS
        cutoffs = self.__cutoff(prequery_file)
        query_file = Queries.HARDY_WEINBERG
        query = self.__prepare_query(query_file, cutoffs)
        return self.bq.run(query, self.__query_name(query_file), save_to_table)

    '''
    Heterozygous Haplotype

    For each variant within the X and Y chromosome, identify heterozygous variants in male genomes.  All positions
    with heterozygous positions outside the pseudo autosomal regions in male genomes are returned.

    By default all output is directed to a BigQuery table as the results may be very large.  Output is directed to
    qc_tables.sex_chromosome_heterozygous_haplotypes
    '''
    def heterozygous_haplotype(self, save_to_table=False):
        print "Running Heterozygous Haplotype"
        query_file = Queries.HETERZYGOUS_HAPLOTYPE
        query = self.__prepare_query(query_file)
        return self.bq.run(query, self.__query_name(query_file), save_to_table)

    '''
    Variant Level Missingness

    For each variant, compute the missingness rate. This query can be used to identify variants with a poor call rate.

    By default all output is directed to a BigQuery table as the results may be very large.  Output is directed to
    qc_tables.variant_level_missingness_fail
    '''
    def missingness_variant_level(self, save_to_table=False):
        print "Running Variant Level Missingness"
        query_file = Queries.MISSINGNESS_VARIANT_LEVEL
        query = self.__prepare_query(query_file)
        return self.bq.run(query, self.__query_name(query_file), save_to_table)

    '''
    Ti/Tv By Alternate Allele Counts

    Check whether the ratio of transitions vs. transversions in SNPs appears to be resonable across the range of rare
    variants to common variants. This query may help to identify problems with rare or common variants.

    By default all output is directed to a BigQuery table as the results may be very large.  Output is directed to
    qc_tables.titv_alternate_alleles
    '''
    def titv_by_alternate_allele_counts(self, save_to_table=False):
        print "Running Ti/Tv By Alternate Allele Counts"
        #todo
        return

    '''
    Ti/Tv By Depth

    Check whether the ratio of transitions vs. transversions in SNPs is within a set range at each sequencing depth for
    each sample.

    By default all output is directed to a BigQuery table as the results may be very large.  Output is directed to
    qc_tables.titv_by_depth_fail
    '''
    def titv_by_depth(self, save_to_table=False):
        print "Running Ti/Tv By Depth"
        query_file = Queries.TITV_DEPTH
        query = self.__prepare_query(query_file)
        return self.bq.run(query, self.__query_name(query_file), save_to_table)

    '''
    Ti/Tv By Genomic Window

    Check whether the ratio of transitions vs. transversions in SNPs appears to be reasonable in each window of genomic
    positions. This query may help identify problematic regions.

    By default all output is directed to a BigQuery table as the results may be very large.  Output is directed to
    qc_tables.titv_by_genomic_window_fail
    '''
    def titv_by_genomic_window(self, save_to_table=False):
        print "Running Ti/Tv By Genomic Window"
        query_file = Queries.TITV_GENOMIC_WINDOW
        return self.bq.run(query, self.__query_name(query_file), save_to_table)

    '''
    Private functions used to prepare queries for execution
    '''
    # Set up the query, read it in, apply substitutions
    def __prepare_query(self, query_file, cutoffs={}):
        logging.debug("Preparing query: %s" % query_file)
        raw_query = self.__get_query(query_file)
        all_subs = {}
        presets = self.__get_preset_cutoffs(query_file)
        all_subs.update(presets)
        all_subs.update(cutoffs)
        if query_file in Queries.MAIN_QUERY:
            main_query = Queries.MAIN_QUERY[query_file]
            main = self.__prepare_query(main_query)
            all_subs['_MAIN_QUERY_'] = main
        query = self.__query_substitutions(raw_query, other=all_subs)
        return query

    # Recursively substitute placeholders in each query.  Recursion is required because in some cases items that
    # substituted also contain placeholders.
    def __query_substitutions(self, query, other=None):
        replacements = {
            "_THE_TABLE_": Config.VARIANT_TABLE,
            "_THE_EXPANDED_TABLE_": Config.EXPANDED_TABLE,
            "_PATIENT_INFO_": Config.PATIENT_INFO,
            "_GENOTYPING_TABLE_": Config.GENOTYPING_TABLE
        }
        replacements.update(other)
        count = 0
        for r in replacements:
            count += query.count(r)
        if (count == 0):
            return query
        for r in replacements:
            query = query.replace(r, replacements[r])
        # Recursively make substitutions
        return self.__query_substitutions(query, other)

    # Check if a query requires a main query substitution
    def __main_query(self, query_file):
        if query_file in Queries.MAIN_QUERY:
            main_query = Queries.MAIN_QUERY[query_file]
            prepped_main = self.prepare_query(main_query)
            return prepped_main
        return None

    # Read raw query in from file
    def __get_query(self, file):
        path = os.path.join(self.query_repo, file)
        query = ''
        with open (path, "r") as f:
            query = f.read()
        return query

    # Get the base name of the query file and make some replacements
    def __query_name(self, query_file):
        query_name = query_file.split('.')[0]
        query_name = query_name.replace("-", "_")
        return query_name

    # Get preset cutoffs from query file
    def __get_preset_cutoffs(self, query):
        cutoffs = {}
        if query in Queries.PRESET_CUTOFFS:
            cutoffs = Queries.PRESET_CUTOFFS[query]
        return cutoffs

    # Determine cutoffs for methods that require cutoffs set at three standard deviations from the mean
    def __three_sigma_cutoffs(self, query_file):
        logging.debug("Getting average and standard deviation")
        query = self.__prepare_query(query_file)
        query_name = self.__query_name(query_file)
        result = self.bq.run(query, query_name=query_name)
        average, stddev = self.__get_average_stddev(result)
        logging.debug("Average: %s, Standard Deviation: %s" % (average, stddev))
        max, min = self.__calculate_max_min(average, stddev)
        logging.debug("Max: %s, Min: %s" % (max, min))
        substitutions = self.__create_max_min_substitutions(max, min)
        return substitutions

    # Get the average and standard deviation from a query result
    def __get_average_stddev(self, result):
        for r in result:
            average = r['average']
            stddev = r['stddev']
            return average, stddev

    # Calculate the 3 sigma rule given an average and standard deviation
    def __calculate_max_min(self, average, stddev):
        max = float(average) + (3 * float(stddev))
        min = float(average) - (3 * float(stddev))
        return max, min

    # Create a dictionary of substitutions given a max and min
    def __create_max_min_substitutions(self, max, min):
        dict = {
            "_MAX_VALUE_": "%s" % max,
            "_MIN_VALUE_": "%s" % min,
        }
        return dict

    # Get cutoff from prequery
    def __cutoff(self, query_file):
        logging.debug("Getting cutoff")
        query = self.__prepare_query(query_file)
        query_name = self.__query_name(query_file)
        result = self.bq.run(query, query_name=query_name)
        for r in result:
            cutoff = r['cutoff']
            substitution = {"_CUTOFF_": cutoff}
            return substitution
        return {}

    # Add failed samples to total
    def __collect_failed_samples(self, result, query):
        logging.debug("Collecting failed samples.")
        if result is None:
            return
        for r in result:
            sample_id = r['sample_id']
            if sample_id in self.failed_samples:
                self.failed_samples[sample_id].append(query)
            else:
                self.failed_samples[sample_id] = [query]
