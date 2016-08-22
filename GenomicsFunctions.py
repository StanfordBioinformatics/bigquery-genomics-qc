from __future__ import division
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


class GenomicsUtils(object):
    def __init__(self, verbose=False, client_secrets=None, project_number=None, dataset=None):
        # Set global variables
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
        #self.date = time.strftime("%Y%m%d-%H%M%S")
        #self.setup_log(verbose)
        # Set up API clients

        self.gg = GoogleGenomicsClient(client_secrets=self.client_secrets_path,
                                       project_number=self.project_number,
                                       dataset=self.dataset)

    def bed_coverage(self, bed_file, sample_id=None):
        # Get read groups
        read_groups = self.gg.get_read_groups(self.dataset, sample_id)
        # For each read group
        for r in read_groups:
            print r
            # read in bed file
            coverages = {}
            with open(bed_file) as f:
                for line in f:
                    line = line.rstrip()
                    chr, start, end, gene = line.split("\t")
                    if not gene in coverages:
                        coverages[gene] = {}
                    coverageBuckets = self.gg.coverage_buckets(r, chr, start, end)
                    if coverageBuckets is None:
                        print("No coverage calculated for %s %s %s" % (chr, start, end))
                    for c in coverageBuckets:
                        r_start = c["range"]["start"]
                        r_end = c["range"]["end"]
                        key = "-".join([chr, r_start, r_end])
                        coverages[gene][key] = {
                            "chr": chr,
                            "start": r_start,
                            "end": r_end,
                            "coverage": c["meanCoverage"]
                        }
            gene_coverage = self.percent_covered_gene(coverages)
            self.percent_covered_buckets(gene_coverage)



    def percent_covered_buckets(self, gene_coverage):
        percents = [100, 99.99, 99.95, 99.9, 99.5, 99, 95, 90, 75, 50]
        gene_count = len(gene_coverage)
        print "total genes: %s" % gene_count
        for p in percents:
            count = 0
            for g in gene_coverage:
                if gene_coverage[g] >= p/100:
                    count += 1
            print "%s,%s" % (p, count)


    def percent_covered_gene(self, coverages, cutoff=10):
        result = {}
        for gene in coverages:
            total_bases = 0
            total_covered = 0
            for k in coverages[gene]:
                range = coverages[gene][k]
                start = range["start"]
                end = range["end"]
                length = int(end) - int(start)
                coverage = range["coverage"]
                total_bases += length
                print ("%s\t%s\t%s\t%s" % (gene, start, total_bases, coverage))
                if coverage > cutoff:
                    total_covered += length
            percent_covered = total_covered / total_bases
            result[gene] = percent_covered
        print result
        return result

