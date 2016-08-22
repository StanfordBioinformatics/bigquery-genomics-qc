#!/usr/bin/python
from __future__ import division
import argparse

# 1: chromosome
# 2: start
# 3: end
# 4: gene
# 5: average coverage
# 6: 100x
# 7: 50x
# 8: 30x
# 9: 20x
# 10: 10x
# 11: 5x
# 12: 1x



def per_base_coverage(input):
    # Read input
    gene_sizes, gene_coverage, gene_metrics = process_file(input)
    coverage_percentages = percent_coverage(gene_sizes, gene_coverage)

    #gene_info(gene_metrics, coverage_percentages)
    ten_x_buckets(coverage_percentages)

def process_file(input):
    gene_sizes = {}
    gene_coverage = {}
    gene_metrics = {}
    with open(input) as f:
        for line in f:
            line = line.rstrip()
            values = line.split(",")
            gene = values[3]
            segment_length = int(values[2]) - int(values[1]) + 1
            if gene not in gene_sizes:
                gene_sizes[gene] = segment_length
                gene_coverage[gene] = [0,0,0,0,0,0,0]
                gene_metrics[gene] = {"min": 100000, # arbitrarily large minimum coverage
                                "max": 0,
                                "average": 0,
                                "length": 0,
                                "zero": 0}
            else:
                gene_sizes[gene] += segment_length
            for i in range(0, 7):
                gene_coverage[gene][i] += int(values[i+7])

            if gene_metrics[gene]["min"] > int(values[5]):
                gene_metrics[gene]["min"] = int(values[5])
            if gene_metrics[gene]["max"] < int(values[6]):
                gene_metrics[gene]["max"] = int(values[6])


            gene_metrics[gene]["average"] = ((gene_metrics[gene]["average"] *  gene_metrics[gene]["length"]) +
                                            (int(values[4]) * segment_length)) / (gene_metrics[gene]["length"] + segment_length)

            gene_metrics[gene]["length"] += segment_length
            gene_metrics[gene]["zero"] += segment_length - int(values[13])





    return gene_sizes, gene_coverage, gene_metrics

def percent_coverage(gene_sizes, gene_coverage):
    for g in gene_coverage:
        size = gene_sizes[g]
        for i in range(0, len(gene_coverage[g])):
            gene_coverage[g][i] = gene_coverage[g][i]/size
    return gene_coverage

def ten_x_buckets(coverage_percentages):
        # print out proportion of genes that have X% bases covered at 10x. X = 100, 99.9, 99.5, 99, 95, 90, 50, 10, 0
    cutoffs = [1, .999, .995, .99, .95, .90, .50, .10, 0]
    for c in cutoffs:
        count = 0
        for g in coverage_percentages:
            if coverage_percentages[g][4] >= c:
                count += 1
        #print "%s: %s" % (c, count)
        print "%s," % count,

def gene_info(gene_metrics, coverage_percentages):
    for g in gene_metrics:
        print "%s,%s,%s,%s,%s,%s,%s, %s, %s" %(g, gene_metrics[g]["min"], gene_metrics[g]["max"],
                                    gene_metrics[g]["average"], gene_metrics[g]["length"],
                                    coverage_percentages[g][4], coverage_percentages[g][0],
                                    coverage_percentages[g][6], gene_metrics[g]["zero"])

def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'This script calculates per base coverage based on output from an in-progress dataflow job')

    parser.add_argument("--input", default=None,
                                help="Input file. --input")

    options = parser.parse_args()
    if options.input is None:
        print "Exiting, specify input. --input."
        exit(0)

    return options

if __name__ == "__main__":
    options = parse_command_line()
    per_base_coverage(input=options.input)