#!/bin/env python

'''
BWA pipeline for Kat Holt.

Authors: Bernie Pope, Kat Holt.

Description:

This program implements a workflow pipeline for next generation
sequencing, predominantly the read mapping task, up until variant
detection.

It uses the Ruffus library to make the description of the pipeline
more declarative.

It supports parallel evaluation of independent pipeline stages,
and can run stages on a cluster environment.

The pipeline is configured by an options file in YAML format,
including the actual commands which are run at each stage.
'''

from ruffus import *
import os.path
import shutil
from utils import (runStage, splitPath, getOptions, initLog, getCommand)
from chrom_info import (chromInfo)
import sys
import glob

# Read the configuation options from file, determine the reference file
# and list of sequence files.
options = getOptions()
reference = options['reference']
#sequences = options['sequences']
sequencePatterns = options['sequences']
sequences = []
if type(sequencePatterns) == list:
    for pattern in sequencePatterns:
        sequences.append(glob.glob(pattern))
else:
    sequences = glob.glob(sequencePatterns)
# Start the logging process.
logger = initLog(options)
# Get information about chromosomes in the reference file
chromosomes = chromInfo(reference)

# Index the reference file.
@files(reference, reference + '.bwt', logger)
def mkRefDataBase(reference, output, logger):
    runStage('mkRefDataBase', logger, options, reference, output)

# Align sequence reads to the reference genome.
@follows(mkRefDataBase)
@transform(sequences, suffix('.fastq'), '.sai', logger)
def alignSequence(sequence, output, logger):
    runStage('alignSequence', logger, options, reference, sequence, output)

input = r'(.+)_1\.sai'
extraInputs = [r'\1_2.sai', r'\1_1.fastq', r'\1_2.fastq']

# Convert alignments to SAM format.
@transform(alignSequence, regex(input), add_inputs(extraInputs), r'\1.sam', logger)
def alignToSam(inputs, output, logger):
   align1, [align2, seq1, seq2] = inputs
   runStage('alignToSam', logger, options, reference, align2, align2, seq1, seq2, output)

# Convert SAM alignments to BAM format.
@transform(alignToSam, suffix('.sam'), '.bam', logger)
def samToBam(samFile, output, logger):
    runStage('samToBam', logger, options, reference, samFile, output)

# Sort BAM alignments by (leftmost?) coordinates.
#@transform(samToBam, suffix('.bam'), '.sorted.bam', logger)
#def sortBam(bamFile, output, logger):
#    (prefix, name, ext) = splitPath(output)
#    outFile = os.path.join(prefix,name)
#    runStage('sortBam', logger, options, bamFile, outFile)

# Sort BAM alignments.
@transform(samToBam, suffix('.bam'), '.sorted.bam', logger)
def sortBam(bamFile, output, logger):
    (prefix, name, ext) = splitPath(output)
    outFile = os.path.join(prefix,name)
    runStage('sortBam', logger, options, bamFile, outFile)

# Convert BAM alignment to pileup format.
#@transform(sortBam, suffix('.bam'), '.full.pileup', logger)
@transform(sortBam, regex(r'(.+)\.sorted\.bam'), r'\1.full.pileup', logger)
def pileupFull(bamFile, output, logger):
    runStage('pileupFull', logger, options, reference, bamFile, output)

# Call SNPs
#@transform(sortBam, suffix('.bam'), '.var.pileup', logger)
@transform(sortBam, regex(r'(.+)\.sorted\.bam'), r'\1.var.pileup', logger)
def callSNPs(bamFile, output, logger):
    runStage('callSNPs', logger, options, reference, bamFile, output)

(prefix,seqName,ext) = splitPath(sequences[0])
fullPileups = glob.glob(os.path.join(prefix, '*.full.pileup'))

def filterInputs():
    for chromName,chromLength in chromosomes:
        for fullPileup in fullPileups:
            (prefix,seqName,ext) = splitPath(fullPileup)
            varPileup = os.path.join(prefix, seqName + '.var.pileup')
            output = os.path.join(prefix, seqName + '.' + chromName + 'var_filtered.pileup')
            print([fullPileup, output, varPileup, chromName, chromLength, logger])
            yield([fullPileup, output, varPileup, chromName, chromLength, logger])

@follows(pileupFull)
@follows(callSNPs)
@files(filterInputs)
def varFilter(fullPileup, output, varPileup, chromName, chromLength, logger):
    runStage('varFilter', logger, options, fullPileup, varPileup, chromName, chromLength, output)

# Define the end-points of the pipeline.
pipeline = [varFilter]

# Invoke the pipeline.
pipelineOptions = options['pipeline']
if pipelineOptions['style'] == 'run':
    # Perform the pipeline steps.
    pipeline_run(pipeline, multiprocess = pipelineOptions['procs'], logger = black_hole_logger)
elif pipelineOptions['style'] == 'flowchart':
    # Draw the pipeline as a diagram.
    pipeline_printout_graph ('flowchart.svg', 'svg', pipeline, no_key_legend = False)
