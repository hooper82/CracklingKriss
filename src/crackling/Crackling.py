'''
https://github.com/bmds-lab/Crackling

Author: Jake Bradford, Dimitri Perrin

Config:
    - See config.ini
'''

import argparse
import ast
import csv
import joblib
import os
import re
import sys
import time
import tempfile
import psutil

from datetime import datetime

from ConfigManager import ConfigManager
from Paginator import Paginator
from Batchinator import Batchinator
from Constants import *
from Helpers import *
import multiprocessing as mp

#########################################
##   Removing targets with leading T   ##
#########################################


def removeT(key, lproxy):
    d = lproxy[key]
    if (key[-2:] == 'GG' and key[0] == 'T') or (key[:2] == 'CC' and key[-1] == 'A'):
        d['passedAvoidLeadingT'] = CODE_REJECTED
    else:
        d['passedAvoidLeadingT'] = CODE_ACCEPTED
    lproxy[key] = d

#########################################
##    AT% ideally is between 20-65%    ##
#########################################


def Ideal_AT(key, lproxy):
    d = lproxy[key]
    AT = AT_percentage(key[0:20])
    if AT < 20 or AT > 65:
        d['passedATPercent'] = CODE_REJECTED
    else:
        d['passedATPercent'] = CODE_ACCEPTED
    d['AT'] = AT
    lproxy[key] = d

############################################
##   Removing targets that contain TTTT   ##
############################################


def remove_TTTT(key, lproxy):
    d = lproxy[key]
    if 'TTTT' in key:
        d['passedTTTT'] = CODE_REJECTED
    else:
        d['passedTTTT'] = CODE_ACCEPTED
    lproxy[key] = d

    #########################################
    ##         sgRNAScorer 2.0 model       ##
    #########################################


def sgRNAScorer(key, lproxy):
    d = lproxy[key]

    # binary encoding
    encoding = {
        'A': '0001',    'C': '0010',    'T': '0100',    'G': '1000',
        'K': '1100',    'M': '0011',    'R': '1001',    'Y': '0110',
        'S': '1010',    'W': '0101',    'B': '1110',    'V': '1011',
        'H': '0111',    'D': '1101',    'N': '1111'
    }

    clfLinear = joblib.load(configMngr['sgrnascorer2']['model'])

    sequence = d.upper()
    entryList = []

    for x in range(0, 20):
        for y in range(0, 4):
            entryList.append(int(encoding[sequence[x]][y]))

        # predict based on the entry
        prediction = clfLinear.predict([entryList])
        score = clfLinear.decision_function([entryList])[0]

        d['sgrnascorer2score'] = score

        if float(score) < float(float(configMngr['sgrnascorer2']['score-threshold'])):
            d['acceptedBySgRnaScorer'] = CODE_REJECTED
        else:
            d['acceptedBySgRnaScorer'] = CODE_ACCEPTED
    lproxy[key] = d


def SecondsConversion(list):
    totalTime = 0
    totalTime += list[0] * 3600
    totalTime += list[1] * 60
    totalTime += list[2]
    totalString = f'{totalTime}.' + str(list[3])
    totalTime = float(totalString)
    print(f'total time in seconds convertor = {totalTime}')
    return totalTime


def Crackling(configMngr):
    totalSizeBytes = configMngr.getDatasetSizeBytes()
    completedSizeBytes = 0

    _stdout = sys.stdout
    _stderr = sys.stderr

    sys.stdout = configMngr.getLogMethod()
    sys.stderr = configMngr.getErrLogMethod()

    lastRunTimeSec = 0
    lastScaffoldSizeBytes = 0
    totalRunTimeSec = 0

    ####################################
    ###     Run-time Optimisation     ##
    ####################################
    def filterCandidateGuides(dictCandidateGuides, module):
        nonlocal configMngr

        module = module.lower()

        optimisation = configMngr['general']['optimisation']

        consensusN = int(configMngr['consensus']['n'])

        for target23 in dictCandidateGuides:
            doAssess = True

            if optimisation == 'ultralow':
                doAssess = True

            if optimisation == 'low':
                # Never assess guides that appear twice
                if (candidateGuides[target23]['seenDuplicate'] == CODE_REJECTED):
                    doAssess = False

            if optimisation == 'medium':
                # Never assess guides that appear twice
                if (candidateGuides[target23]['seenDuplicate'] == CODE_REJECTED):
                    doAssess = False

                # For mm10db:
                if (module == MODULE_MM10DB):
                    # if any of the mm10db tests have failed, then fail them all
                    if (CODE_REJECTED in [
                        candidateGuides[target23]['passedAvoidLeadingT'],
                        candidateGuides[target23]['passedATPercent'],
                        candidateGuides[target23]['passedTTTT'],
                        candidateGuides[target23]['passedSecondaryStructure'],
                        candidateGuides[target23]['acceptedByMm10db'],
                    ]):
                        doAssess = False

                # For CHOPCHOP:
                # Always assess, unless the guide is seen multiple times

                # For sgRNAScorer2:
                # Always assess, unless the guide is seen multiple times

                # For specificity:
                if (module == MODULE_SPECIFICITY):
                    # don't assess if they failed consensus
                    if (int(candidateGuides[target23]['consensusCount']) < consensusN):
                        doAssess = False

                    # don't assess if they failed Bowtie
                    if (candidateGuides[target23]['passedBowtie'] == CODE_REJECTED):
                        doAssess = False

            if optimisation == 'high':
                # Never assess guides that appear twice
                if (candidateGuides[target23]['seenDuplicate'] == CODE_REJECTED):
                    doAssess = False

                # For mm10db:
                if (module == MODULE_MM10DB):
                    # if any of the mm10db tests have failed, then fail them all
                    if (CODE_REJECTED in [
                        candidateGuides[target23]['passedAvoidLeadingT'],
                        candidateGuides[target23]['passedATPercent'],
                        candidateGuides[target23]['passedTTTT'],
                        candidateGuides[target23]['passedSecondaryStructure'],
                        candidateGuides[target23]['acceptedByMm10db'],
                    ]):
                        doAssess = False

                # For CHOPCHOP:
                if (module == MODULE_CHOPCHOP):
                    if consensusN == 1 and candidateGuides[target23]['acceptedByMm10db'] == CODE_ACCEPTED:
                        doAssess = False

                # For sgRNAScorer2:
                if (module == MODULE_SGRNASCORER2):
                    currentConsensus = ((int)(candidateGuides[target23]['acceptedByMm10db'] == CODE_ACCEPTED) +    # mm10db accepted
                                        (int)(candidateGuides[target23]['passedG20'] == CODE_ACCEPTED))                            # chopchop-g20 accepted

                    # if the guide is further than one test away from passing
                    # the consensus approach, then skip it.
                    if (currentConsensus < (consensusN - 1)):
                        doAssess = False

                # For specificity:
                if (module == MODULE_SPECIFICITY):
                    # don't assess if they failed consensus
                    if (int(candidateGuides[target23]['consensusCount']) < consensusN):
                        doAssess = False

                    # don't assess if they failed Bowtie
                    if (candidateGuides[target23]['passedBowtie'] == CODE_REJECTED):
                        doAssess = False

            if doAssess:
                yield target23

    def processSequence(sequence):
        # Patterns for guide matching
        pattern_forward = r'(?=([ATCG]{21}GG))'
        pattern_reverse = r'(?=(CC[ACGT]{21}))'

        # New sequence deteced, process sequence
        # once for forward, once for reverse
        for pattern, strand, seqModifier in [
            [pattern_forward, '+', lambda x: x],
            [pattern_reverse, '-', lambda x: rc(x)]
        ]:
            p = re.compile(pattern)
            for m in p.finditer(sequence):
                target23 = seqModifier(seq[m.start(): m.start() + 23])
                yield [target23, seqHeader,  m.start(),  m.start() + 23, strand]

    ###################################
    ##   Processing the input file   ##
    ###################################

    printer('Analysing files...')

    # Sets to keep track of Guides and sequences seen before
    candidateGuides = set()
    duplicateGuides = set()
    recordedSequences = set()

    for seqFilePath in configMngr.getIterFilesToProcess():
        # Run start time
        start_time = time.time()

        printer(f'Identifying possible target sites in: {seqFilePath}')

        completedPercent = round(
            (float(completedSizeBytes) / float(totalSizeBytes) * 100.0), 3)
        printer(
            f'{completedSizeBytes} of {totalSizeBytes} bytes processed ({completedPercent}%)')

        lastScaffoldSizeBytes = os.path.getsize(seqFilePath)

        completedSizeBytes += lastScaffoldSizeBytes

        # We first remove all the line breaks within a given sequence (FASTA format)
        with open(seqFilePath, 'r') as inFile, tempfile.NamedTemporaryFile(mode='w', delete=False) as parsedFile:
            for line in inFile:
                line = line.strip()
                if line[0] == '>':
                    # this is the header line for a new sequence, so we break the previous line and write the header as a new line
                    parsedFile.write('\n'+line+'\n')
                else:
                    # this is (part of) the sequence; we write it without line break
                    parsedFile.write(line.strip())

        guideBatchinator = Batchinator(int(configMngr['input']['batch-size']))

        with open(parsedFile.name, 'r') as inFile:
            seqHeader = ''
            seq = ''
            for line in inFile:
                # Remove garbage from line
                line = line.strip()
                # Some lines (e.g., first line in file) can be just a line break, move to next line
                if line == '':
                    continue
                # Header line, start of a new sequence
                elif line[0] == '>':
                    # If we haven't seen the sequence OR we have found a sequence without header
                    if (seqHeader not in recordedSequences) or (seqHeader == '' and seq != ''):
                        # Record header
                        recordedSequences.add(seqHeader)
                        # Process the sequence
                        for guide in processSequence(seq):
                            # Check if guide has been seen before
                            if guide[0] not in candidateGuides:
                                # Record guide
                                candidateGuides.add(guide[0])
                                # Record candidate guide to temp file
                                guideBatchinator.recordEntry(guide)
                            else:
                                # Record duplicate guide
                                duplicateGuides.add(guide[0])
                    # Update sequence and sequence header
                    seqHeader = line[1:]
                    seq = ''
                # Sequence line, section of existing sequence
                else:
                    # Append section to total sequence
                    seq += line.strip()

            # Process the last sequence
            for guide in processSequence(seq):
                # Check if guide has been seen before
                if guide[0] not in candidateGuides:
                    # Record guide
                    candidateGuides.add(guide[0])
                    # Record candidate guide to temp file
                    guideBatchinator.recordEntry(guide)
                else:
                    # Record duplicate guide
                    duplicateGuides.add(guide[0])

        printer(f'Identified {len(candidateGuides)} possible target sites.')

        printer(
            f'\t{len(duplicateGuides)} of {len(candidateGuides)} were seen more than once.')

        # Update total time
        preprocessingTime = time.time() - start_time
        totalRunTimeSec += preprocessingTime

    # Write header line for output file
    with open(configMngr['output']['file'], 'a+') as fOpen:
        csvWriter = csv.writer(fOpen, delimiter=configMngr['output']['delimiter'],
                               quotechar='"', dialect='unix', quoting=csv.QUOTE_MINIMAL)

        csvWriter.writerow(DEFAULT_GUIDE_PROPERTIES_ORDER)

    # Clean up unused variables
    os.unlink(parsedFile.name)
    del candidateGuides
    del recordedSequences

    for batchFile in guideBatchinator:
        # Run start time
        start_time = time.time()

        printer('Processing batch file...')

        # Create new candidate guide dictionary
        candidateGuides = {}
        # Load guides from temp file
        with open(batchFile, 'r') as inputFp:
            # Create csv reader to parse temp file
            csvReader = csv.reader(inputFp, delimiter=configMngr['output']['delimiter'],
                                   quotechar='"', dialect='unix', quoting=csv.QUOTE_MINIMAL)
            # Rebuild dictonary from temp file
            for row in csvReader:
                candidateGuides[row[0]] = DEFAULT_GUIDE_PROPERTIES.copy()
                candidateGuides[row[0]]['seq'] = row[0]
                if row[0] in duplicateGuides:
                    candidateGuides[row[0]]['header'] = CODE_AMBIGUOUS
                    candidateGuides[row[0]]['start'] = CODE_AMBIGUOUS
                    candidateGuides[row[0]]['end'] = CODE_AMBIGUOUS
                    candidateGuides[row[0]]['strand'] = CODE_AMBIGUOUS
                    candidateGuides[row[0]]['seenDuplicate'] = CODE_REJECTED
                else:
                    candidateGuides[row[0]]['header'] = row[1]
                    candidateGuides[row[0]]['start'] = row[2]
                    candidateGuides[row[0]]['end'] = row[3]
                    candidateGuides[row[0]]['strand'] = row[4]

        printer(
            f'Loaded batch {guideBatchinator.currentBatch} of {len(guideBatchinator.batchFiles)}')

        ###################################
        ##        Multiprocessing        ##
        ###################################

        manager = mp.Manager()
        mpd = manager.dict(candidateGuides)
        pool = mp.Pool(mp.cpu_count())

        # Start timing
        start = datetime.now()
        printer('Starting process sequence')

        # Call function using starmap and shared result dict
        # this is where the multiprocessing is implemented, however it tends to add to runtime
        pool.starmap(sgRNAScorer, [(key, mpd)
                     for key in mpd])
        pool.close()

        print({mpd[c]['passedATPercent']
              for c in mpd})

        # Stop timing
        end = datetime.now()
        print('Finished program')

        # convert times into seconds
        timeTaken = (end - start).total_seconds()
        print(f'Time taken = {timeTaken:.4f} seconds')

        #########################################
        ##           Begin output              ##
        #########################################
        printer('Writing results to file.')

        # Write guides to file. Include scores etc.
        with open(configMngr['output']['file'], 'a+') as fOpen:
            csvWriter = csv.writer(fOpen, delimiter=configMngr['output']['delimiter'],
                                   quotechar='"', dialect='unix', quoting=csv.QUOTE_MINIMAL)

            for target23 in mpd:
                output = [mpd[target23][x]
                          for x in DEFAULT_GUIDE_PROPERTIES_ORDER]

                csvWriter.writerow(output)

        #########################################
        ##              Clean up               ##
        #########################################
        printer('Cleaning auxiliary files')
        for f in [
            configMngr['rnafold']['input'],
            configMngr['rnafold']['output'],
            configMngr['offtargetscore']['input'],
            configMngr['offtargetscore']['output'],
            configMngr['bowtie2']['input'],
            configMngr['bowtie2']['output'],
        ]:
            try:
                os.remove(f)
            except:
                pass

        #########################################
        ##               Done                  ##
        #########################################
        printer('Done.')

        printer(f'{len(candidateGuides)} guides evaluated.')

        printer('Ran in {} (dd hh:mm:ss) or {} seconds'.format(
            time.strftime('%d %H:%M:%S', time.gmtime(
                (time.time() - start_time))),
            (time.time() - start_time)
        ))

        lastRunTimeSec = time.time() - start_time
        totalRunTimeSec += lastRunTimeSec

    printer('Total run time (dd hh:mm:ss) {} or {} seconds'.format(
        time.strftime('%d %H:%M:%S', time.gmtime(totalRunTimeSec)),
        totalRunTimeSec
    ))

    sys.stdout = _stdout
    sys.stderr = _stderr


if __name__ == '__main__':
    # load in config
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', help='Configuration file',
                        default=None, required=True)
    args = parser.parse_args()

    configMngr = ConfigManager(
        args.c, lambda x: print(f'configMngr says: {x}'))

    if not configMngr.isConfigured():
        print('Something went wrong with reading the configuration.')
        exit()
    else:
        printer('Crackling is starting...')

        Crackling(configMngr)

    print('Goodbye.')
