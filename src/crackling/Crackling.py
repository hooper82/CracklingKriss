'''
https://github.com/bmds-lab/Crackling

Author: Jake Bradford, Dimitri Perrin

Config:
    - See config.ini
'''

import ast, csv, joblib, os, re, sys, time, tempfile

from datetime import datetime
import multiprocessing as mp

from crackling.Paginator import Paginator
from crackling.Batchinator import Batchinator
from crackling.Constants import *
from crackling.Helpers import *
from crackling.FileProcessor import find_candidates_in_file


def sgRNAScorer(key, lproxy, sgrnascorer_model):
    d = lproxy[key]

    # binary encoding
    encoding = {
        'A': '0001',    'C': '0010',    'T': '0100',    'G': '1000',
        'K': '1100',    'M': '0011',    'R': '1001',    'Y': '0110',
        'S': '1010',    'W': '0101',    'B': '1110',    'V': '1011',
        'H': '0111',    'D': '1101',    'N': '1111'
    }

    clfLinear = joblib.load(sgrnascorer_model)

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



def Crackling(configMngr):
    totalSizeBytes = configMngr.getDatasetSizeBytes()
    completedSizeBytes = 0

    _stdout = sys.stdout
    _stderr = sys.stderr

    sys.stdout = configMngr.getLogMethod()
    sys.stderr = configMngr.getErrLogMethod()

    lastRunTimeSec = 0
    seqFileSize = 0
    totalRunTimeSec = 0

    startTime = time.time()

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
                if (candidateGuides[target23]['isUnique'] == CODE_REJECTED):
                    doAssess = False

            if optimisation == 'medium':
                # Never assess guides that appear twice
                if (candidateGuides[target23]['isUnique'] == CODE_REJECTED):
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
                if (candidateGuides[target23]['isUnique'] == CODE_REJECTED):
                    doAssess = False

                # For efficiency
                if module in [MODULE_CHOPCHOP, MODULE_MM10DB, MODULE_SGRNASCORER2]:

                    # `consensusCount` cannot be used yet as it may not have been
                    # calculated. instead, calculate it on-the-go.
                    countAlreadyAccepted = sum([
                        candidateGuides[target23]['acceptedByMm10db'] == CODE_ACCEPTED,
                        candidateGuides[target23]['passedG20'] == CODE_ACCEPTED,
                        candidateGuides[target23]['acceptedBySgRnaScorer'] == CODE_ACCEPTED,
                    ])

                    countAlreadyAssessed = sum([
                        candidateGuides[target23]['acceptedByMm10db'] in [CODE_ACCEPTED, CODE_REJECTED],
                        candidateGuides[target23]['passedG20'] in [CODE_ACCEPTED, CODE_REJECTED],
                        candidateGuides[target23]['acceptedBySgRnaScorer'] in [CODE_ACCEPTED, CODE_REJECTED],
                    ])

                    countToolsInConsensus = sum([
                        configMngr['consensus'].getboolean('mm10db'),
                        configMngr['consensus'].getboolean('chopchop'),
                        configMngr['consensus'].getboolean('sgRNAScorer2'),
                    ])

                    # Do not assess if passed consensus already
                    if countAlreadyAccepted >= consensusN:
                        doAssess = False

                    # Do not assess if there are not enough remaining tests to pass consensus
                    #   i.e. if the number of remaining tests is less than what is needed to pass then do not assess
                    if countToolsInConsensus - countAlreadyAssessed < consensusN - countAlreadyAccepted:
                        doAssess = False

                    # For mm10db:
                    if module == MODULE_MM10DB:
                        # if any of the mm10db tests have failed, then fail them all
                        if (CODE_REJECTED in [
                            candidateGuides[target23]['passedAvoidLeadingT'],
                            candidateGuides[target23]['passedATPercent'],
                            candidateGuides[target23]['passedTTTT'],
                            candidateGuides[target23]['passedSecondaryStructure'],
                            candidateGuides[target23]['acceptedByMm10db'],
                        ]):
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


    ###################################
    ##   Processing the input file   ##
    ###################################

    printer('Analysing files...')

    # Sets to keep track of Guides and sequences seen before
    candidateGuides = set()
    duplicateGuides = set()
    recordedSequences = set()

    guideBatchinator = Batchinator(int(configMngr['input']['batch-size']))

    printer(f'Batchinator is writing to: {guideBatchinator.workingDir.name}')

    for seqFilePath in configMngr.getIterFilesToProcess():

        start_time = time.time()

        candidateGuides, duplicateGuides, recordedSequences, fileSize, numIdentifiedGuides, numDuplicateGuides = find_candidates_in_file(guideBatchinator, seqFilePath, candidateGuides, duplicateGuides, recordedSequences)
        completedSizeBytes += fileSize

        duplicatePercent = round(numDuplicateGuides / numIdentifiedGuides * 100.0, 3)
        printer(f'\tIdentified {numIdentifiedGuides:,} possible target sites in this file.')
        printer(f'\tOf these, {len(duplicateGuides):,} are not unique. These sites occur a total of {numDuplicateGuides} times.')
        printer(f'\tRemoving {numDuplicateGuides:,} of {numIdentifiedGuides:,} ({duplicatePercent}%) guides.')
        printer(f'\t{len(candidateGuides):,} distinct guides have been discovered so far.')

        completedPercent = round(completedSizeBytes / totalSizeBytes * 100.0, 3)
        printer(f'\tExtracted from {completedPercent}% of input')

    # Write header line for output file
    with open(configMngr['output']['file'], 'a+') as fOpen:
        csvWriter = csv.writer(fOpen, delimiter=configMngr['output']['delimiter'],
                        quotechar='"',dialect='unix', quoting=csv.QUOTE_MINIMAL)

        csvWriter.writerow(DEFAULT_GUIDE_PROPERTIES_ORDER)

    # Clean up unused variables
    del candidateGuides
    del recordedSequences


    batchFileId = 0
    for batchFile in guideBatchinator:
        batchStartTime = time.time()

        printer(f'Processing batch file {(batchFileId+1):,} of {len(guideBatchinator)}')

        # Create new candidate guide dictionary
        candidateGuides = {}
        # Load guides from temp file
        with open(batchFile, 'r') as inputFp:
            # Create csv reader to parse temp file
            csvReader = csv.reader(inputFp, delimiter=configMngr['output']['delimiter'],
                quotechar='"',dialect='unix', quoting=csv.QUOTE_MINIMAL)
            # Rebuild dictonary from temp file
            for row in csvReader:
                candidateGuides[row[0]] = DEFAULT_GUIDE_PROPERTIES.copy()
                candidateGuides[row[0]]['seq'] = row[0]
                if row[0] in duplicateGuides:
                    candidateGuides[row[0]]['header'] = CODE_AMBIGUOUS
                    candidateGuides[row[0]]['start'] = CODE_AMBIGUOUS
                    candidateGuides[row[0]]['end'] = CODE_AMBIGUOUS
                    candidateGuides[row[0]]['strand'] = CODE_AMBIGUOUS
                    candidateGuides[row[0]]['isUnique'] = CODE_REJECTED
                else:
                    candidateGuides[row[0]]['header'] = row[1]
                    candidateGuides[row[0]]['start'] = row[2]
                    candidateGuides[row[0]]['end'] = row[3]
                    candidateGuides[row[0]]['strand'] = row[4]

        printer(f'\tLoaded {len(candidateGuides):,} guides')

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
        pool.starmap(sgRNAScorer, [(key, mpd, configMngr['sgrnascorer2']['model'])
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
    print('This file is not callable')
    exit(1)
