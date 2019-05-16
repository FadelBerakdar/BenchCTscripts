#!/usr/bin/env python3
"""
    **
    * @file         compareTP.py
    * @brief        Python Script to Compare True/False Positives BenchCT Log Files
    * @copyright    © 2019 Novocraft Technologies Sdn Bhd. All rights reserved.
    * @author       Fadel Berakdar.
    * @license      This script is released under MIT License
    * @date         12/05/2019
    **
"""
    
import re
import argparse

spJ2readID = dict()

    
def getSpliceJunctions(_startPos, _CIGAR):
    """
    Finds splice junctions in a CIGAR string
    :param _startPos: start Postion in the reference
    :param _CIGAR:    CIGAR string
    :return: list of splice junctions [start, end, lengthN]
    """
    assert(_startPos >= 0)

    spJnList = []
    startPos = _startPos
    pattern = re.compile('\d+[A-Z]')
    pairs    = pattern.findall(_CIGAR)

    for pair in pairs:
        operation = pair[-1]
        length    = int(pair[0:-1])

        if operation == 'N':
            spJnList.append([startPos, startPos + length, pair])
            startPos += length

        elif operation in "MDX=":
            startPos += length

        elif operation in "IS":
            pass

        else:
            print("Error\t" + pair, '\t', length, '\t', operation)
            assert(True == False)

    return spJnList


def parseTruePosFile(_path2TP, _readIDColumn):
    """
    Parse true positives log file produced by benchCT
    :param _path2TP: path to true positives log file
    :param _readIDColumn: readID column is 0 for false  positives, 1 for true positives
    :return: spliceJunctions set contains hashes of splice junctions, hash = chromosome_start_end_lengthN
    """
    file = open(_path2TP, "r")
    spliceJunctions = set()

    for line in file.readlines():
        readID   = line.split('\t')[_readIDColumn]
        fields   = readID.split(':')
        pairInfo = fields[1].split(';')

        readOneInfo = pairInfo[0].split(',')
        readTwoInfo = pairInfo[1].split(',')

        chromsome1 = readOneInfo[0]
        chromsome2 = readTwoInfo[0]

        strand1 = readOneInfo[2]
        strand2 = readTwoInfo[2]

        startPos1 = int(readOneInfo[1]) + 1
        startPos2 = int(readTwoInfo[1]) + 1


        CIGAR1 = readOneInfo[3].split(':')[0]
        CIGAR2 = readTwoInfo[3]

        if CIGAR2[len(CIGAR2) - 2] == '/':
            CIGAR2 = CIGAR2[0:-2]

        for spJ in  getSpliceJunctions(startPos1, CIGAR1):
            # spliceJunctions.append([readID, spJ])
            spJHash = chromsome1  + "_" + str(spJ[2]) + "_" + str(spJ[0]) + "_" + str(spJ[1])
            spliceJunctions.add(spJHash)
            if spJHash in spJ2readID.keys():
                spJ2readID[spJHash].append(readID)
            else:
                spJ2readID[spJHash] = [readID]

            if readID == "77871310:X,130340251,+,100M;X,130340337,-,16M4828N84M:AAAAAAAAAAAAAAAoJ/1":
                print("\t\t" + spJHash)



        for spJ in  getSpliceJunctions(startPos2, CIGAR2):
            # spliceJunctions.append([readID, spJ])
            spJHash = chromsome2 + "_" + str(spJ[2]) + "_" + str(spJ[0]) + "_" + str(spJ[1])
            spliceJunctions.add(spJHash)
            
            if spJHash in spJ2readID.keys():
                spJ2readID[spJHash].append(readID)
            else:
                spJ2readID[spJHash] = [readID]

            if readID == "77871310:X,130340251,+,100M;X,130340337,-,16M4828N84M:AAAAAAAAAAAAAAAoJ/1":
                print("\t" + spJHash)

    file.close()

    return spliceJunctions


parser = argparse.ArgumentParser(description='Compare true positives BenchCT log files')

parser.add_argument('-1',
                    metavar='',
                    type=str,
                    # nargs=1,
                    dest="input1",
                    action="store",
                    required=True,
                    help="Path to first true positives log file")

parser.add_argument('-2',
                    metavar='',
                    type=str,
                    # nargs=1,
                    dest="input2",
                    action="store",
                    required=True,
                    help="Path to second true positives log file")

parser.add_argument('-c',
                    metavar='',
                    type=int,
                    # nargs=1,
                    dest="column",
                    action="store",
                    required=True,
                    help="0 if false positives, 1 if true positive. This is because true positives files contains extra first column")

if __name__ == '__main__':
    args = parser.parse_args()
    set1 = parseTruePosFile(args.input1, args.column)
    set2 = parseTruePosFile(args.input2, args.column)

    one_two = sorted(list(set1.intersection(set2)))
    one_only = sorted(list(set1.difference(set2)))
    two_only = sorted(list(set2.difference(set1)))

    for sp in one_two:
        print("1_2\t{}\t{}".format(sp, spJ2readID[sp]))

    for sp in one_only:
        print("1_only\t{}\t{}".format(sp, spJ2readID[sp]))

    for sp in two_only:
        print("2_only\t{}\t{}".format(sp, spJ2readID[sp]))
