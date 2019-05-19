#!/usr/bin/env python3
"""
    **
    * @file         convertFASTQ2SAM.py
    * @brief        Convert SimCT FASTQ to SAM
    * @copyright    © 2019 Novocraft Technologies Sdn Bhd. All rights reserved.
    * @author       Fadel Berakdar.
    * @license      This script is released under MIT License
    * @date         26/04/2019
    **
"""

import argparse
import gzip

def convert(_path2Reads1, _path2Reads2, _path2output):
    """
    Parse input SimCT FASTQ files ling by line, and write output SAM
    :param _path2Reads1: path to first mates FASTQ file
    :param _path2Reads2: path to second mates FASTQ file
    :param _path2output: path to output SAM file
    :return: non
    """
    inputFile1 = gzip.open(_path2Reads1, "r")
    inputFile2 = gzip.open(_path2Reads2, "r")
    outputFile = open(_path2output, "w")

    while True:
        read1 = [inputFile1.readline().decode('ascii'),
                 inputFile1.readline().decode('ascii'),
                 inputFile1.readline(),
                 inputFile1.readline().decode('ascii')]

        read2 = [inputFile2.readline().decode('ascii'),
                 inputFile2.readline().decode('ascii'),
                 inputFile2.readline(),
                 inputFile2.readline().decode('ascii')]

        if read1[0] == "":
            break

        fields = read1[0].split(':')
        pairInfo = fields[1].split(';')

        readOneInfo = pairInfo[0].split(',')
        readTwoInfo = pairInfo[1].split(',')

        chromsome1 = readOneInfo[0]
        chromsome2 = readTwoInfo[0]

        strand1 = readOneInfo[2]
        strand2 = readTwoInfo[2]

        startPos1 = int(readOneInfo[1]) + 1
        startPos2 = int(readTwoInfo[1]) + 1
        distance  = startPos2 + len(read2[1]) - startPos1

        CIGAR1 = readOneInfo[3].split(':')[0]
        CIGAR2 = readTwoInfo[3]

        flag1 = ""
        flag2 = ""

        if CIGAR2[len(CIGAR2) - 3] == '/':
            CIGAR2 = CIGAR2[0:-3]

        if strand1 == "+" and strand2 == "-":
            flag1 = "99"
            flag2 = "147"

        elif strand1 == "-" and strand2 == "+":
            flag1 = "83"
            flag2 = "163"

        else:
            assert(True == False)

        outputFile.write(read1[0][1:-1] +
                         "\t" +
                         flag1 +
                         "\t" +
                         chromsome1 +
                         "\t" +
                         str(startPos1) +
                         "\t60\t" +
                         CIGAR1 +
                         "\t=\t" +
                         str(startPos2) +
                         "\t" +
                         str(distance) +
                         "\t" +
                         read1[1][:-1] +
                         "\t" +
                         read1[3])

        outputFile.write(read2[0][1:-1] +
                         "\t" +
                         flag2 +
                         "\t" +
                         chromsome2 +
                         "\t" +
                         str(startPos2) +
                         "\t60\t" +
                         CIGAR2 +
                         "\t=\t" +
                         str(startPos1) +
                         "\t" +
                         str(-distance) +
                         "\t" +
                         read2[1][:-1] +
                         "\t" +
                         read2[3])

    inputFile2.close()
    inputFile1.close()
    outputFile.close()

parser = argparse.ArgumentParser(description='Convert SimCT FASTQ to SAM')

parser.add_argument('-1',
                    metavar='',
                    type=str,
                    # nargs=1,
                    dest="input1",
                    action="store",
                    required=True,
                    help="Path to first mates FASTQ file")

parser.add_argument('-2',
                    metavar='',
                    type=str,
                    # nargs=1,
                    dest="input2",
                    action="store",
                    required=True,
                    help="Path to second mates FASTQ file")

parser.add_argument('-o',
                    metavar='',
                    type=str,
                    # nargs=1,
                    dest="output",
                    action="store",
                    required=True,
                    help="Path to output file")

if __name__ == '__main__':
    args = parser.parse_args()
    convert(args.input1, args.input2, args.output)