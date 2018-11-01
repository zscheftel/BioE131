#!/usr/bin/python

"""

Needleman-Wunsch Aligner
Bioengineering 131/231, Fall 2018

Command-line script to read in the contents of a multi-FASTA containing two
sequences and a score matrix, perform a global alignment, and print the
resulting alignment matrix and optimal alignment to STDOUT.

"""

import os
import sys
from Bio import SeqIO

class NWAligner:
    def __init__(self, score_matrix_fname):
        self.score_matrix, self.gap_penalty = self.load_score_matrix(score_matrix_fname)

    @staticmethod
    def load_score_matrix(fname):
        """
        Input: (String) A path to a scoring matrix file.
        Output: (Dictionary) A nested dictionary where keys are strings
                and elements are scores as integers.
                (Integer) The gap penalty score

        Example:

        >>> matrix, gap_penalty = NWAligner.load_score_matrix('/home/bioe131/BLOSUM62')
        >>> matrix['A']['A']
        4
        >>> matrix['W']['W']
        11
        >>> gap_penalty
        -4

        """

        score_matrix = {}
        gap_penalty = None
        columns = []

        with open(fname) as fp:
            for line_num, line in enumerate(fp):
                # ignore comments in matrix file
                line = line.rstrip('\n').split()
                if len(line) == 1:
                    gap_penalty = int(line[0])
                elif line[0] == '#':
                    continue
                elif line[0] == 'A':
                    columns = line
                    score_matrix = {column:{} for column in columns}
                else:
                    pairing = zip(columns, line)
                    score_matrix[columns[line_num-1]] = {column:row for column, row in pairing}
                    continue
        return score_matrix, gap_penalty

    @staticmethod
    def load_FASTA(fname):
        """
        Input: (String) A path to a FASTA file containing exactly two sequences.
        Output: (List) A list containing two strings: one for each sequence.

        Example:

        >>> seqs = NWAligner.load_FASTA('example.fa')
        >>> seqs[0]
        'YAADSKATPGNPAFHQDEIFLARIAFIYQMWDGGQLKLIDYAPHHVMCEE'
        >>> seqs[1]
        'WVGQPNMKVQHWSNMKACCVKFITWTFIAPEKHACKWTETAYQADCDIIW'
        >>> len(seqs)
        2

        """

        fasta_sequences = list(SeqIO.parse(open(fname),'fasta'))
        if len(fasta_sequences) > 2:
            raise Exception('Number of sequences should not exceed 2. The number of sequences was: {}'.format(len(fasta_sequences)))
        sequences = []
        for i in fasta_sequences:
            sequences.append(str(i.seq))
        return sequences


    def align(self, seq_x, seq_y, print_matrix = False):
        """
        Input: (Strings) Two sequences to be aligned (seq_x and seq_y).
               (Boolean) If print_matrix is True, print the dynamic programming
                         matrix before traceback.
        Output: (Tuple of strings) Two sequences, aligned.

        Example:

        >>> aligner = NWAligner('BLOSUM62')
        >>> seqs = aligner.load_FASTA('example.fa')
        >>> aligner.align(seqs[0], seqs[1])
        ('YAAD-SKATPGNPAF---HQDEIF--L-AR--IA-FIYQM-WDGGQLK-LIDYAPH-HVM-C---E-------E---',
         'W---VGQ--P-N--MKVQH----WSNMKA-CCV-KFI---TW------TFI--APEKH--ACKWTETAYQADCDIIW')

        """

        ###
        ### INITIALIZATION
        ###

        # create two empty matrices with sizes based on the input sequences.
        # one contains the dynamic programming matrix, the other contains
        # pointers we'll use during traceback
        matrix = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]
        pointers = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]

        ### TODO ###
        # Fill the top row of the matrix with scores for gaps.
        # Fill the first column of the matrix with scores for gaps.

        matrix[0][0] = 0
        for i in range(1, len(matrix[0])):
            matrix[0][i] = self.gap_penalty
            matrix[i][0] = self.gap_penalty

        ###
        ### RECURSION
        ###

        # fill the dynamic programming and pointer matrices
        for x in range(1, len(seq_x) + 1):
            for y in range(1, len(seq_y) + 1):
                match_score = self.score_matrix[seq_x[x - 1]][seq_y[y - 1]]

                ### TODO ###
                # Take the maximum score of three possibilities:
                #   1) The element in the matrix diagonal from this one
                #      plus the score of an exact match
                #   2) The element to the left plus a gap penalty
                #   3) The element above plus a gap penalty
                # ... and set the current element (matrix[x][y]) equal to that
                #
                # Keep track of which of these choices you made by setting
                #   the same element (i.e., pointers[x][y]) to some value that
                #   has meaning to you.

                upLeft = matrix[x - 1][y - 1] + int(match_score)
                left = matrix[x - 1][y] + self.gap_penalty
                up = matrix[x][y - 1] + self.gap_penalty
                #Finding the maximum
                largestElement = max(upLeft, left, up)
                matrix[x][y] = largestElement
                #For the pointers matrix, Up is 1, upLeft is 0, Left is -1
                if up == largestElement:
                    pointers[x][y] = 1
                elif upLeft == largestElement:
                    pointers[x][y] = 0
                else:
                    pointers[x][y] = -1

        # print the dynamic programming matrix
        if print_matrix:
            for x in range(len(seq_x) + 1):
                print (" ".join(map(lambda i: str(int(i)), matrix[x])))

        ###
        ### TRACEBACK
        ###

        # starting from the bottom right corner, follow the pointers back
        x, y = len(seq_x)-1, len(seq_y)-1

        # fill these lists with the aligned sequences
        align_x = []
        align_y = []

        while x > 0 or y > 0:
            move = pointers[x][y]

            ### TODO ###
            # Follow pointers back through the matrix to the origin.
            # Depending on which "move" you made at each element in the
            #   matrix, you'll either align seq_x to seq_y, seq_x to a gap, or
            #   seq_y to a gap.

            if move == 1:
                align_x = ['-'] + align_x
                align_y = [seq_y[y]] + align_y
                x = x - 1
            elif move == 0:
                align_x = [seq_x[x]] + align_x
                align_y = [seq_y[y]] + align_y
                x = x - 1
                y = y - 1
            else:
                align_x = [seq_x[x]] + align_x
                align_y = ["-"] + align_y
                y = y - 1

        # flip the alignments, as they're reversed
        return ("".join(align_x[::]), "".join(align_y[::]))

###                                      ###
### NO NEED TO EDIT CODE BELOW THIS LINE ###
###                                      ###

if __name__ == '__main__':
    def usage():
        print('usage: %s matrixfilename stringfilename')
        sys.exit(1)

    if len(sys.argv) != 3:
        usage()

    for fname in sys.argv[1:]:
        if not os.path.isfile(fname):
            print ('Can not open %s' % (fname,))
            usage()

    aligner = NWAligner(sys.argv[1])
    seqs = aligner.load_FASTA(sys.argv[2])
    result = aligner.align(seqs[0], seqs[1])

    print('>seq1\n%s\n>seq2\n%s' % (result[0], result[1]))
