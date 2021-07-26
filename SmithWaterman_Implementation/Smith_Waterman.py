#!/usr/bin/env python

"""
Zerbini Andrea - MAT. 224099
UniTN - Quantitative and Computational Biology (QCB)
Bioinformatics - Algorithms for Bioinformatics
"""

import numpy as np
import pandas as pd
import sys

# if I need to get only a single character
# EXAMPLE:
# in interactive mode, press F to choose a file and press S to insert a string

# getch = _Getch()
# class _GetchWindows:
#     def __init__(self):
#         import msvcrt

#     def __call__(self):
#         import msvcrt
#         return msvcrt.getch()

def create_matrix(seq1, seq2):
    """
    Creates a matrix using a 2D numpy array
    having as dimensions the length of the two input sequences + 1
    """

    length1 = len(seq1)
    length2 = len(seq2)

    # length+1 considering the first column and the first row
    # of the matrix that should contain only zeros
    matrix = np.zeros((length2+1,length1+1))

    return matrix

def smith_waterman(seq1, seq2, match_score=3, mismatch_score=-3, gap_penalty=-2):
    """
    Implementation of the Smith Waterman algorithm;
    given two sequences 
    and, optionally, the values match_score, mismatch_score and gap_penalty
    returns the scoring matrix as a 2D numpy array
    """
    # SCORES
    #     match_score = 3
    #     mismatch_score = -3
    #     gap_penalty = -2
    
    # create the matrix
    matrix = create_matrix(seq1,seq2)

    n_rows, n_columns = matrix.shape

    for i in range(1,n_rows):
        for j in range(1,n_columns):
            # DIAGONAL SCORE
            if seq1[j-1]==seq2[i-1]: # j is for the 1st string, i for the 2nd string; just because of how I want it to be displayed
                diagonal_score = matrix[i-1,j-1]+match_score
            else:
                diagonal_score = matrix[i-1,j-1]+mismatch_score

            # VERTICAL SCORE
            vertical_score = matrix[i,j-1]+gap_penalty
            #HORIZONTAL SCORE
            horizontal_score = matrix[i-1,j]+gap_penalty

            matrix[i,j] = max(0,diagonal_score,vertical_score,horizontal_score)

    return matrix

def smith_waterman_traceback(matrix, seq1, seq2, match_score=3, mismatch_score=-3, gap_penalty=-2):
    """
    Traceback wrapper.
    For every couple of coordinates containing the maximum value, 
    calls the recursive function to find every alignment.
    Returns all the alignments found as a list of tuples (align_seq1, align_seq2)
    """
    # not using max_score at the moment
    max_score, list_of_coordinates = find_max(matrix)

    res = []

    # for every pair of coordinates containing the maximum value
    for row,column in list_of_coordinates:
        tmp_res = smith_waterman_traceback_rec(matrix, row, column, [], [], seq1, seq2, match_score, mismatch_score, gap_penalty)

        # asserting that the traceback_rec function returned a list of tuples
        assert len(tmp_res)%2 == 0, "Something went wrong with the size of the result in the traceback function!"

        # this is needed to return a list of tuples, each containing the two alignment sequences
        # because of the structure of the traceback_rec function
        # the value returned is a list of multiple alignments, 
        # containing 2 values for each "path"
        for i in range(0,len(tmp_res),2):
            res.append((tmp_res[i],tmp_res[i+1]))

    return res

def smith_waterman_traceback_rec(matrix, row, column, align_seq1, align_seq2, seq1, seq2, match_score, mismatch_score, gap_penalty):
    """
    Traceback implementation of the Smith Waterman algorithm.
    Given a matrix and the two sequences that created the matrix,
    return a list of tuples containing all the local alignments with the maximum score
    """
    # starting from the score in the position of the maximum score
    # calculate the score that the "parent" cells should contain 
    # and check which one is actually correct

    # could have also implemented the other way, 
    # starting from all the 3 "parent" cells (upper, left-most and diagonal)
    # and check which one gets the score of the current cell

    res = []

    # stopping when finding a zero, meaning that the algorithm reached the start of the alignment
    if matrix[row][column] == 0:
        # return the two sequences as a tuple
        # reverting the sequences because of them being built starting from the last character of the string
        return (align_seq1[::-1], align_seq2[::-1])
    else:
        #basically the same scores in smith_waterman but in reverse
        if seq1[column-1]==seq2[row-1]: # -1 because the 1st row and the 1st column of the matrix are filled with zeros
            diagonal_score = matrix[row][column] - match_score
        else:
            diagonal_score = matrix[row][column] - mismatch_score

        horizontal_score = matrix[row][column] - gap_penalty
        vertical_score = matrix[row][column] - gap_penalty

        # check if the values obtained are the values of the "parents"
        if matrix[row-1][column-1] == diagonal_score:
            # need to create copies of the original sequences to append the correct characters
            # row-=1
            # column-=1
            # align_seq1.append(seq1[column])
            # align_seq2.append(seq2[row])
            align_tmp1 = align_seq1.copy()
            align_tmp2 = align_seq2.copy()
            align_tmp1.append(seq1[column-1])
            align_tmp2.append(seq2[row-1])

            # recursion considering the diagonal element
            res+=smith_waterman_traceback_rec(matrix,row-1,column-1,align_tmp1,align_tmp2,seq1,seq2,match_score,mismatch_score,gap_penalty)

        if matrix[row-1][column] == vertical_score:
            # need to create copies of the original sequences to append the correct characters
            # row-=1
            # align_seq1.append("-")
            # align_seq2.append(seq2[row])
            align_tmp1 = align_seq1.copy()
            align_tmp2 = align_seq2.copy()
            align_tmp1.append("-")
            align_tmp2.append(seq2[row-1])

            # recursion considering the vertical element
            res+=smith_waterman_traceback_rec(matrix,row-1,column,align_tmp1,align_tmp2,seq1,seq2,match_score,mismatch_score,gap_penalty)

        if matrix[row][column-1] == horizontal_score:
            # need to create copies of the original sequences to append the correct characters
            # column-=1
            # align_seq1.append(seq1[column])
            # align_seq2.append("-")
            align_tmp1 = align_seq1.copy()
            align_tmp2 = align_seq2.copy()
            align_tmp1.append(seq1[column-1])
            align_tmp2.append("-")

            # recursion considering the horizontal element
            res+=smith_waterman_traceback_rec(matrix,row,column-1,align_tmp1,align_tmp2,seq1,seq2,match_score,mismatch_score,gap_penalty)

        return res

def print_matrix(matrix, seq1, seq2):
    """
    Given a matrix as a 2D numpy array and two sequences,
    prints and returns a DataFrame 
    having the characters of the two input strings as the names of rows and columns.
    """
    # adding a space before every sequence because
    # the first row and the first column are filled with zeros
    # and aren't characters of the sequences
    matrix = pd.DataFrame(matrix, columns=list(" "+seq1), index=list(" "+seq2))
    matrix=matrix.astype(int)
    print("Dynamic programming matrix visualization:\n")
    print(matrix)

    return matrix

def find_max(matrix):
    """
    Given a matrix as a 2D numpy array,
    finds the coordinates of the maximum value and returns them as a list, 
    in case the matrix contains multiple times the same maximum value.
    Also prints the maximum value and its coordinates.
    """
    # find the maximum value contained in the matrix using numpy
    max_score = int(np.amax(matrix))
    # search for the coordinates of the cells containig the max value
    result = np.where(matrix == max_score)
    # zip the 2 arrays to get the exact coordinates
    list_of_coordinates = list(zip(result[0], result[1]))

    # exit the program if no alignment was found (max=0)
    if max_score == 0:
        print("\nNo alignment was found!")
        sys.exit(0)

    print("\nThe maximum local alignment score in the matrix is",max_score,"found at these coordinates:",list_of_coordinates,"\n")

    return max_score, list_of_coordinates

# to be used instead of find_max if I want to find
# not just the best alignment but also the others
# having a score higher than a treshold
def find_highest_scores(matrix, min_score=3):
    """
    Given a matrix and, optionally, a minimum score,
    returns a list of tuples having a score higher than the minimum score,
    each containing the value of the item and its coordinates in the matrix (X,Y)
    """
    n_rows, n_columns = matrix.shape

    highest_scores = []

    for i in range(n_rows):
        for j in range(n_columns):
            if matrix[i,j] > min_score:
                highest_scores.append((int(matrix[i,j]),i,j))

    return highest_scores

def print_res(align_res):
    """
    Prints the alignments found in the matrix
    """
    print("Alignments found:")
    # joining the characters of the two lists contained in every tuple
    for align in align_res:
        print(" ".join(align[0]))
        print(" ".join(align[1]),"\n")

def read_fasta_seq(file):
    """
    Read a file, assuming it's in a fasta format
    and returns the two sequences
    """
    seq_dict = {}
    header = ""

    with open(file,"r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line
                seq_dict[header] = ""
            else:
                line = line.upper()
                seq_dict[header] += line

    # if the file doesn't contain EXACTLY 2 sequences
    if not len(seq_dict.values()) == 2:
        print("INVALID FILE: the fasta file doesn't contain EXACTLY 2 sequences")
        sys.exit(-1)

    return seq_dict

def help_function():
    """
    Print the help function, explaining the input arguments
    """
    print("|=================================================================|")
    print("| This is a Python implementation of the Smith Waterman algorithm.|")
    print("|=================================================================|")

    print("\n- To use the program in an interactive mode, please use -i or --interactive as the only parameter.")
    print("\tExample: \n\t\tpython Smith_Waterman.py -i")

    print("\n- If you want to input two sequences from a file, use -f or --file as the first parameter, \n"+
        "  followed by the name of the file containing the sequences. \n"+
        "  The file is expected to be formatted as a fasta file containing two sequences.")
    print("\tExample: \n\t\tpython Smith_Waterman.py -f file.fasta")

    print("\n\nDEFAULT SCORES:\n"+
        "\t- match_score = 3\n"+
        "\t- mismatch_score = -3\n"+
        "\t- gap_penalty = -2\n")
    sys.exit(1)

def input_check():
    """
    Checks the input parameters and returns the two input sequences
    """
    # no arguments
    if len(sys.argv) < 2:
        help_function()
    else:
        if sys.argv[1] == "-h" or sys.argv[1] == "--help":
            help_function()
        elif sys.argv[1] == "-i" or sys.argv[1] == "--interactive":
            print("| ================================= |")
            print("|    Entering interactive mode...   |")
            print("| ================================= |")
            print("|")
            seq1 = input("| Please insert the first sequence: \n| ")
            seq2 = input("| Please insert the second sequence: \n| ")
            print("|")
            print("| ================================= |\n")
        elif sys.argv[1] == "-f" or sys.argv[1] == "--file":
            if len(sys.argv) < 3:
                print("Not enough parameters!")
                sys.exit(-1)
            else:
                seq_dict = read_fasta_seq(sys.argv[2])
                
                seq1 = list(seq_dict.values())[0]
                seq2 = list(seq_dict.values())[1]

    # exit the program if one or more strings are empty
    if not seq1 or not seq2:
        print("\n\nERROR:")
        print("\tOne or more input strings are empty!")
        sys.exit(-1)

    # converting the two sequences to upper case
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    return seq1, seq2

if __name__ == '__main__':
    # SCORES
    #   match_score = 3
    #   mismatch_score = -3
    #   gap_penalty = -2
    # EXAMPLE SEQUENCES
    #   seq1 = "TGTTATG"
    #   seq2 = "TAGGTGTGTG"

    # get the two input sequences
    # either from a file or from the interactive mode
    seq1, seq2 = input_check()

    # computing the scoring matrix and saving it
    matrix = smith_waterman(seq1,seq2,5,-3,-2)

    # just to print the matrix with the characters of the sequences along the rows and columns
    print_matrix(matrix,seq1,seq2)

    # traceback function, to find the alignments
    align_res = smith_waterman_traceback(matrix,seq1,seq2,5,-3,-2)

    # print the alignments found
    print_res(align_res)