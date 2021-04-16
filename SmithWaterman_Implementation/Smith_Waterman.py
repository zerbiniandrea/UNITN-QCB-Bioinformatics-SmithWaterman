#!/usr/bin/env python

# TODO
# - Traceback with multiple max values
# - Sequences as INPUT from USER

import numpy as np
import pandas as pd

# needed if I want an Object instead of a list of tuples in find_highest_scores
# class matrixElement:
#     def __init__(self,score,x,y):
#         self.score = score
#         self.x = x
#         self.y = y

def create_matrix(seq1, seq2):
    """
    Creates a matrix using a 2D numpy array
    having as dimensions the length of the two input sequences + 1
    """

    length1 = len(seq1)
    length2 = len(seq2)

    matrix = np.zeros((length2+1,length1+1))

    # SET THE TWO SEQUENCES AS NAMES OF ROWS AND COLUMNS IN A DATAFRAME
    #matrix = pd.DataFrame(matrix, columns=list(seq1)+1, index=list(seq2)+1)
    
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

    #https://stackoverflow.com/questions/12666494/how-do-i-decide-which-way-to-backtrack-in-the-smith-waterman-algorithm

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

# TODO
# atm works only for 1 max value
# need a loop for every max value + check multiple paths (take every option instead of elif)
def smith_waterman_traceback(matrix, seq1, seq2, match_score=3, mismatch_score=-3, gap_penalty=-2):
    list_of_coordinates = find_max(matrix)

    row = list_of_coordinates[0][0]
    column = list_of_coordinates[0][1]

    align_seq1 = []
    align_seq2 = []

    #for x,y in list_of_coordinates:

    # starting from the score in the position of the maximum score
    # calculate the score that the "parent" cells should contain 
    # and check which one is actually correct

    # could have also implemented the other way, 
    # starting from all the 3 "parent" cells (upper, left-most and diagonal)
    # and check which one gets the score of the current cell

    while not matrix[row][column] == 0:
        #basically the same scores in smith_waterman but in reverse
        if seq1[column-1]==seq2[row-1]: # -1 because the 1st row and the 1st column of the matrix are filled with zeros
            diagonal_score = matrix[row][column] - match_score
        else:
            diagonal_score = matrix[row][column] - mismatch_score

        horizontal_score = matrix[row][column] - gap_penalty
        vertical_score = matrix[row][column] - gap_penalty

        # check which cell contains the values obtained
        # and add the corresponding characters of the sequences to a list
        if matrix[row-1][column-1] == diagonal_score:
            row-=1
            column-=1
            align_seq1.append(seq1[column])
            align_seq2.append(seq2[row])
        elif matrix[row-1][column] == vertical_score:
            row-=1
            align_seq1.append("-")
            align_seq2.append(seq2[row])
        elif matrix[row][column-1] == horizontal_score:
            column-=1
            align_seq1.append(seq1[column])
            align_seq2.append("-")

    # revert the two sequences
    align_seq1 = align_seq1[::-1]
    align_seq2 = align_seq2[::-1]

    print(align_seq1)
    print(align_seq2)

def print_matrix(matrix, seq1, seq2):
    """
    Given a matrix as a 2D numpy array and two sequences,
    prints and returns a DataFrame 
    having the characters of the two input strings as the names of rows and columns.
    """
    matrix = pd.DataFrame(matrix, columns=list(" "+seq1), index=list(" "+seq2))
    matrix=matrix.astype(int)
    print(matrix)

    return matrix

def find_max(matrix):
    """
    Given a matrix as a 2D numpy array,
    finds the coordinates of the maximum value and returns them as a list, 
    in case the matrix contains multiple times the same maximum value.
    Also prints the maximum value and its coordinates.
    """
    max_score = int(np.amax(matrix))
    result = np.where(matrix == max_score)
    # zip the 2 arrays to get the exact coordinates
    list_of_coordinates = list(zip(result[0], result[1]))
    
    print("\nThe maximum score in the matrix is",max_score,"found at these coordinates:",list_of_coordinates)

    #for coord in list_of_coordinates:
    #    print(coord)

    return list_of_coordinates

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

if __name__ == '__main__':
    seq1 = "TGACTGTTGT"
    seq2 = "TAGGTGTGTG"

    matrix = smith_waterman(seq1,seq2)

    print_matrix(matrix,seq1,seq2)
    smith_waterman_traceback(matrix,seq1,seq2)
