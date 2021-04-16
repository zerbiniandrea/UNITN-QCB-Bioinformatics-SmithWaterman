#!/usr/bin/env python

import numpy as np
import pandas as pd

class matrixElement:
    def __init__(self,score,x,y):
        self.score = score
        self.x = x
        self.y = y

# CREATES A MATRIX WITH THE SHAPE OF THE TWO SEQUENCES
def create_matrix(seq1, seq2):
    length1 = len(seq1)
    length2 = len(seq2)

    matrix = np.zeros((length2+1,length1+1))

    # SET THE TWO SEQUENCES AS NAMES OF ROWS AND COLUMNS IN A DATAFRAME
    #matrix = pd.DataFrame(matrix, columns=list(seq1)+1, index=list(seq2)+1)
    
    return matrix

# EXECUTE SMITH WATERMAN ALGORITHM, RETURNING THE SCORING MATRIX 
def smith_waterman(seq1,seq2, match_score=3, mismatch_score=-3, gap_penalty=-2):
    # SCORES
    # match_score = 3
    # mismatch_score = -3
    # gap_penalty = -2

    # creating the matrix
    matrix = create_matrix(seq1,seq2)

    nRows, nColumns = matrix.shape

    #https://stackoverflow.com/questions/12666494/how-do-i-decide-which-way-to-backtrack-in-the-smith-waterman-algorithm

    for i in range(1,nRows):
        for j in range(1,nColumns):
            # DIAGONAL SCORE
            if seq1[j-1]==seq2[i-1]: # j is for the 1st string, i for the 2nd string; just because of how I want it to be displayed
                diagonalScore = matrix[i-1,j-1]+match_score
            else:
                diagonalScore = matrix[i-1,j-1]+mismatch_score

            # VERTICAL SCORE
            verticalScore = matrix[i,j-1]+gap_penalty
            #HORIZONTAL SCORE
            horizontalScore = matrix[i-1,j]+gap_penalty

            matrix[i,j] = max(0,diagonalScore,verticalScore,horizontalScore)

    return matrix

def print_matrix(matrix,seq1,seq2):
    matrix = pd.DataFrame(matrix, columns=list(" "+seq1), index=list(" "+seq2))
    matrix=matrix.astype(int)
    print(matrix)

    return matrix

def find_max(matrix):
    result = np.where(matrix == np.amax(matrix))
    # zip the 2 arrays to get the exact coordinates
    listOfCordinates = list(zip(result[0], result[1]))
    
    for cord in listOfCordinates:
        print(cord)

    return listOfCordinates

def find_best_scores(matrix, min_score=3):
    nRows, nColumns = matrix.shape

    best_scores = []

    for i in range(nRows):
        for j in range(nColumns):
            if matrix[i,j] > min_score:
                best_scores.append(matrixElement(int(matrix[i,j]),i,j))

    return best_scores

if __name__ == '__main__':
    seq1 = "TGTTACGG"
    seq2 = "GGTTGACTA"

    matrix = smith_waterman(seq1,seq2)

    print_matrix(matrix,seq1,seq2)
