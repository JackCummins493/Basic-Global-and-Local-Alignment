import pandas as pd
import numpy as np

blosum62 = pd.read_csv("blosum62.csv")
blosum62.set_index("Unnamed: 0", inplace=True)
advanced_nucleotide = pd.read_csv("pyrimidene_purine_nucleotide_scoring.csv")
advanced_nucleotide.set_index("Unnamed: 0", inplace=True)
basic_amino_acid = pd.read_csv("basic_amino_acid_scoring.csv")
basic_amino_acid.set_index("Unnamed: 0", inplace=True)
basic_nucleotide = pd.read_csv("basic_nucleotide_scoring.csv")
basic_nucleotide.set_index("Unnamed: 0", inplace=True)
pam250 = pd.read_csv("PAM250.csv")
pam250.set_index("Unnamed: 0", inplace=True)
pam250.index = pam250.index.str.strip()
pam250.columns.str.strip()
# gap penalty should be negative, and if using match_score and mismatch_penalty, match_score should be positive, mismatch_penalty should be negative
# table choices: 'simple_nucleotide', 'purine_pyrimidine_nucleotide', 'simple_aminoacid', 'BLOSUM_aminoacid', 'PAM_aminoacid', or don't specify and enter match score and mismatch penalty
def align(sequence1, sequence2, gap_penalty, table=None, match_score = None, mismatch_penalty = None):
    # set scoring based off of value entered for table
    if table == "simple_nucleotide":
        global basic_nucleotide
        scoring = basic_nucleotide
    elif table == "purine_pyrimidine_nucleotide":
        global advanced_nucleotide
        scoring = advanced_nucleotide
    elif table == "simple_aminoacid":
        global basic_amino_acid
        scoring = basic_amino_acid
    elif table == "BLOSUM_aminoacid":
        global blosum62
        scoring = blosum62
    elif table == "PAM_aminoacid":
        global pam250
        scoring = pam250
    elif table == None:
        if mismatch_penalty == None or match_score == None:
            return "If you are not using a score table, please specify the match score and the mismatch penalty."
    else:
        return "Please enter a valid table."
    
    sequence1_list = list(sequence1.upper())
    sequence2_list = list(sequence2.upper())
    
    score_matrix = np.zeros((len(sequence1_list) + 1, len(sequence2_list) + 1)) # define score matrix with correct dimensions
    score = 0
    for i, value in np.ndenumerate(score_matrix[0]): # set default value for first row of score matrix
        score_matrix[0, i] = score
        score += gap_penalty
    score = 0
    for i, value in np.ndenumerate(score_matrix.T[0]): # set default value for first column of score matrix
        score_matrix[i, 0] = score
        score += gap_penalty
    

    traceback_matrix = np.zeros((len(sequence1_list) + 1, len(sequence2_list) + 1), 'U4') # define traceback matrix
    # default values for first row and column of traceback matrix
    traceback_matrix[0, :] = "left"
    traceback_matrix[:, 0] = "up"
    traceback_matrix[0,0] = "done"
    highest_row = 0
    highest_column = 0
    highest_score = 0
    if 'scoring' not in locals(): # if the code is using match and mismatch scores
        for i, value in np.ndenumerate(score_matrix): # iterate through every value and record index
            if i[0] == 0 or i[1] == 0:
                continue
            if sequence1_list[i[0]-1] == sequence2_list[i[1]-1]: # if there is a match
                diag_score = score_matrix[i[0]-1, i[1]-1] + match_score
            else: # if there is not a match
                diag_score = score_matrix[i[0]-1, i[1]-1] + mismatch_penalty 
            left_score = score_matrix[i[0], i[1]-1] + gap_penalty
            up_score = score_matrix[i[0]-1, i[1]] + gap_penalty
            if left_score >= diag_score and left_score >= up_score and left_score >= 0: # if the left gap score is the greatest
                score_matrix[i[0], i[1]] = left_score
                traceback_matrix[i[0], i[1]] = 'left'
            elif up_score >= diag_score and up_score >= 0: # if the up gap score is the greatest
                score_matrix[i[0], i[1]] = up_score
                traceback_matrix[i[0], i[1]] = 'up'
            elif diag_score >= 0: # if the no gap score is the greatest
                score_matrix[i[0], i[1]] = diag_score
                traceback_matrix[i[0], i[1]] = 'diag'
            else:
                score_matrix[i[0], i[1]] = 0
                traceback_matrix[i[0], i[1]] = 'done'
            if score_matrix[i[0], i[1]] > highest_score:
                highest_row = i[0]
                highest_column = i[1]
                highest_score = score_matrix[i[0], i[1]]
    else: # if there is a scoring table
        for i, rows in np.ndenumerate(score_matrix):
            if i[0] == 0 or i[1] == 0:
                continue
            try:
                match_points = scoring.loc[sequence1_list[i[0]-1], sequence2_list[i[1]-1]] # find the score of the match with the two letters using the table
            except: # if one of the letters in the sequence is not in the table 
                return "Make sure that your sequences match the input table"
            diag_score = score_matrix[i[0]-1, i[1]-1] + match_points
            left_score = score_matrix[i[0], i[1]-1] + gap_penalty
            up_score = score_matrix[i[0]-1, i[1]] + gap_penalty
            # same as before
            if left_score >= diag_score and left_score >= up_score and left_score >= 0:
                score_matrix[i[0], i[1]] = left_score
                traceback_matrix[i[0], i[1]] = 'left'
            elif up_score >= diag_score and up_score >= 0:
                score_matrix[i[0], i[1]] = up_score
                traceback_matrix[i[0], i[1]] = 'up'
            elif diag_score >= 0: 
                score_matrix[i[0], i[1]] = diag_score
                traceback_matrix[i[0], i[1]] = 'diag'
            else:
                score_matrix[i[0], i[1]] = 0
                traceback_matrix[i[0], i[1]] = 'done'
                
            if score_matrix[i[0], i[1]] > highest_score:
                highest_row = i[0]
                highest_column = i[1]
                highest_score = score_matrix[i[0], i[1]]
    i = highest_row
    j = highest_column
    curr = traceback_matrix[i, j] # start at the bottom right of the matrix
    final_alignment_sequence1 = ""
    final_alignment_sequence2 = ""
    sequence1_list = sequence1_list[:highest_row]
    sequence2_list = sequence2_list[:highest_column]
    # record the final sequence alignment based on where the traceback matrix is pointing by taking the last value of the sequence list and popping it to the front of the final alignment
    while curr != 'done':
        if curr == 'diag':
            final_alignment_sequence1 = sequence1_list.pop() + final_alignment_sequence1
            final_alignment_sequence2 = sequence2_list.pop() + final_alignment_sequence2
            i -= 1
            j -= 1
        elif curr == 'up':
            final_alignment_sequence1 = sequence1_list.pop() + final_alignment_sequence1
            final_alignment_sequence2 = "-" + final_alignment_sequence2
            i -= 1
        elif curr == 'left':
            final_alignment_sequence2 = sequence2_list.pop() + final_alignment_sequence2
            final_alignment_sequence1 = "-" + final_alignment_sequence1
            j -= 1
        curr = traceback_matrix[i, j]
    return final_alignment_sequence1 + "\n" + final_alignment_sequence2 + "\nScore: " + str(highest_score)