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
    
    if 'scoring' not in locals(): # if the code is using match and mismatch penalties
        for i, value in np.ndenumerate(score_matrix): # iterate through every value and record index
            if i[0] == 0 or i[1] == 0:
                continue
            if sequence1_list[i[0]-1] == sequence2_list[i[1]-1]: # if there is a match
                diag_score = score_matrix[i[0]-1, i[1]-1] + match_score
            else: # if there is not a match
                diag_score = score_matrix[i[0]-1, i[1]-1] + mismatch_penalty 
            left_score = score_matrix[i[0], i[1]-1] + gap_penalty
            up_score = score_matrix[i[0]-1, i[1]] + gap_penalty
            if left_score >= diag_score and left_score >= up_score: # if the left gap score is the greatest
                score_matrix[i[0], i[1]] = left_score
                traceback_matrix[i[0], i[1]] = 'left'
            elif up_score >= diag_score: # if the up gap score is the greatest
                score_matrix[i[0], i[1]] = up_score
                traceback_matrix[i[0], i[1]] = 'up'
            else: # if the no gap score is the greatest
                score_matrix[i[0], i[1]] = diag_score
                traceback_matrix[i[0], i[1]] = 'diag'
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
            if left_score >= diag_score and left_score >= up_score:
                score_matrix[i[0], i[1]] = left_score
                traceback_matrix[i[0], i[1]] = 'left'
            elif up_score >= diag_score:
                score_matrix[i[0], i[1]] = up_score
                traceback_matrix[i[0], i[1]] = 'up'
            else: 
                score_matrix[i[0], i[1]] = diag_score
                traceback_matrix[i[0], i[1]] = 'diag'
    
    i = score_matrix.shape[0] - 1
    j = score_matrix.shape[1] - 1
    curr = traceback_matrix[i, j] # start at the bottom right of the matrix
    final_alignment_sequence1 = ""
    final_alignment_sequence2 = ""
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
    return final_alignment_sequence1 + "\n" + final_alignment_sequence2 + "\nScore: " + str(score_matrix[-1, -1])

# adding gap extention penalty
def align_with_affine_gaps(sequence1, sequence2, open_gap_penalty, extend_gap_penalty, table=None, match_score=None, mismatch_penalty=None):
    # scoring based on table or match and mismatch score
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
    
    score_matrix = np.zeros((len(sequence1_list) + 1, len(sequence2_list) + 1, 3)) # DIMENSION 1: rows, DIMENSION 2: columns, DIMENSION 3: confined to 3 values, one value for each "Level of Manhattan"
    # set default values for score matrix
    score_matrix[0,:,0] = np.NINF
    score_matrix[:,0,0] = np.NINF
    score_matrix[0,0,0] = 0
    score_matrix[0,:,1] = np.NINF
    score = open_gap_penalty
    for i, value in np.ndenumerate(score_matrix[:,0,1]):
        score_matrix[i[0],0,1] = score
        score += extend_gap_penalty
    score_matrix[:,0,2] = np.NINF
    score = open_gap_penalty
    for i, value in np.ndenumerate(score_matrix[0,:,2]):
        score_matrix[0,i[0],2] = score
        score += extend_gap_penalty

    # set default values for traceback matrix
    traceback_matrix = np.zeros((len(sequence1_list) + 1, len(sequence2_list) + 1, 3), 'U1')
    traceback_matrix[:,0,1] = 'x'
    traceback_matrix[0,:,1] = 'y'
    traceback_matrix[0,0,:] = 'd'
    
    if 'scoring' not in locals(): # if using mismatch and match penalties
        for i, value in np.ndenumerate(score_matrix):
            if i[0] == 0 or i[1] == 0:
                continue
            if i[2] == 0: # alignment score
                if sequence1_list[i[0]-1] == sequence2_list[i[1]-1]: # if there is a match
                    m_score = score_matrix[i[0]-1, i[1]-1, 0] + match_score
                    x_score = score_matrix[i[0]-1, i[1]-1, 1] + match_score
                    y_score = score_matrix[i[0]-1, i[1]-1, 2] + match_score
                else: # if there is a mismatch
                    m_score = score_matrix[i[0]-1, i[1]-1, 0] + mismatch_penalty
                    x_score = score_matrix[i[0]-1, i[1]-1, 1] + mismatch_penalty
                    y_score = score_matrix[i[0]-1, i[1]-1, 2] + mismatch_penalty
                if y_score >= x_score and y_score >= m_score: # if the score coming from y matrix is greatest
                    score_matrix[i[0], i[1], 0] = y_score
                    traceback_matrix[i[0], i[1], 0] = "y"
                elif x_score >= m_score: # if the score coming from x matrix is greatest
                    score_matrix[i[0], i[1], 0] = x_score
                    traceback_matrix[i[0], i[1], 0] = "x"
                else: # if the score coming from m matrix is greatest
                    score_matrix[i[0], i[1], 0] = m_score
                    traceback_matrix[i[0], i[1], 0] = "m"
                    
            elif i[2] == 1: # gap on x axis score
                m_score = score_matrix[i[0] - 1, i[1], 0] + open_gap_penalty + extend_gap_penalty # open gap penalty
                x_score = score_matrix[i[0] - 1, i[1], 1] + extend_gap_penalty # extend gap penalty
                if x_score >= m_score:
                    score_matrix[i[0], i[1], 1] = x_score
                    traceback_matrix[i[0], i[1], 1] = "x"
                else:
                    score_matrix[i[0], i[1], 1] = m_score
                    traceback_matrix[i[0], i[1], 1] = "m"
            
            elif i[2] == 2: # gap on y axis score
                m_score = score_matrix[i[0], i[1]-1, 0] + open_gap_penalty + extend_gap_penalty
                y_score = score_matrix[i[0], i[1]-1, 2] + extend_gap_penalty
                if y_score >= m_score:
                    score_matrix[i[0], i[1], 2] = y_score
                    traceback_matrix[i[0], i[1], 2] = "y"
                else:
                    score_matrix[i[0], i[1], 2] = m_score
                    traceback_matrix[i[0], i[1], 2] = "m"
                    
    else: # if there is a scoring table
    # same as without except for alignment score
        for i, value in np.ndenumerate(score_matrix):
            if i[0] == 0 or i[1] == 0:
                continue
            
            if i[2] == 0:
                try:
                    table_score = scoring.loc[sequence1_list[i[0]-1], sequence2_list[i[1]-1]] # get score from table
                except:
                    return "Make sure that your sequences match the input table"
                m_score = score_matrix[i[0]-1, i[1]-1, 0] + table_score # add score from table to score
                x_score = score_matrix[i[0]-1, i[1]-1, 1] + table_score
                y_score = score_matrix[i[0]-1, i[1]-1, 2] + table_score
                if y_score >= x_score and y_score >= m_score:
                    score_matrix[i[0], i[1], 0] = y_score
                    traceback_matrix[i[0], i[1], 0] = "y"
                elif x_score >= m_score:
                    score_matrix[i[0], i[1], 0] = x_score
                    traceback_matrix[i[0], i[1], 0] = "x"
                else:
                    score_matrix[i[0], i[1], 0] = m_score
                    traceback_matrix[i[0], i[1], 0] = "m"
                    
            elif i[2] == 1:
                m_score = score_matrix[i[0] - 1, i[1], 0] + open_gap_penalty + extend_gap_penalty
                x_score = score_matrix[i[0] - 1, i[1], 1] + extend_gap_penalty
                if x_score >= m_score:
                    score_matrix[i[0], i[1], 1] = x_score
                    traceback_matrix[i[0], i[1], 1] = "x"
                else:
                    score_matrix[i[0], i[1], 1] = m_score
                    traceback_matrix[i[0], i[1], 1] = "m"
            
            elif i[2] == 2:
                m_score = score_matrix[i[0], i[1]-1, 0] + open_gap_penalty + extend_gap_penalty
                y_score = score_matrix[i[0], i[1]-1, 2] + extend_gap_penalty
                if y_score >= m_score:
                    score_matrix[i[0], i[1], 2] = y_score
                    traceback_matrix[i[0], i[1], 2] = "y"
                else:
                    score_matrix[i[0], i[1], 2] = m_score
                    traceback_matrix[i[0], i[1], 2] = "m"
                    
    if score_matrix[-1, -1, 2] >= score_matrix[-1, -1, 1] and score_matrix[-1, -1, 2] >= score_matrix[-1, -1, 0]: # if y is the best matrix to start on
        k = 2
    elif score_matrix[-1, -1, 1] >= score_matrix[-1, -1, 0]: # if x is the best matrix to start on
        k = 1
    else: # if m is the best matrix to start on
        k = 0
    i = score_matrix.shape[0] - 1
    j = score_matrix.shape[1] - 1
    curr = traceback_matrix[i, j, k]
    score = str(score_matrix[i, j, k]) # alignment score
    final_alignment_sequence1 = ""
    final_alignment_sequence2 = ""
    while curr != 'd':
        if k == 0: # if you are aligning
            i -= 1
            j -= 1
            final_alignment_sequence1 = sequence1_list.pop() + final_alignment_sequence1 # continue first sequence
            final_alignment_sequence2 = sequence2_list.pop() + final_alignment_sequence2 # gap on second sequence
        elif k == 1: # if you are creating a gap on the x-axis
            i -= 1
            final_alignment_sequence1 = sequence1_list.pop() + final_alignment_sequence1 # continue first sequence
            final_alignment_sequence2 = "-" + final_alignment_sequence2 # gap on second sequence
        elif k == 2: # if you are creating a gap on the y-axis
            j -= 1
            final_alignment_sequence1 = "-" + final_alignment_sequence1
            final_alignment_sequence2 = sequence2_list.pop() + final_alignment_sequence2
        if curr == 'm': # if the traceback matrix is pointing to m
            k = 0
            curr = traceback_matrix[i, j, 0] # move to m matrix
        elif curr == 'x': # if traceback matrix is pointing to x
            k = 1
            curr = traceback_matrix[i, j, 1]
        elif curr == 'y': # if traceback matrix is pointing to y
            k = 2
            curr = traceback_matrix[i, j, 2]    
    return final_alignment_sequence1 + "\n" + final_alignment_sequence2 + "\nScore: " + score # return final alignment and score

print(align_with_affine_gaps("CGGGGTATTTAACGTA", "AATTACGTATTTT", open_gap_penalty=-2, extend_gap_penalty=-0.5, table = "purine_pyrimidine_nucleotide"))