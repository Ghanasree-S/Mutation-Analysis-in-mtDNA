import numpy as np  # Import NumPy for matrix operations

# Scoring scheme
MATCH_SCORE = 2
MISMATCH_PENALTY = -1
GAP_PENALTY = -2

def local_alignment(seq1, seq2):
    """Smith-Waterman Algorithm for Local Alignment"""
    
    m, n = len(seq1), len(seq2)

    # Initialize scoring matrix
    score_matrix = np.zeros((m+1, n+1))

    # Track the max score for traceback
    max_score = 0
    max_pos = None

    # Fill the scoring matrix using Smith-Waterman logic
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = score_matrix[i-1][j-1] + (MATCH_SCORE if seq1[i-1] == seq2[j-1] else MISMATCH_PENALTY)
            delete = score_matrix[i-1][j] + GAP_PENALTY
            insert = score_matrix[i][j-1] + GAP_PENALTY

            score_matrix[i][j] = max(0, match, delete, insert)  # Local alignment allows scores to reset at 0

            # Track the maximum score and position for traceback
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    # If no alignment was found, return empty strings
    if max_pos is None:
        return "", "", 0

    # Traceback to reconstruct the aligned sequences
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = max_pos

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        current_score = score_matrix[i][j]

        if current_score == score_matrix[i-1][j-1] + (MATCH_SCORE if seq1[i-1] == seq2[j-1] else MISMATCH_PENALTY):
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif current_score == score_matrix[i-1][j] + GAP_PENALTY:
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append('-')
            i -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j-1])
            j -= 1

    # Reverse the aligned sequences to restore original order
    aligned_seq1 = ''.join(aligned_seq1[::-1])
    aligned_seq2 = ''.join(aligned_seq2[::-1])

    return aligned_seq1, aligned_seq2, int(max_score)  # Ensure max_score is returned as an integer
