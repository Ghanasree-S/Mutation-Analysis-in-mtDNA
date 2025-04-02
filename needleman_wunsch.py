import numpy as np  # Import NumPy for matrix operations (O(1))

# Scoring scheme (constants, O(1))
match_score = 1
mismatch_penalty = -1
gap_penalty = -2

def global_alignment(seq1, seq2):
    m, n = len(seq1), len(seq2)  # Get lengths of both sequences (O(1))

    # Initialize the scoring matrix (O(mn))
    score_matrix = np.zeros((m+1, n+1))  # Create a (m+1) x (n+1) matrix (O(mn))

    # Initialize gap penalties in first row and column (O(m+n))
    for i in range(m+1):  # O(m)
        score_matrix[i][0] = i * gap_penalty
    for j in range(n+1):  # O(n)
        score_matrix[0][j] = j * gap_penalty

    # Fill the scoring matrix (O(mn))
    for i in range(1, m+1):  # O(m)
        for j in range(1, n+1):  # O(n)
            # Compute possible values for the cell
            match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)  # O(1)
            delete = score_matrix[i-1][j] + gap_penalty  # O(1)
            insert = score_matrix[i][j-1] + gap_penalty  # O(1)

            # Take the maximum value among match, delete, and insert
            score_matrix[i][j] = max(match, delete, insert)  # O(1)

    # Traceback to find the best alignment (O(m+n))
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = m, n

    while i > 0 or j > 0:  # O(m + n)
        current_score = score_matrix[i][j]  # O(1)

        if i > 0 and j > 0 and (current_score == score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)):  
            aligned_seq1.append(seq1[i-1])  # O(1)
            aligned_seq2.append(seq2[j-1])  # O(1)
            i -= 1
            j -= 1
        elif i > 0 and (current_score == score_matrix[i-1][j] + gap_penalty):  
            aligned_seq1.append(seq1[i-1])  # O(1)
            aligned_seq2.append('-')  # O(1)
            i -= 1
        else:  
            aligned_seq1.append('-')  # O(1)
            aligned_seq2.append(seq2[j-1])  # O(1)
            j -= 1

    # Reverse the aligned sequences (O(m + n))
    aligned_seq1 = aligned_seq1[::-1]  # O(m)
    aligned_seq2 = aligned_seq2[::-1]  # O(n)

    # final alignment score (O(1))
    final_score = score_matrix[m][n]

    # Return results (O(1))
    return ''.join(aligned_seq1), ''.join(aligned_seq2), final_score

