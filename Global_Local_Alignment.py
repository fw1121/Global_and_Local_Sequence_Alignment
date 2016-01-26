# initial matrix all is 0
def zeros(shape):
    return [[0 for i in range(shape[1])] for i in range(shape[0])]

# define scores matrix
alphabet = ['A', 'C', 'G', 'T']
score = [[5, -4, -2, -4, -8],
         [-4, 5, -4, -2, -8],
         [-2, -4, 5, -4, -8],
         [-4, -2, -4, 5, -8],
         [-8, -8, -8, -8, -8]]
# match: 5
# synonym mismatch: -2
# no-synonymous mismatch: -4
# gap: -8

def finalize(seq1, seq2):
    align1 = seq1[::-1]    #reverse sequence 1
    align2 = seq2[::-1]    #reverse sequence 2
    
    #calcuate identity, score and aligned sequeces
    symbol = ''
    scoreAll = 0
    identity = 0
    for i in range(0,len(align1)):
        # if two AAs are the same, then output the letter
        if align1[i] == align2[i]:                
            symbol = symbol + "|"
            identity += 1
            scoreAll += score[alphabet.index(align1[i])][alphabet.index(align2[i])]
    
        # if they are not identical and none of them is gap
        elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-': 
            scoreAll += score[alphabet.index(align1[i])][alphabet.index(align2[i])]
            symbol += ' '
        
        #if one of them is a gap, output a space
        elif align1[i] == '-' or align2[i] == '-':          
            symbol += ' '
            scoreAll += score[-1][-1]
    
    identity = float(identity) / len(align1) * 100
    
    print 'Identity =', "%3.3f" % identity, 'percent'
    print 'Score =', scoreAll
    print "Seq1: " + align1
    print "      " + symbol
    print "Seq2: " + align2
    
def globalAlignment(seq1, seq2):
    """
    global alignment
    """
    # Generate Dynamic Programming (DP) table
    m, n = len(seq1), len(seq2)
    scoreRecord = zeros((m+1, n+1))
    # initial DP of first column and first row
    for i in range(0, m + 1):
        scoreRecord[i][0] = score[alphabet.index(seq1[i-1])][-1] * i
    for j in range(0, n + 1):
        scoreRecord[0][j] = score[-1][alphabet.index(seq2[j-1])] * j

    # alignment
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match  = scoreRecord[i-1][j-1] + score[alphabet.index(seq1[i-1])][alphabet.index(seq2[j-1])]
            delete = scoreRecord[i-1][j] + score[alphabet.index(seq1[i-1])][-1]
            insert = scoreRecord[i][j-1] + score[-1][alphabet.index(seq2[j-1])]
            scoreRecord[i][j] = max(match, delete, insert)

    # Traceback and compute the alignment
    align1, align2 = '', ''
    i, j = m, n
    while i > 0 and j > 0:
        score_current  = scoreRecord[i][j]
        score_diagonal = scoreRecord[i-1][j-1]
        score_up = scoreRecord[i][j-1]
        score_left = scoreRecord[i-1][j]
        
        if score_current == score_diagonal + score[alphabet.index(seq1[i-1])][alphabet.index(seq2[j-1])]:
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif score_current == score_left + score[-1][alphabet.index(seq2[j-1])]:
            align1 += seq1[i-1]
            align2 += '-'
            i -= 1
        elif score_current == score_up + score[alphabet.index(seq1[i-1])][-1]:
            align1 += '-'
            align2 += seq2[j-1]
            j -= 1
    # Finish trace up to the top left cell
    while i > 0:
        align1 += seq1[i-1]
        align2 += '-'
        i -= 1
    while j > 0:
        align1 += '-'
        align2 += seq2[j-1]
        j -= 1
    finalize(align1, align2)
    
def localAlignment(seq1, seq2):
    m, n = len(seq1), len(seq2)
    scoreRecord = zeros((m+1, n+1)) # the DP table  
    
    # Initialize first column
    for i in range(1, len(x)+1):
        scoreRecord[i][0] = scoreRecord[i-1][0] + score[alphabet.index(x[i-1])][-1]
    # Calculate DP table and mark pointers
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match  = scoreRecord[i-1][j-1] + score[alphabet.index(seq1[i-1])][alphabet.index(seq2[j-1])]
            delete = scoreRecord[i-1][j] + score[alphabet.index(seq1[i-1])][-1]
            insert = scoreRecord[i][j-1] + score[-1][alphabet.index(seq2[j-1])]
            scoreRecord[i][j] = max(0, match, delete, insert)

    # find max
    i_index = 0
    j_index = 0
    max_score = 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if scoreRecord[i][j] >= max_score:
                max_score = scoreRecord[i][j]
                i_index = i
                j_index = j

    align1, align2 = '', ''    # initial sequences
    #traceback, follow pointers
    while i_index > 0 and j_index > 0:
        score_current  = scoreRecord[i_index][j_index]
        score_diagonal = scoreRecord[i_index-1][j_index-1]
        score_up = scoreRecord[i_index][j_index-1]
        score_left = scoreRecord[i_index-1][j_index]
        if score_current == 0:
            break 
        elif score_current == score_left + score[-1][alphabet.index(seq2[j-1])]:
            align1 += seq1[i_index-1]
            align2 += '-'
            i_index -= 1
        elif score_current == score_up + score[alphabet.index(seq1[i-1])][-1]:
            align1 += '-'
            align2 += seq2[j_index-1]
            j_index -= 1
        else:
            align1 += seq1[i_index-1]
            align2 += seq2[j_index-1]
            i_index -= 1
            j_index -= 1
    finalize(align1, align2)
