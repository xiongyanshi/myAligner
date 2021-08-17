import sys

def theta(a, b):
    if a == '-' or b == '-' or a != b:   # gap or mismatch
        return -1
    elif a == b:                         # match
        return 1

def make_score_matrix(seq1, seq2):
    """
    return score matrix and map(each score from which direction)
    0: diagnosis
    1: up
    2: left
    """
    seq1 = '-' + seq1
    seq2 = '-' + seq2
    score_mat = {}
    trace_mat = {}

    for i,p in enumerate(seq1):
        score_mat[i] = {}
        trace_mat[i] = {}
        for j,q in enumerate(seq2):
            if i == 0:                    # first row, gap in seq1
                score_mat[i][j] = -j
                trace_mat[i][j] = 1
                continue
            if j == 0:                    # first column, gap in seq2
                score_mat[i][j] = -i
                trace_mat[i][j] = 2
                continue
            ul = score_mat[i-1][j-1] + theta(p, q)     # from up-left, mark 0
            l  = score_mat[i][j-1]   + theta('-', q)   # from left, mark 1, gap in seq1
            u  = score_mat[i-1][j]   + theta(p, '-')   # from up, mark 2, gap in seq2
            picked = max([ul,l,u])
            score_mat[i][j] = picked
            trace_mat[i][j] = [ul, l, u].index(picked)   # record which direction
    return score_mat, trace_mat

def traceback(seq1, seq2, trace_mat):
    '''
    find one optimal traceback path from trace matrix, return path code
    -!- CAUTIOUS: if multiple equally possible path exits, only return one of them -!-
    '''
    seq1, seq2 = '-' + seq1, '-' + seq2
    i, j = len(seq1) - 1, len(seq2) - 1
    path_code = ''
    while i > 0 or j > 0:
        direction = trace_mat[i][j]
        if direction == 0:                    # from up-left direction
            i = i-1
            j = j-1
            path_code = '0' + path_code
        elif direction == 1:                  # from left
            j = j-1
            path_code = '1' + path_code
        elif direction == 2:                  # from up
            i = i-1
            path_code = '2' + path_code
    return path_code

def print_m(seq1, seq2, m):
    """print score matrix or trace matrix"""
    seq1 = '-' + seq1; seq2 = '-' + seq2
    print()
    print(' '.join(['%3s' % i for i in ' '+seq2]))
    for i, p in enumerate(seq1):
        line = [p] + [m[i][j] for j in range(len(seq2))]
        print(' '.join(['%3s' % i for i in line]))
    print()
    return

def pretty_print_align(seq1, seq2, path_code):
    '''
    return pair alignment result string from
    path code: 0 for match, 1 for gap in seq1, 2 for gap in seq2
    '''
    align1 = ''
    middle = ''
    align2 = ''
    for p in path_code:
        if p == '0':
            align1 = align1 + seq1[0]
            align2 = align2 + seq2[0]
            if seq1[0] == seq2[0]:
                middle = middle + '|'
            else:
                middle = middle + ' '
            seq1 = seq1[1:]
            seq2 = seq2[1:]
        elif p == '1':
            align1 = align1 + '-'
            align2 = align2 + seq2[0]
            middle = middle + ' '
            seq2 = seq2[1:]
        elif p == '2':
            align1 = align1 + seq1[0]
            align2 = align2 + '-'
            middle = middle + ' '
            seq1 = seq1[1:]

    #print('\n  ' + align1 + '\n  ' + middle + '\n  ' + align2 + '\n')
    print('\n  %s\n  %s  \n  %s\n' % (align1, middle, align2))
    return

def usage():
    print('Usage:\n\tpython nwAligner.py seq1 seq2\n')
    return

def main():
    try:
        seq1, seq2 = map(str.upper, sys.argv[1:3])
    except:
        seq1, seq2 = 'TCATC','TCATGGC'
        usage()
        print('--------Demo:-------\n')

    score_mat, trace_mat = make_score_matrix(seq1, seq2)
    #print_m(seq1, seq2, score_mat)
    #print_m(seq1, seq2, trace_mat)

    path_code = traceback(seq1, seq2, trace_mat)
    pretty_print_align(seq1, seq2, path_code)
    #print(path_code)

if __name__ == '__main__':
    main()

