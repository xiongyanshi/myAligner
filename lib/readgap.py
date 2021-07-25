import re

"""
open=-2
extend=-2
"""

def readgap(filename):
    gap_open, gap_extend = re.search(r'open=(.*)\nextend=(.*)', open(filename).read()).groups()
    gap_open, gap_extend = map(lambda x: int(str.strip(x)), [gap_open, gap_extend])
    return {'gap_open':gap_open, 'gap_extend':gap_extend}

def main():
    filename = 'dna_gap.default.txt'
    print(readgap(filename))

if __name__ == '__main__':
    main()
