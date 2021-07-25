import re, os

class Fasta:
    def __init__(self, data=None):
        if not data:
            pass
        elif os.path.isfile(data):
            Fasta.from_file(self, data)
        else:
            Fasta.from_string(self, data)

    def from_string(self, stringContent):
        """parse sequence file, get first two reads"""
        try:
            # more than 2 > in fasta file?
            name1, seq1, name2, seq2 = re.match(
                                 r'^>(.+?)\n([\w\W]+?)\n>(.+?)\n([\w\W]+)\n>.+',
                                 stringContent).groups()
        except AttributeError:
            # no just 2 > records in fasta file, safe to parse until file end
            name1, seq1, name2, seq2 = re.match(
                                 r'^>(.+?)\n([\w\W]+?)\n>(.+?)\n([\w\W]+)',
                                 stringContent).groups()
        seq1 = ''.join(seq1.splitlines()).upper()
        seq2 = ''.join(seq2.splitlines()).upper()

        self.name1 = name1
        self.seq1  = seq1
        self.name2 = name2
        self.seq2  = seq2

        self.read1 = '>%s\n%s' % (name1, seq1)
        self.read2 = '>%s\n%s' % (name2, seq2)

        return self

    def from_file(self, filename):
        return Fasta.from_string(self, open(filename).read())

    def __repr__(self):
        return '>%s\n%s\n>%s\n%s' % (self.name1, self.seq1,
                                     self.name2, self.seq2)

    def to_file(self, outfilename):
        with open(outfilename, 'w') as handle:
            handle.write('>%s\n%s\n>%s\n%s\n' % (self.name1, self.seq1,
                                               self.name2, self.seq2))

def main():
    fasta1 = Fasta('../demo/input1.fa')
    print(fasta1)
    fasta2 = Fasta('../demo/input2.fa')
    print(fasta2)

if __name__ == '__main__':
    main()
