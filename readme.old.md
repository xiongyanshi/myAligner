myAligner
> my Alignment program.

# Usage
## Getting start
```
# pairwise alignment with needleman-wunch algorithm:
$ python align_nw.py actt attc

  ACTT-
  | ||   
  A-TTC

# pairwise alignment with smith-waterman aligorithm:
$ python align_sw.py actt attc

  ACTT-
    || 
  -ATTC


```
## Advanced usage
### use your own substitution matrix
substitution matrix is a tab delimited table, first character of first row doesn't matter, only int will be accepted (default, lib/dna_sub.default.mat):
```
\	a	t	c	g
a	3	-3	-3	-3
t	-3	3	-3	-3
c	-3	-3	3	-3
g	-3	-3	-3	3
```
For instance, if you prefer more similarity between purines or pyrimidines, you may wanna use this one (lib/dna_sub.1.mat) :
```
\	a	t	c	g
a	3	-3	-3	1
t	-3	3	1	-3
c	-3	1	3	-3
g	1	-3	-3	3
```
### Another gap penalty scheme
Default gap penalty configuration file is "lib/dna_gap.default.txt" :
```
open=0
extend=-2
```
This means opening a gap will cost nothing, and extending a gap will -2. In other words, there is no difference between opening a gap and extending a gap. As such, "-" match to "A" means -2, "---" match to "ACG" means -2*3 = -6.

But, you may prefer opening a gap cost much more than extending a gap, this will make alignment program trying to find continuous body, for example in lib/dna_gap.1.txt:
```
open=-5
extend=-2
```
In this way, "-" match to "A" means -5 + -2*1 =-7, "---" match to "ACG" means -5 + -2*3 = -11.



