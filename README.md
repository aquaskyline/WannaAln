## WannaAln
WannaAln takes in a pair of FASTQ file and an arbitary length pattern, calculates the edit distance of 4 pairs including "R1 x pattern", "R2 x pattern", "R1 x ReverseComplement(pattern)" and "R2 x ReverseComplement(pattern)", and output the reads to a pair of new FASTQ file if any of the 4 edit distances go below a threshold.

## Parameters
```
wannaAln
Ruibang Luo
rluo5@jhu.edu

-a [STR]   Input R1
-b [STR]   Input R2
-x [STR]   Output R1, gzipped
-y [STR]   Output R2, gzipped
-p [STR]   Matching pattern
-m [INT]   Edit distance allowed, default 1
-h         Show help information
```

## Getting Started
```
https://github.com/aquaskyline/WannaAln.git
cd WannaAln
make
cd test
sh test.sh
```
