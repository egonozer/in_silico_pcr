# IN\_SILICO_PCR
***
## CONTENTS

1. **Introduction**
2. **Requirements**
3. **Installation**
4. **Usage**
5. **Output**
6. **License**
7. **Contact**

***
## 1. INTRODUCTION

This script searches nucleotide sequences using a set of primer sequences and attempts to anticipate the results of a PCR reaction using these same sequences. Results are based on sequence identity. Real-world factors such as melting temperature, annealing time, salt concentration, or secondary structure are not accounted for in this script.

This program was adapted and extended from a php script ([link](http://www.biophp.org/minitools/pcr_amplification/)) authored by Joseba Bikandi which was covered by the [GNU GPL v2](http://www.biophp.org/minitools/pcr_amplification/) license. 

## 2. REQUIREMENTS

Perl

## 3. INSTALLATION

Download and save.

## 4. USAGE
Basic command: `perl in_silico_PCR.pl -s query.fasta -a (forward primer sequence) -b (reverse primer sequence) > results.txt 2> amplicons.fasta`

For list of options, call the script without any inputs: `perl in_silico_PCR.pl`

### 4.1 Required Inputs:
`-s`  
Path to sequence file in fasta format. Can be a multi-fasta file with multiple sequence records.

**And Either**  
`-a`  
Sequence of forward primer, 5' → 3'  
`-b`  
Sequence of reverse primer, 5' → 3'  
**Or**  
`-p`  
Path to file containing one or more primer sequence pairs.  
File format should be:  
_forward primer sequence (5' → 3') \<tab> reverse primer sequence (5' → 3') \<tab> primer set ID (no spaces)_  
One primer set per line.  
Example:  

```
ATGACATAAGA	TACTACTGAA	primer_set_1
GATATACTGAG	CTAGACTACT	dnaA
```

Primer sequences containing [IUPAC ambiguity codes](https://www.bioinformatics.org/sms/iupac.html) are accepted.

### 4.2 Optional Inputs:

* `-l`  
Maximum "amplicon" length, in bp.  
(default: 3000)

* `-m`  
Allow up to one base mismatch per primer sequence.  
(default: no mismatches tolerated)

* `-i`  
Allow up to one base insertion or deletion per primer sequence.  
(default: no indels tolerated) 

* `-e`  
Exclude primer sequences from amplicon output sequence, position, and length.  
(default: primer sequences included in result)

* `-r`  
Ensure that output amplicons are in hte orientation given by the order of the primers, i.e. forward primer at 5' end, reverse primer at 3' end.  
(default: amplicons will be output in the orientation they are found in the query sequence)

* `-c`  
If no amplicons are found in the sequence, will instead output and and all primer hits followed by sequence either to the end of the contig or to the maximum amplicon length (given by `-l`), whichever is shorter.   
(default: only sequences where both primers hit in the correct orientation and within the maximum amplicon length are output) 

## 5. OUTPUT

A summary of each _in silico_ "PCR" reaction result will be output to STDOUT. The output can be directed to a file instead of the screen by using `> results_file.txt` on the command line. There are four columns in the results:

1. "AmpID": Name of the amplicon. It will include the name of the primer set (third column in the file given to `-p`) followed by "amp" for full amplicon (both primers in correct orientation within maximum distance), followed by a number to represent the number of amplicons found in the the input sequence. If the `-c` option is given above, then if one or both primer sequences is found in the input sequence, but not in the correct orientation or within the maximum amplicon size, then the prefix will be followed by "p1" to represent the first primer sequence or "p2" to represent the second primer sequence.
2. "SequenceId": The name of the input sequence in which the primer sequences were found. Will be "No amplification" if no amplicon was identified.
3. "Position in sequence": The 5'-most position (1-based) at which the amplicon or primer sequence starts in the input sequence. This value will vary depending on whether `-e` is given or not (see Options above)
4. "Length": The total length of the amplified sequence. This value will vary depending on whether `-e` is given or not (see Options above)

"Amplicon" sequences in fasta format will be output to STDERR. The output can be directed to a file instead of the screen by using `2> amplicons.fasta` on the command line.  


## 6. LICENSE:

in\_silico_pcr  
Copyright (C) 2017 Egon A. Ozer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  See LICENSE.txt

## 7. CONTACT:

Contact [Egon Ozer](e-ozer@northwestern.edu) with questions or comments.



> Written with [MacDown](https://macdown.uranusjr.com/).
 