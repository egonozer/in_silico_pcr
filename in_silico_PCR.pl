#!/usr/bin/env perl

#in_silico_pcr.pl
#Copyright (C) 2017-2025 Egon A. Ozer
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see [http://www.gnu.org/licenses/].

my $version = "0.6";

## changes from v0.5.1
## added option for more than 1 mismatch with the -m option
## added Needleman-Wunsch alignment of the primer sequences to the matched sequence. Now outputs the numbers of primer mismatches and indels
## changed to use GetOpts::Long for getting options so as to allow the -m option given on its own to correspond to 1 mismatch (for backwards compatibility) 
## added option -f to print amplicon sequences to a file instead of to STDERR
## added inosine (I) as possible base in primer sequences

## changes from v0.5
## fix bug resulting in an error when forward primer sequence was found at the end of the input sequence. This only resulted in an error when -c was not given.
## include bzip2 support

## changes from v0.4
## fixed bug where if amplicons started at the same position on multiple different contigs, only one would be reported

## changes from v0.3
## can take a list of primer sequences. This can save time on repeatedly reading the sequence file.
## changed output format somewhat
## when only one primer is bound and the -c option is given, only output sequence from the primer up to the maximum band length (-l)
## fixed some bugs that could have resulted in incorrect reverse complementation if there was in indel in the first primer sequence

## changes from v0.2
## when using option -c, will now output all primer binding results, regardless of whether the distance from the primer to the end of the sequence is less than the maximum amplicon length or not. 
## added ability to find primer sequences with up to one indel

use strict;
use Cwd;
use warnings;

## Usage
my $usage = "
in_silico_PCR.pl (version $version)

Adapted from Joseba Bikandi's php script
(http://www.biophp.org/minitools/pcr_amplification/)

Required:
  -s <string>   Sequence file in fasta format
                (can be gzipped or bzipped, must have .gz or .bz2 suffix)
  
Options:
  -a <string>   Forward primer
  -b <string>   Reverse primer
    OR
  -p            file with multiple primer sets
                Format should be:
                    forward_primer<tab>reverse_primer<tab>output_prefix
                with one primer set per line.
                If this option is given, values given to -a and -b will be
                ignored
  -l <number>   Maximum length of the resulting \"band\" (default: 3000 bases)
  -m            Number of mismatches per primer sequence to allow
                (default, no mismatches)
  -i            Allow up to one insertion or deletion per primer sequence
                (default, no indels)
  -e            Exclude primer sequences from amplicon sequence, position, and
                length (default: primer sequences included)
  -r            Ensure that amplicons are in the orientation given by the order
                of the primers, i.e. forward primer at 5' end, reverse primer
                at 3' end (default: amplicons will be output in the orientation
                found in the input sequence)
  -c            If no amplicons are found in the sequence, will instead output
                all primer hits followed by sequence to the end of the sequence
                unit or the maximum \"band\" length, whichever is shorter. Could
                be helpful when trying to amplify across a sequencing contig
                break.
  -f            Filename to output amplicon sequences (default: amplicon sequences
                output to STDERR)

";

use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
my ($opt_s, $opt_a, $opt_b, $opt_p, $opt_l, $opt_m, $opt_i, $opt_e, $opt_r, $opt_c, $opt_f);
GetOptions(
    's=s'   => \$opt_s,
    'a=s'   => \$opt_a,
    'b=s'   => \$opt_b,
    'p=s'   => \$opt_p,
    'l=i'   => \$opt_l,
    'm:s'   => \$opt_m,
    'i'     => \$opt_i,
    'e'     => \$opt_e,
    'r'     => \$opt_r,
    'c'     => \$opt_c,
    'f=s'   => \$opt_f,
);

die $usage unless ($opt_s);
die $usage unless (($opt_a and $opt_b) or $opt_p);

my $seqfile     = $opt_s if $opt_s;
my $primer1     = $opt_a if $opt_a;
my $primer2     = $opt_b if $opt_b;
my $pfile       = $opt_p if $opt_p;
my $maxlength   = $opt_l ? $opt_l : 3000;

my $mm = 0;
if (defined $opt_m){
    if ($opt_m =~ m/^\s*$/){ ## backwards compatibility with prior versions. If -m is given without arguments, it will default to 1
        $mm = 1;
    } else {
        $opt_m =~ s/\s//g;
        if ($opt_m =~ m/^\d+$/){
            $mm = $opt_m;
        } else {
            die "ERROR: Number of potential mismatches (-m) must be 0 or a positive integer. Give '-m' without arguments to allow for up to 1 mismatch\n";
        }
    }
}

my $in;
if ($seqfile =~ m/\.gz$/){
    open ($in, "gzip -cd '$seqfile' |") or die "Can't open $seqfile: $!\n";
} elsif ($seqfile =~ m/\.bz2$/){
    open ($in, "bzip2 -cdk '$seqfile' |") or die "Can't open $seqfile: $!\n";
} else {
    open ($in, "<", $seqfile) or die "Can't open $seqfile\n";
}
my %sequences;
my ($id, $seq);
while (my $line = <$in>){
    chomp $line;
    if ($line =~ m/^>/){
        if ($id){
            $sequences{$id} = $seq;
            $seq = "";
        }
        $id = substr($line,1);
        $id =~ s/\s.*$//; #delete everything after the first space
        next;
    }
    $line = uc($line);
    $line =~ s/\s|\W|\d//g;
    $seq .= $line;
}
close ($in);
$sequences{$id} = $seq;
$seq = "";

my @iupac = (["R", "[AG]"],
             ["Y", "[CT]"],
             ["S", "[GC]"],
             ["W", "[AT]"],
             ["K", "[GT]"],
             ["M", "[AC]"],
             ["B", "[CGT]"],
             ["D", "[AGT]"],
             ["H", "[ACT]"],
             ["V", "[ACG]"],
             ["I", "."],
             ["N", "."],
             ["U", "T"]);

my @prims;
if ($pfile){
    open (my $pin, "<", $pfile) or die "Can't open $pfile: $!\n";
    while (my $fline = <$pin>){
        $fline =~ s/\R/\012/g; #converts to UNIX-style line endings
        my @lines = split("\n", $fline); #need to split lines by line-ending character in the case of Mac-formatted files which only have CR line terminators, not both CR and LF like DOS
        while (@lines){
            my $line = shift @lines;
            my ($for, $rev, $ipref) = split("\t", $line);
            $for =~ s/\s//g;
            $rev =~ s/\s//g;
            push @prims, ([$for, $rev, $ipref]);
        }
    }
    close ($pin);
} else {
    push @prims, ([$primer1, $primer2, ""]);
}

my $ampfile;
if ($opt_f){
    open ($ampfile, ">$opt_f") or die "ERROR: Can't open $opt_f for writing: $!\n";
}

print "AmpId\tSequenceId\tPositionInSequence\tLength\tP1Mismatch\tP1Indel\tP2Mismatch\tP2Indel\n";

foreach my $slice (@prims){
    my $pref;
    ($primer1, $primer2, $pref) = @{$slice};
    
    ## All non-word characters (\W) and digits(\d) are removed from primers and from sequence file
    ## I'll need to modify a bit to accept multi-contig files
    
    $primer1 = uc($primer1);
    $primer1 =~ s/\W|\d//g;
    $primer2 = uc($primer2);
    $primer2 =~ s/\W|\d//g;
    
    ##SET PATTERNS FROM PRIMERS
    my $pattern1 = $primer1;
    my $pattern2 = $primer2;
    ## If missmatches are allowed, create new pattern
    ## example: pattern="ACGT"; to allow one missmatch pattern=".CGT|A.GT|AC.T|ACG."
    if ($mm > 0){
        $pattern1 = includeN($primer1);
        $pattern2 = includeN($primer2);
    }
    ## If one indel is allowed, create new pattern
    ## example: pattern="ACGT"; to allow one indel pattern="A.CGT|AC.GT|ACG.T|AGT|ACT"
    ## will only add indels to original sequence, not sequences with mismatches
    my ($indel1, $indel2);
    if ($opt_i){
        $indel1 = indel($primer1);
        $indel2 = indel($primer2);
    }
    
    ## Change non-standard nucleotides in primers
    for my $i (0 .. $#iupac){
        my ($from, $to) = @{$iupac[$i]};
        $pattern1 =~ s/$from/$to/g;
        $pattern2 =~ s/$from/$to/g;
        if ($opt_i){
            $indel1 =~ s/$from/$to/g;
            $indel2 =~ s/$from/$to/g;
        }
    }
    
    ## SET PATTERN
    my $start_pattern = "$pattern1|$pattern2";
    my $end_pattern = reverse ($start_pattern);
    $end_pattern =~ tr/ACTG[]/TGAC][/;
    my ($full_start, $full_end) = ($start_pattern, $end_pattern);
    if ($opt_i){
        my $start_indel = "$indel1|$indel2";
        my $end_indel = reverse($start_indel);
        $end_indel =~ tr/ACTG[]/TGAC][/;
        $full_start .= "|$start_indel";
        $full_end .= "|$end_indel";
    }
    ## For determining strand
    my ($for_pattern, $rev_pattern) = ($pattern1, $pattern2);
    if ($opt_i){
        $for_pattern .= "|$indel1";
        $rev_pattern .= "|$indel2";
    }
    
    my %results_hash = Amplify($full_start, $full_end, $maxlength, $for_pattern, $primer1, $primer2);
    
    ## PRINT RESULTS
    if (scalar keys %results_hash > 0){
        my $count = 0;
        foreach my $id (sort {$a cmp $b} keys %results_hash){
            foreach my $key (sort {$a<=>$b} keys %{$results_hash{$id}}){
                $count++;
                my ($val, $rc, $fmm, $fid, $rmm, $rid) = @{$results_hash{$id}{$key}};
                my $sequence = $sequences{$id};
                my $amp = substr($sequence, $key, $val);
                if ($opt_i){
                    #check to make sure that an extra base wasn't added to the end of the sequence due to an indel pattern
                    my @check = split(/($start_pattern)/, $amp);
                    if (length($check[0]) == 1){
                        $amp = substr($amp, 1);
                        $val--;
                        $key++;
                    } elsif (length($check[$#check]) == 1){
                        $amp = substr($amp, 0, length($amp) - 1);
                        $val--;
                    }
                    my @check2 = split(/$end_pattern/, $amp);
                    if (length($check2[0]) == 1){
                        $amp = substr($amp, 1);
                        $val--;
                        $key++;
                    } elsif (length($check2[$#check2]) == 1){
                        $amp = substr($amp, 0, length($amp) - 1);
                        $val--;
                    }
                }
                my $ampID;
                $ampID = "$pref\_" if $pref;
                $ampID .= "amp_$count";
                print "$ampID\t$id\t". ($key + 1) ."\t$val\t$fmm\t$fid\t$rmm\t$rid";
                if ($opt_r){
                    print "\tcomplement" if $rc eq "y";
                }
                print "\n";
                if ($opt_r and $rc eq "y"){
                    $amp = reverse($amp);
                    $amp =~ tr/ACTGRYKMBVDHactgrykmbvdh/TGACYRMKVBHDtgacyrmkvbhd/;
                }
                my $outamp = ">$ampID\n$amp\n";
                if ($opt_f){
                    print $ampfile $outamp;
                } else {
                    print STDERR $outamp;
                }
            }
        } 
    } else {
        if ($opt_c){
            my @c_array = C_amplify($for_pattern, $rev_pattern, $maxlength, $primer1, $primer2);
            if (scalar @c_array > 0){
                @c_array = sort{$a->[0] cmp $b->[0]} @c_array;
                foreach my $slice (@c_array){
                    my ($type, $id, $pos, $seq, $fmm, $fid, $rmm, $rid) = @{$slice};
                    if ($opt_i){
                        #check to make sure that an extra base wasn't added to the end of the sequence due to an indel pattern
                        my @check = split(/($start_pattern)/, $seq);
                        if (length($check[0]) == 1){
                            $seq = substr($seq, 1);
                            $pos++;
                        } elsif (length($check[$#check]) == 1){
                            $seq = substr($seq, 0, length($seq) - 1);
                        }
                        my @check2 = split(/$end_pattern/, $seq);
                        if (@check2){
                            if (length($check2[0]) == 1){
                                $seq = substr($seq, 1);
                                $pos++;
                            } elsif (length($check2[$#check2]) == 1){
                                $seq = substr($seq, 0, length($seq) - 1);
                            }
                        }
                    }
                    $type = "$pref\_$type" if $pref;
                    my $seqleng = length($seq);
                    if ($seqleng > $maxlength){
                        $seq = substr($seq, 0, $maxlength);
                        $seqleng = $maxlength;
                    }
                    print "$type\t$id\t$pos\t$maxlength\t$fmm\t$fid\t$rmm\t$rid\n";
                    my $outamp = ">$type\n$seq\n";
                    if ($opt_f){
                        print $ampfile $outamp;
                    } else {
                        print STDERR $outamp;
                    }
                }
            } else {
                print "$pref\t" if $pref;
                print "No amplification\n";
            }
        } else {
            print "$pref\t" if $pref;
            print "No amplification\n";
        }
    }
}
close ($ampfile) if $opt_f;

#----------------------------------------------------------------------
sub includeN {
    my $pattern = shift;
    my $leng = length($pattern);
    my $this_mm = $mm;
    $this_mm = $leng if ($mm > $leng);

    my %patterns;
    $patterns{$pattern} = 1;
    for my $i (1 .. $this_mm){
        my %newpatterns;
        foreach my $patternkey (keys %patterns){
            for my $j (0 .. ($leng-1)){
                my $newpattern = $patternkey;
                next if substr($newpattern, $j, 1) eq ".";
                substr($newpattern, $j, 1, ".");
                unless ($patterns{$newpattern}){
                    $newpatterns{$newpattern} = 1;
                }
            }
        }
        %patterns = %newpatterns;
    }
    my $new_pattern = join("|", sort keys %patterns);

    return ($new_pattern);
}

sub indel {
    my $pattern = shift;
    my ($ins, $del, $indel);
    if (length($pattern) > 2){
        $ins = substr($pattern, 0, 1) . "." . substr($pattern, 1);
        my $pos = 1;
        while ($pos < length($pattern) - 1){
            $ins .= "|" . substr($pattern, 0, 1 + $pos) . "." . substr($pattern, 1 + $pos);
            $pos++;
        }
        
        $del = substr($pattern, 0, 1) . substr($pattern, 2);
        $pos = 2;
        while ($pos < length($pattern) - 1){
            $del .= "|" . substr($pattern, 0, $pos) . substr($pattern, 1 + $pos);
            $pos++;
        }
        
        $indel = $ins . "|" . $del;
    }
    return ($indel);
}

sub Amplify {
    my $start_pattern = shift;
    my $end_pattern = shift;
    my $maxlength = shift;
    my $for_pattern = shift;
    my $p1 = shift;
    my $p2 = shift;

    my %results_hash;
    ## SPLIT SEQUENCE BASED IN $start_pattern (start positions of amplicons)
    
    foreach my $id (sort keys %sequences){
        my $sequence = $sequences{$id};
        my @fragments = split(/($start_pattern)/, $sequence);
        ####### print STDERR "fragments: ", join(" - ", @fragments), "\n";
        my $maxfragments = scalar @fragments;
        my $position = length($fragments[0]);
        if ($maxfragments > 1){
            for (my $m = 1; $m < $maxfragments ; $m+= 2){
                next unless $fragments[$m + 1];
                my $subfragment_to_maximum = substr($fragments[$m + 1], 0, $maxlength);
                my @fragments2 = split(/($end_pattern)/, $subfragment_to_maximum);
                ########### print STDERR "fragments2: ", join(" - ", @fragments2), "\n";
                
                if (scalar @fragments2 > 1){
                    my $lenfragment = length($fragments[$m].$fragments2[0].$fragments2[1]);
                    my $rc = "y";
                    my ($fmm, $fid, $rmm, $rid) = (0) x 4; ## mismatches and indels
                    my ($fprim, $rprim);
                    ############ print STDERR "fragments[m]:$fragments[$m] fragments2[1]:$fragments2[1]\n";
                    if ($fragments[$m] =~ m/$for_pattern/){
                        ############ print STDERR "Forward\n";
                        $rc = "n";
                        $fprim = $fragments[$m];
                        $rprim = reverse($fragments2[1]);
                        $rprim =~ tr/ACTGRYKMBVDHactgrykmbvdh/TGACYRMKVBHDtgacyrmkvbhd/;
                    } else {
                        $fprim = reverse($fragments2[1]);
                        $fprim =~ tr/ACTGRYKMBVDHactgrykmbvdh/TGACYRMKVBHDtgacyrmkvbhd/;
                        $rprim = $fragments[$m];
                    }
                    ($fmm, $fid) = nw($fprim, $p1) unless $fprim eq $p1;
                    ($rmm, $rid) = nw($rprim, $p2) unless $rprim eq $p2;
                    if (!$opt_e){
                        $results_hash{$id}{$position} = ([$lenfragment, $rc, $fmm, $fid, $rmm, $rid]);
                    } else {
                        my $new_pos = $position + length $fragments[$m];
                        my $new_len = length $fragments2[0];
                        $results_hash{$id}{$new_pos} = ([$new_len, $rc, $fmm, $fid, $rmm, $rid]);
                    }
                }
                $position += length($fragments[$m]) + length($fragments[$m+1]);
            }
        }
    }
    return(%results_hash);
}

sub C_amplify {
    my $pattern1 = shift;
    my $pattern2 = shift;
    my $maxlength = shift;
    my $p1 = shift;
    my $p2 = shift;
    
    my $pattern1_rc = reverse($pattern1);
    $pattern1_rc =~ tr/ACTG[]/TGAC][/;
    my $pattern2_rc = reverse($pattern2);
    $pattern2_rc =~ tr/ACTG[]/TGAC][/;
    
    my @results;
    my @counts = (0) x 4;
    
    foreach my $id (sort keys %sequences){
        my $sequence = $sequences{$id};
        my @fragments = split(/($pattern1)/, $sequence);
        my $maxfragments = scalar @fragments;
        my $position = 0;
        #$position = length($fragments[0]) if ($fragments[0] !~ m/$pattern1/);
        if (@fragments > 1){
            while (@fragments){
                my $frag = shift @fragments;
                my $newpos = $position + length($frag);
                if ($frag =~ m/$pattern1/){
                    my $outpos = $position;
                    $outpos = $newpos if $opt_e;
                    my $subfragment = substr($sequence, $outpos);
                    $counts[0]++;
                    my ($pmm, $pid) = (0) x 2;
                    ($pmm, $pid) = nw($frag, $p1) unless $frag eq $p1;
                    push @results, (["p1\_$counts[0]", $id, $outpos, $subfragment, $pmm, $pid, "-", "-"]);
                }
                $position = $newpos;
            }
        }
        
        @fragments = split(/($pattern1_rc)/, $sequence);
        $maxfragments = scalar @fragments;
        $position = 0;
        #$position = length($fragments[0]) if ($fragments[0] !~ m/$pattern1_rc/);
        if (@fragments > 1){
            while (@fragments){
                my $frag = shift @fragments;
                my $newpos = $position + length($frag);
                if ($frag =~ m/$pattern1_rc/){
                    my $outpos = $newpos;
                    $outpos = $position if $opt_e;
                    my $subfragment = substr($sequence, 0, $outpos);
                    $subfragment = reverse($subfragment);
                    $subfragment =~ tr/ACTGRYKMBVDHactgrykmbvdh/TGACYRMKVBHDtgacyrmkvbhd/;
                    $counts[1]++;
                    $frag = reverse($frag);
                    $frag =~ tr/ACTGRYKMBVDHactgrykmbvdh/TGACYRMKVBHDtgacyrmkvbhd/;
                    my ($pmm, $pid) = (0) x 2;
                    ($pmm, $pid) = nw($frag, $p1) unless $frag eq $p1;
                    push @results, (["p1rc_$counts[1]", $id, $outpos, $subfragment, $pmm, $pid, "-", "-"]);
                }
                $position = $newpos;
            }
        }
        
        @fragments = split(/($pattern2)/, $sequence);
        $maxfragments = scalar @fragments;
        $position = 0;
        #$position = length($fragments[0]) if ($fragments[0] !~ m/$pattern2/);
        if (@fragments > 1){
            while (@fragments){
                my $frag = shift @fragments;
                my $newpos = $position + length($frag);
                if ($frag =~ m/$pattern2/){
                    my $outpos = $position;
                    $outpos = $newpos if $opt_e;
                    my $subfragment = substr($sequence, $outpos);
                    $subfragment = reverse($subfragment);
                    $subfragment =~ tr/ACTGRYKMBVDHactgrykmbvdh/TGACYRMKVBHDtgacyrmkvbhd/;
                    $counts[2]++;
                    my ($pmm, $pid) = (0) x 2;
                    ($pmm, $pid) = nw($frag, $p2) unless $frag eq $p2;
                    push @results, (["p2_$counts[2]", $id, $outpos, $subfragment, "-", "-", $pmm, $pid]);
                }
                $position = $newpos;
            }
        }
        
        @fragments = split(/($pattern2_rc)/, $sequence);
        $maxfragments = scalar @fragments;
        $position = 0;
        #$position = length($fragments[0]) if ($fragments[0] !~ m/$pattern2_rc/);
        if (@fragments > 1){
            while (@fragments){
                my $frag = shift @fragments;
                my $newpos = $position + length($frag);
                if ($frag =~ m/$pattern2_rc/){
                    my $outpos = $newpos;
                    $outpos = $position if $opt_e;
                    my $subfragment = substr($sequence, 0, $outpos);
                    $counts[3]++;
                    $frag = reverse($frag);
                    $frag =~ tr/ACTGRYKMBVDHactgrykmbvdh/TGACYRMKVBHDtgacyrmkvbhd/;
                    my ($pmm, $pid) = (0) x 2;
                    ($pmm, $pid) = nw($frag, $p2) unless $frag eq $p2;
                    push @results, (["p2rc_$counts[3]", $id, $outpos, $subfragment, "-", "-", $pmm, $pid]);
                }
                $position = $newpos;
            }
        }
    }
    return(@results);
}

sub nw{
    #needleman-wunsch alignment
    #adapted from http://www.tcoffee.org/Courses/Exercises/Aln/Perl_examples/aln_algorithms.html
    my $seq1 = uc(shift);
    my $seq2 = uc(shift);

    #Parameters:
    my $match=10;
    my $mismatch=-10;
    my $gep=-10;

    $seq1 =~ s/\s//g;
    $seq1 =~ s/-//g;
    $seq2 =~ s/\s//g;
    $seq2 =~ s/-//g;
    my $len1 = length($seq1);
    my $len2 = length($seq2);

    my @smat;
    my @tb;
    for my $i (0 .. $len1){
        $smat[$i][0] = $i * $gep;
        $tb[$i][0] = 1;
    }
    for my $i (0 .. $len2){
        $smat[0][$i] = $i * $gep;
        $tb[0][$i] = -1;
    }

    for my $i (1 .. $len1){
        for my $j (1 .. $len2){
            #calculate the score
            my $score = $mismatch;
            $score = $match if (substr($seq1, $i-1, 1) eq substr($seq2, $j-1, 1));
            my $sub = $smat[$i-1][$j-1]+$score;
            my $del = $smat[$i][$j-1]+$gep;
            my $ins = $smat[$i-1][$j]+$gep;
            if ($sub > $del and $sub > $ins){
                $smat[$i][$j]=$sub;
                $tb[$i][$j]=0;
            } elsif ($del >$ins){
                $smat[$i][$j]=$del;
                $tb[$i][$j]=-1;
            } else {
                $smat[$i][$j]=$ins;
                $tb[$i][$j]=1;
            }
        }
    }
    my $i = $len1;
    my $j = $len2;
    my ($aln1, $aln2);
    my $mm = 0;
    my $indel = 0;
    while (!($i == 0 and $j == 0)){
        if ($tb[$i][$j] == 0){
            $aln1 .= substr($seq1, --$i, 1);
            $aln2 .= substr($seq2, --$j, 1);
            $mm++ if substr($aln1, -1, 1) ne substr($aln2, -1, 1);
        } elsif ($tb[$i][$j] == -1){
            $aln1 .= "-";
            $aln2 .= substr($seq2, --$j, 1);
            $indel++;
        } elsif ($tb[$i][$j] == 1){
            $aln1 .= substr($seq1, --$i, 1);
            $aln2 .= "-";
            $indel++;
        }
    }
    $aln1 = reverse($aln1);
    $aln2 = reverse($aln2);

    ## check for extra bases added by the indel pattern
    if (length($aln2) > length($seq2)) {
        if ($aln2 =~ m/^(-+)/){
            for my $i (1 .. length($1)){
                $indel--;
            }
        }
        if ($aln2 =~ m/(-+)$/){
            for my $i (1 .. length($1)) {
                $indel--;
            }
        }
    }
    #return($aln1, $aln2, $mm, $indel);
    return($mm, $indel);
}
