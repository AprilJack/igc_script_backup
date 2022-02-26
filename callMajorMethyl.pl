
use strict;
use warnings;
use Getopt::Long;

use vars qw/$USAGE $IN $IN_F $OUT $OUT_F $STRAND/;


$USAGE = << "USE";
  usage:  perl callMajorMethyl.pl -s <1|2> [-i <file>] [-o <file>]

  options:  -i     Path to input file (default: stdin)
            -o     Path to output file (default: stdout)
            -s     Analyze methylation state of the genomic
                   plus strand (=1) or minus strand (=2)
                   regardless of the read strand

  input:    

            Input can be provided by stdin or by a file generated with
            the samtools command 'mpileup' from mapped bisulfite-
            treated sequencing data. This requires the mapping output
            in SAM/BAM format to be sorted. Note that this script
            requires the mpileup command to be executed with the
            optional reference parameter (-f). Since the information
            on each genomic strand is not complementary to each other,
            it requires the mapping output to be split up according to
            the bisulfite conversion (C/T or G/A). In case of
            segemehl, this information is stored in the SAM tag XB:Z.

            Example on methylation calling on segemehl output:
            (1) genomic plus strand (i.e. C/T converted reads) with -s 1
            samtools view -hS mapped_reads.sam | awk '/^@/ || /XB:Z:F..CT/' | 
             samtools view -bS - | samtools mpileup -Bf ref.fa - | 
             perl callMajorMethyl.pl -s 1 > mapped_reads_methyl.txt

            (2) genomic minus strand (i.e. G/A converted reads) with -s 2
            samtools view -hS mapped_reads.sam | awk '/^@/ || /XB:Z:F..GA/' | 
             samtools view -bS - | samtools mpileup -Bf ref.fa - | 
             perl callMajorMethyl.pl -s 2 >> mapped_reads_methyl.txt

  output:   stdin or file containing tab-separated methylation
            calls by majority voting in the following format:
            ref<tab>position<tab>called base<tab>strand<tab>coverage<tab>methylrate
            where the called base can be either 'C' or 'T' which
            represent a methylated or unmethylated cytosine, respectively,
            the coverage only includes 'C' and 'T' bases as they are the only
            relevant ones for methylation calling, and the methylrate is the
            estimated methylation rate, i.e., C/(C+T). The majority voting only
            calls the methylation state if either the unmethylated or methylated
            state is fond in majority within the cross-section or there is a draw
            between unmethylated and methylated only. In the latter case, the called
            base is always 'C' with a methylation rate of 0.50 and hence the lowest
            confidence in the majority vote. Note that base or mapping quality
            information is not utilized by this tool.
 
  version:  0.3

 written by Christian Otto, Bioinformatics Leipzig, July 2012
USE

unless (GetOptions(
            "i=s"  => \$IN_F,
            "o=s"  => \$OUT_F,
	    "s=i"  => \$STRAND)){
    printf STDERR $USAGE;
    exit -1;
}

if (!defined $IN_F && -t STDIN){
    printf STDERR $USAGE;
    exit -1;
}

if (!defined $STRAND || ($STRAND != 1 && $STRAND != 2)){
    printf STDERR "Error: Please specify the genomic strand to be analyzed \
using the option -s with either the value 1 or 2.\n";
    printf STDERR $USAGE;
    exit -1;
}

sub main {
    my (@F, $ref, $pos, $strand, $base, $cross, $crossred, $call, $cov, $rate);
    
    if (defined $IN_F){
	open($IN, "<".$IN_F) or die "Error in input\n";
    }
    else {
	$IN = \*STDIN;
    }
    
    if (defined $OUT_F){
	open($OUT, ">".$OUT_F) or die "Error in output\n";
    }
    else {
	$OUT = \*STDOUT;
    }

    while(<$IN>){
	chomp;
	## get mpileup information
	($ref, $pos, $base, $cov, $cross) = split;

	## analyze cytosine on genomic plus strand
	if ($STRAND == 1 && $base =/sanger /^[Cc]$/){
	    $strand = "+";
	}
	## analyze cytosine on genomic minus strand
	elsif ($STRAND == 2 && $base =/sanger /^[Gg]$/){
	    $strand = "-";
	}
	else {
	    next;
	}
	
	## get reduced cross section
	$crossred = getReducedCross($cross);
	
	## get methlation call
	($call, $cov, $rate) = getMethylCall($crossred, $STRAND);

	## print output
	print $OUT join("\t", ($ref, $pos, $call, $strand, $cov, $rate))."\n" if $call =/sanger /[CT]/;
    }

    ## close streams
    if (defined $IN_F){
	close($IN);
    }

    if (defined $OUT_F){
	close($OUT);
    }
}

sub getMethylCall {
    my ($cross) = @_;
    my (%occ, $ch, $cnt, $call, $cov, $rate);

    $cov = length($cross);
    
    ## disregard read strand
    $cross =/sanger tr/a-z,/A-Z./;

    ## count char occurences
    %occ = ();
    $occ{"."} = 0; $occ{"T"} = 0; $occ{"A"} = 0;
    for (my $i = 0; $i < length($cross); $i++){
	$ch = substr($cross, $i, 1);
	$occ{$ch} = 0 if !defined $occ{$ch};
	$occ{$ch}++;
    }

    ## majority voting
    $cnt = 0;
    $ch = "";
    foreach(keys(%occ)){
	## update most-frequently
	## occuring character
	if ($occ{$_} > $cnt){
	    $ch = $_;
	    $cnt = $occ{$_};
	}
	elsif ($occ{$_} == $cnt){
	    ## clear majority symbol
	    ## at draws except if only
	    ## methylated and unmethylated
	    ## bases are both in majority
	    ## ==> methylation rate of 0.5
	    if ($ch !/sanger /[.AT]/ || $_ !/sanger /[.AT]/){
		$ch = "";
	    }
	}
    }

    ## methylation calling
    ## C = methylated on strand
    ## T = unmethylated on strand
    ## V = variant at this position
    ## N = not covered
    ## A = ambigious calling
    
    if ($ch =/sanger /[.]/ || ($STRAND == 1 && $ch =/sanger /T/) || 
	($STRAND == 2 && $ch =/sanger /A/)){	
	$cov = $occ{"."};
	$cov += $occ{"T"} if $STRAND == 1;
	$cov += $occ{"A"} if $STRAND == 2;
	$call = ($occ{"."}/$cov >= 0.5) ? "C" : "T";
	$rate = sprintf "%.2f", $occ{"."}/$cov;
    }
    elsif ($ch eq ""){
	if ($cross eq ""){
	    $call = "N";
	    $rate = 0;
	}
	else {
	    $call = "A";
	    $rate = 0;
	}	    
    }
    else {
	$call = "V";
	$rate = 0;
    }
    
    return ($call, $cov, $rate);
}

sub getReducedCross {
    my ($cross) = @_;
    ## exclude read start (and additional char
    ## for mapping qual) and read end
    ## e.g. "^/sanger" ==> "" and "$" ==> ""
    $cross =/sanger s/\^.//g;
    $cross =/sanger s/\$//g;

    ## exclude indels
    ## e.g. "[+-]2AT" ==> ""
    $cross =/sanger s/[+-](\d+)([a-zA-Z]+)/substr($2,$1)/eg;
    
    ## exclude potential reference skips
    ## or deleted bases
    $cross =/sanger s/\*//g;
    $cross =/sanger s/[<>]//g;

    return $cross;
}

main();
