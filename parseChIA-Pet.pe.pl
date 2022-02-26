#!/usr/bin/perl -w

#print matches, mode=1 for good data (i.e. interactions that are between DNA from same sample)
##              mode=2 for bad data (i.e. interactions that are between different samples)
my $mode=1;
my $minLen = 15;

if (@ARGV < 2) {
	print STDERR "<fq 1> <fq 2>\n";
	exit;
}

my %seq = ();

open $file, $ARGV[0] or die "Could not open $ARGV[0]\n";
my $T = 0;
my $t = 0;
my $L = 0;
($fq,$name,$code,$L) = getRead($file);
while ($name ne '') {
	$T++;
	if ($code ne '' && $L >= $minLen) {
		$seq{$name} = {s=>$fq,c=>$code};
		$t++;
	}
	($fq,$name,$code,$L) = getRead($file);
}
close $file;
my $r = $t/$T;
print STDERR "\tTotal=$T, $t good($r)\n";


open OUT, ">" . $ARGV[0] . ".$mode.fq";
open OUT2, ">" . $ARGV[1] . ".$mode.fq";


open $file, $ARGV[1];
$t=0;
my $tt=0;
($fq,$name,$code,$L) = getRead($file);
while ($name ne '') {
	if ($code ne '' && $L >= $minLen) {
		$t++;
	}
	if (exists($seq{$name}) && $code ne '' && $L >= $minLen) {
		$tt++;
		if ($code eq $seq{$name}->{'c'} && $mode eq '1') {
			print OUT2 $fq;
			print OUT $seq{$name}->{'s'};
		}
		if ($code ne $seq{$name}->{'c'} && $mode eq '2') {
			print OUT2 $fq;
			print OUT $seq{$name}->{'s'};
		}
	}
	($fq,$name,$code,$L) = getRead($file);
}
close $file;
$r = $t/$T;
print STDERR "\tTotal=$T, $t good($r)\n";
$r = $tt/$T;
print STDERR "\tBoth $tt good($r)\n";


close OUT;
close OUT2;


sub getRead {
	my ($f) = @_;
	my $fq = '';
	my $name = '';
	my $code = '';
	my $seq = '';
	my $c = 0;
	my $L = 0;
	while (<$file>) {
		$c++;
		last if ($c >4);
		if ($c==1) {
			$fq .= $_;
			$name = $_;
			chomp $name;
			$name =/sanger s/ .*$//;
		} elsif ($c == 2) {
			my $x = $_;
			if ($x =/sanger /^(.*).GT[AT]GGA(....)ATATCGCGGCC/) {
				#print STDERR "Good = $_";
				$seq = $1;
				$code = $2;
				$fq .= $seq . "\n";
			} else {
				#print STDERR "bad  = $_";
				$fq .= $_;
			}
		} elsif ($c == 3) {
			$fq .= $_;
		} elsif ($c == 4) {
			if ($seq ne '') {
				$L = length($seq);
				$fq .= substr($_,0,$L) . "\n";
			} else {
				$fq .= $_;
			}
		}
	}
	return ($fq, $name, $code,$L);
}
