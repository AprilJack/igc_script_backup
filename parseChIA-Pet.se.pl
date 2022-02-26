#!/usr/bin/perl -w

#print matches
my $mode='match';
my $minLen = 15;

if (@ARGV < 1) {
	print STDERR "<fq 1>\n";
	exit;
}

my %seq = ();

open $file, $ARGV[0] or die "Could not open $ARGV[0]\n";
my $T = 0;
my $t = 0;
my $tshort = 0;
my $tmis = 0;
my $tnocode = 0;
my $L1 = 0;

open OUT, ">" . $ARGV[0] . ".R1.match.fq";
open OUT2, ">" . $ARGV[0] . ".R2.match.fq";
open OUT3, ">" . $ARGV[0] . ".R1.mismatch.fq";
open OUT4, ">" . $ARGV[0] . ".R2.mismatch.fq";
($name,$fq1,$fq2,$code1,$code2,$L1,$L2) = getReadSE($file);
while ($name ne '') {
	$T++;
	if ($code1 eq '' || $code2 eq '') {
		$tnocode++;
		($name,$fq1,$fq2,$code1,$code2,$L1,$L2) = getReadSE($file);
		next;
	}
	if ($L1 < $minLen || $L2 < $minLen) {
		$tshort++;
		($name,$fq1,$fq2,$code1,$code2,$L1,$L2) = getReadSE($file);
		next;
	}
	if ($code1 eq $code2) {
		$t++;
		print OUT $fq1;
		print OUT2 $fq2;
	} else {
		$tmis++;
		print OUT3 $fq1;
		print OUT4 $fq2;
	}
	($name,$fq1,$fq2,$code1,$code2,$L1,$L2) = getReadSE($file);
}
close $file;
print STDERR "\tTotal=$T:\n";
my $r = sprintf("%.3f",$t/$T);
print STDERR "\t\t$t matching ($r)\n";
$r = sprintf("%.3f",$tmis/$T);
print STDERR "\t\t$tmis not matching ($r)\n";
$r = sprintf("%.3f",$tshort/$T);
print STDERR "\t\t$tshort with short sequences ($r)\n";
$r = sprintf("%.3f",$tnocode/$T);
print STDERR "\t\t$tnocode don't have correct adapter sequence ($r)\n";

close OUT;
close OUT2;

exit;

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

sub getReadSE {
	my ($f) = @_;
	my $fq1 = '';
	my $fq2 = '';
	my $name = '';
	my $code1 = '';
	my $code2 = '';
	my $seq1 = '';
	my $seq2 = '';
	my $c = 0;
	my $L1 = 0;
	my $L2 = 0;
	while (<$file>) {
		$c++;
		#print STDERR "$c $_";
		if ($c==1) {
			$fq1 .= $_;
			$fq2 .= $_;
			$name = $_;
			chomp $name;
			$name =/sanger s/ .*$//;
		} elsif ($c == 2) {
			chomp;
			my $x = $_;
			if ($x =/sanger /^(.*).GTTGGA(....)ATATCGCGGCCGCGATAT(....)TCCAAC.(.*)/) {
				#print STDERR "Good = $_";
				$seq1 = $1;
				$code1 = $2;
				$code2 = revopp($3);
				$seq2 = revopp($4);
				$fq1 .= $seq1 . "\n";
				$fq2 .= $seq2 . "\n";
			} else {
		#		print STDERR "bad $c = $_\n";
				$fq1 .= $_;
				$fq2 .= $_;
			}
		} elsif ($c == 3) {
			$fq1 .= $_;
			$fq2 .= $_;
		} elsif ($c == 4) {
			if ($seq1 ne '') {
				$L1 = length($seq1);
				$L2 = length($seq2);
				$fq1 .= substr($_,0,$L1) . "\n";
				$fq2 .= reverse(substr($_,$L1+6+4+18+4+6,$L2)) . "\n";
			} else {
				$fq1 .= $_;
				$fq2 .= $_;
			}
		}
		last if ($c >= 4);
	}
	return ($name,$fq1, $fq2, $code1, $code2, $L1, $L2);
}
sub revopp {
	my ($s) = @_;
	$s = reverse($s);
	$s =/sanger s/A/X/g;
	$s =/sanger s/T/A/g;
	$s =/sanger s/X/T/g;
	$s =/sanger s/C/X/g;
	$s =/sanger s/G/C/g;
	$s =/sanger s/X/G/g;
	return $s;
}
