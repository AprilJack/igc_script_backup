#!/usr/bin/perl -w
#

my $minDist = 0;
my $maxDist = 1e6;
my $inter=0;

my $color = "255,0,0";

foreach(@ARGV) {
	my $file = $_;
	open IN, $file;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		next if (@line < 10);
		my $c1 = $line[1];
		my $p1 = $line[2];
		my $d1 = $line[3];
		my $l1 = $line[5];
		my $c2 = $line[6];
		my $p2 = $line[7];
		my $d2 = $line[8];
		my $l2 = $line[9];
		if ($c1 eq $c2) {
			next if ($p2 <= $p1);
			my $diff = $p2-$p1;
			next if ($diff < $minDist || $diff > $maxDist);
			$p2+=$l2;
			next if ($p2 < $p1+$l1);
			my $name= "$c1:$p1-$p2";
			my $v = $p2-$p1-$l2;
			print "$c1\t$p1\t$p2\t$name\t1\t+\t$p1\t$p2\t$color\t2\t$l1,$l2,\t0,$v,\n";
		} else {
			next if ($inter ==0);
			my $v = $p2+$l2;
			my $name= "$c2:$p2-$v";
			$v = $p1+$l1;
			print "$c1\t$p1\t$v\t$name\t1\t+\t$p1\t$v\t$color\t1\t$l1,\t0\n";
		}
	}
	close IN;
}
