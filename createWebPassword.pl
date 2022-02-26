#!/usr/bin/perl -w
#
#use POSIX;
use Cwd;
use Cwd 'abs_path';

#my $baseDir = "/home/illumina/fastq/";
#my $baseURL = "http://igc2.snl.salk.edu/illumina/runs/";

my $externalSymDir = "/genomics/www/ngs/data/";

my $baseDir = "/srv/www/igc1.salk.edu/public_html/";
my $baseURL = "http://igc1.salk.edu";

sub cmd {
	print STDERR "\n\tUsage: createWebPassword.pl [options] -u <username> -p <password>\n";
	print STDERR "\n\tRun this program in the directory that you want to be password protected.\n";
	print STDERR "\n\tExample: createWebPassword.pl -u GageLab -p Just4Rusty\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-u <username> (user name for access, required)\n";
	print STDERR "\t\t-p <password> (password for access, required)\n";
	print STDERR "\t\t-n <name at loginscreen> (default: protected data)\n";
	print STDERR "\t\t-d <directory> (Directory if different than current)\n";
	print STDERR "\t\t-baseUrl <url> (Base URL for output,\n";
	print STDERR "\t\t\t\tdefault: $baseURL)\n";
	print STDERR "\t\t-baseDir <directory> (Webserver root directory/subdirectory\n";
	print STDERR "\t\t\t\tdefault: $baseDir)\n";
	#print STDERR "\t\t-igc1 (reset base url/dir for igc1, only works in /genomics/www/ directory for outside access)\n";
	print STDERR "\t\t-external (setup access for external users)\n";
	print STDERR "\n";
	exit;
}

my $password = '';
my $username = '';
my $externalFlag = 0;
my $folderName = 'Protected_Data';
my $directory= getcwd();

if (@ARGV < 4) {
	cmd();
}


for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-u') {
		$username = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-p') {
		$password = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-d') {
		$directory = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-n') {
		$folderName = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-baseUrl') {
		$baseURL = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-external') {
		$externalFlag = 1;
	} elsif ($ARGV[$i] eq '-baseDir') {
		$baseDir = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-igc1') {
		$baseDir = "/genomics/www/";
		$baseURL = "http://homer.salk.edu/";
	} else {
		cmd();
	}
}

if ($password eq '' || $username eq '') {
	print STDERR "!!! Username (-u) and password (-p) are required!\n";
	exit;
}

my $absDir = abs_path($directory);


my $url = '';
if ($externalFlag) {
	my $symName = $absDir;
	$symName =~ s/^\/kinison//;
	$symName =~ s/\/genomics\/illumina\/fastq\///;
	$symName =~ s/\/boll\/cellar\/illumina\/fastq\///;
	$symName =~ s/\/home\/illumina\/fastq\///;
	$symName =~ s/\//\-/g;
	print STDERR "$directory -> $symName\n";
		
	my $sympath = $externalSymDir . $symName;
	$url = "http://homer.salk.edu/ngs/data/$symName/";

	print STDERR "ln -s $absDir $sympath\n";
	print STDERR "URL: $url\n";

	`ln -s $absDir $sympath`;
}
	

print STDERR "\tDirectory: $directory\n";	
print STDERR "\tAbs path:  $absDir\n";	


`htpasswd -bc $absDir/.htpasswd $username $password`;

open OUT, ">$absDir/.htaccess" or die "Could not open $absDir/.htaccess, probably don't have write access\n";
print OUT "AuthName \"$folderName\"\n";
print OUT "AuthType Basic\n";
print OUT "AuthUserFile $absDir/.htpasswd\n";
print OUT "require valid-user\n";
close OUT;

if ($externalFlag == 0) {
	my $dirDiff = $absDir;
	$dirDiff =~ s/^$baseDir//;
	$url = $baseURL. "/". $dirDiff;
}

print STDERR "\n\t$url\n";
print STDERR "\tuser: $username\n";
print STDERR "\tpass: $password\n";
print STDERR "\n";
exit;

