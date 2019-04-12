#######################################################################################
#	AUTHOR 	: Vaibhavi Pathak
#	PURPOSE	: excel feature of vlookup using perl
#	USAGE 	: perl perl_vlookup.pl file1 file2
#	OUTPUT 	: file1_combined.txt
#######################################################################################	

#!/usr/bin/perl -w 
use strict;

my %hash;
my $num_args = $#ARGV + 1;
if ($num_args != 2) {
    print "\nUsage: perl_vlookup.pl vcffile annovarfile\n";
    exit;
}
my $inputfile1 = $ARGV[0];
my $inputfile2 = $ARGV[1];

open(FILE,$inputfile1) or die "cant open file\n";
my @vcf = <FILE>;
close(FILE);

open(FILE,$inputfile2) or die "cant open file\n";
my @annovar=<FILE>;
close(FILE);
open(OUTFILE,">$inputfile1"."_combined.txt") or die "cant open file\n";

foreach my $line (@vcf)
{
	chomp($line);
	if($line !~ /^#/)
	{
		my @data = split("\t",$line);
		my ($chr,$start,$rd) = @data[0,1,9];
		my $match_str = $chr."_".$start;
		$hash{$match_str}=$rd;
	}
}
	foreach my $line1 (@annovar)
	{
		chomp($line1);
		my @data1 = split("\t",$line1);
		my ($chr1,$start1) = @data1[1,2];
		my $match_str1 = $chr1."_".$start1;
		if(defined($hash{$match_str1}))
		{	
			my @rd_str = split(":",$hash{$match_str1});
			my ($geno,$AD,$DP) = @rd_str[0,1,2]; 
			print OUTFILE "$line1\t$geno\t$AD\t$DP\n";
			next;
		}
	}
close(OUTFILE);	

