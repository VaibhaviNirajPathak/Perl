#!/usr/bin/perl -w
use strict;
#######################################################################
### Author :  Vaibhavi Pathak
### PURPOSE: to take first 3 hit_def of blast output
### USAGE : perl read_blast_xml_reduce_multispecies.pl blast_output.xml
#######################################################################

my $xml = $ARGV[0];
my $len_str = length($xml);

my $t = substr($xml,0,$len_str-8);
print "$t\n";
my $out_filename = join("_",$t,"new");
print "$out_filename\n";
open (IN, $xml) || die "Cannot open blast output file\n";

open(OUT,">$out_filename.xml") or die "cant open output file\n";
my $new_line;
while (chomp(my $line = <IN>))
{
=if
	if($line =~ /<Iteration_query-def>/)
	{
		print "$line\t";
	}
=cut
	if($line =~ /<Hit_def>/)
	{
		if($line =~ /gt\;/)
		{
		my @data = split(';',$line);
		#print "$#data\n";
		if($#data >= 3)
		{
		$data[2] =~ /(.*)&gt/;
		my $trd_hit = $1;
		#print "$trd_hit\n";
		my $n_line = join(";",$data[0],$data[1],$trd_hit);
		$new_line = join("",$n_line,"</Hit_def>");
		}
		else
		{
			$new_line = $line;
		}
		#print "$line\n";
		}
		else
		{
			$new_line = $line;
		}
		print OUT "$new_line\n";
	}
	else
	{
		print OUT "$line\n";

	}
}
