#####################################################
#Author: Vaibhavi Niraj Pathak
#Usage : perl grouping_and_matrix.pl
#input files : a intersect file of SNP, gene and QTL information.
#output : a Matrix representing unique SNPs and gene count with respect to each QTL and each chromosome 
###########################################################################################


#!/usr/bin/perl -w
use strict;
use Cwd qw(cwd);

my $dir = $ENV{'PWD'};
open(FILE,$dir."/SNP_gene_QTL_intersect_50k.bed") or die "cant open input file\n";
my @file = <FILE>;
close(FILE);
open(OUTFILE,">temp_output.txt") or die "cant open output file";
my %SNP_hash;
my %Gene_hash;
my %found;
my @QTL;
my @CHR;
my %seen;
foreach my $line (@file)
{
	chomp($line);
	my @data = split('\t',$line);
	my $chr = $data[0];
	my $SNP = $data[1];
	my $gene = $data[3];
	my $QTL = $data[9];

	push(@QTL,$QTL);
	push(@CHR,$chr);

	 $SNP_hash{$chr}{$QTL}{$SNP};
	if(exists($SNP_hash{$chr}{$QTL}{$SNP}))
	{
		$SNP_hash{$chr}{$QTL}{$SNP}+=1;	
	}	
	else
	{
		$SNP_hash{$chr}{$QTL}{$SNP}=1;
	}
	 $Gene_hash{$chr}{$QTL}{$gene};
	if(exists($Gene_hash{$chr}{$QTL}{$gene}))
	{
	
		$Gene_hash{$chr}{$QTL}{$gene}+=1;
	}
	else
	{
		$Gene_hash{$chr}{$QTL}{$gene}=1;
	}
}
my @unique_QTL = grep { !$seen{$_}++} @QTL;
my @unique_chr = grep { !$seen{$_}++} @CHR;
my $LINE = join ("\t\t", @unique_QTL);
print OUTFILE "\t"."\t".$LINE."\n";
my $len_QTL=$#unique_QTL;
print OUTFILE "Chromosome\t";
for(my $i=0;$i<=$len_QTL;$i++)
{
	print OUTFILE "SNP_count\tGene_count\t";
}
print OUTFILE "\n";
foreach my $C (@unique_chr)
{
	print OUTFILE "$C\t";

	foreach my $Q (@unique_QTL)
	{
		my $SNP_count=0;my $GENE_count=0;		
			foreach my $SNP (keys %{$SNP_hash{$C}{$Q} })
			{
				$SNP_count++;
			}
			foreach my $GENE (keys %{$Gene_hash{$C}{$Q} } )
			{
				$GENE_count++;							
			}
			print OUTFILE "$SNP_count\t$GENE_count\t";		
	}
	print OUTFILE "\n";
}




