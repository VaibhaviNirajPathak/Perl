#!/usr/bin/perl -w
#####################################################
#Author: Vaibhavi Niraj Pathak
#Usage : perl QTLwise_SNP_genecount.pl
#input files : SNPs coordinates, Gene IDs/Genename and QTL name(a combined intersect file)
###########################################################################################

open(FILE,"SNP_gene_QTL_intersect_50k.bed") or die "cant open file\n";
@file=<FILE>;
close(FILE);

open(OUTFILE,">QTLwise_output_50k.txt") or die "cant open file";
my %SNP_hash;
my %Gene_hash;
my %found;

foreach my $line (@file)
{
	chomp($line);
	my @data = split('\t',$line);
	my $chr = $data[0];
	my $SNP = $data[1];
	my $gene = $data[3];
	my $QTL = $data[9];
	 $SNP_hash{$QTL}{$chr}{$SNP};
	if(exists($SNP_hash{$QTL}{$chr}{$SNP}))
	{
		$SNP_hash{$QTL}{$chr}{$SNP}+=1;	
	}	
	else
	{
		$SNP_hash{$QTL}{$chr}{$SNP}=1;
	}
	 $Gene_hash{$QTL}{$chr}{$gene};
	if(exists($Gene_hash{$QTL}{$chr}{$gene}))
	{
	
		$Gene_hash{$QTL}{$chr}{$gene}+=1;
	}
	else
	{
		$Gene_hash{$QTL}{$chr}{$gene}=1;
	}
}
	foreach $QTL (sort keys %SNP_hash)
	 {
		my%SNPcount;
		my $GENEcount;
		print OUTFILE "\n$QTL\n";
		foreach $chr (sort keys %{ $SNP_hash{$QTL} }) 
		{
			$SNP_count=0;$GENE_count=0;
			@SNP=();@GENE=();
			foreach $SNP (keys %{$SNP_hash{$QTL}{$chr} })
			{
				$SNP_count++;
				push(@SNP,$SNP);
					
			}
			foreach $GENE (keys %{$Gene_hash{$QTL}{$chr} } )
			{
				$GENE_count++;
				push(@GENE,$GENE);
				
			}
			$SNPcount{$chr}=$SNP_count;	
			$GENEcount{$chr}=$GENE_count;
			print OUTFILE "$chr\t$SNPcount{$chr}\t$GENEcount{$chr}\n";
		}
	}
	



