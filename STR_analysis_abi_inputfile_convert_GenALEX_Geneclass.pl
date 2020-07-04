=if
		AUTHOR : VAIBHAVI PATHAK
		USAGE : perl STR_bioinfo.pl STR_input_filename
		INPUT : ABI 3730XL DNA sequencer file, STR(Thermo fisher in this case) control file
		REMARKS : This file takes ABI 3730 XL DNA sequencer input file,normalizes it with respect to control, adjust the peak values into odd/even number and convert to GenALex and Geneclass2 fromat.
=cut

#!/usr/bin/perl -w
use strict;
use Cwd qw(cwd);
use POSIX;

my $num_args = $#ARGV + 1;
if ($num_args != 1) {
    print "\nUsage: STR_bioinfo.pl STR_input_filename\n";
    exit;
}
my $input_filename=$ARGV[0];
my %ctrl;
my %addnorm;
my %subnorm;
my $SAM={};
my $final_hash={};
my @Allele;
my @filtered_file;
my $Dir = $ENV{'PWD'};
my $ctrl_sam_name = "08_CORP_35";
my $infile4 = $Dir."/control.txt";
open(FILE,$infile4) or die "cant open file\n";
while(chomp(my $line4 =<FILE>))
{
	my($lo,$a1,$a2)  = split("\t",$line4);
	$ctrl{$lo}{'p1'} = $a1;
	$ctrl{$lo}{'p2'} = $a2; 
	push(@Allele,$lo);
}
open (FILE, $input_filename) || die "Cannot open the file2\n";
chomp(my $line2 = <FILE>);
while (chomp(my $line2 = <FILE>))
{
	my @data = split("\t",$line2);
	my $sam_name = $data[1];
	if (!exists($SAM->{$sam_name}))
	{
		if($sam_name ne $ctrl_sam_name)
		{ $SAM ->{$sam_name}++;}
	}
	if($sam_name eq $ctrl_sam_name)
	{
	$data[3] =~ /(.*)\((.*)\)/g;
	my $loc = $1;
	my $p1 = $data[4];
	my $p2 = $data[6];
	if($p2 == " ") {$p2 = $p1};
	if(($ctrl{$loc}{'p1'}) >= $p1)
		{	
			my $add = $ctrl{$loc}{'p1'} - $p1;
			$addnorm{$loc}{'p1'} = $add;
		} 
		else
		{	
			my $sub = $p1 - $ctrl{$loc}{'p1'};
			$subnorm{$loc}{'p1'} = $sub;
		}
		if(($ctrl{$loc}{'p2'}) >= $p2)
		{
			my $add = $ctrl{$loc}{'p2'} - $p2;
			$addnorm{$loc}{'p2'} = $add;
		} 
		else
		{
			my $sub = $p2 - $ctrl{$loc}{'p2'};
			$subnorm{$loc}{'p2'} = $sub;
		}
	}
	else
	{
		push(@filtered_file,$line2);
	}
} 
foreach my $line3 (@filtered_file)
{
	my @data1 = split("\t",$line3);
	my $samname = $data1[1];
	$data1[3] =~ /(.*)\((.*)\)/g;
	my $l= $1;
	my $peak1 = $data1[4];
	my $peak2 = $data1[6];
	if($peak2 == " ") {$peak2 = $peak1};
	if(exists($subnorm{$l}{'p1'}))
		{
			$peak1 = $peak1 - $subnorm{$l}{'p1'};
			print "$peak1\n";	
		}
		elsif(exists($addnorm{$l}{'p1'}))
		{
			$peak1 = $peak1 + $addnorm{$l}{'p1'};	}
		else
		{$peak1 = $peak1;}
		if(exists($subnorm{$l}{'p2'}))
		{
			$peak2 = $peak2 - $subnorm{$l}{'p2'};}
		elsif(exists($addnorm{$l}{'p2'}))
		{
			$peak2 = $peak2 + $addnorm{$l}{'p2'};	}
		else {$peak2 = $peak2}	

		if(($ctrl{$l}{'p1'}) % 2 == "0")
		{	
			if($peak1 % 2 == "0"){$peak1=$peak1;}
			else{$peak1=$peak1+1;}
		}
		else
		{	
			if($peak1 % 2 == "0"){$peak1=$peak1+1;}
			else{$peak1=$peak1;}
		}				
		if(($ctrl{$l}{'p2'}) % 2 == "0")
		{	
			if($peak2 % 2 == "0"){$peak2=$peak2;}
			else{$peak2=$peak2+1;}
		}
		else
		{
			if($peak2 % 2 == "0"){$peak2=$peak2+1;}
			else{$peak2=$peak2;}
		}
		print "$samname\t$l\t$peak1\t$peak2\n";
		$final_hash->{$samname}->{$l}->{'1'}= $peak1;
		$final_hash->{$samname}->{$l}->{'2'}=$peak2;
}


open (OUT, ">GenALex_format_output.txt") || die;
my @HEADER = ();
foreach my $Allele(@Allele)
{
	push (@HEADER, $Allele, $Allele."_2");
}
my $LINE = join ("\t", @HEADER);
print OUT "Sample"."\t".$LINE."\n";
foreach my $k (sort keys %{$SAM})
{
	print OUT $k;
	foreach my $A (@Allele)
	{
		if (exists($final_hash->{$k}->{$A}))
		{
	print OUT "\t".floor($final_hash->{$k}->{$A}->{'1'})."\t".floor($final_hash->{$k}->{$A}->{'2'});
		}
		else {
			print OUT "\t"."0"."\t"."0";
			}
	}
	print OUT "\n";
}

open(OUT1,">Geneclass2_format_output.txt") || die;
print OUT1 "Title line:\"DATA\"\n";
foreach my $Allele (@Allele)
{
	print OUT1 "$Allele\n";
}
my $i=1;
foreach my $k1 (sort keys %{$SAM})
{	
	print OUT1 "Pop\n p$i, ";
	foreach my $A1 (@Allele)
	{
	print OUT1 floor($final_hash->{$k1}->{$A1}->{'1'}).floor($final_hash->{$k1}->{$A1}->{'2'})." " ;
	}
	print OUT1 "\n";
	$i++;
}
