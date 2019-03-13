#!/usr/bin/perl -w
#####################################################
#Author: Vaibhavi Niraj Pathak
#Usage : perl bin_known_snps.pl
#input files : directory of bams,panel runs list, manifest file
###########################################################################################


$dirin =  "/StrandxNAS/cases/exported_BAMs/2015";
open(OUTPUT1,">binned_SNPs.txt") or die "cant open output file1\n";
opendir(DIRIN,$dirin);
#@dirs = grep {-d "$dirin/$_" && ! /^\.{1,2}$/ && ! /^[U|V]/} readdir DIRIN;
open(FILE1,"clinical_exom_runs_list.txt") or die "cant open file\n";
@dirs = <FILE1>;
close(FILE1);

$manifest = "Trusight_one_Manifest_modified_zerobased_14_11_14_20bps_ext.bed";

%bin_hash=();%file_hash=();%tot_hash=();%mani_hash=();


foreach $dir (@dirs)
{
	chomp($dir);
	$dir =~ s/\r//g;
	opendir(DIR,"$dirin/$dir") or die "cant open directory\n";
	@files = grep { /\.vcf$/ && ! /_snp\.vcf|_CNV\.vcf|_SV\.vcf$/ } readdir DIR;
	closedir(DIR);
	foreach $file (@files)
	{
	chomp($file);
	open(FILE,"$dirin/$dir/$file") or die "cant open file\n";
	@F = <FILE>;
	close(FILE);
	@data=();@data1=();$SR='';$tot_snp=0;
	$counter10=0;$counter20=0;$counter30=0;$counter40=0;$counter50=0;$counter60=0;$counter70=0;$counter80=0;$counter90=0;$counter100=0;	
	$norm10=0;$norm20=0;$norm30=0;$norm40=0;$norm50=0;$norm60=0;$norm70=0;$norm80=0;$norm90=0;$norm100=0;
	@filtered_vcf=();$myfilter='';
	$myfilter = `intersectBed -a $dirin/$dir/$file -b $manifest`;
	@filtered_vcf = split('\n',$myfilter);
	foreach $line (@filtered_vcf)
	{
		chomp($line);
		if($line !~ /^[##|#]/g)
		{
			@data = split('\t',$line);
			$ref_length = length($data[3]);
			$alt_length = length($data[4]);
			if(($data[2] =~ /^rs(.*)$/) && ($ref_length==1) && ($alt_length==1))
			{
				@data1=split(':',$data[9]); 
				$SR = $data1[3];
					if($SR>0 && $SR<10)
					{
						$counter10++;					
					}
					elsif($SR>=10 && $SR<20)
					{
						$counter20++;
					}
					elsif($SR>=20 && $SR<30)
					{
						$counter30++;
					}
					elsif($SR>=30 && $SR<40)
					{
						$counter40++;
					}
					elsif($SR>=40 && $SR<50)
					{
						$counter50++;
					}
					elsif($SR>=50 && $SR<60)
					{
						$counter60++;
					}
					elsif($SR>=60 && $SR<70)
					{
						$counter70++;
					}
					elsif($SR>=70 && $SR<80)
					{
						$counter80++;
					}
					elsif($SR>=80 && $SR<90)
					{
						$counter90++;
					}
					elsif($SR>=90 && $SR<=100)
					{
						$counter100++;
					}
					else
					{
						#print "$dir\n";
					}
				}			
		}
	}
	$tot_snp = $counter10+$counter20+$counter30+$counter40+$counter50+$counter60+$counter70+$counter80+$counter90+$counter100;
	if($counter10 != '0'){	$norm10 = sprintf("%.2f", ($counter10/$tot_snp)*100); } else { $norm10=0; }
	if($counter20 != '0'){ $norm20 =  sprintf("%.2f", ($counter20/$tot_snp)*100); } else { $norm20=0; }
	if($counter30 != '0'){ $norm30 =  sprintf("%.2f", ($counter30/$tot_snp)*100); } else { $norm30=0; }
	if($counter40 != '0'){ $norm40 =  sprintf("%.2f", ($counter40/$tot_snp)*100); } else { $norm40=0; }
	if($counter50 != '0'){ $norm50 =  sprintf("%.2f", ($counter50/$tot_snp)*100); } else { $norm50=0; }
	if($counter60 != '0'){ $norm60 =  sprintf("%.2f", ($counter60/$tot_snp)*100); } else { $norm60=0; }
	if($counter70 != '0'){ $norm70 =  sprintf("%.2f", ($counter70/$tot_snp)*100); } else { $norm70=0; }
	if($counter80 != '0'){ $norm80 =  sprintf("%.2f", ($counter80/$tot_snp)*100); } else { $norm80=0; }
	if($counter90 != '0'){ $norm90 =  sprintf("%.2f", ($counter90/$tot_snp)*100); } else { $norm90=0; }
	if($counter100 != '0'){ $norm100 =  sprintf("%.2f", ($counter100/$tot_snp)*100); } else { $norm100=0; }
	$bin_hash{$file}{10}=$norm10;
	$bin_hash{$file}{20}=$norm20;
	$bin_hash{$file}{30}=$norm30;
	$bin_hash{$file}{40}=$norm40;
	$bin_hash{$file}{50}=$norm50;
	$bin_hash{$file}{60}=$norm60;
	$bin_hash{$file}{70}=$norm70;
	$bin_hash{$file}{80}=$norm80;
	$bin_hash{$file}{90}=$norm90;
	$bin_hash{$file}{100}=$norm100;
	@fl=(10,20,30,40,50,60,70,80,90,100);
	$tot_hash{$file}=$tot_snp;
	}
}
closedir(DIRIN);

	print OUTPUT1 join("\t", q{NAME}, @fl, q{total_SNP}), qq{\n};
	foreach my $k (sort keys %bin_hash) 
	{
   	 my @data_final;
   	 foreach my $l (@fl)
	 {
        	if (defined $bin_hash{$k}{$l})
	       	{
           	 push @data_final, $bin_hash{$k}{$l};
        	}
	       	else
	       	{
            	push @data_final, q{ };
        	}
   	}
	 if(exists $tot_hash{$k})
	 {
	 	push @data_final,$tot_hash{$k};
	 }
	print OUTPUT1 join("\t", $k, @data_final), qq{\n};
 	}
