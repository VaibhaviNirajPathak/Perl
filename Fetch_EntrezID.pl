####################################################################################################
# AUTHOR : Vaibhavi Niraj Pathak
# USAGE : perl Fetch_EntrezID.pl
# Requirements : HGNC_file, ANY manifest file in bed format(Chr,start,end,genename).
# it fetches EntrezID from HGNC file for the respective genename form Manifest file
#####################################################################################################


open(FILE,"HGNC_file.txt") or die "cant open file\n";
@HGNC=<FILE>;
close(FILE);

open(FILE1,"InhT1_CNV-manifest_12-may-16.bed") or die "cant open file\n";
@manifest=<FILE1>;
close(FILE);

open(OUTFILE,"InhT1_CNV-manifest_12-may-16_new.bed") or die "cant open\n";

foreach $line (@manifest)
{
	chomp($line);
	@data = split('\t',$line);
	foreach $line1 (@HGNC)
	{
		chomp($line1);
		@data1 = split('\t',$line1);
		if($data[3] eq $data1[1])
		{
			print OUTFILE "$data[0]\t$data[1]\t$data[2]\t$data[6]\t$data[7]\t$data[3]\t$data1[10]\n";
		}
	}
}
