use List::Util qw( min max );

open(OUTPUT1,">colon_positive_sample_MSI_profile.txt") or die "cant open output file1\n";

$infile1 = "MSI_coordinates_mono_152_probe_footprints_new1.bed";

$dirout="intermediate_file";

$sample_file = "colon_positive_sample_list.txt";
open(FILE,"$sample_file") or die "cant open sample file\n";
@list = <FILE>;
close(FILE);

%final_hash=();%MS_hash=();
foreach $f_line (@list)
{
	chomp($f_line);
	$dirin='';$file='';
	($dirin,$file)=$f_line=~m/(.*)\/(.*)$/;
	chomp($dirin);
	chomp($file);
	$file=~s/\r//g;
	
	$infile2 = $file;
	$outfile1='';$outfile2='';$outfile3='';$outfile4='';$outfile5='';$outfile6='';$outfile7='';
	$outfile1 = $file."output.txt";
	#$outfile2 = $file."temp.txt";
	$outfile3 = $file."temp1.txt";
	#$outfile4 = $file."temp2.txt";
	$outfile5 =  $file."bam_filtered.txt";
	#$outfile6 = $file."t.txt";
	#$outfile7 = "printhash.txt";

	$outfile1=~s/\r//g;$outfile5=~s/\r//g;$outfile3=~s/\r//g;
	#$outfile2=~s/\r//g;$outfile4=~s/\r//g;$outfile6=~s/\r//g;$outfile7=~s/\r//g;

	open(OUTFILE1,">$dirout/$outfile1") or die "cant open outfile1\n";
	
	open(OUTFILE3,">$dirout/$outfile3") or die "cant open outfile3\n";
	#open(OUTFILE4,">$dirout/$outfile4") or die "cant open outfile4\n";
	open(OUTFILE5,">$dirout/$outfile5") or die "cant open outfile5\n";
	#open(OUTFILE6,">$dirout/$outfile6") or die "cant open outfile6\n";
	#open(OUTFILE7,">$dirout/$outfile7") or die "cant open outfile7\n";

	system("intersectBed -abam $dirin/$file -b $infile1 -bed > $dirout/$outfile1");

	$create_fasta = "samtools view $dirin/$file | awk '{OFS=\"\\t\"; print \$1\"\\t\"\$5\"\\t\"\$7\"\\t\"\$11}' -> $dirout/$outfile5";
	system($create_fasta);

 	 %myhash = ();
	%N_num = ();
  	$id = '';
  	$count=0;
    open F,"$dirout/$outfile5",or die $!;
    while(<F>)
    {
        chomp;
	$count++;
	($read_id,$read_start,$read_value) = split('\t', $_, 3);
	$myhash{$read_id}{$read_start} = $read_value;
    }
=if
	foreach $reaid (keys %myhash)
	{
		foreach $reastart (keys %{$myhash{$reaid}})
		{
			print OUTFILE7 "$reaid\t$reastart\t$myhash{$reaid}{$reastart}\n";
		}
	}
=cut
	open(FILE2,"$dirout/$outfile5") or die "cant open infile2\n";
	@infile2 = <FILE2>;
	close(FILE2);

	open(FILE1,"$infile1") or die "cant open infile1\n";
	@infile1 = <FILE1>;
	close(FILE1);

	open(FILE3,"$dirout/$outfile1") or die "cant open infile3\n";
	@infile3 = <FILE3>;
	close(FILE3);

	foreach $line1 (@infile1)
	{
		chomp($line1);
		@data1 = split('\t',$line1);
		$N_count=0;$unique_MSid=0;
	
		foreach $line2 (@infile3)
		{
			chomp($line2);
			@data2 = split('\t',$line2);
			$data2[3] =~ /(.+)\s(.+)/;
			$r_id = $1;
			if(($data1[0] eq $data2[0]) && ($data1[1] eq $data2[1]) && ($data1[2] eq $data2[2]))
			{
				$unique_MSid=$data1[0].":".$data1[1].":".$data1[2];
					
				if(($data2[1] >= ($data2[6]-3)) && ($data2[2] <= ($data2[7]+3)))
				{
				if(exists($myhash{$r_id}{($data2[6]+1)}))
				{
					print OUTFILE3 "$line2\t$myhash{$r_id}{($data2[6]+1)}\n";
				}		
				$N_count++;
				}
			}
		}
	
	if(!exists $N_num{$unique_MSid})
	{
			$N_num{$unique_MSid}=$N_count;
	}

       }

	foreach $id (keys %N_num)
	{
		chomp($id);
		print OUTFILE2 "$id\t$N_num{$id}\n";
		
	}

	open(FILE4,"$dirout/$outfile3") or die "cant open infile3\n";
	@infile4 = <FILE4>;
	close(FILE4);
	%new_hash =();	$value=0;

	foreach $line3 (@infile4)
	{
		chomp($line3);
		@data3 = split('\t',$line3);
		$start=$data3[6];
		$chr=$data3[0];

		$start_MS_inread = $data3[1]-$data3[6];
		$cigar=$data3[12];
		@cigar_len = split(/\D+/,$cigar);
		@cigar_char = split(/\d+/,$cigar);
		shift @cigar_char;
		$pos=0;
		if($cigar_char[0] eq 'S')
		{
			$start_MS_inread = $cigar_len[0] + $start_MS_inread;
		}

		for($i=0;$i<=$#cigar_len;$i++)
		{
			
			if( $cigar_char[$i] eq 'M')
			{
				$pos = $cigar_len[$i] + $pos;
			}
			elsif($cigar_char[$i] eq 'D')
			{
				$pos = $pos - $cigar_len[$i];
			}
			elsif($cigar_char[$i] eq 'I')
			{
				$pos = $pos + $cigar_len[$i];
			}
			if($pos >= $start_MS_inread)
			{
				$offset = ($start_MS_inread+1);
				last;
			}
		}
		

	$seq = $data3[13];
	$base =	substr($seq,$offset,1);
	$homo_len=1;
	for($i=1;$i<=length($seq);$i++)
	{
		
		$base2 = substr($seq,($offset+$i),1);
		if($base eq $base2)
		{
			$homo_len = $homo_len+1;
		}
		else
		{
			last;
		}
	}
	print OUTFILE4 "$line3\t$offset\t$start_MS_inread\t$pos\t$homo_len\t$cigar\n";

	$key1 = $data3[0].":".$data3[1].":".$data3[2];
	$key2=$homo_len;
	if(!exists $new_hash{$key1}{$key2})
	{
		$new_hash{$key1}{$key2}=1;
	}
	else
	{
		$new_hash{$key1}{$key2}++;
	}
      }

	foreach  $MS (sort keys %new_hash)
	{
	
		print OUTFILE6 "$MS: \n";
		$C=0;@sorted_array=();$NormC=0;$r_count=0;
	 	foreach $H_len (sort keys%{$new_hash{$MS}})
	 	{
		 	$C++;
			push(@sorted_array,$new_hash{$MS}{$H_len});
		
   	 	}
		$max = max(@sorted_array);
	
		for($i=0;$i<=$#sorted_array;$i++)
		{
			$NormC = ($sorted_array[$i]*100)/$max;
			if($NormC>5)
			{
				$r_count++;
			}
		}
		$k1=$file;
		if($N_num{$MS} <= '50')
		{
			#print "$MS,N-count=$N_num{$MS}\n";
			$r_count=0;
		}
		$final_hash{$k1}{$MS} = $r_count;
		$MS_hash{$MS}++;
 	}
	@fl=();$dirin='';$file='';
	
}

	my @sets = sort keys %MS_hash;
	print OUTPUT1 join("\t", q{NAME}, @sets), qq{\n};
	foreach my $k (sort keys %final_hash) 
	{
   	 my @data_final;
   	 foreach my $l (@sets)
	 {
        	if (defined $final_hash{$k}{$l})
	       	{
           	 push @data_final, $final_hash{$k}{$l};
        	}
	       	else
	       	{
            	push @data_final, q{ };
        	}
    	}
    print OUTPUT1 join("\t", $k, @data_final), qq{\n};
	}



