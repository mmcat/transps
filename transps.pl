#!/usr/bin/perl -w 

############################################
# Distributed under the GNC Public License #
# see the accompanying file LICENSE.txt    #
############################################

use strict;
use Getopt::Long;
use String::LCSS_XS qw(lcss lcss_all);
use List::Util qw(min max);
use POSIX qw(ceil floor);

our ($transcriptome,  $blastFile,%parseResults,$mapFile);
our ($overlappercent,$extentrate,$contigscafdist,$evalue);


processOpt();
parseBlast();
getMap();
scaffold();

sub processOpt{
	GetOptions('t=s' => \$transcriptome, 
#		'r=s' => \$ref,
		'b=s' => \$blastFile,
		'per=f' => \$overlappercent,
		'rate=f' => \$extentrate,
		'dist=i'=> \$contigscafdist,
		'evalue=f' => \$evalue) or usage("Invalid commmand line options.");
	usage("The transcriptome must be specified!\n\n") unless defined $transcriptome;
#	usage("The reference must be specified!\n") unless defined $ref;
	usage("The blastFile must be specified!\n\n") unless defined $blastFile;

	if(!defined $overlappercent){
		$overlappercent = 0.4;
	}
	if(!defined $extentrate){
		$extentrate = 0.1;
	}
	if(!defined $contigscafdist){
		$contigscafdist = 100;
	}
	if(!defined $evalue){
		$evalue = 1e-05;
	}
#$evalue = 1e-05;
}

sub usage{
	my $message = $_[0];
	print STDERR (
		"\n",
		$message,
#		"Usage: perl transScanffold.pl\n-t <transcriptome file>\n-r <reference file>\n-b <blastx output file>\n"
		"Usage: ./transps.pl -t <transcriptome file> -b <blastx output file> [options]\n
Options:
--per <percentage of overlapping>
--rate <rate of extention>
--dist <distance between two scaffoding contigs>
--evalue <evalue cutoff>\n");

	die("\n");

}

sub parseBlast{
#	$parseOut = $blastFile.".parse";
	open(IN,$blastFile) or die("Cannot open file for reading $blastFile:$!\n");
#	my %results;
	while(<IN>){
		chomp;
		my @line = split('\t',$_);
		my $qseqid = $line[0];
		if(!exists($parseResults{$qseqid})){

			$parseResults{$qseqid} = {'qseqid' => $line[0],
				'sseqid' => $line[1],
				'pident' => $line[2],
				'length' => $line[3],
				'mismatch' => $line[4],
				'gapopen' => $line[5],
				'qstart' => $line[6],
				'qend' => $line[7],
				'sstart' => $line[8],
				'send' => $line[9],
				'evalue' => $line[10],
				'bitscore' => $line[11]
			};

		}
		elsif($parseResults{$qseqid}->{'bitscore'}<$line[11] || 
			($parseResults{$qseqid}->{'bitscore'}==$line[11] && $parseResults{$qseqid}->{'evalue'}>$line[10])){

			$parseResults{$qseqid} = {'qseqid' => $line[0],
				'sseqid' => $line[1],
				'pident' => $line[2],
				'length' => $line[3],
				'mismatch' => $line[4],
				'gapopen' => $line[5],
				'qstart' => $line[6],
				'qend' => $line[7],
				'sstart' => $line[8],
				'send' => $line[9],
				'evalue' => $line[10],
				'bitscore' => $line[11]
			};

		}
	}
	close(IN);
=pod
	open(OUT,"> $parseOut") or die("Cannot open parseOut for writing $parseOut:$!\n");
	foreach (keys %results){
		print OUT ">",join("\t",$results{$_}->{'qseqid'},$results{$_}->{'sseqid'},"Score=".$results{$_}->{'bitscore'}),"\n";
		print OUT join("\t",$results{$_}->{'qstart'},$results{$_}->{'qend'},$results{$_}->{'sstart'},$results{$_}->{'send'},$results{$_}->{'pident'},$results{$_}->{'length'},$results{$_}->{'evalue'},$results{$_}->{'mismatch'},$results{$_}->{'gapopen'}),"\n";
	}
	close($parseOut);
=cut

}

sub getMap{
	$mapFile = $transcriptome.".map";
        my %pep;
	my ($pkey,$tkey,$tmapref);
	foreach (keys %parseResults){
		
		$pkey = $parseResults{$_}->{'sseqid'};
		$tkey = $parseResults{$_}->{'qseqid'};

		my $info = join("\t",$parseResults{$_}->{'qstart'},$parseResults{$_}->{'qend'},$parseResults{$_}->{'sstart'},$parseResults{$_}->{'send'},$parseResults{$_}->{'pident'},$parseResults{$_}->{'length'},$parseResults{$_}->{'evalue'},$parseResults{$_}->{'mismatch'},$parseResults{$_}->{'gapopen'},$parseResults{$_}->{'bitscore'});

		if(!exists $pep{$pkey}){
			$pep{$pkey}={$tkey => $info};
		}
		else{
			$pep{$pkey}->{$tkey} = $info;
		}
	}

	open(OUT,"> $mapFile") or die("Cannot open file for writing $mapFile:$!\n");

	foreach (keys %pep){
		my %tmap = %{$pep{$_}};
		print OUT ">$_\n";
		my $tkey;
		foreach $tkey (keys %tmap){
			print OUT $tkey,"\t",$tmap{$tkey},"\n";
		}
		print OUT "--\n";

	}
	close(OUT);	

}

sub readSeqFromFasta{
	my $file = $_[0];
	my (%map,$id,$seq);
	open(IN,$file) or die("Cannot open file for reading $file:$!\n");
	while(<IN>){
		chomp;
		if($_ =~ /^>/){
			if(defined $id && defined $seq){
				$map{$id} = $seq;
				$seq="";
			}
			($id) = split(' ',$_);
			$id=~s/>//;
		}
		else {
			$seq=$seq.$_;
		}
	}
	$map{$id} = $seq;
	close(IN);
	return \%map;

}

sub printSeq{
	my ($fh,$seq) = @_;
	my @chunks = split /(.{80})/, $seq;
	foreach (@chunks) {
		if($_ ne ""){
			print $fh "$_\n";
		}
	  }
}

sub scaffold{
	my ($SCA,$UNUSE,$ACC);
	open($SCA,"> $transcriptome.scaffold") or die("Cannot open file for writting:$transcriptome.scaffold:$!\n");
	open($UNUSE,"> $transcriptome.unused") or die("Cannot open file for writting:$transcriptome.unused:$!\n");
	open($ACC,"> $transcriptome.accept") or die ("Cannot open file for writing:$transcriptome.accepted:$!\n");

	my($pmapref,$tmapref,%pmap,%tmap);
#	$pmapref = readSeqFromFasta($ref);
	$tmapref = readSeqFromFasta($transcriptome);

#	%pmap = %{$pmapref};
	%tmap = %{$tmapref};

	open(IN,$mapFile) or die("Cannot open file for reading $mapFile:$!\n");

	my @match = <IN>;
	close(IN);
	my $i = 0;
	my ($pepname,@intervals);
	my %nullhash=("hitleft" => -1);

#parse each 1-to-1 or 1-to-many relationships
	while($i<@match){
		chomp $match[$i];
        	unshift(@intervals,\%nullhash);
		if($match[$i]=~/^>/){
			my $pepname=substr($match[$i],1);
			
			print "scanffolding contigs matching $pepname...\n";
			$i=$i+1;
#put each albopictus sequences (object) to a hash
		while($match[$i]!~/--/){
			chomp $match[$i];
			my @array=split('\t',$match[$i]);
#			my $qleft = min($array[1],$array[2]);
#			my $qright = max($array[1],$array[2]); 
			my %segment=(
					"name"=>$array[0],
					"queryleft"=>$array[1],
					"queryright"=>$array[2],
					"hitleft"=>$array[3],
					"hitright"=>$array[4],
					"perid"=>$array[5],
					"evalue"=>$array[7],
					"score"=>$array[10],
					"seq"=>$tmap{$array[0]});
			if($segment{queryleft} !~ /^[0-9]/ || $segment{evalue} > $evalue || $segment{queryleft}>$segment{queryright}){
#				if($segment{name} eq ""){
#					print "check6!\n";
#				}
#				print OUT3 ">",$segment{name},"|",$segment{queryleft},"\n";
                                print $UNUSE ">", $segment{name},"@",$pepname,"\n";
				printSeq($UNUSE,$segment{seq});
#				print $UNUSE $segment{seq},"\n";
				$i=$i+1;
				next;
			}
			
			
			my $j=1;
			my $n=scalar(@intervals);
# put the reference of each hash to an array as a group that match with the same agypti sequence
			if($n==1){push(@intervals,\%segment);}
			else{
				while ($j<$n && $segment{hitleft} > $intervals[$j]->{hitleft}){  # sort the intervals according to the left coordinate.
					$intervals[$j-1]=$intervals[$j];
					$j=$j+1;
			}
			
			$intervals[$j-1]=\%segment;
			unshift(@intervals,\%nullhash);	}		
			$i=$i+1;

		}
	       $i=$i+1;
	       
	       shift(@intervals);

# Regarded the oocyte contig as accepted if it is a 1-to-1 relationship

		if(scalar(@intervals)==1){
			print $ACC ">",$intervals[0]->{name},"&",$intervals[0]->{queryleft},"..",$intervals[0]->{queryright},"@",$pepname,"\n";
			printSeq($ACC,$intervals[0]->{seq});
#			print $ACC $intervals[0]->{seq},"\n";
#			print OUT2 ">",$pepname,"\n";
#			print OUT2 $map2{$pepname},"\n";
		}
# Filter some oocyte contigs if they meet the constraints(If the overlapping rate of two contigs is greater than 0.4, remove the contig with a lower score to 'unused')
		else{
#			sort { $intervals[$a]->{hitleft} <=> $intervals[$b]->{hitleft} } @intervals;
#			my $pname=substr($pepname,1);
			my $pname = $pepname;
			for (my $j=0; $j<@intervals-1; $j++){
				my $len1=$intervals[$j]->{hitright}-$intervals[$j]->{hitleft}+1;
				my $len2=$intervals[$j+1]->{hitright}-$intervals[$j+1]->{hitleft}+1;
				my $len;
				$len = min($len1,$len2);
				my $overlap_len=0;
				if($intervals[$j]->{hitright}<=$intervals[$j+1]->{hitright}){
					$overlap_len=$intervals[$j]->{hitright}-$intervals[$j+1]->{hitleft}+1;
				}
				else {$overlap_len=$intervals[$j+1]->{hitright}-$intervals[$j+1]->{hitleft}+1;}
				my $overlap_per=$overlap_len/$len;
				if($overlap_per>=$overlappercent){
					if ($intervals[$j]->{score}>=$intervals[$j+1]->{score}){
#						if($intervals[$j+1]->{name} eq ""){
#							print "check2\n";
#						}
						 
						 print $UNUSE ">", $intervals[$j+1]->{name},"@",$pname,"\n";
						 printSeq($UNUSE,$intervals[$j+1]->{seq});
#						 print $UNUSE $intervals[$j+1]->{seq},"\n";
                                                 splice(@intervals,$j+1,1);
						 $j=$j-1;
					}
					else {
#						if($intervals[$j]->{name} eq ""){
#							print "check1\n";
#						}
						
						print $UNUSE ">",$intervals[$j]->{name},"@",$pname,"\n";
						printSeq($UNUSE,$intervals[$j]->{seq});
#						print $UNUSE $intervals[$j]->{seq},"\n";	
                                                splice(@intervals,$j,1);					
						$j=$j-1;}

				}
			}
# Put the contig to 'accepted' if only one contig left after filtering
			if(scalar(@intervals)==1){
				         print $ACC ">",$intervals[0]->{name},"&",$intervals[0]->{queryleft},"..",$intervals[0]->{queryright},"@",$pname,"\n";
					 printSeq($ACC,$intervals[0]->{seq});
#					 print $ACC $intervals[0]->{seq},"\n";
					 $#intervals=-1;
					 next;
                                       }
# Begin to scaffolding
			my $cur_interval=$intervals[0];
		
			my $cur_right=$cur_interval->{hitright};
			my $flag=0;
			
			for ( my $j=1;$j<@intervals; $j++){
				if(length($intervals[$j])==0){
					print "WARNING:Cannot find sequence for $intervals[$j]->{name}\n";
				}
				my $sequence;
				my %merged_seg;
				my $scafdist = $contigscafdist;

				if($cur_right-$intervals[$j]->{hitleft}<=$scafdist && $cur_right <= $intervals[$j]->{hitright}){
					my $scaf_n=(abs($cur_right-$intervals[$j]->{hitleft})+1)*3;
					my $tmpstr="";
                                     
					my ($str1,$str2,$addn1,$addn2);
					$str1=substr($cur_interval->{seq},0,$cur_interval->{queryright});
					$str2=substr($intervals[$j]->{seq},$intervals[$j]->{queryleft}-1);
				        $addn1=0;$addn2=0;	
					my $diff=$cur_right-$intervals[$j]->{hitleft};
					if($diff>0){
						
					        my ($str11,$str22);
					        $str11 = substr($str1,-$scaf_n);
					        $str22 = substr($str2,0,$scaf_n);
						$tmpstr = jointStr1($str11,$str22);
						$sequence=substr($str1,0,-$scaf_n).$tmpstr.substr($str2,$scaf_n);
					}
					else{
						my $rate=$extentrate;
						my $n1=length($cur_interval->{seq})-$cur_interval->{queryright};
						my $n2=$intervals[$j]->{queryleft}-1;
						if(min($n1,$n2)<ceil($scaf_n*$rate) && (ceil($scaf_n*$rate)-min($n1,$n2))<min($cur_interval->{queryright}-1,length($intervals[$j]->{seq})-$intervals[$j]->{queryleft})){
							$addn1=ceil($scaf_n*$rate)-min($n1,$n2);
							$addn2=$addn1;
						}
						
						my($str11,$str22);
						
						$str11 = substr($cur_interval->{seq},$cur_interval->{queryright},$scaf_n+$addn1);
						$str11=substr($cur_interval->{seq},$cur_interval->{queryright}-$addn1,$addn1).$str11."N"x($scaf_n+$addn1-length($str11));
						my $off = $intervals[$j]->{queryleft}-$scaf_n-$addn2;
						if($off<=0){
							$off=0;$str22=substr($intervals[$j]->{seq},$off,$intervals[$j]->{queryleft}-1);
							my $gap="N"x($scaf_n+$addn2-length($str22));
							$str22=$gap.$str22.substr($intervals[$j]->{seq},$intervals[$j]->{queryleft}-1,$addn2);
						
						}
						else{
						$str22 = substr($intervals[$j]->{seq},$off-1,$scaf_n+$addn2);
						$str22 = $str22.substr($intervals[$j]->{seq},$intervals[$j]->{queryleft}-1,$addn2);
						}


						$tmpstr = jointStr2($str11,$str22,$addn1);
					        if($cur_interval->{evalue}<1e-05 && $intervals[$j]->{evalue}<1e-05) {
						           $sequence=substr($str1,0,length($str1)-$addn1).$tmpstr.substr($str2,$addn2);
					
					        }
					       else {
#						       if($intervals[$j]->{name} eq ""){
#							       print "check3\n";
#						       }
						       $sequence=$cur_interval->{seq};
						       print $UNUSE ">",$intervals[$j]->{name},"@",$pname,"\n";
						       printSeq($UNUSE,$intervals[$j]->{seq});
#						       print $UNUSE $intervals[$j]->{seq},"\n";
						       next;
					
					       }


					}
					

					my $tmp_query_right;

					if($diff>=0){

					        $tmp_query_right=$cur_interval->{queryright}+$intervals[$j]->{queryright}-$intervals[$j]->{queryleft}-length($tmpstr)+1;
					}
					else{
						$tmp_query_right=$cur_interval->{queryright}-$addn1+$intervals[$j]->{queryright}-$intervals[$j]->{queryleft}-$addn2+1+length($tmpstr);
					}
	
					%merged_seg = (
							"name"=> $cur_interval->{name}."&".$intervals[$j]->{name},
							"queryleft"=> $cur_interval->{queryleft},
							"queryright"=> $tmp_query_right,
							"hitleft"=> $cur_interval->{hitleft},
							"hitright"=> $intervals[$j]->{hitright},
							"perid"=> ($cur_interval->{perid}+$intervals[$j]->{perid})/2,
							"evalue"=> ($cur_interval->{evalue}+$intervals[$j]->{evalue})/2,
							"score"=> ($cur_interval->{score}+$intervals[$j]->{score})/2,
							"seq"=> $sequence
							);

					$cur_interval=\%merged_seg;
					$cur_right=$cur_interval->{hitright};
					$flag=1;

				}
				else{
#					if($intervals[$j]->{name} eq ""){
#						print "check4!\n";
#					}
				       print $UNUSE ">",$intervals[$j]->{name},"@",$pname,"\n";
				       printSeq($UNUSE,$intervals[$j]->{seq});
#				       print $UNUSE $intervals[$j]->{seq},"\n";	
				
				}				
				
				}#for
			if($flag==1){
				print $SCA ">",$cur_interval->{name},"&",$cur_interval->{queryleft},"..",$cur_interval->{queryright},"@", $pname,"\n";
				printSeq($SCA,$cur_interval->{seq});
#				print $SCA $cur_interval->{seq},"\n";
#				print OUT2 ">",$pepname,"\n";
#				print OUT2 $pmap{substr($pepname,1)}, "\n";
			}
			else {
				if(exists $cur_interval->{name} && $cur_interval->{name} ne ""){
					print $UNUSE ">",$cur_interval->{name},"@",$pname,"\n";
					printSeq($UNUSE,$cur_interval->{seq});
#					print $UNUSE $cur_interval->{seq},"\n";
				}
			}

			
			}#else
                       

		}#if($ref[$i]=~/>/)
		$#intervals=-1;        

	}#while($i<@ref)
	close($SCA);
	close($ACC);
	close($UNUSE);

}

sub jointStr1{

	my ($str1,$str2)=@_;
	my $lcs = lcss($str1,$str2);
        my $ps1 = index($str1,$lcs);
	my $pe1 = $ps1 +length($lcs)-1;
	my $ps2 = index($str2,$lcs);								
	my $pe2 = $ps2 +length($lcs)-1;
	my $flag=0;
	if($ps2>$ps1){
		$flag=1;
		my $str = $str1;
		$str1 = $str2;
		$str2 = $str; 	
		my $ps = $ps1;
		$ps1 = $ps2;
		$ps2 = $ps;
		my $pe = $pe1;
		$pe1 = $pe2;
		$pe2 = $pe;
	}
	my $tmpstr;
	if($flag==0){
	     
		$tmpstr=substr($str1,0,$pe1+1);
	
                $tmpstr = $tmpstr.substr($str2,$pe2+1);
	}
        else{
		$tmpstr=substr($str2,0,$pe2+1);
		$tmpstr=$tmpstr.substr($str1,$pe1+1);
	
	}	
	if(length($tmpstr)%3!=0){
		my $nn=3-length($tmpstr)%3;
		$tmpstr=$tmpstr."N"x$nn;
	}
	return $tmpstr;
}

sub jointStr2{

	my ($str1,$str2,$add)=@_;
	my $lcs = lcss($str1,$str2);
	if($lcs eq "N" x length($lcs)){
	$str1=substr($str1,0,length($str1)-$add);
	$str2=substr($str2,$add);
	$lcs=lcss($str1,$str2);
	}
        my $ps1 = index($str1,$lcs);
	my $pe1 = $ps1 +length($lcs)-1;
	my $ps2 = index($str2,$lcs);								
	my $pe2 = $ps2 +length($lcs)-1;
	my $flag=0;
	if($ps2>$ps1){
		my $str = $str1;
		$str1 = $str2;
		$str2 = $str; 	
		my $ps = $ps1;
		$ps1 = $ps2;
		$ps2 = $ps;
		my $pe = $pe1;
		$pe1 = $pe2;
		$pe2 = $pe;
		$flag=1;
	}
	my $sn = $ps1 - $ps2;
        my $tmpstr="";
        for(my $i=0;$i<$sn;$i++){
		if(substr($str1,$i,1) eq "N"){next;}
		elsif(substr($str2,0,1) eq "N" ||$flag==0) {$tmpstr=$tmpstr.substr($str1,$i,1);}
		else {next;}	
                
	}	
	for(my $i=$sn;$i<$ps1;$i++){
		if(substr($str1,$i,1) eq substr($str2,$i-$sn,1)){
			$tmpstr=$tmpstr.substr($str1,$i,1);
		}
		elsif(substr($str1,$i,1) eq "N" && substr($str2,$i-$sn,1) ne "N"){$tmpstr=$tmpstr.substr($str2,$i-$sn,1);}
		elsif(substr($str1,$i,1) ne "N" && substr($str2,$i-$sn,1) eq "N"){$tmpstr=$tmpstr.substr($str1,$i,1);}
		elsif($flag==0){
			$tmpstr=$tmpstr.substr($str1,$i,1);
		}
		else{$tmpstr=$tmpstr.substr($str2,$i-$sn,1);}
	}
	$tmpstr = $tmpstr.$lcs;
	for (my $i=$pe1+1;$i<length($str1);$i++){
		if($pe2+1<length($str2) && substr($str1,$i,1) eq substr($str2,$pe2+$i-$pe1,1)){
			$tmpstr=$tmpstr.substr($str1,$i,1);
		}
		elsif($pe2+1<length($str2) && substr($str1,$i,1) eq "N" && substr($str2,$pe2+$i-$pe1,1) ne "N"){$tmpstr=$tmpstr.substr($str2,$pe2+$i-$pe1,1);}
		elsif($pe2+1<length($str2) && substr($str1,$i,1) ne "N" && substr($str2,$pe2+$i-$pe1,1) eq "N"){$tmpstr=$tmpstr.substr($str1,$i,1);}
		elsif($flag==0){
			$tmpstr=$tmpstr.substr($str2,$pe2+$i-$pe1,1);
		
		}
		else{$tmpstr=$tmpstr.substr($str1,$i,1);}
        }
	for(my $i=$pe2+length($str1)-$pe1;$i<length($str2);$i++){
	        if(substr($str2,$i,1) eq "N"){next;}
		elsif(substr($str1,-1) eq "N" || $flag==0) {$tmpstr=$tmpstr.substr($str2,$i,1);}
		else{next;}
             	
        }
	if((length($tmpstr)-2*$add)%3!=0){
		my $nn = 3-(length($tmpstr)-2*$add)%3;
		$tmpstr=substr($tmpstr,0,length($tmpstr)-$add).("N"x$nn).substr($tmpstr,length($tmpstr)-$add);
	}
	
	return $tmpstr;

}

