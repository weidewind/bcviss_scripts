#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use File::Find;
use File::Basename;
use File::Path 2.09 qw(make_path);
use File::Copy;
use XML::Simple;
use Class::Struct;
use File::Temp qw/ tempfile tempdir /;
use ABI;
#use Bio::SeqIO;
use BCVIndel::Project;

#######################
#These variables contain shell commands that must be set in the environment
my $TraceTunerCmd=$ENV{'TTUNERCMD'};
die "\nTTUNERCMD environment variable is not set!" unless defined($TraceTunerCmd);
my $PolyScanCMD=$ENV{'POLYSCANCMD'};
die "\nPOLYSCANCMD environment variable is not set!" unless defined($PolyScanCMD);
my $BCVCMD=$ENV{'BCVCMD'};
die "\nBCVCMD environment variable is not set!" unless defined($BCVCMD);
my $BcvHome=$ENV{'BCV_HOME'};
#######################

my %sources;
my $InFileExt=".ab1"; #default value

struct ReadInfo => {
	source => '$',
	bcv_cfg => '$',
	SampleName => '$',
	PrimerName => '$',
	Direction => '$',
	LeftTrim => '$',
	RightTrim => '$'
};

my %seen_path;

sub get_sources{
	if (-f && /$InFileExt$/) {
		print "found $File::Find::name\n";
		my $name=$_;
		$name=~s/$InFileExt$//;
		$name=~s/\s/_/g;
		my $p=ReadInfo->new();
		$p->source($File::Find::name);
		$sources{$name}=$p;
		++$seen_path{$File::Find::dir};
	}elsif (-d && $seen_path{$File::Find::dir}) {
		$File::Find::prune = 1;
	}
};

sub parse_read{
	my $p=shift;
	my $rprimer_directions=shift;
	my $delim=shift;
	$delim="_" unless defined $delim;
	return 1 if(defined($p->Direction()) &&  defined($p->SampleName()) && defined($p->PrimerName()));
	if(!defined $p->SampleName()){
		#my $in = Bio::SeqIO->new(-file => $p->source(), -format => 'abi');
		#while ( my $seq = $in->next_seq() ) {
		#	$p->SampleName($seq->display_id());
		#	last;
		#};
		my $in=ABI->new($p->source());
		my $tmp=$in->get_sample_name();
		$p->SampleName($tmp);
	};
	if(!defined $p->PrimerName()){
		return 0 if(!defined($p->SampleName()) && !defined($p->Direction()));
		my $tmp=$p->SampleName();
		my $b=undef;
		foreach my $prm(keys %{$rprimer_directions}){
			if($tmp=~m/$prm/i){
				$b=$rprimer_directions->{$prm};
				$tmp=~s/$prm//i;
				$tmp=~s/$delim{2}//;
				$tmp=~s/$delim\s*$//;
				$tmp=~s/^\s*$delim//;
				$p->PrimerName($prm);
				$p->SampleName($tmp);
				$p->Direction($b) unless defined $p->Direction();
				last;
			};
		};
	};
	return 0 unless defined $p->Direction();
	return 1;
};

sub parse_expfile{
	my $p=shift;
	my $file_name=shift;
	return 1 if(defined($p->LeftTrim)&&defined($p->RightTrim));
	open INF, "<$file_name" or return 0;
	my ($ql,$qr,$sl,$sr);
	my @cs;
	my $seq="";
	while(<INF>){
		if(/^QL\s+(\d+)/){
			#$ql=$1;
		}elsif(/^QR\s+(\d+)/){
			#$qr=$1;
		}elsif(/^SL\s+(\d+)/){
			$sl=$1;
			$p->LeftTrim($sl);
		}elsif(/^SR\s+(\d+)/){
			$sr=$1;
			$p->RightTrim($sr);
		}elsif(/^CS\s+(\d+)\.\.(\d+)/){
			@cs=($1,$2);
		}elsif(/^SQ\s*$/){
			while(<INF>){
				if(/\/\//){
					last;
				}else{
					$seq.=$_;
				}
			}
		};
	};
	close INF;
	$seq=~s/\s//g;
	$p->LeftTrim($ql<$sl?$sl+1:$ql+1) if(defined($ql)&&defined($sl));
	$p->RightTrim($qr<$sr?$qr-1:$sr-1) if(defined($qr)&&defined($sr));
	if(@cs==2){
		#remove cloning vector sequence
		my @I=($p->LeftTrim,$p->RightTrim);
		$I[0]=0 unless defined $I[0];
		$I[1]=length($seq) unless defined $I[1];
		my @J=($I[0]<$cs[0]?$cs[0]:$I[0],$I[1]<$cs[1]?$I[1]:$cs[1]);
		if($J[0]<$J[1]){
			my @A=($I[0],$J[0]-1);
			my @B=($J[1]+1,$I[1]);
			my $la=$A[1]-$A[0]+1;
			my $lb=$B[1]-$B[0]+1;
			@A=@B if $la<$lb;
			if($A[0]<$A[1]){
				$p->LeftTrim($A[0]);
				$p->RightTrim($A[1]);
			}
		}
	}	
	return defined($p->LeftTrim)||defined($p->RightTrim);
};


sub trim_bqs{
	my $p=shift;
	my $infile=shift;
	my $max_length= shift;
	my ($name,$dir,$ext) = fileparse($infile,'\.?[^\.]*$');
	my $outfile=$dir.$name.".trim.bqs";
	open OPF, ">$outfile" or die "\nUnable to open output file: $outfile!";
	open INF, "<$infile" or die "\nUnable to open input file: $infile!";
	while(<INF>){
		if(/^RD:/){
			print OPF;
			last;
		};
	};
	die "\nWrong BQS file $infile!" if eof(INF);
	my $ltrim=0;
	$ltrim=$p->LeftTrim if(defined($p->LeftTrim));
	while(<INF>){
		if(/^(\w)\t(\d+)/){
			my $sym=$1;
			my $pos=$2;
			my $str=$';
			last if(defined($p->RightTrim)&&($pos>=$p->RightTrim));
			if($pos>$ltrim){
				$pos-=$ltrim;
				print OPF "$sym\t$pos\t$str";
			};
			last if(defined($max_length)&&($pos==$max_length));
		}
	};
	close INF;
	close OPF;
	return $outfile;
};

sub trim_fpoly{
	my $p=shift;
	my $infile=shift;
	my $max_length= shift;
	my ($name,$dir,$ext) = fileparse($infile,'\.?[^\.]*$');
	my $outfile=$dir.$name.".trim.fpoly";
	open OPF, ">$outfile" or die "\nUnable to open output file: $outfile!";
	open INF, "<$infile" or die "\nUnable to open input file: $infile!";
	#skip first line
	while(<INF>){
		if(/[^\s]/){
			print OPF;
			last;
		};
	};
	my $ltrim=0;
	my $pos=1;
	$ltrim=$p->LeftTrim if(defined($p->LeftTrim));
	while(<INF>){
		if(/[^\s]/){
			last if(defined($p->RightTrim)&&($pos>=$p->RightTrim));
			if($pos>$ltrim){
				print OPF;
			};
			last if(defined($max_length)&&($pos-$ltrim==$max_length));
			$pos++;
		}
	};
	close INF;
	close OPF;
	return $outfile;
};


#Build BCV project that is a set of configuration (CFG, bcvindel,prj.xml), input (BQS, FPOLY) 
#and optionally output files files in specified foulder
#Args:
#project_file - name of a XML project file
#returns 0 if succeeded.
sub build_project{
	my $in_project_file=shift;
	print "\nBCV Project File. XML input: $in_project_file ";
	
	my ($name,$dir,$ext) = fileparse($in_project_file,'\.?[^\.]*$');
	my $LogFile=$dir.$name.".log";
	open LOGF, ">$LogFile" or die "\nUnable to open the file: $LogFile!";
	
	#read a project file
	my $bcv_prj = XML::Simple->new();
	
	my $data   = $bcv_prj->XMLin($in_project_file, 
		ForceArray => ['Flag','Read','Primer','PrimerDirection'],
		SuppressEmpty => 1,
		ContentKey => '-content' );
	#use Data::Dumper;
	#print Dumper($data);
	#exit;
	$InFileExt=$data->{InFileExt};
	$InFileExt=~s/^\s*//;
	#directory with source files
	my $InPath=$data->{InPath};
	#read ABIF files
	find (\&get_sources,$InPath);
	die "\nNo source files with type: $InFileExt in the foulder: $InPath or read permission error occured!" if scalar(keys %sources)==0;
	#read user specified read info
	foreach my $rname(keys %{$data->{READS}->{Read}}){
		my $name=$rname;
		$name=~s/\.?[^\.]*$//;
		$name=~s/\s/_/g;
		if(defined $sources{$name}){
			foreach my $rkey(keys %{$data->{READS}->{Read}->{$rname}}){
				if($rkey eq 'SampleName'){
					$sources{$name}->SampleName($data->{READS}->{Read}->{$rname}->{$rkey});
				}elsif($rkey eq 'Direction'){
					$_=$data->{READS}->{Read}->{$rname}->{$rkey};
					if(m/reverse/i){
						$sources{$name}->Direction(0);
					}else{
						$sources{$name}->Direction(1);
					};
				}elsif($rkey eq 'PrimerName'){
					$sources{$name}->PrimerName($data->{READS}->{Read}->{$rname}->{$rkey});
				}elsif($rkey eq 'LeftTrim'){
					$sources{$name}->LeftTrim($data->{READS}->{Read}->{$rname}->{$rkey});
				}elsif($rkey eq 'RightTrim'){
					$sources{$name}->RightTrim($data->{READS}->{Read}->{$rname}->{$rkey});
				}
			};
		}else{
			print STDERR "\nThe read: $rname is absent at the specified location: $InPath. Ignored!";
		};
	};
	
	#Get usecase
	my $UseCase=$data->{UseCase};
	#Possible use cases
	my @usecases=('basecall','indelcall','deconvolution'); # 'deconvolition' is default
	if(!defined $UseCase){
		$UseCase=2;
	}else{
		$_=$UseCase;
		if(m/deconvolution/i){
			$UseCase=2;
		}elsif(m/indelcall/i){
			$UseCase=1;
		}elsif(m/basecall/i){
			$UseCase=0;
		}else{
			die "\nI don't know the UseCase: $UseCase. Terminated!";
		};
	};

	

	#Get flags
	#set default values
	my %Flags=(
		'ParseSampleNames' => 1,
		'RunBCV' => 1);
	foreach my $flag(keys %Flags){
		my $f_lag=$data->{FLAGS}->{Flag}->{$flag};
		if(defined $f_lag){
			$Flags{$flag}=0 if(m/false/i);
		}; 
	};
	
	my %primer_directions=%{$data->{PRIMERS}->{PrimerDirection}} if $Flags{ParseSampleNames} && defined $data->{PRIMERS}->{PrimerDirection};
	foreach my $primer(keys %primer_directions){
		$_=$primer_directions{$primer};
		if(m/reverse/i){
			$primer_directions{$primer}=0;
		}else{
			$primer_directions{$primer}=1;
		};
	};
	my $primer_delim=$data->{PRIMERS}->{PrimerDelim};
	#read basic parameters:
	#BCV vocabulary in FASTA or GFAS format
	my $Vocabulary=$data->{Vocabulary};
	if(!defined $Vocabulary){
		print STDERR "\nWarning: Missing value for the Vocabulary! It will be taken from the BCV configuration file." unless $UseCase==0;
	}else{
		if(defined $BcvHome){
			$Vocabulary=~s/^\${BCV_HOME}/$BcvHome/;
		}
		open RFILE, "<$Vocabulary" or die "\nUnable to open the file: $Vocabulary! Terminated: $!";
		close RFILE;
	};
	#BCV configuration file with basic parameters
	my $BCV_CFG=$data->{BCV_CFG};
	if(!defined $BCV_CFG){
		die "\nMissing value for the BCV_CFG (configuration file for BCV program)!";
	}else{
		if(defined $BcvHome){
			$BCV_CFG=~s/^\${BCV_HOME}/$BcvHome/;
		}
		
	$BCV_CFG =~ s/\r//; #fisa
		open RFILE, "<$BCV_CFG" or die "\nUnable to open the file: $BCV_CFG! Terminated!";
		close RFILE;
	};
	#The directory where the input and configuration files for BCV program have to be stored
	my $BCV_Input=$data->{BCV_Input};
	if(!defined $BCV_Input){
		die "\nMissing value for the BCV_Input directory!";
	}else{
		$BCV_Input=~s/\/\s*$//;
		$BCV_Input.='/';
		die "\nThe directory doesn't exists: $BCV_Input!" unless -d $BCV_Input;
	};
	#The directory where the BCV program data will be outputed
	my $BCV_Output=$data->{BCV_Output};
	if(!defined $BCV_Output){
		die "\nMissing value for the BCV_Output directory!";
	}else{
		$BCV_Output=~s/\/\s*$//;
		$BCV_Output.='/';
		die "\nThe directory doesn't exists: $BCV_Output!" unless -d $BCV_Output;
	};
	#The directory where the staden::pregap output files located
	my $pregap4_dir=$data->{Pregap_Dir};
	if(defined $pregap4_dir){
		$pregap4_dir=~s/\/\s*$//;
		$pregap4_dir.='/';
		die "\nThe directory doesn't exists: $pregap4_dir!" unless -d $pregap4_dir;
	};
	my $max_seqlength=$data->{MaxSequenceLength};
	if(defined $max_seqlength){
		$max_seqlength=~s/[^\d]//;
		die "\nError in sequence length parameter: $data->{MaxSequenceLength}!" if $max_seqlength eq "";
	};

	#parse read info from sample names
	#Hash of read info objects for which the BCV run could be performed
	#Hash keys '$name' will be used for naming of BCV output foulders
	my %sources_ok;
	foreach my $name(keys %sources){
		print "\n$name";
		parse_read($sources{$name},\%primer_directions,$primer_delim);
		if(defined $sources{$name}->Direction()){
			$sources_ok{$name}=$sources{$name};
			if(defined $pregap4_dir){
				parse_expfile($sources_ok{$name},$pregap4_dir.$name.".exp");
			}
		};
	};
	my $nOK=keys %sources_ok;
	die "\nUnable to obtain the read directions from their source files!" unless $nOK;
	print LOGF  "[Configuration]\nTotal reads OK=$nOK\n";
	
	#building CFG files 	
	foreach my $src(keys %sources_ok){
		my $out=$BCV_Output.$src."/";
		my $in=$BCV_Input.$src."/";
		my @ps_prj=(
			$in."ps_prj/chromat_dir/",
			$in."ps_prj/edit_dir/",
			$in."ps_prj/phd_dir/",
			$in."ps_prj/poly_dir/"
	    );
		#build PolyScan project directories
		make_path(@ps_prj);
		#build BCV output directory
		make_path($out);
		#run the TT
		my $source_file=$sources_ok{$src}->source();
		#coping a source file
		my $ofile_name=$ps_prj[0].$src.$InFileExt;
		copy($source_file, $ofile_name);
		my $cmd=$TraceTunerCmd." -mix -pd $ps_prj[2] -cd $ps_prj[0] -dd $ps_prj[3] $ofile_name";
		$cmd =~ s/\r//; #fisa
		print "\n$cmd\n";
		print LOGF  "\n$cmd\n";
		system($cmd);
		if ($? == -1) {
			print "failed to execute: $!\n";
			exit;
		}
		elsif ($? & 127) {
			printf  LOGF  "\tchild died with signal %d, %s coredump\n",
				($? & 127),  ($? & 128) ? 'with' : 'without';
		}
		else {
			if($? >> 8!=0){
				printf LOGF  "\tchild exited with value %d\n", $? >> 8;
			};
		};
		#run the PS
		$cmd=$PolyScanCMD." -tQA -pd $in"."ps_prj/";
		$cmd =~ s/\r//; #fisa
		print "\n$cmd\n";
		print LOGF  "\n$cmd\n";
		system($cmd);
		if ($? == -1) {
			print "failed to execute: $!\n";
			exit;
		}
		elsif ($? & 127) {
			printf  LOGF  "\tchild died with signal %d, %s coredump\n",
				($? & 127),  ($? & 128) ? 'with' : 'without';
		}
		else {
				if($? >> 8!=0){
				printf LOGF  "\tchild exited with value %d\n", $? >> 8;
			};
		};
		#make CFG file
		my $fpoly=$ps_prj[3].$src.$InFileExt.".fpoly";
		my $bqs=$ps_prj[3].$src.$InFileExt.".bqs";
		if(defined($sources_ok{$src}->LeftTrim)||defined($sources_ok{$src}->RightTrim)||defined($max_seqlength)){
			$fpoly=trim_fpoly($sources_ok{$src},$fpoly,$max_seqlength);
			$bqs=trim_bqs($sources_ok{$src},$bqs,$max_seqlength);
		};
		$ofile_name=$BCV_Input.$src.$InFileExt.".cfg";
		open OFILE, ">$ofile_name" or die "\nUnable to open $ofile_name file. Terminated!";
		print OFILE  "#[Configuration]\n\n#configuration file\ncfg=\"$BCV_CFG\"\n\n";
		print OFILE  "#[Output]\n\n#output directory\nout=\"$out\"\n\n";
		if(defined $Vocabulary){
			print OFILE  "#[Dictionary_Alignment]\n\ndict=\"$Vocabulary\"\n\n";
		};
		
		print OFILE  "#[Chromatogram]\n\n#complementary strand\ncompl=";
		if($sources_ok{$src}->Direction()){
			print OFILE  "no";
		}else{
			print OFILE  "yes";
		};
		print OFILE  "\n\n#fpoly file path\nfpoly=\"$fpoly\"\n\n#bqs file path\nbqs=\"$bqs\"";
		if($UseCase==0){
			print OFILE  "\n\n#[Dictionary_Alignment]\n\n#output Viterbi graph?\nvit_graph=yes";
			print OFILE  "\n\n#[Sampling]\n\n#take samples? (no => uses Viterbi path as sample)\nsample=no";
			print OFILE  "\n\n#[Decomposition]\n\n#make decomposition?\ndecomp=no";
		}elsif($UseCase==1){
			print OFILE  "\n\n#[Dictionary_Alignment]\n\n#output Viterbi graph?\nvit_graph=yes";
			print OFILE  "\n\n#[Decomposition]\n\n#make decomposition?\ndecomp=yes";
			print OFILE  "\n\n#log decomposition?\nlog_decomp=yes";
			print OFILE  "\n\n#[Clustering]\n\n#do clustering\ncluster=no";
		}elsif($UseCase==2){
			print OFILE  "\n\n#[Decomposition]\n\n#make decomposition?\ndecomp=yes";
			print OFILE  "\n\n#log decomposition?\nlog_decomp=yes";
			print OFILE  "\n\n#[Clustering]\n\n#do clustering\ncluster=yes";
		};
		close OFILE;
		$sources_ok{$src}->bcv_cfg($ofile_name);
		print "\nBuilt: $ofile_name\n";
	};
  print "\nBuilding BCV configurations done!";
  if($Flags{RunBCV}){
  	#Temp file for redirecting the STDERR from BCV
  	my $template="bcv_stderr.tmpXXXXX"; 
		my ($fh, $filename) = tempfile($template, TMPDIR => 1);
		die "\nUnable to create tmp file $filename: $template!" unless defined $fh;
		close $fh;
  	print "\nRunning BCV main program:";
  	print LOGF  "\n[BCV]\n";
  	my $n=0;
  	foreach my $src(keys %sources_ok){
  		my $cfg=$sources_ok{$src}->bcv_cfg();
  		my $cmd=$BCVCMD." --cfg $cfg 2>&1 1>$filename";
		$cmd =~ s/\r//; #fisa
  		print "\n$cmd\n";
    	print LOGF  "\n$cmd\n";
			system($cmd);
			if ($? == -1) {
				print "failed to execute: $!\n";
				exit;
    	}
    	elsif ($? & 127) {
        printf  LOGF  "\tchild died with signal %d, %s coredump\n",
            ($? & 127),  ($? & 128) ? 'with' : 'without';
    	}
    	else {
    		if($? >> 8!=0){
        	printf LOGF  "\tchild exited with value %d\n", $? >> 8;
        	open TMPF, "<$filename";
        	while(<TMPF>){
        		print LOGF "\t$_" if /[^\s]+/;
        	};
        	close TMPF;
        }else{
        	$n++;
        };
    	};
    }
    print "\nSuccessfully processed $n of $nOK configurations. Done!";
    print LOGF  "\nSuccessfully processed $n of $nOK configurations!\n";
 		unlink $filename;
 		
 		if($UseCase==1){
 			print LOGF  "\n[BCVIndels]\n";
 			#Build a project file for the Indellcalling 
 			#set the path to BCV vocabulary file
 			if(!defined $Vocabulary){
 				my @cfgs=$BCV_CFG;
 				while(@cfgs){
					my $cfg=shift @cfgs;
					open RD, "< $cfg" or die "\nUnable to open BCV configuration file - $cfg";
					while(<RD>){
						last if defined $Vocabulary;
						if(/^\s*cfg\s*=\s*\"(.*)\"/){
							push @cfgs, $1;
						}elsif(/^\s*dict\s*=\s*\"(.*)\"/){
							$Vocabulary=$1;
							if(defined $BcvHome){
								$Vocabulary=~s/^\${BCV_HOME}/$BcvHome/;
							}
						};         
					};      
					close RD;          
				};
			};
 			#Build configuration files for bcvindels.pl script:
 			my %samples;
 			foreach my $src(keys %sources_ok){
 				my $smpl=$sources_ok{$src}->SampleName();
				if(!defined($smpl)){
					print LOGF "\nThe name of the sample: $src was not defined!";
					next;
				}
 				if(!defined $samples{$smpl}){
 					$samples{$smpl}={};
 					$samples{$smpl}->{out}=[];
 					$samples{$smpl}->{directions}=[];
 				}
 				my $out=$BCV_Output.$src."/";
 				push @{$samples{$smpl}->{out}},$out;
 				push @{$samples{$smpl}->{directions}}	,$sources_ok{$src}->Direction()
 			};
 			$n=0;
 			$nOK=keys %samples;
 			print LOGF  "\nTotal samples detected=$nOK\n";
 			foreach my $smpl(keys %samples){
 				my $ofile_name=$BCV_Input.$smpl.".gfas";
 				print LOGF  "\n$ofile_name";
 				if(BCVIndel::Project::build_indel_project($smpl,$Vocabulary,$ofile_name,\@{$samples{$smpl}->{out}},\@{$samples{$smpl}->{directions}})){
 					print LOGF  "\t - ok";
 					$n++;
 				}else{
 					warn "\nUnable to build indel project for the sample: $smpl!";
 					print LOGF  "\t - error";
 				}
			};
			print "\nBCVIndel configuration files successfully built for $n of $nOK samples! Done!";
			print LOGF  "\nBCVIndel configuration files successfully built for $n of $nOK samples!";
 		};
 	};
 	close LOGF;
};

build_project @ARGV;
