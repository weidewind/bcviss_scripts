#!/usr/bin/perl
use IPC::System::Simple qw(capture);
	use Data::Dumper;
	use XML::Simple;


## INPUT: $dir - directory which contains .bcvrun.prj.xml file and a folder with chromatograms
##        $prj - name of the .bcvrun.prj.xml file

my ($dir, $prj) = @ARGV;

open LOG, ">>$dir/pipelinelog.txt";
open STDERR, ">>$dir/pipelineErrlog.txt";

## check if input exist 
if (!-e $dir || !defined $dir) {
	print STDERR "No such directory: ".$dir.". Exiting\n";
	close STDERR;
	close LOG;
	exit 2;
}
   
my $path_to_prj = $dir."/".$prj; 

if (!-e $path_to_prj) {
	print STDERR "No such file: ". $path_to_prj.". Exiting\n";
	close STDERR;
	close LOG;
	exit 2;
}

## read config file
my $bcv_prj = XML::Simple->new();
		my $data   = $bcv_prj->XMLin($path_to_prj, 
		ForceArray => ['Flag','Read','Primer','PrimerDirection'],
		SuppressEmpty => 1,
		ContentKey => '-content' );
my $mode = $data->{Mode};
## run bcv


if ($mode eq "PIPELINE" || $mode eq "BCV"){
	eval {
		my $logs = capture($^X, "/home/bcviss/BCV/bcvrun-0.2.3/script/bcv_run.pl", $path_to_prj);
		print LOG $logs;
	};
	if ($@){
		print STDERR "BCV error: $@\n";
		close STDERR;
		close LOG;
		exit 3;
	}
}

my $bcv_output = $data->{BCV_Output}; # update 8_10

## check if output folder exists
if (!-e $bcv_output){
	print STDERR "No such directory: ".$bcv_output.". It seems that your configuration file is malformed.\n Exiting\n";
	close STDERR;
	close LOG;
	exit 2;
}	




if ($mode eq "PIPELINE" || $mode eq "STAP") {
    	
       my @chrom_array; # holds chromatogram names	
    	
       opendir DIR, $bcv_output or do {
		print STDERR "Can't open directory $bcv_output: $!\n Exiting\n";	
		close STDERR;
		close LOG;
		exit 2;
	   };
       while (my $folder = readdir(DIR)) {
       	print LOG $folder."\n";
       	if (-d $bcv_output."/".$folder){
    		if (($folder ne ".") & ($folder ne "..")){	
      			push @chrom_array, $folder;
     			print $folder."\t";
      		}
    	}
       }
       closedir DIR;
       
     ## run STAP "-g", "$STAPDatabase", 
       foreach my $chrom_name(@chrom_array){
       	my $path_to_chrom_folder = $bcv_output."/".$chrom_name;
       	my $path_to_chrom = $bcv_output."/".$chrom_name."/".$chrom_name.".ab1.cluster.fasta";
		if (!-e $path_to_chrom) {
			print STDERR "Was going to run STAP, but couldn't find ab1.cluster.fasta file: ".$path_to_chrom.". Looks like BCV hasn't produced it.\n";
			print STDERR "Skipping ".$chrom_name."\n";
			next;
		}
    	eval {
    		$logs = capture($^X, "/home/bcviss/STAP/scripts/rRNA_pipeline.pl",  "-d", "P", "-x", "$path_to_prj", "-i", "$path_to_chrom", "-o",  "$path_to_chrom_folder");
    		print LOG $logs;
    	};
    	if ($@){
    		print STDERR "STAP error: $@\n";
    		close STDERR;
			close LOG;
    		exit 3;
    	}
      }
    
    ## write summary  
       eval {
    		$logs = capture($^X, "/home/bcviss/pipelineFiles/stap_results.pl", "-i", "$bcv_output", "-x", "$path_to_prj" );
    		print LOG $logs;
    	};
    	if ($@){
    		print STDERR "STAP results error: $@\n";
    		close STDERR;
			close LOG;
    		exit 3;
    	}
    	
    	}
   
   
   close STDERR;
   close LOG;
    

