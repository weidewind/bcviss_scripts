#!/usr/bin/perl
use lib "/home/bcviss/STAP/scripts" ;
## Usage: rRNA_pipeline_for_one.pl -i input_fasta (one sequence)  -n 5 (number of Unclassified to skip)  
##                                 -t 1 (number of step for tree2 zoom-in) 
##                                 -o output_dir (default current) -d domain (P or E)
##                                 -l 0.0<=number<=1.0 default - 0.5 (minimal fraqtion of sequence length of the longest one to be processed by pipeline)
##                                 -f 0.0<=number<=1.0 default - 0.05 (minimal fraction in mixture)
## fisa                            -g all - whole greengenes database, otu - 99 otu gg db, named - only named isolates from gg db (default).
## fisa							   -c clusterization
## fisa							   -x configuration file


use strict;
use vars qw($opt_i $opt_n $opt_t $opt_o $opt_d $opt_l $opt_f $opt_g $opt_c $opt_x);
use Getopt::Std;
use DWU_TOOLBOX;
use Time::HiRes; #fisa 23-06
use Bio::AlignIO; #fisa 14-08
use File::Copy; #fisa 25-08

my $start_time = [Time::HiRes::gettimeofday()]; #fisa 23-06

###neva
use OSPath;
my $OS_name=$^O;

getopt('intodlfgcx');
my $confidence_threshold = 0.8;
my $length_threshold = 200;



my $tree_type="ml"; ##only build ml tree
#my $boot=0; ## no bootstraping
my $boot=-4; ## Approximate LRT for branches
#neva 12.03.11
my $min_confidence=0.8;

#fisa
$| = 1;

if ($opt_x){
	if ($opt_g){
		print STDERR "\nPlease choose either -g or -x option, but not both";
		die;
	}
}
elsif (!$opt_g) {$opt_g="named";} ## if no config and no opt_g 
my $tool=DWU_TOOLBOX->new($opt_x, $opt_g); ## this step setups all the db and bin paths

if(!(defined $opt_t)) {$opt_t=1;} ## number of steps to trace back for the zoomin tree2 (always leave 2 steps from the top) 
if(!($opt_t>=0)) {$opt_t=1;}

if(!$opt_n) {$opt_n=5;} ## number of Unclassified records to skip

if(!$opt_l) {$opt_l=0.5;} #neva
if(!($opt_l<=1.0)) {$opt_l=1.0;};
if(!($opt_l>=0)) {$opt_l=0;};

if(!$opt_f) {$opt_f=0.05;} #neva
if(!($opt_f<=1.0)) {$opt_f=1.0;};
if(!($opt_f>=0)) {$opt_f=0;};

if(!($opt_c)) {$opt_c=1;};


if(!$opt_o) {
    $opt_o=$ENV{PWD};
}
$opt_o=~s/\/$//;
###neva
if($OS_name=~m/cygwin/i){
	$opt_o=convert_cygwin2perl_win_path $opt_o;
};

$opt_o.="/";  ## output directory, global


my $in_acc;
my @in_acc; #neva
my $in_seq;
my @in_seq; #neva
my $in_seq_count=0;
#sequence fraqtion in mixture;
my @in_fraq;
my $in_fraq;

my $tree2_taxonomy_level; #neva
my %taxonomy_hash; #neva
my $errlog; #neva

open(IN,"$opt_i") || die "cannot open input file $opt_i \n";
my $max_seq_len=0;

while(<IN>){
    if($_=~/^>/) {
        $in_seq_count++;
        #neva if($in_seq_count >= 2) {last;}
        #neva
        if( defined $in_seq){
        	$in_seq=uc $in_seq;
					$in_seq=~s/[^A-Z]//g;
        	push @in_seq, $in_seq ;
        	my $l=length $in_seq;
        	$max_seq_len=$l if($l>$max_seq_len);
        };
        if($_=~m/^.+?\s+(0\.\d+)/){
        	$in_fraq=$1;
        };
        ($in_acc)=split(/\s+/,$_);
        if(defined $in_fraq){
        	my $str=sprintf("%.3f", $in_fraq);
        	$in_acc.="_$str" ;
        	push @in_fraq,$in_fraq;
        };
				$in_acc=~s/^>//;
        $in_acc=~s/[^a-zA-Z0-9_\.]+/_/g;
        push @in_acc, $in_acc; #neva
        $in_seq=undef; #neva
        $in_fraq=undef;
    }elsif(/^\s*;/){
    	next;
    }else{
    		$_=~s/;(.+)//;
        $_=~s/[\-\s]+//g; #neva
				$in_seq.=$_;
				$in_seq.=$_ if $_ ne "";
    }
}
close IN;

#neva
if( defined $in_seq){
	$in_seq=uc $in_seq;
	$in_seq=~s/[^A-Z]//g;
	push @in_seq, $in_seq; 
	my $l=length $in_seq;
  $max_seq_len=$l if($l>$max_seq_len);
};
if(scalar @in_seq != scalar @in_acc){
	die "\nWrong input FASTA file: $opt_i";
};
if(scalar @in_seq != scalar @in_fraq){
	print STDERR "\nThe file $opt_i isn't a valid mixture! Not all sequences have estimated frequencies!";
	@in_fraq=();
};
###

for(my $idx=0;$idx<@in_acc;$idx++){ #neva
	$in_acc=$in_acc[$idx]; #neva
	$in_seq=$in_seq[$idx]; #neva
	$errlog=$opt_o.$in_acc.".errlog";   ##update 10/16/06
	my $l=length $in_seq;
	if($l/$max_seq_len<$opt_l){
		open(ERRLOG,">>$errlog");      ##update 10/16/06
	  print ERRLOG "length of sequence $in_acc is less than $opt_l fraction of the longest one - ignored!\n"; ##update 10/16/06
	  close ERRLOG;   ##update 10/16/06
		next ; #neva
	};
	if(@in_fraq>0 && $in_fraq[$idx]<$opt_f){
		open(ERRLOG,">>$errlog");      ##update 10/16/06
	  print ERRLOG "sequence $in_acc has a muxture fraction less than $opt_f - ignored!\n"; ##update 10/16/06
	  close ERRLOG;   ##update 10/16/06
		next ; #neva
	};
#### $in_seq and $in_acc is the input seq and acc#########
	

	########################################################################################################
	###step one : decide the orietation of the input sequences
	my $seq_file=&create_seq_file($tool,$in_acc,$in_seq);
	next if !defined $seq_file;
	
	$in_seq='';
	open(IN,$seq_file) || die "cannot open file $seq_file";
	while(<IN>){
	    if(/>/){next;}
	    $_=~s/\s+//g;
	    $in_seq.=$_;
	}
	
	close IN;
	### the above steps are neccessary in case the sequence is the minus strand 
	
	###########decide which database to chose###########
	
	my $domain_depth;
	
	if($opt_d=~/^[PABE]/i){
	    $opt_d=substr($opt_d,0,1);
	    $opt_d=uc $opt_d;
	    $domain_depth="INPUT";
	}
	else{
	($opt_d,$domain_depth)=&get_domain_info($tool,$in_acc,$in_seq);
	}
	
	$domain_depth=$opt_d."|".$domain_depth;
	
	
	my $db=&set_db($tool,$opt_d);
	
	##############build tree one
	
	############blast to get a list
	### blast seqdb and repdb
	my $tmp_tag=$tool->{tmp_path}.$in_acc.".tmp.".int(rand(99999)).".".substr(time(),5);
	my $blast_list=$tmp_tag.".list";
	my $seq_blast_list=$blast_list.".seqdb";
	my $rep_blast_list=$blast_list.".repdb";
	print "blastn against ".$db->{seq}."   ..............\n";
	
	
		#my $top_blast_hit_number=&get_list_by_blast($tool,$seq_file,$db->{seq},$seq_blast_list,50,10,30,0.05); ## top 50 hits, low 10 hits at 30 hit interval,e cutoff 0.05)
	my @top_blast_hit_number_and_best_hit_info= @{ &get_list_by_blast($tool,$seq_file,$db->{seq},$seq_blast_list,50,10,30,0.05,0) }; ## top 50 hits, low 10 hits at 30 hit interval,e cutoff 0.05, 0 relatives) #fisa 24_06_14
	my $top_blast_hit_number = $top_blast_hit_number_and_best_hit_info[0]; #fisa 24_06_14
	my @best_hit_info = @{ $top_blast_hit_number_and_best_hit_info[1] }; #fisa 24_06_14
	
	
	if($top_blast_hit_number <= 3){  ##update 10/16/06
	  open(ERRLOG,">>$errlog");      ##update 10/16/06
	  print ERRLOG "$in_acc have less than 3 blastn hits in the database ".$db->{seq}."\n"; ##update 10/16/06
	  close ERRLOG;   ##update 10/16/06
	  unlink  $blast_list; ##update 10/16/06
	  unlink  $seq_blast_list; ##update 10/16/06 #fisa commented out
	  unlink  $rep_blast_list; ##update 10/16/06
	  #exit; ##update 10/16/06
	  next;
	}
	
	
	print "blastn against ".$db->{rep_seq}."   ..............\n";
	

	&get_list_by_blast($tool,$seq_file,$db->{rep_seq},$rep_blast_list,10,0,1,0.05,0); ## top 10 hits only at e cutoff of 0.05
	
	
	my $res_cl=&combine_list($seq_blast_list,$rep_blast_list,$blast_list); ## combine list and get rid of redundency
	next if(!$res_cl);
	
	##fisa 24-06-14 commented out 
#	my $best_blast; 
#	open(BLASTIN, $seq_blast_list);
#	
#	my $from_top=0;                ##11/1/06
#	while(<BLASTIN>){              ##11/1/06
#	    my @t=split(/\t/,$_);      ##11/1/06
#	if($from_top == $opt_n){       ##11/1/06
#	    $best_blast=$_;            ##11/1/06
#	    last;                      ##11/1/06
#	 }                             ##11/1/06
#	  if($t[0]=~/PROK3/){          ##11/1/06
#	      $from_top++;             ##11/1/06
#	      next;                    ##11/1/06
#	  }                            ##11/1/06
#	else{                          ##11/1/06 
#	    $best_blast=$_;            ##11/1/06
#	    last;                      ##11/1/06
#	}                              ##11/1/06
#	}                              ##11/1/06
#	close BLASTIN;
	
	my $best_blast = $best_hit_info[1]; ##fisa 24-06-14
	
	
	$best_blast=~s/\s+$//;
	unlink $seq_blast_list;
	unlink $rep_blast_list;
	
	### get taxonomy info from $blast_list;
	#my %taxonomy_hash; #neva
	
	
	&get_taxa_for_list($tool,$db,$blast_list,\%taxonomy_hash);
		
	
	my $best_blast_taxa=&get_taxa_from_hash($best_blast,\%taxonomy_hash);

	
	##my $best_blast=$best_blast."|".$best_blast_taxa; ##fisa 24-06-14 commented out
	my $best_blast=$best_blast."|".$best_blast_taxa."\tIdentity=".$best_hit_info[2]."\tCovered=".($best_hit_info[7]-$best_hit_info[6]+1)
	."\tE-value=".$best_hit_info[10]."\tBit-score=".$best_hit_info[11]; ##fisa 24-06-14 
	
	# fisa 8-08-14
	if (($best_hit_info[7]-$best_hit_info[6]+1) < $length_threshold){
		open(ERRLOG,">>$errlog");     
		print ERRLOG "length of sequence $in_acc covered by best blast hit is less than $length_threshold - ignored!\n"; 
		close ERRLOG;   
		next ; 
	}
	##
	
	###########build alignments
	my $alignment_file=$tmp_tag.".ali";
	
	&build_alignment($tool,$seq_file,$in_acc,$blast_list,$db,$alignment_file);
	unlink $blast_list;
	
	##fisa 14-08
	if ($opt_c){
	    my $clustered_alignment_file = $tmp_tag.".clustered.ali";
		$tool->clusterize($alignment_file, $clustered_alignment_file);
		$alignment_file = $clustered_alignment_file;
	}
	#
	
	##########build tree 1
	
	my $tree1_file=&build_tree($tool, $alignment_file, $opt_o.$in_acc, 1, $tree_type, $boot, 1);
	
	unlink $alignment_file; 
	$alignment_file=$tmp_tag.".ali";

	######## processing tree1
	my $tree_tab=$tool->njtab($in_acc,$tree1_file);
	
	my $treetab1_file=$tree1_file.".tab";
	
	&output_treetab_to_file($tree_tab, $treetab1_file,\%taxonomy_hash);
	
	
	if(!(-s $treetab1_file)){
	    &generate_report($in_acc,$opt_o,$best_blast,$domain_depth); ##update 10/16/06
	}
	else{  ##update 10/16/06
	    &add_taxa_to_tree($tree1_file,$treetab1_file); ##update 10/16/06
	} ##update 10/16/06
	
	$tree2_taxonomy_level=&tree2_taxonomy_level_bcv($tree_tab,\%taxonomy_hash,$opt_t,$opt_n); #neva
	
	if(!$tree2_taxonomy_level){
	    print STDERR "Cannot define taxonomy level for database creation\n";
	    &generate_report($in_acc,$opt_o,$best_blast,$domain_depth,$treetab1_file); ###update 10/16/2006
		next; ##fisa 08-08-14
	}
	else {
		print STDOUT "Using $tree2_taxonomy_level taxonomy level to create database\n";
	}
	

	
	my $new_blastn_db=$tmp_tag.".blastdb";
	&build_blastn_db_from_taxa_id($tree2_taxonomy_level,$db, $tool, $new_blastn_db);
	
	
	
	print "blastn against new_blastn_db for tree2 building ..............\n";
	&get_list_by_blast($tool,$seq_file,$new_blastn_db,$seq_blast_list,40,10,30,0.05,10); ## top 40 hits, low 10 hots at 30 hit interval,e cutoff 0.05, 10 more distant relatives (maximum)) 
	
	
	unlink  $new_blastn_db;
	$tool->clean_formatdb($new_blastn_db);
	
	
	
	my $outgroup_acc=&chose_outgroup_acc($tool,$seq_file,$db,$tree2_taxonomy_level,$rep_blast_list);
	if (!$outgroup_acc) {
	print STDERR "Cannot find outgroup\n";
	 &generate_report($in_acc,$opt_o,$best_blast,$domain_depth,$treetab1_file); ##update 10/16/2006
	}
	
	
	$res_cl=&combine_list($seq_blast_list,$rep_blast_list,$blast_list); ## combine list and get rid of redundency

	
	unlink $seq_blast_list; 
	unlink $rep_blast_list;
	next if(!$res_cl);
	&get_taxa_for_list($tool,$db,$blast_list,\%taxonomy_hash);

	&build_alignment($tool,$seq_file,$in_acc,$blast_list,$db,$alignment_file);
	unlink $blast_list;
	
	##fisa 14-08
	if ($opt_c){
	    my $clustered_alignment_file = $tmp_tag.".clustered.ali";
		$tool->clusterize($alignment_file, $clustered_alignment_file, $opt_o);
		$alignment_file = $clustered_alignment_file;
	}
	#
	
	##########build tree 2
	my $tree2_file=&build_tree($tool, $alignment_file, $opt_o.$in_acc ,2, $tree_type,$boot,1);
	copy ($tree2_file, $tree2_file.".native") or warn "CANNOT COPY ".$tree2_file;
	
	unlink $alignment_file;
	
	if(!(-s $tree2_file)){      ####update 10/16/06
	&generate_report($in_acc,$opt_o,$best_blast,$domain_depth,$treetab1_file); ####update 10/16/06
	}                           ###update 10/16/06
	
	
	######## processing tree2
	$tree_tab=$tool->njtab($in_acc,$tree2_file,$outgroup_acc);
	
	my $treetab2_file=$tree2_file.".tab";
	
	&output_treetab_to_file($tree_tab, $treetab2_file,\%taxonomy_hash);
	
	if(-s $treetab2_file){ ##update 10/16/06
	   &add_taxa_to_tree($tree2_file,$treetab2_file); ##update 10/16/06
	&generate_report($in_acc,$opt_o,$best_blast,$domain_depth,$treetab1_file,$treetab2_file);  ##update 10/16/06
	
	   
	
	}                               ##update 10/16/06
	else{                                                                      ##update 10/16/06
	&generate_report($in_acc,$opt_o,$best_blast,$domain_depth,$treetab1_file); ##update 10/16/06
	} ##update 10/16/06

}

print STDOUT "Time used ", Time::HiRes::tv_interval($start_time), "\n", ##fisa 23-06-14
##exit 1; #neva 
exit 0; #fisa 


sub add_taxa_to_tree{
    my ($tree_file,$tab)=@_;
    my $new_file=$tree_file.".with_taxonomy";
    my $tree;

    open(TXTREIN, $tree_file);
    while(<TXTREIN>){
     $tree.=$_;
     }
     close TXTREIN;

     open(TXTREIN,$tab);
while(<TXTREIN>){
    my @t=split(/\t/);
    $t[-1]=~s/\|/_/g;
    $t[1]=~s/\./\\\./g;
    my ($acc)=split(/_/,$t[1]);
    my $old_acc=$t[1];
    my $new_acc=$acc."_".$t[-1];

   if($old_acc && $new_acc){
   $tree=~s/\b$old_acc\b/$new_acc/;
   }
   }
close TXTREIN;
   
open(TXTREOUT,">$new_file");
print TXTREOUT $tree;
close TXTREOUT;

return;

}


#neva 12.03.11
sub generate_report{
    my ($in_acc,$opt_o,$best_blast,$domain_depth,$treetab1_file,$treetab2_file)=@_;
    my $result_file=$opt_o.$in_acc.".results";


    open (RESULTS,">$result_file") || die "cannot output results to file $result_file\n";
    print RESULTS $in_acc;   ###update 10/16/06

  	if($best_blast){ 
  		print RESULTS "\tBLAST=".$best_blast; ##Update 10/16/06
   	}

		
    if(open(TABIN,$treetab1_file)){
				my $line=<TABIN>;
        $line=~s/\s+$//; 
        my @t=split(/\t/,$line);
        my $a1=$t[1];
        my $a2=$t[-1];
        if($boot!=0){
        	my $conf=$t[-2];
        	if($conf<$min_confidence){
        		my $tmp=$a1;
        		$a1=~s/^\w+?_//;
        		while(<TABIN>){
        			s/\s+$//;
        			$line=$_;
        			@t=split(/\t/,$line);
        			my ($b1,$b2)=($t[1],$t[-1]);     
        			$conf=$t[-2]; 		
        			if($conf>=$min_confidence){
        				$b1=~s/^\w+?_//;
        				last if($a1 eq $b1);
        				if("$a1 $b1" =~ /^([A-Za-z]+\d+(\.\d+)+)((\.\d+)*) \1((\.\d+)*)$/){
        					$a1=$1;
        				};
								$a2=~s/\s//g;
								$b2=~s/\s//g;
        				if("$a2 $b2" =~ /^(\w+(\|\w+)*)((\|\w+)*) \1((\|\w+)*)$/){
        					$a2=$1;
        				};    			
        				last;
        			};
        		};
        	};
        };
        print RESULTS "\tTREE1=".$a1."|".$a2;
				close TABIN;
    }
    if(open(TABIN,$treetab2_file)){
				my $line=<TABIN>;
        $line=~s/\s+$//; 
        my @t=split(/\t/,$line);
        my $a1=$t[1];
        my $a2=$t[-1];
        if($boot!=0){
        	my $conf=$t[-2];
        	if($conf<$min_confidence){
        		my $tmp=$a1;
        		$a1=~s/^\w+?_//;
        		while(<TABIN>){
        			s/\s+$//;
        			$line=$_;
        			@t=split(/\t/,$line);
        			my ($b1,$b2)=($t[1],$t[-1]);     
        			$conf=$t[-2]; 		
        			if($conf>=$min_confidence){
        				$b1=~s/^\w+?_//;
        				last if($a1 eq $b1);
        				if("$a1 $b1" =~ /^([A-Za-z]+\d+(\.\d+)+)((\.\d+)*) \1((\.\d+)*)$/){
        					$a1=$1;
        				};
								$a2=~s/\s//g;
								$b2=~s/\s//g;
        				if("$a2 $b2" =~ /^(\w+(\|\w+)*)((\|\w+)*) \1((\|\w+)*)$/){
        					$a2=$1;
        				};    			
        				last;
        			};
        		};
        	};
        };
        print RESULTS "\tTREE2=".$a1."|".$a2;
				close TABIN;
    }

    if($domain_depth){
       print RESULTS "\tDOMAIN=".$domain_depth;   
    }
    print RESULTS "\n";
    close RESULTS;
}

sub build_tree{
    my ($tool,$alignment_file,$tree_file_head,$tree_file_tail,$tree_type,$boot,$optimize)=@_;
    my $tree_file; 

    $tree_file=$tree_file_head.".mltree".$tree_file_tail;
    $tool->rRNA_gde2mltree_linux($alignment_file,$tree_file,$boot,$optimize);

    return $tree_file;
}

sub get_taxa_for_list{
    my ($tool,$db,$blast_list,$hash_ref)=@_;

    open(LIST, $blast_list) || return 0;
while(<LIST>){
    my ($acc)=split(/\s+/);
    my @t=split(/_/,$acc);
    my $taxa_id=$t[-1];
my $taxa=$tool->get_taxa_by_id($taxa_id,$db->{xml_index});
$taxonomy_hash{$acc}=$taxa;
}
close LIST;

    return 1;
}


sub chose_outgroup_acc{
    my ($tool,$seq_file,$db,$tree2_taxonomy_level,$rep_blast_list)=@_;

    my @rep_out_group;

     &get_list_by_blast($tool,$seq_file,$db->{rep_seq},$rep_blast_list,100,0,1,0.05, 0); ## top 50 hits only at e cutoff of 0.05


    open (REPIN,$rep_blast_list) || return 0;
while(<REPIN>){
    my ($this_acc)=split(/\s+/,$_);
    my @t=split(/_/,$this_acc);
    my @rep_level=split(/\./,$t[1]);
    my @tree2_level=split(/\./,$tree2_taxonomy_level);

## if rep acc have more levels than tree2 level
## shrink rep acc to the same level and compare
    while(@rep_level > @tree2_level){
      pop @rep_level;
    }
## if rep acc have fewer levels than tree2 level
## shrink tree2 level to compare
    while(@rep_level < @tree2_level){
      pop @tree2_level;
    }
## at the same level, tree2 and rep acc cannot be the same

    my $rep_level=join("_",@rep_level);
    my $tree2_level=join("_",@tree2_level);
 
    if($rep_level eq $tree2_level){next;}

    push(@rep_out_group, $this_acc);
    if(@rep_out_group >=2 ){last;} 
    
}
close REPIN;

    unlink $rep_blast_list;

    if(@rep_out_group){
    open(REPOUT,">$rep_blast_list") || return 0;
    print REPOUT join("\n",@rep_out_group)."\n";
    close REPOUT;    
    return $rep_out_group[-1];
    }
    else{return 0;}
}






sub build_blastn_db_from_taxa_id{
my ($taxa_id, $db, $tool, $blast_db_file)=@_;
my $acc_list=$tool->get_acc_by_id($tree2_taxonomy_level,$db->{xml_index});
my $acc_list_file=$blast_db_file.".acclist.from.tree1";
open(LISTOUT,">$acc_list_file") || die ; 
print LISTOUT $acc_list; 
close LISTOUT;
$tool->fetch_seq_by_accession($acc_list_file,$db->{seq_index},$blast_db_file);
$tool->run_formatdb($blast_db_file);
unlink $acc_list_file;
return;
}

#########start to decide which level to use to draw tree2
## go up $opt_t step(s)
## if I hit unclassified or clade less than 50 records go to the next step untill $opt_n records have been screened
## if my current level is phylum level or lower, I will draw tree2 or copy tree1 and exit
sub tree2_taxonomy_level{
    my ($tree_tab, $taxonomy_hash, $opt_t, $opt_n)=@_;
    my $line_count=0;
    my $taxa_id;
    my @tab=split(/\n/,$tree_tab);
    
    foreach my $line(@tab) {
        $line_count++;
        if($line_count > $opt_n){$taxa_id="";last;}
		my @t=split(/\t/,$line);
        my $taxa=$taxonomy_hash->{$t[1]};
		
        my @taxa=split(/\t/,$taxa);
        my $i=0;


        while($i < $opt_t){
            pop @taxa; 
            if(@taxa <= 2) {last;}
            $i++;
           }
      
        while(@taxa >= 2){
            my $this_tab=pop(@taxa);
            my ($id,$name,$number)=split(/:/,$this_tab);
            if($number < 50) {next;}
            if($name=~/Unclassif/i) {
		next;
	    }
            $taxa_id=$id;
            last;
        }
            
	if($taxa_id){last;}
    }

    return $taxa_id;
}

sub tree2_taxonomy_level_bcv{
    my ($tree_tab, $taxonomy_hash, $opt_t, $opt_n)=@_;
    my $line_count=0;
    my $my_id = "";
    my @tab=split(/\n/,$tree_tab);
	my $nearest = (split(/\t/,$tab[0]))[1];
	my $nearest_taxonomy = (split(/K/,$nearest))[2];
	my $nearest_rank = scalar split (/\./, $nearest_taxonomy);
	
    my @taxonomies;
	my $confidence_level;
	my $limiting_node;
    foreach my $line(@tab) {
		my @t1=split(/\t/,$line);
        my @t2=split(/K/,$t1[1]);
		my $taxonomy = $t2[2];
		
		if ($limiting_node){
			if ($limiting_node != $t1[5]){
				last;
			}	
		}
		else {
			if ($t1[6] >= $confidence_threshold){
				$limiting_node = $t1[5]
			}
		}
		push @taxonomies, $taxonomy;
		$line_count++;
		}
		
		
		my $common_taxonomy = longest_common_prefix(@taxonomies);
		if ((substr $common_taxonomy,-1,1) eq '.') {
			chop $common_taxonomy; 
		}

		my $rank = scalar split(/\./, $common_taxonomy);
		## if ($rank > 2){ ##fisa commented this out - cheating 
			my $taxa_tab = (split /\t/, $taxonomy_hash->{$nearest})[$rank-1];
			my ($id,$name,$number)=split(/:/,$taxa_tab);
			$my_id = $id;
		
		##fisa commented this out at 15_08, hoping it won't be necessary since the rank is rised 
		#	while($line_count < $opt_n & ($name =~ m/Unclassif/i || $number < 50)){
		#		my $line = $tab[$line_count];
		#		my @t1=split(/\t/,$line);
		#		my @t2=split(/K/,$t1[1]);
	#			my $taxonomy = $t2[2];
	#			$common_taxonomy = longest_common_prefix(($common_taxonomy,$taxonomy));
#				print STDERR "this is common_taxonomy again ".$common_taxonomy."\n";
	#			my $rank = scalar split(/\./, $common_taxonomy);
	#			if ($rank <= 2){last;}
	#			$taxa_tab = (split /\t/, $taxonomy_hash->{$nearest})[$rank-1];
#				($my_id,$name,$number)=split(/:/,$taxa_tab);#
#				$line_count++;
#			}

		#} ##fisa cheating 
		
		my @arr_my_id = split /\./, $my_id;
		my $counter = 0;
		
		if (scalar @arr_my_id == $nearest_rank){ ## fisa 15_08
			while ($counter < $opt_t){
				if (scalar @arr_my_id == 1) {last;} ## fisa cheating, isn't she? If you can't find a proper subtree on this step, you should quit. But this line guarantees you won't.
				pop @arr_my_id;
				$counter++;
			}
		}
		return join ".", @arr_my_id;
		
       
}




sub output_treetab_to_file{
    my ($treetab,$file,$taxonomy_hash)=@_;

open(TREETAB, ">$file") || die "cannot output treetab $file\n";
my @tree_tab=split(/\n/,$treetab);
#########output tab for tree1
foreach my $line(@tree_tab){
    my @t=split(/\t/,$line);
my $taxa_info=&get_taxa_from_hash($t[1],$taxonomy_hash);
$line.="\t".$taxa_info."\n";
print TREETAB $line;
}
close TREETAB;
#######finish output
    return;
}



sub get_taxa_from_hash{
    my ($acc,$taxonomy_hash)=@_;
    my $taxa_info=$taxonomy_hash->{$acc};
    my @taxa_info=split(/\t/,$taxa_info);
    my @tmp_taxa;
foreach my $tmp_taxa(@taxa_info){
    my @tmp=split(/:/,$tmp_taxa);
    push(@tmp_taxa,$tmp[1]);
}
return join("|",@tmp_taxa);

}



sub build_alignment{
my ($tool,$seq_file,$in_acc,$blast_list,$db,$output)=@_;
my $gde_file=$output.".tmpgde";
my $gde;

$tool->fetch_seq_by_accession($blast_list, $db->{ali_index},$gde_file);

open(GDEIN,$gde_file) || die "cannot open alignment file $gde_file \n";
while(<GDEIN>){
    $gde.=$_;
}
close GDEIN;

$gde=$tool->trim_all_gap($gde);

open(GDEOUT,">$gde_file") || die "cannot output alignment file $gde_file \n";
print GDEOUT $gde;
close GDEOUT;

my $tmp_one_seq_file=$gde_file.".one_seq";

$tool->align_to_profile($seq_file,$gde_file,$tmp_one_seq_file);


system "cat $tmp_one_seq_file >> $gde_file";

unlink $tmp_one_seq_file;


$gde="";
open(GDEIN,$gde_file) || die "cannot open alignment file $gde_file \n";
while(<GDEIN>){
    $gde.=$_;
}
close GDEIN;
unlink $gde_file;

$gde=$tool->trim_gap_to_reference($gde, $in_acc, 10, 0); ### >=10 gaps in query will be removed
$gde=$tool->trim_gde_by_alignment_score($gde,'IUB',30,80); ### IUB is the matrix, 20 is the score cutoff, 80 is the percentage gap cutoff



open(GDEOUT,">$output") || die "cannot output alignment file $output \n";
print GDEOUT $gde;
close GDEOUT;

}



sub combine_list{
    my ($list1,$list2,$list_out)=@_;
    my %list;
#neva
  if(!open (LISTIN,$list1)){ 
  	open(ERRLOG,">>$errlog");      ##update 10/16/06
  	print ERRLOG "cannot open list fiile $list1\n"; 
  	close ERRLOG;   ##update 10/16/06
  	return 0;
  };
    while(<LISTIN>){
	my ($acc)=split(/\s+/);
	$list{$acc}=1;
    }
    close LISTIN;
  if(!open (LISTIN,$list2)){ 
  	open(ERRLOG,">>$errlog");      ##update 10/16/06
  	print ERRLOG "cannot open list fiile $list2\n"; 
  	close ERRLOG;   ##update 10/16/06
  	return 0;
  };
    while(<LISTIN>){
	my ($acc)=split(/\s+/);
	$list{$acc}=1;
    }
    close LISTIN;
    
   if(!open (LISTOUT,">$list_out")){ 
   	open(ERRLOG,">>$errlog");      ##update 10/16/06
   	print STDERR "cannot open list file $list_out\n"; 
   	close ERRLOG;   ##update 10/16/06
   	return 0;
   };
   foreach my $acc(keys %list){
       print LISTOUT $acc."\n";
   } 
    close LISTOUT;
   1;
}




sub get_list_by_blast{
    my ($tool,$seq_file,$db,$output,$top_hit_number, $low_hit_number, $low_hit_interval,$e_cutoff, $relatives_max_number)=@_;
    my $tmp_blast_report=$output.".".int(rand(999)).".".substr(time(),5).".blastn";
    if (!(length($e_cutoff) >= 1)){
	$e_cutoff=1e-2;
    }
    my $v=1500;
    my $b=1500;

    if(!(length($top_hit_number) >=1)){
      $top_hit_number=50; 
       }
    if(!(length($low_hit_number) >=1)){
      $low_hit_number=0; 
       }
     if(!($low_hit_interval >=1)){
	 $low_hit_interval=10;
       } 
     if(!($relatives_max_number >=1)){
	 $relatives_max_number=0;
       }  	   

    $tool->run_blastn($seq_file, $db, $tmp_blast_report, $e_cutoff, $v, $b);

    my %in;

my    $top_count=0;
my    $low_count=0;
my    $low_all_count=0;   
    my $id=0;

    my @list;
	my @additional_list; # fisa 11-08, updating relatives
    open(TMPBLAST,$tmp_blast_report) || die "cannot open blast report $tmp_blast_report\n";
    my $is_first = 1; #fisa 24_06_14
	my @best_hit_info; #fisa 24_06_14
	my $best_hit_tax;
	my $relatives_count = 0;
	my %relatives;
	
	while(<TMPBLAST>){
	#my ($q_acc,$h_acc)=split(/\s+/);
	
	my ($q_acc,$h_acc)=split(/\s+/);  #fisa 24_06_14
	if ($is_first == 1) {
	@best_hit_info =  split(/\s+/); #fisa 24_06_14
	my @temp = split /K/, $best_hit_info[1];
	$best_hit_tax = $temp[2];
	} #fisa 24_06_14
	$is_first = 0; #fisa 24_06_14
        if($in{$h_acc}){next;}
        if($q_acc eq $h_acc) {next;}       

    ##    $in{$h_acc}=1; commented out by fisa (06-08, adding relatives)
        if($top_count < $top_hit_number){
	    $in{$h_acc}=1; ##fisa
		push(@list,$h_acc);
            $top_count++;
            next;            
		}
		##fisa 11-08
		elsif ($top_count < $top_hit_number+$relatives_max_number){
			push(@additional_list,$h_acc);
			$in{$h_acc}=1;
            $top_count++;
            next;            

		}
		#
		elsif ($relatives_count < $relatives_max_number || $low_count < $low_hit_number){
		if ($relatives_count < $relatives_max_number){
			my @temp = split /K/, $h_acc;
			my $hit_tax = $temp[2];
			if ($hit_tax ne $best_hit_tax){
				if(!exists $relatives{$hit_tax} || $relatives{$hit_tax} == 1){
					$in{$h_acc}=1;
					$relatives_count++;
					push(@list,$h_acc);
					if (exists $relatives{$hit_tax}) {$relatives{$hit_tax} = 2;}
					else {$relatives{$hit_tax} = 1;}
				}
			}
		}
		if($low_count < $low_hit_number && !(exists $in{$h_acc})){
			$low_all_count++;
            if((($low_all_count)%($low_hit_interval)) == 0){
				$in{$h_acc}=1; ##fisa
				$low_count++;
                push(@list,$h_acc);
			}       
		}
		}
    
        else{
	    last;
	}
    }
	
	while ($relatives_count < $relatives_max_number && @additional_list){	
		my $h_acc = shift @additional_list;
		push(@list, $h_acc);
		$in{$h_acc}=1;
		$relatives_count++; 
	}
    close TMPBLAST;          

    #unlink $tmp_blast_report;
  
    open(TMPLIST,">$output") || die "cannot output temp blast list file $output \n";
    foreach my $acc(@list){
	print TMPLIST $acc."\n";
    }
    close TMPLIST;

    undef %in;
    my @array = ($top_count, \@best_hit_info); #fisa 24_06_14
	return \@array; #fisa 24_06_14
   # return $top_count; ##update 10/16/06

}



sub set_db {
    my ($tool,$opt_d)=@_;
    if(!($opt_d=~/[PABE]/i)){
	die "cannot decide which database to use\n";
    }

    my $db;
    if($opt_d=~/[PAB]/){
	$db->{ali}=$tool->{db_path}."prok.ali";
        $db->{ali_index}=$tool->{db_path}."prok.ali.index";
        $db->{rep_seq}=$tool->{db_path}."prok.rep.seq";
	$db->{seq}=$tool->{db_path}."prok.seq";
        $db->{seq_index}=$tool->{db_path}."prok.seq.index";
        $db->{xml_index}=$tool->{db_path}."prok.xml.index";
    }    
    else{
	$db->{ali}=$tool->{db_path}."euk.ali";
        $db->{ali_index}=$tool->{db_path}."euk.ali.index";
        $db->{rep_seq}=$tool->{db_path}."euk.rep.seq";
	$db->{seq}=$tool->{db_path}."euk.seq";
        $db->{seq_index}=$tool->{db_path}."euk.seq.index";
        $db->{xml_index}=$tool->{db_path}."euk.xml.index";
    }
    return $db;
}

sub get_domain_info{
    my ($tool,$in_acc,$in_seq)=@_;

### make alignment


##$tool->align_to_profile($seq, $profile,$output);

my $tmp_seq=$tool->{tmp_path}."TMP.".$in_acc.".".int(rand(9999999))."D".substr(time(),5).".seqfile";
open(TMPOUT,">$tmp_seq");
print TMPOUT ">getdomain $in_acc \n";
print TMPOUT $in_seq."\n";
close TMPOUT;

my $ali_file=$tmp_seq.".ali";
my $profile=$tool->{db_path}."euro.rep.trim.ali";
$tool->align_to_profile($tmp_seq, $profile,$ali_file);

my $c="cat $profile >> $ali_file";
system $c;
my $tmp_tree=$tmp_seq.".tmp_tree";

    $tool->rRNA_gde2mltree_linux($ali_file,$tmp_tree,0,0);

    my $marker=$tool->{db_path}."euro.rep.trim.marker";
    my $domain_line=$tool->dist_3point_tree('getdomain', $tmp_tree, $marker);

    unlink $tmp_seq;
    unlink $ali_file;
    unlink $tmp_tree;
    
    my @t=split(/\s+/,$domain_line);
    
    return (substr($t[1],0,1),$t[2]);
}


sub create_seq_file {  ## get strand info based on blast 
    my ($tool,$acc,$seq)=@_;
    my $tmp_seqfile=$tool->{tmp_path}."TMP.".$acc.".".int(rand(9999999))."D".substr(time(),5).".seqfile";

    open(TMPSEQ,">$tmp_seqfile") || die "Cannot create tmp file\n";
    print TMPSEQ ">".$acc."\n".$seq."\n";
    close TMPSEQ;
    my $tmp_out=$tmp_seqfile.".out";
    my $e=1;
    my $v=1;
    my $b=1;
    my $db=$tool->{db_path}."combo.rep.seq"; 
   
    $tool->run_blastn($tmp_seqfile, $db, $tmp_out, $e, $v, $b);
 
   open(TMPTAB,$tmp_out) || die "cannot open tmp file $tmp_out\n";
    my $line=<TMPTAB>;
    close TMPTAB;
    unlink $tmp_seqfile;
    unlink $tmp_out;
    my $in_seq_strand=0;
    my @t=split(/\t/,$line);
    if((($t[7]-$t[6])*($t[9]-$t[8])) > 0){$in_seq_strand=1;}
    elsif((($t[7]-$t[6])*($t[9]-$t[8])) < 0) {$in_seq_strand=-1;}

    my $seq_file=$opt_o.$in_acc.".seq";
 
if($in_seq_strand > 0) {
    open(OUT,">$seq_file");
    print OUT ">".$in_acc."\n";
    print OUT $in_seq."\n";
    close OUT;
}
elsif($in_seq_strand <0){
    $in_seq=reverse $in_seq;
    $in_seq=~tr/ATGC/TACG/;
    open(OUT,">$seq_file");
    print OUT ">".$in_acc." complementary strand \n";
    print OUT $in_seq."\n";
    close OUT;
}
else{
    
    open(OUT,">>$errlog"); ##update 10/16/06
    print OUT "Sequence $in_acc may not be SSU rRNA sequences\n";
    close OUT;  
    $seq_file=undef;
    }
    return $seq_file;
}

sub longest_common_prefix {
    my $prefix = shift;
    for (@_) {
	chop $prefix while (! /^\Q$prefix\E/);
    }
    return $prefix;
}


