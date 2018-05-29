#!/usr/bin/perl
use lib "C:/Users/weidewind/Programs/STAP/stap/stap_native/" ;
use strict;
use DWU_TOOLBOX;
## align_to_rRNA_profile.pl takes a ss-rRNA multifasta format, it will output alignments of the sequences in
## the multifasta file using clustalw profile aligning methods

my $desc="
    Usage: align_to_rRNA_profile.pl -i input_seq -d domain -o output_file
           -i: sequence input in fasta format file (only 1 sequence in one file)
           -d: domain (E for eukaryotes, P for prokaryotes:default)
           -o: output filename
         ";

################################################
my %db;

my %opt=@ARGV;
my $opt_i=$opt{'-i'};
my $opt_d=$opt{'-d'};
my $opt_o=$opt{'-o'};

if($opt_d=~/^(A|B|P)/i){$opt_d='P';}
if($opt_d=~/^E/i){$opt_d='E';}


if(!($opt_d=~/^(P|E)$/)){die "Please specify domain: align_to_rRNA_profile.pl -i input_fasta -o output -d P(or E)\n";   
}
my $db_seq_number=20;

my $tool=DWU_TOOLBOX->new(); ## this step setup all the db and bin paths

my $this_seq;
my $tmp_seq;
my $this_acc;
my $tmp_acc;

open(OUT,">$opt_o") || die "cannot create file $opt_o\n";
open(IN,"$opt_i") || die "cannot open input file $opt_i \n";
while(<IN>){
    
    if($_=~/^>/) {
         $this_acc=$tmp_acc;
         $this_seq=$tmp_seq;
         $tmp_seq='';
        ($tmp_acc)=split(/\s+/,$_);
         $tmp_acc=~s/^>//;
         $tmp_acc=~s/[^a-zA-Z0-9_\.]+/_/g;
        if($this_acc && $this_seq){
        my $gde=run_one_seq($tool,$this_acc,$this_seq,$opt_o);
         if($gde){
         print OUT $gde;
          } 
        }
    }    
    else{
        $_=~s/\s+//g;
        $tmp_seq.=$_;
    }
  }
close IN;

       
       $this_acc=$tmp_acc;
       $this_seq=$tmp_seq;
      
      if($this_acc && $this_seq){
      
        my $gde=run_one_seq($tool,$this_acc,$this_seq,$opt_o);
        if($gde){
         print OUT $gde;
          }
        }

close OUT;





sub run_one_seq{
my ($tool,$in_acc,$in_seq,$opt_o)=@_;
    $in_seq=uc $in_seq;
my $errlog=$opt_o."_".$in_acc.".profile_ali.errlog";   ##update 10/16/06
###step one : decide the orietation of the input sequences
my $seq_file=&create_seq_file($tool,$in_acc,$in_seq,$opt_o,$errlog);

if(!(-s $seq_file)){return;}
$in_seq='';
open(SEQIN,$seq_file) || die "cannot open file $seq_file";
while(<SEQIN>){
    if(/>/){next;}
    $_=~s/\s+//g;
    $in_seq.=$_;
}
close SEQIN;
### the above steps are neccessary in case the sequence is the minus strand 

my $db=&set_db($tool,$opt_d);


############blast to get a list
### blast seqdb and repdb
my $tmp_tag=$tool->{tmp_path}.$in_acc.".tmp.".int(rand(99999)).".".substr(time(),5);
my $seq_blast_list=$tmp_tag.".list";

print "blastn against ".$db->{seq}."   ..............\n";
my $top_blast_hit_number=&get_list_by_blast($tool,$seq_file,$db->{seq},$seq_blast_list,20,0,0,1e-4); 
## top 20 hits, low 0 hits at 0 hit interval,e cutoff 1e-4

if($top_blast_hit_number <= 3){  
  open(ERRLOG,">>$errlog");      
  print ERRLOG "$in_acc have less than 3 blastn hits in the database ".$db->{seq}."\n";
  close ERRLOG;   
  unlink  $seq_blast_list;
  unlink  $seq_file;
  return; 
}


my $gde=&build_alignment($tool,$seq_file,$in_acc,$seq_blast_list,$db,$opt_o);

unlink $seq_blast_list;
unlink $seq_file;

if($gde=~/^>[^\n]+\n(-|[A-Z])/){
return $gde;
}

}




sub make_all_gap_mask{
my $gde=shift;
my %seq;
my @line=split(/\n/,$gde);
my $acc;
foreach my $line(@line){
if($line=~/^(>|\%|\#)/){
($acc)=split(/\s+/,$line);
  }
else{
$line=~s/\s+//g;
$line=uc $line;
$seq{$acc}.=$line;
    }
 }

my @mask;

foreach my $acc(keys %seq){ 
my @seq=split(//,$seq{$acc});
my $len=scalar @seq;

if(!((scalar @mask) >= 1)){
my $i;
for $i(0..$len-1){
$mask[$i]=0;
    }
       }

my $i;
for $i(0..$len-1){
if($seq[$i]=~/[A-Za-z]/){
  $mask[$i]=1;
  }
    }
        }
my $mask=join("",@mask);
return $mask;

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

my $mask=&make_all_gap_mask($gde); ###write this one 
my @mask=split(//,$mask);


$gde=$tool->trim_all_gap($gde);

open(GDEOUT,">$gde_file") || die "cannot output alignment file $gde_file \n";
print GDEOUT $gde;
close GDEOUT;

my $tmp_one_seq_file=$gde_file.".one_seq";

$tool->align_to_profile($seq_file,$gde_file,$tmp_one_seq_file);

$gde='';

my $acc;
my $tmp_seq;
my @tmp_seq;
my $new_seq;

open(GDEIN, $tmp_one_seq_file) || die "cannot open tmp file $tmp_one_seq_file \n";
while(<GDEIN>){
if(/>/){
($acc)=split(/\s+/);
  }
else{
$_=~s/\s+//g;
$tmp_seq.=$_;
    }
     }
close GDEIN;

@tmp_seq=split(//,$tmp_seq);

while(@mask){
my $ok=shift @mask;
if($ok=~/1/){
my $nt=shift(@tmp_seq);
$new_seq.=$nt;
  }
else{
$new_seq.="-";
    }
      }

$gde=$acc."\n".$new_seq."\n";
unlink $tmp_one_seq_file;
unlink $gde_file;

return $gde;

}

####################################



sub create_seq_file {  ## get strand info based on blast 
    my ($tool,$in_acc,$in_seq,$opt_o,$errlog)=@_;
    my $tmp_seqfile=$tool->{tmp_path}."TMP.".$in_acc.".".int(rand(9999999))."D".substr(time(),5).".seqfile";

    open(TMPSEQ,">$tmp_seqfile") || die "Cannot create tmp file\n";
    print TMPSEQ ">".$in_acc."\n".$in_seq."\n";
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

    my $seq_file=$opt_o."_".$in_acc.".seq";
 
if($in_seq_strand > 0) {
    open(SEQOUT,">$seq_file");
    print SEQOUT ">".$in_acc."\n";
    print SEQOUT $in_seq."\n";
    close SEQOUT;
}
elsif($in_seq_strand <0){
    $in_seq=reverse $in_seq;
    $in_seq=~tr/ATGC/TACG/;
    open(SEQOUT ,">$seq_file");
    print SEQOUT ">".$in_acc." complementary strand \n";
    print SEQOUT $in_seq."\n";
    close SEQOUT;
}
else{   
    open(ERRLOG,">>$errlog"); ##update 10/16/06
    print ERRLOG "Sequence $in_acc may not be SSU rRNA sequences\n";
    close ERRLOG;  
    return;
    }
    return $seq_file;
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

sub get_list_by_blast{
    my ($tool,$seq_file,$db,$output,$top_hit_number, $low_hit_number, $low_hit_interval,$e_cutoff)=@_;
    my $tmp_blast_report=$output.".".int(rand(999)).".".substr(time(),5).".blastn";
    if (!(length($e_cutoff) >= 1)){
        $e_cutoff=1e-2;
    }
    

    if(!(length($top_hit_number) >=1)){
      $top_hit_number=50; 
       }
    if(!(length($low_hit_number) >=1)){
      $low_hit_number=0; 
       }
     if(!($low_hit_interval >=1)){
         $low_hit_interval=10;
       }   

   my $v=$top_hit_number;
   my $b=$top_hit_number;


    $tool->run_blastn($seq_file, $db, $tmp_blast_report, $e_cutoff, $v, $b);

    my %in;

my    $top_count=0;
my    $low_count=0;
my    $low_all_count=0;   
    my $id=0;

    my @list;
    open(TMPBLAST,$tmp_blast_report) || die "cannot open blast report $tmp_blast_report\n";
    while(<TMPBLAST>){
        my ($q_acc,$h_acc)=split(/\s+/);
        if($in{$h_acc}){next;}
        if($q_acc eq $h_acc) {next;}       

        $in{$h_acc}=1;
        if($top_count < $top_hit_number){
            push(@list,$h_acc);
            $top_count++;
            next;            
        }
        elsif($low_count < $low_hit_number){
            $low_all_count++;
            if((($low_all_count)%($low_hit_interval)) == 0){
                $low_count++;
                push(@list,$h_acc);
        }             
        }
        else{
            last;
        }
    }
    close TMPBLAST;          

    unlink $tmp_blast_report;
  
    open(TMPLIST,">$output") || die "cannot output temp blast list file $output \n";
    foreach my $acc(@list){
        print TMPLIST $acc."\n";
    }
    close TMPLIST;

    undef %in;

    return $top_count; ##update 10/16/06

}




