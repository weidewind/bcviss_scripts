
package DWU_TOOLBOX;
use strict;
use Bio::Index::Fasta;
use XML::DOM;
use Bio::Phylo::IO;
use Bio::Phylo::Forest::NodeRole;
use Bio::Phylo::Forest::TreeRole;
use Bio::Phylo::Taxa::Taxon;
use Bio::Phylo::Factory;
use Bio::Phylo::Listable;
use Bio::TreeIO;
use IO::String;
use Bio::SeqIO;
use File::Copy;
use XML::Simple;
use IO::Handle;


###neva
use OSPath;
my $OS_name=$^O;
sub new {  ## set all the path here
    my $class=shift;
	my $cfg=shift;
	my $db=shift;
    my $self={};
	
	my $db_name;
	if (!$cfg){
		if ($db eq "named") {
			$db_name = "db_named_dir";
		}
			elsif ($db eq "otu") {
			$db_name = "db_otu_dir";
		}
			elsif ($db eq "all") {
			$db_name = "db_dir";
		}
		else {
		print STDERR "Invalid database name. Choose from 'named', 'otu' or 'all'.";
		}
		$self->{blast_path}="";
		$self->{phyml_path}="";
		$self->{clustalw_path}=""; ##C:/Users/weidewind/Programs/clustalw/src/
		$self->{cdhit_path}="/home/bcviss/pipelinePrograms/cdhit/";
		$self->{db_path}="/home/bcviss/STAP/".$db_name."/";
		$self->{db_acclist_file_path}="/home/bcviss/STAP/".$db_name."/index_dir/";
		$self->{tmp_path}="/home/bcviss/STAP/tmp/";
	}
	else {
		my $bcv_prj = XML::Simple->new();
		my $data   = $bcv_prj->XMLin($cfg, 
		ForceArray => ['Flag','Read','Primer','PrimerDirection'],
		SuppressEmpty => 1,
		ContentKey => '-content' );
		$db=$data->{Database};
		$db_name = option_to_db_name($db);

		$self->{blast_path}=$data->{BLASTPath};
		$self->{phyml_path}=$data->{PhyMLPath};
		$self->{clustalw_path}=$data->{ClustalWPath}; ##C:/Users/weidewind/Programs/clustalw/src/
		$self->{cdhit_path}=$data->{CDHitPath};
		$self->{db_path}=$data->{DBPath}.$db_name."/";
		$self->{db_acclist_file_path}=$data->{DBPath}.$db_name."/index_dir/";
		$self->{tmp_path}=$data->{tmpPath};
		
	}
    bless $self, $class;
    return $self;
    }
    
sub option_to_db_name{
my $db = shift;
my $db_name;
if ($db eq "named") {
			$db_name = "db_named_dir";
		}
elsif ($db eq "otu") {
			$db_name = "db_otu_dir";
		}
elsif ($db eq "all") {
			$db_name = "db_dir";
		}
else {
		print STDERR "Invalid database name. Choose from 'named', 'otu' or 'all'.";
		}
return $db_name;
}
################  Don't alter any codes below this line  #####################


#############Usage: $newgde=DWU_TOOLBOX->trim_gap_to_reference($class, $gde0, $opt_t, $opt_G, $opt_e);
#############  $gde0: input gde (or fasta)
############   $opt_t: the accession of a sequence in the alignments, columns contains a strech of gaps for the accession will be trimed
############   $opt_G:  gap cutoffs for -t switch (default 10 for a valid -t input)
############   $opt_e: 1 trim the gaps in the reference sequences from both end, the middle part will be left alone, opt_G won't matter
sub trim_all_gap{
    my ($class,$gde0)=@_;
    my $gde1=DNATRIM->trim_all_gap($gde0);
    return $gde1;
}

sub trim_gap_to_reference{
    my ($class,$gde0,$opt_t,$opt_G,$opt_e)=@_;
    my $gde1=DNATRIM->trim_gap_to_reference($gde0, $opt_t, $opt_G,$opt_e);
    return $gde1; 
}
sub trim_gde_by_alignment_score  {
    my ($class,$gde,$matrix_type,$cutoff,$gap_cutoff)=@_; 
    my $gde1=DNATRIM->trim_gde_by_alignment_score( $gde,$matrix_type,$cutoff,$gap_cutoff);
    return $gde1;
}

#############Usage: DWU_TOOLBOX->run_blast($seq_file, $db, $output,$e_value,$v_cutoff,$b_cutoff)
sub run_blastn{ 
    my ($class,$seq,$db,$output,$e,$v,$b)=@_;
    if(!(length($e) >=1 )){$e=1e-4;}
    if(!($v >=1 )){$v=20;}
    if(!($b >=1 )){$b=20;}
  
  if((-e $seq) && (-e $db)){
    my $blast=$class->{blast_path}."blastall -p blastn -i $seq -d $db -o $output -m 8 -e $e -v $v -b $b";
    print $blast."\n";
    system $blast;
    }
    return;
}
##############################





#############Usage: DWU_TOOLBOX->run_formatdb($seq_file)
sub run_formatdb{ 
    my ($class,$seq)=@_;
    my $formatdb=$class->{blast_path}."formatdb -i $seq -p F -o T";
    print $formatdb."\n";
    system $formatdb;
    return;
}
##############################

#######Usage: DWU_TOOLBOX->clean_formatdb($seq_file)

sub clean_formatdb {
    my ($class,$seq)=@_;
    my $nhr=$seq.".nhr";
    my $nin=$seq.".nin";
    my $nsd=$seq.".nsd";
    my $nsi=$seq.".nsi";
    my $nsq=$seq.".nsq";

    if(-e $seq) {unlink $seq;}
    if(-e $nhr) {unlink $nhr;}
    if(-e $nin) {unlink $nin;}
    if(-e $nsd) {unlink $nsd;}
    if(-e $nsi) {unlink $nsi;}
    if(-e $nsq) {unlink $nsq;}

    return;
}

################################


#############Usage: DWU_TOOLBOX->njtab($taxon_query,$treefile,$out_group)
############ output format:
##   tab[0] query taxon
##   tab[1] one taxon from the tree
##   tab[2] distance between the query taxon and the taxon in tab[1]
##   tab[3] distance between the nodes that the query taxon and tab[1] taxon 
##          attached to
##   tab[4] number of nodes between query taxon and tab[1] taxon
##   tab[5] degree of separation between the query taxon and tab[1] taxon. if
##          a valid outgroup has been defined by -o switch, the degree of separation
##          will be based on the user defined rooted tree, otherwise, the tree 
##          will be rooted by mid-point rooting by default
##         (the tab report is sort by tab[5] and tab[3])
##########################################################################################



sub njtab {
my ($class,$taxon_input,$treefile,$out_group)=@_;
if (!(-e $treefile)) {die "cannot find treefile $treefile\n";}

my $obj=DWU_PARSETREE->new($treefile); #treefile instead of tree
if (!$obj){
    die "EXIT: Please check the format of the tree file $treefile \n";
}

if (!($obj->check_taxon($taxon_input) eq 'T')) 
{die "Taxon $taxon_input is not included in the tree $treefile \n";}

if ((($out_group) && (!($obj->check_taxon($out_group)))) || ($out_group eq $taxon_input)) 
{die "Invalid out group input $out_group\n";}

my $list=$obj->taxon_dist($taxon_input,$out_group);
$obj->end();

return join("\n",@$list);
}

###################################################


sub clusterize{
	my ($class, $alignment_file, $clustered_alignment_file, $output_dir) = @_;
	
	my $rand = int(rand(999));
	my $temp_file = $alignment_file.$rand.".temp";
	my $cluster_temp_file = $alignment_file.$rand.".clustered.temp";
	
    my $in  = Bio::SeqIO->new(-file => $alignment_file , '-format' => 'Fasta');
   # my $out = Bio::SeqIO->newFh(-file => ">$temp_file" , '-format' => 'Fasta');
    open NOGAP, ">$temp_file";
    while ( my $seq = $in->next_seq() ) {
	print NOGAP ">".$seq->id."\n";
	my $sequence = $seq->seq;
	$sequence =~ s/-//g;
	print NOGAP $sequence."\n";
    }
	close NOGAP;
	
    my $command=$class->{cdhit_path}."cd-hit-est";
    $command .=" -i $temp_file -o $cluster_temp_file -c 1 -n 9 -d 45";
    system $command;
	
	my %ids;
	my $cluster_file  = Bio::SeqIO->new(-file => $cluster_temp_file, '-format' => 'Fasta');
	while ( my $seq = $cluster_file->next_seq() ) {
		$ids{$seq->id}  = 1;
	}
	
	open CLUSTERED, ">$clustered_alignment_file";
	$in  = Bio::SeqIO->new(-file => $alignment_file , '-format' => 'Fasta');
	while ( my $seq = $in->next_seq() ) {
		if (exists $ids{$seq->id}){
			print CLUSTERED ">".$seq->id."\n";
			print CLUSTERED $seq->seq."\n";
		}
	}
	close CLUSTERED;
	
	if ($output_dir){
		my @bcv_cluster_arr = split /[\\\/]/, $cluster_temp_file;
		my $bcv_cluster = pop @bcv_cluster_arr;
		@bcv_cluster_arr = split /[\.]/, $bcv_cluster;
		$bcv_cluster = $bcv_cluster_arr[0].".".$bcv_cluster_arr[1];
		print STDOUT $cluster_temp_file.".clstr"."\n";
		print STDOUT $output_dir.$bcv_cluster.".temp.clstr"."\n";
		copy($cluster_temp_file.".clstr", $output_dir.$bcv_cluster.".temp.clstr") or die "Copy failed: $!";
	}
	
	unlink $temp_file;
	unlink $cluster_temp_file;
	
}


sub dist_3point_tree {  ### Given 3 taxa, a tree can be splited at the join point
                        ### this function report a taxonomy's distance (in term of nodes) 
                        ### to the join-point 
###Usage: DWU_TOOLBOX->dist_3point_tree(input_taxomomy, tree_file, 3_taxa_maker_file)
### return: a string: input_taxon, the marker's side, distance_in_terms_of_node_number

    my ($class,$input_taxon, $treefile, $markerfile)=@_;

    my $marker_input;

open(I,"$markerfile") || die "Cannot open three point marker file $markerfile\n";
while(<I>){
    $marker_input.=$_;
}
close I;
$marker_input=~s/^\s+//;
$marker_input=~s/\s+$//;
my @temp=split(/\s+/,$marker_input);
my @taxa_array;
my %marker;
    foreach my $marker(@temp){
	if($marker{$marker}){next;}
        $marker{$marker}=$marker;
        push(@taxa_array,$marker);
    }

   
if(!(@taxa_array == 3)){
    die "Error: the script only take 3 taxa as markers, no more, no less. \n";
}

my $tree;

open (IN, "$treefile")|| die "cannot open tree file $treefile \n";

while(<IN>){
    $tree.=$_;
}
close IN;

$tree=~s/\s+//g;
$tree=~s/;$//;

my $pa_input=$input_taxon;
$pa_input=~s/\./\\\./g;
if($pa_input && (!($tree=~/$pa_input/))) {die "Please input a valid taxon\n";}

foreach my $marker(keys %marker){
    my $m_pa=$marker;
    $m_pa=~s/\./\\\./g;
if(!($tree=~/$m_pa/)) {die "Please input valid taxa as markers\n";}
}

my $obj=DWU_PARSETREE->new($tree);
if (!$obj){
    die "EXIT: Please check the format of the tree file $treefile \n";
}

my $ref=$obj->get_intersect_nodes(\@taxa_array);

my @output;

if($ref) {
foreach my $node(@$ref) {
my $node1=$obj->get_immediate_nodes_for_a_taxon($input_taxon);
my $path=$obj->path_between_two_taxa($node,$node1);
my $txt=$obj->break_node($node);
my @t=split(/\t/,$path);
my $node_numer_output=@t-1;
my @group=split(/\n/,$txt);

foreach my $group(@group){
    if($group=~/\b$pa_input\b/){
      $txt=$group;
      last; 
    }
}

foreach my $marker(keys %marker){
    my $m_pa=$marker;
    $m_pa=~s/\./\\\./g;
    if($txt=~/\b$m_pa\b/){
    $txt=$marker;
    last; 
    }
}

push(@output,$input_taxon."\t".$txt."\t".$node_numer_output);
}
}

$obj->end();

return join("\n",@output)."\n";

}



######### get_taxa_by_id################
########Usage: DWU_TOOLBOX->get_taxa_by_id($query_id,$db)
######  Usage: DWU_TOOLBOX->get_acc_by_id($query_id,$db)
####### return:taxa hierachy or a list of sequence accession under the taxomomy id
################################

sub get_taxa_by_id {
    my ($class,$query,$db)=@_;  ##the db should be fullpath
my $parser = new XML::DOM::Parser;
my $doc = $parser->parsefile ("$db");
my $nodes = $doc->getElementsByTagName ("$query");
if (!$nodes || !($nodes->item(0))){
    $doc->dispose();
    die "cannot fetch taxa info for $query in $db \n";
     }

my $this_node=$nodes->item (0);
my @taxa;
my $taxa_line=&private_get_taxa_list($this_node,\@taxa);
 
 my $taxa_line=join("\t",@taxa);

    $doc->dispose();
    return  $taxa_line;

}


sub get_acc_by_id {
    my ($class,$query,$db)=@_;  ##the db should be fullpath
my $parser = new XML::DOM::Parser;
my $doc = $parser->parsefile ("$db");
my $nodes = $doc->getElementsByTagName ("$query");
if (!$nodes || !($nodes->item(0))){
    $doc->dispose();
    die "cannot fetch taxa info for $query in $db \n";
     }

my $this_node=$nodes->item (0);

my @file_list;

&private_get_file_list($this_node,\@file_list);
 
    my $acc_list;
 my $file_line=join("\t",@file_list);
    $doc->dispose();
    foreach my $file(@file_list){
	my $full_file=$class->{db_acclist_file_path}.$file;
        if(open(FILEACC,$full_file)){
	    while(<FILEACC>){
                my ($this_acc)=split(/\s+/);
		$acc_list.=$this_acc."\n";
	    }
	    close FILEACC;
	}
        else {
	    print "cannot open $full_file\n";
	}
    }

    return $acc_list;

}



sub private_get_taxa_list{
 
    my ($up,$ref_taxa)=@_;
    
    push(@$ref_taxa, $up->getNodeName().":".$up->getAttribute("name").":".$up->getAttribute("number"));
while(1){
    my $parent_node=$up->getParentNode();
    my $id=$parent_node->getNodeName();
    if($id eq 'root'){last;}
    my $name=$parent_node->getAttribute("name");
    my $number=$parent_node->getAttribute("number");
    push(@$ref_taxa,$id.":".$name.":".$number);
    $up=$parent_node;   
 }
@$ref_taxa=reverse @$ref_taxa;
    
    
  }

sub private_get_file_list {
    my ($down,$ref_file_list)=@_;
   
    my @children = $down->getChildNodes();
   if ($down->getNodeName eq "file"){ 
       my $file_name=$down->getFirstChild()->getData;
       $file_name=~s/^\s+//;
       $file_name=~s/\s+$//;
   push(@$ref_file_list,$file_name);
  }
   my @nodes=$down->getChildNodes();
    foreach my $kid(@nodes){
	&private_get_file_list($kid,$ref_file_list);
    }
   
  }


############## fetch_seq_by_accession############
############# Usage: DWU_TOOLBOX->fetch_seq_by_accession($acc_list_file, $index_file,$output);
###############################
sub  fetch_seq_by_accession {

my ($class,$list_file,$index_file,$out)=@_;
my $inx = Bio::Index::Fasta->new($index_file);
my @acc;


open (FETCHIN,$list_file) || die "cannot open accession list $list_file\n";
while(<FETCHIN>){
   
    $_=~s/^\s+//g;
    $_=~s/\s+$//g;
     push(@acc,$_);
}
close FETCHIN;       

open(FETCHOUT,">$out") || die "cannot output to file $out \n";
foreach my $acc(@acc){
    
      my $seq = $inx->fetch($acc);  # Returns Bio::Seq object

      if($seq){
print FETCHOUT ">".$seq->primary_id()." ";        
print FETCHOUT $seq->desc."\n";
print FETCHOUT $seq->seq()."\n";
}
  }
close FETCHOUT;
}

###############end of fetch_seq_by_accession##############



################ align_to_profile###########
############### Usage: DWU_TOOLBOX->align_to_profile($seq, $profile,$output);
###########################################

sub align_to_profile {
    my ($class, $opt_i, $opt_p, $opt_o)=@_;
    my $tmp_dir=$class->{tmp_path};
 
  
my @in=split(/\//,$opt_i);
my $temp_tag=$tmp_dir."DWU".$in[-1];
my $acc;

    

if ((!(-e $opt_i))||(!(-e $opt_p))||(!$opt_o)) {
    die "Invalid input and output for align_to_profile\n";
}
open(INSEQ,"$opt_i")||  die "Invalid input for align_to_profile\n";
open(OUTPUT,">$opt_o")||  die "Invalid output for align_to_profile\n";
my $name;
my $seq;
my $tempname;
my $tempseq;
my $temp_tag=$tmp_dir."T".int(rand(100)).(substr(time(),5,4));

while(<INSEQ>){
    if(/^>/){
	$name=$tempname;
	$seq=$tempseq;
	$tempseq='';
        ($tempname)=split(/\s+/,$_);
	$tempname=~s/>//g;
        if($seq && $name){
      	    my $seq_ali=&align_to_profile_get_ali($class,$name,$seq,$opt_p,$temp_tag);
            if($seq_ali){
            print OUTPUT ">".$name."\n".$seq_ali."\n";
	}
	}
  
    }
    else{       
        $_=uc $_;
	$tempseq.=$_;
    }
              }
close INSEQ;
	$name=$tempname;
	$seq=$tempseq;
   if($seq && $name){
      	    my $seq_ali=&align_to_profile_get_ali($class,$name,$seq,$opt_p,$temp_tag);
	    if($seq_ali) {
            print OUTPUT ">".$name."\n".$seq_ali."\n";
	}
	}
close OUTPUT;

    return 1;

}



sub align_to_profile_get_ali{
    my ($class,$acc,$seq,$opt_p,$temp_tag)=@_;
    $acc=~s/[^A-Za-z0-9]//g;
    $acc=substr($acc,0,10);
    my $seq_id="TMPSEQIDYADA";
   
  
    my %seq_hash;
    my $seq_file=$temp_tag.$acc."s.seq";
    my $ali_file=$temp_tag.$acc."g.gde";
    my $profile_file=$temp_tag.$acc."p.pro";
        my $dnd_file=$temp_tag.$acc."p.dnd";

    
    open (TMP,">$seq_file") || return;
    print TMP ">".$seq_id."\n".$seq."\n";
    close TMP;


    
   
  ##  system "cp $opt_p $profile_file";
    my $tmp_ali_profile;
    my $i=0;
    open(OPT_O,"$opt_p") || die "cannot open profile $opt_p \n";
    while(<OPT_O>){
	if(/>/){
	    $_=">ID$i\n";
            $i++;
	}
        $tmp_ali_profile.=$_;
    }
    close OPT_O;

    open(OPT_O,">$profile_file") || die "cannot create temp profile file\n";

    print OPT_O $tmp_ali_profile;
    close OPT_O;

    my $clustalw=$class->{clustalw_path}."clustalw"; #fisa: clustalw2 on vgnki, clustalw on fbb

    my $c="$clustalw -profile1=$profile_file -profile2=$seq_file -output=gde -case=upper -outfile=$ali_file  -QUICKTREE=FAST";
#    if($OS_name=~m/cygwin|win32/i){
#    	$c="$clustalw /profile1=$profile_file /profile2=$seq_file /output=gde /case=upper /outfile=$ali_file  /QUICKTREE=FAST";
#    };
    print $c;
    system $c;
 
    


    $seq='';
    $acc='';

    unlink $seq_file;
    unlink $profile_file; #fisa commented out 13-08
    unlink $dnd_file;

    open(TMPIN,"$ali_file") || die "cannot open alignment file $ali_file\n";
    while(<TMPIN>){
	if(/^(\#|\%|>)/){           
	    ($acc)=split(/\s+/);
            $acc=substr($acc,1); 
	}
        else {
        $_=~s/\s+//g;
        $_=uc $_; 
        if($acc eq $seq_id){
	    $seq.=$_;
	}
        else{
           $seq_hash{$acc}.=$_;                      
        }
	}
    }
    close TMPIN;
 
   # unlink $ali_file; #fisa 13_08
  
   
    my $mask=&align_to_profile_make_mask(\%seq_hash);
   
    my $return_str=&align_to_profile_trim_by_mask($seq,$mask);

   

   
    foreach $acc(keys %seq_hash){
	delete $seq_hash{$acc};
    }
    return $return_str;
}




sub align_to_profile_make_mask{
    my $hash_ref=shift;
    my $acc;
  
    my $length=0;

    my %mask;
    foreach my $acc(keys %$hash_ref){
	my $seq=$hash_ref->{$acc};
        $length=length($seq);
        last;
    }
    my $i;
    for $i(0..$length-1){
	$mask{$i}=0;
    }

    foreach my $acc(keys %$hash_ref){
	my $seq=$hash_ref->{$acc};
        my @temp=split(//,$seq);
        for $i(0..$length-1){
	    if($temp[$i]=~/[A-Za-z]/){$mask{$i}=1;}
	}
    }

    my $mask;
    for $i(0..$length-1){
	$mask.=$mask{$i};
    }

    return $mask;

}

sub align_to_profile_trim_by_mask{
    my($seq,$mask)=@_;

    my $i;
    my $return_str;

    my @seq=split(//,$seq);
    my @mask=split(//,$mask);

    for $i(0..length($seq)-1) {
    if($mask[$i]){
    $return_str.=$seq[$i];
    } 

    }

    return $return_str;
}

############################end of align to profile


sub rRNA_gde2mltree_linux {
my ($class,$gde_file,$opt_o,$boot,$optimize)=@_;
my $gde;
my $tree;
     
   open (GDEIN, "$gde_file") || die "cannot open file $gde_file \n";
   while(<GDEIN>){
    $gde.=$_;
   }
   close GDEIN;

$gde_file=~s/\/$//;
my @t=split(/\//,$gde_file);
my $tmp_tag=$class->{tmp_path}.$t[-1].".".int(rand(1000)).substr(time(),5).".phyml";


my $tree_obj=GDE2PHY->mltree($gde,$class->{phyml_path},$tmp_tag,$boot,$optimize);

if($tree_obj){
$tree=$tree_obj->get_tree();
}

open (MLOUT,">$opt_o") || die "cannot create output file $opt_o \n";
print MLOUT $tree;
close MLOUT;

}



###private package

package GDE2PHY;


sub mltree{

    my ($class,$gde,$phyml_dir,$tmp_tag,$boot,$optimize)=@_;
    my $self={}; 
    if (!($gde=~/^(%|\#|>)/)){print "ERROR: Invalid GDE input for PHYML\n"; return;}
    $self->{seqinput}=$gde;
    $self->{phyml}=$phyml_dir.'phyml'; ##fisa 
  
    $self=&gde2input($self);
 
    if(!$self) {return;}
   
    my $tmp_seq_file=$tmp_tag;
    
    

    open(PMLOUT,">$tmp_seq_file") || die "cannot create temp file \n";
             
   print PMLOUT $self->{seqinput};
         
   close PMLOUT;
 
  
    

    my $stat=$tmp_seq_file."_phyml_stats.txt";
    my $tree=$tmp_seq_file."_phyml_tree.txt";
    my $lk=$tmp_seq_file."_phyml_lk.txt";
    my $rates=$tmp_seq_file."_phyml_rates.txt";

  
    my $c_phyml=$self->{phyml};
 
    #neva 12.01.11
    #Next line deprecated because PHYML's interface changes in version 3.0
    #$c_phyml.=" $tmp_seq_file 0 i 1 $boot HKY 4.0 e 1 1.0 BIONJ";
	

	$c_phyml.=" -i $tmp_seq_file -d nt -n 1 -b $boot -m HKY85 -t 4.0 -v e -c 1 -a 1.0";
	
		#fisa 19-08: spr instead of default nni
   ## $c_phyml.=" -i $tmp_seq_file -d nt -n 1 -b $boot -m HKY85 -t 4.0 -v e -c 1 -a 1.0 -s SPR";
    if($optimize){
    $c_phyml.=" -o tlr";
    }     
    else{
    $c_phyml.=" -o lr";
    }
    print $c_phyml."\n";

    system $c_phyml;

    unlink $stat;
    unlink $lk;
    unlink $tmp_seq_file; #fisa  commented out 13-08
    unlink $rates;

   
    open(TMPTREE,$tree) || die "cannot open temp tree file \n";
    while(<TMPTREE>){
	$self->{tree}.=$_;
    }
    close TMPTREE;
    unlink $tree;

    $self->{tree}=~s/\)0:/\):/g;
   
    bless $self,$class;

    return $self;

}






sub gde2input{
    my $class=shift;
    my $gde=$class->{seqinput};
    my @temp=split(/%|\#|>/,$gde);
    my %seq;
    my %name;
    my @order;
    my $id=0;
    my $sp_count;
    my $col_count;
    my $phy;

foreach my $temp(@temp){
    if (!($temp=~/\w/)) {next;}
    my @ln=split(/\n/,$temp);
    $id++;
    my $this_id='PD'.$id.'WU';
    $name{$this_id}=shift(@ln);
    $seq{$this_id}=join("",@ln);
    $seq{$this_id}=~s/\s+//g;
    if ($seq{$this_id}=~/^[^A-Za-z]+$/) {
	delete $name{$this_id};
        delete $seq{$this_id};
        next;     
##print "EXIT because of all-gap sequences\n"; 
##return;
}  
    if (!$col_count) {$col_count=length($seq{$this_id});}
    $sp_count++;
    push(@order,$this_id);
}

    if ($id <= 2) {print "EXIT: Less than two entries, cannot build trees\n"; return;}

$phy.=$sp_count." ".$col_count."\n";

    my $pos=0;

    foreach my $id(@order){
	my $this_id=$id."         ";
           $this_id=substr($this_id,0,10);
           $phy.=$this_id;
           $phy.=substr($seq{$id},$pos,50);
           $phy.="\n";  
    }
    $pos+=50;
    $phy.="\n";
    while ($pos <= $col_count){

    foreach my $id(@order){   
           $phy.=substr($seq{$id},$pos,50);
           $phy.="\n";  
	   }
    $pos+=50;
    $phy.="\n";
    }

    $phy=~s/\s+$//;

  

    $class->{seqinput}=$phy; 
    $class->{idhash}=\%name;

    return $class;
}



sub get_seqinput{
    my $class=shift;
    if ($class->{seqinput}) {return $class->{seqinput};}
    return;
	     }




sub get_tree{
    my $class=shift;
    my %name=%{$class->{idhash}};
    my $tree=$class->{tree};
     
    if (!$tree) {return $tree;}   

    while($tree=~m/(PD\d+WU)/g) {
     my $pos=pos $tree;
     my $this_id=$1;
     my $name=$name{$this_id};   
     $name=~s/\(//g;
     $name=~s/\)//g;
     $name=~s/;//g;
     $name=~s/,//g;
     $name=~s/://g;
    
 
    $pos-=length($this_id);
    
    substr($tree, $pos, length($this_id))=$name;


    }
    return $tree;
}











package DWU_PARSETREE;


sub new {
  my $tree_file = $_[1];
  my $tree_string = treeString($tree_file);
  my $self = {};
  $self->{tree_file} = $tree_file;
  $self->{tree} = parseTree ($tree_file);
 # $self->{bioperl-tree} = parseBioPerlTree($tree_string);  #fisa 11-07-14 used for rerooting
  bless $self;
  return $self;
}


sub treeString{
					my $tree_file = $_[0];
					open TREE, "< $tree_file" or die "Невозможно открыть файл дерева ".$tree_file."\n";
					# get a newick string from some source
 					my $tree_string;
 					my $t;
 					while($t = <TREE>){
 						$t =~ s/\n$//;
 						$tree_string = $tree_string.$t;
 					}
 					 close TREE;
					 
					 return $tree_string;
}


sub parseTree {
					my $tree_file = $_[0];
					my $tree_string = treeString($tree_file);

 					# Call class method parse from Bio::Phylo::IO
 					# note: newick parser returns 'Bio::Phylo::Forest'
                    # Call ->first to retrieve the first tree of the forest.
 					my $tree = Bio::Phylo::IO->parse(
  					  -string => $tree_string,
   					  -format => 'newick'
 					)->first;

 					return $tree;
	
}   


sub parseBioPerlTree{
	my $tree_string = $_[0];
	my $io = IO::String->new($tree_string);
	my $treeio = Bio::TreeIO->new(	-fh => $io,
									-format => 'newick');
	my $tree = $treeio->next_tree;
	return $tree;
}

sub end
{
    my $class=shift;
    undef $class->{tree};
		return 1;
}

sub check_taxon{
    my ($class,$taxon)=@_;
    my $tax = $class->{tree}->get_by_name($taxon);
   
    if ($tax) {return 'T';}
    else {return 0;}
}



### $obj->taxon_dist("taxon1")
### output a reference to an array of taxons with different distance to the inputtaxon
### [0] taxon_input
### [1] taxon of interest
### [2] number of nodes that seperate the two taxons
### [3] the distance between taxon_input and taxon of interest
### [4] the distance bwtween the nodes that taxon_input and taxon of interest directlly linked to
### [5] ID of the taxon groups that defined by taxon_input in a midpoint rooting tree
### [6] confidence of taxon [5]
  

#neva 11.03.11
sub taxon_dist {
	##
    my ($class, $query, $outgroup)=@_;
    my $tree = $class->{tree};
    my $query_node=$tree->get_by_name($query);

	#fisa 11-07-14 reroot
	my $new_root;
	my $bioperl_tree;
	my @bioperl_new_root;
	#my $debug = $tree->to_newick('-nodelabels' => 1);


	if (!defined $outgroup){
	my $query_node_length = $query_node->get_branch_length;
	$query_node->set_branch_length (0);
	$new_root= get_midpoint_bcv($tree); ##node closest to the midpoint of the longest path in the tree
	
	my $new_root_name=$new_root->get_name;
	$new_root->set_name("midpoint_".$new_root_name);
	$query_node->set_branch_length ($query_node_length);
	print STDOUT " new root name ".$new_root-> get_name."\n";
	my $string = $tree->to_newick('-nodelabels' => 1);
	print STDOUT " almost nothing changed ".$string."\n";
	
	$bioperl_tree = parseBioPerlTree($string);
	@bioperl_new_root = $bioperl_tree->find_node(-id => "midpoint_".$new_root_name);
	print STDOUT " midpoint \n";

	}
	else {
		my $string = $tree->to_newick('-nodelabels' => 1);
		$bioperl_tree = parseBioPerlTree($string);
		@bioperl_new_root = $bioperl_tree->find_node(-id => $outgroup);
		}

	if ($bioperl_new_root[0] -> depth != 0){

	$bioperl_tree->get_root_node->bootstrap("0.00000");
	$bioperl_tree->get_root_node->id("0.00000");
	my $is_zero = 0;
	if($bioperl_new_root[0]-> branch_length() == 0){
		$is_zero = 1;
		$bioperl_new_root[0]-> branch_length(0.2);
	}
	
	
	my $rooty = $bioperl_tree->reroot_at_midpoint($bioperl_new_root[0]);
	if ($is_zero == 1) {
	$rooty->branch_length(0);
	$bioperl_new_root[0] ->branch_length(0);
	}
	
	my $tree_file = $class -> {tree_file};
	#my $tmp_file = $tree_file.".tmp";
    my $out = Bio::TreeIO->new(-format => 'newick',
                           -file   => ">$tree_file");
	$out->write_tree($bioperl_tree);					   
	#$class ->{tree} = parseTree($tmp_file);	
    #$tree = $class ->{tree}; 
	$tree =  parseTree($tree_file);
	$query_node=$tree->get_by_name($query);
    }	
	#
    
	#set unique names for internal nodes
    my @internals = @{ $tree ->get_internals };
    my $id = 1;
    foreach my $internal_node(@internals){
    	$internal_node -> set_name($id."_".($internal_node ->get_name));
    	$id++;
    }
      
   
##get path from query to root (without query): nodes and distances
    my @path_query_to_root;
    my $temp_node = $query_node;
    my $distance_from_query;
    
    while(!($temp_node -> is_root)){
    	$distance_from_query += $temp_node -> get_branch_length;
    	$temp_node = $temp_node -> get_parent;
    	my @node_info = ($temp_node, $distance_from_query);
    	push @path_query_to_root, [@node_info];
    	
    }   
    
##
    
    my %dist_hash; 
    my $pivotal_node; ## a node on the path from query to root

    PIVOTAL: for (my $i = 0; $i < scalar @path_query_to_root; $i++ ){

    ## take info about nodes in path from query to root
    my @node_info = @{ $path_query_to_root[$i] };
    my $pivotal_node = $node_info[0];
   # print "node_name ".($pivotal_node -> get_name);
   # print "\t";
   # if ($pivotal_node -> is_root) {
   # 		print "root";
   # 	}
    my $pivotal_node_dist = $node_info[1];
    
    ## $distances[0] - numeric distance between this node and query leaf; $distances[1] - number of nodes between this node and query leaf.
    my @distances = ($pivotal_node_dist, $i, 0, 0); ## 
    $dist_hash{$pivotal_node} = \@distances;
    
  #   print "pivotal_dist ".$pivotal_node_dist."\n";
    
    ## take branch that doesn't lead to query
    my $new_pivotal_node = $pivotal_node ->get_child(1);
    if ($new_pivotal_node eq $path_query_to_root[$i-1][0] || $new_pivotal_node eq $query_node){
    	$new_pivotal_node = $pivotal_node ->get_child(0);
    	if ($new_pivotal_node eq $path_query_to_root[$i-1][0] || $new_pivotal_node eq $query_node){
    		$new_pivotal_node = $pivotal_node ->get_child(2);
    	}
    }
    
    $pivotal_node_dist += $new_pivotal_node-> get_branch_length;

    
    
    $new_pivotal_node->visit_depth_first( -pre => sub{ 
    	
    	 my $node = shift;
    	
    	 if ($pivotal_node ->is_root){
    	 	last;
    	 }
    	 
    	 my $numeric_dist = $node -> get_branch_length; # numeric distance between this node and the query leaf
    	 $numeric_dist += $dist_hash{$node ->get_parent}[0];
    	 my $node_dist = $dist_hash{$node -> get_parent}[1] + 1; # number of nodes between this node and the query leaf.
    	 my $numeric_dist_between_nodes = $numeric_dist - ($query_node -> get_branch_length) - ($node -> get_branch_length); # Distance between the nodes that the query and this node are attached to; no sense for internal nodes
    	

    	 my @distances = ($numeric_dist, $node_dist, $numeric_dist_between_nodes, $i+1); # $i+1 - The branching position of this node along the path that connects the query and the outgroup 
    	 $dist_hash{$node} = \@distances;

    	 },
    	
		#a little bit of magic to stop bio::phylo walking away from the subtree
    -pre_sister     => sub {
    	 my $node2 = $_[1]; ##current_ sister
    	 my $node1 = ${$new_pivotal_node->get_sisters}[1]; ##sister
    	 if ($node2->get_name eq $node1->get_name){ ## || $pivotal_node ->is_root was here
    	 	
    	 	next;
    	 	
    	 }
    },
    -with_relatives => 1 
    	 );
    	 

    }

    #calculate distances for the outgroup(s)
       my $root = $tree->get_root;
    my @outsisters = @{ $root ->get_child(0) -> get_sisters };
	foreach my $outsister(@outsisters){
	if (!exists $dist_hash{$outsister}){
		my $new_pivotal_node = $outsister;
		 $new_pivotal_node->visit_depth_first( -pre => sub{ 
    	
    	 my $node = shift;
    	 

    	 my $numeric_dist = $node -> get_branch_length; # numeric distance between this node and the query leaf
    	 $numeric_dist += $dist_hash{$node ->get_parent}[0];
    	 my $node_dist = $dist_hash{$node -> get_parent}[1] + 1; # number of nodes between this node and the query leaf.
    	 my $numeric_dist_between_nodes = $numeric_dist - ($query_node -> get_branch_length) - ($node -> get_branch_length); # Distance between the nodes that the query and this node are attached to; no sense for internal nodes
    	

    	 my @distances = ($numeric_dist, $node_dist, $numeric_dist_between_nodes, (scalar @path_query_to_root)); # scalar $path_query_to_root - The branching position of this node along the path that connects the query and the outgroup 

    	 $dist_hash{$node} = \@distances;
    	
    	 },
    	
    -pre_sister     => sub {
    	 my $node2 = $_[1]; ##current_ sister
    	 if (exists $dist_hash{$node2}){
    	 	next;	
    	 }
    },
    -with_relatives => 1 
    	 );
	}
	
	
	
	#if ($outsister -> is_terminal){
   # my @distances = ($dist_hash{$root}[0]+($outsister -> get_branch_length), 
    #                 $dist_hash{$root}[1]+1,
   #                  $dist_hash{$root}[0] - ($query_node -> get_branch_length),
    #                 scalar @path_query_to_root);
    #$dist_hash{$outsister} = \@distances;  
	# }
	 
    }	

    #sort 
    my $dist_list;	
    my @leaves = @{ $tree -> get_terminals };
    my %hash_sort_branching_position;
    my %hash_sort_node_dist;
    my %hash_sort_taxon_dist;


    foreach my $leaf(@leaves){
    if ($leaf -> get_name ne $query_node -> get_name){
    	my $str  = ($query_node -> get_name)."\t";
    	   $str .= $leaf -> get_name."\t";
    	   my @distances = @{$dist_hash{$leaf}};

    	   $str .= $distances[0]."\t"; # Numeric distance between this node and the query leaf
    	   $str .= $distances[2]."\t"; # Numeric distance between the nodes that the query and this node are attached to
    	   $str .= $distances[1]."\t"; # Number of nodes between this node and the query leaf
    	   $str .= $distances[3]."\t"; # The branching position of this node along the path that connects the query and the outgroup   
		  my @confidence = split /_/, $path_query_to_root[$distances[3]-1][0] ->get_name;
		   my $conf = $confidence[1];
		   if (!$conf =~ /^[\d\.]+$/){ #т.е. он находится по другую сторону от корня, либо был переименован в процессе поиска корня.
			 my $debugger = 0;

			my $tnode = $leaf;
			while(!$tnode->is_root){
				$conf = $tnode->get_name;
				$tnode = $tnode ->get_parent;

			}

			if ($conf =~ /_/){
				my @a = split /_/, $conf;
				$conf = $a[-1];

			}

		   }
		   $str .= $conf;
    	   
    	$hash_sort_branching_position{$str} = $distances[3];
        $hash_sort_node_dist{$str}=$distances[2];
        $hash_sort_taxon_dist{$str}=$distances[0]; 

       push @$dist_list,$str;

    }
    	
    }

    @$dist_list=sort {$hash_sort_branching_position{$a} <=> $hash_sort_branching_position{$b} 
    	               or $hash_sort_node_dist{$a} <=> $hash_sort_node_dist{$b} 
    	               or $hash_sort_taxon_dist{$a} <=> $hash_sort_taxon_dist{$b}        } @$dist_list;

return $dist_list;
    
}



###########start of sub get_intersect_nodes##############


### $obj->get_intersect_nodes(\@array)
### output a reference to an array
### $return_ref="node1"."\t"."node2";
 
  
sub get_intersect_nodes{
    my ($class,$in_ref)=@_;
    my @id;
    my $id;
    my %id;
    my $in_taxa;

    my $return_txt;

### get all the ids for the taxa
    foreach $in_taxa (@$in_ref){
	if($class->{tn_id}->{$in_taxa}){
	    push(@id,$class->{tn_id}->{$in_taxa});
	}
    }

    if(!(@id)){return;}
    if (!(@id >=2 )){return;}



  

####################

    
 
   
  
    my $original_matrix=&copy_matrix($class->{matrix});

#########first clean the NaN branches
    my $matrix_OK=0;
    my $all_clusters=&matrix_to_cluster($class->{matrix});
    foreach my $this_cluster(@$all_clusters){
	my $lineup=join(" ",@$this_cluster);
        my $this_count=0;
        foreach $id(@id){
        if($lineup=~/\b$id\b/){
	    $this_count++;
         } 
	}

	if($this_count == @id){
	    my $this_matrix=&cluster2matrix($class->{matrix},$this_cluster);
	    $class->{matrix}=$this_matrix;
            $matrix_OK=1;
	    last;
        }

    }
#########end

###########???? must cover all the path redo the path thing  (a-b a-c b-c)

    my %path_start_end;
 
    for my $i(@id){
	for my $j(@id){
	    if($i eq $j) {next;}
            my ($s,$e)=sort($i,$j);
            my $s_e=$s."\t".$e;
            $path_start_end{$s_e}=1;
	}
    }


## find the nodes in common among the pathes

my %node_in_between;
my $path_number; 
    foreach my $key(keys %path_start_end){
	my($taxon_a,$taxon_b)=split(/\t/,$key);
        my $path=&path_between_two_taxa($class,$taxon_a,$taxon_b);
  
        my @t=split(/\t/,$path);
        $path_number++;
        foreach my $t(@t){
       if($t=~/^N/){$node_in_between{$t}++;} 
       }  
    }





    my @node_to_be;
    foreach my $node(keys %node_in_between){
   
	if($node_in_between{$node} == $path_number){
        push(@node_to_be,$node); 
	}
    }

    return \@node_to_be;

}


################end of sub break_to_clusters######################


###### private #################


sub get_midpoint_dwu{ ##fisa renamed get_midpoint into get_midpoint_dwu to prevent confusion (we use get_midpoint from bio::phylo)
    my $class=shift;
  
    my $max_dist=0;
    my $max_path='';
    
    foreach my $i(keys %{$class->{matrix}}){
	my $dist_ref=&taxon_path($class,$i);
	foreach my $j(keys %$dist_ref){
        if ($dist_ref->{$j}->{taxon_dist} > $max_dist){
            $max_dist=$dist_ref->{$j}->{taxon_dist};
            $max_path=$j;
          }  
	}
    }

    my @mark=split(/\t/,$max_path);
    my $start=shift(@mark);
    my %mark;
   
    my $first=$start;
    $mark{$first}=0;
    my $length=0;
    foreach my $i(@mark){
	$mark{$i}=$mark{$first}+$class->{matrix}->{$first}->{$i};
        $length=$mark{$i};
        $first=$i;   
    }
    unshift(@mark,$start);
    my $id=0;
    my $mid_len=sprintf("%.5f",$length/2);
    
    foreach my $i(@mark){
	if ($mark{$i} > $mid_len){
        last;  
	}
	$id++;
    }  

    my($mid_a,$mid_b)=($mark[$id-1],$mark[$id]);

    return ($mid_a,$mid_b);

}

## fisa 12_08 Slightly modified method from bio:: phylo
    sub get_midpoint_bcv {
        my $tree     = shift;
        my $root     = $tree->get_root;
        my $nA       = $tree->get_tallest_tip;
        my $nB       = $nA->get_farthest_node(1); #fisa - optional true value => based on patristic distance instead of nodal (as in bio phylo get_midpoint)
        my $A2B_dist = $nA->calc_path_to_root + $nB->calc_path_to_root;
        my $outgroup = $nA;
        my $middist  = $A2B_dist / 2;
        my $cdist    = 0;
        my $current  = $nA;
        while ($current) {

            if ( $cdist > $middist ) {
                last;
            }
            else {
                if ( my $parent = $current->get_parent ) {
                    $cdist += $current->get_branch_length;
                    $current = $parent;
                }
                else {
                    last;
                }
            }
        }
        return $current;
    }


sub taxon_path {
    my ($class,$taxon)=@_;

    my %include;
    my %current;
    my %path;

    my $go_on=1;  

    $include{$taxon}=1;
    $current{$taxon}=1;    
    $path{$taxon}=$taxon;

    while($go_on){

        $go_on=0;

	foreach my $this_tn(keys %current){
	    if ($current{$this_tn}){

		foreach my $i(keys %{$class->{matrix}->{$this_tn}}){
		    if ($class->{matrix}->{$this_tn}->{$i} != -99999){
			if ($include{$i}) {next;}
                        $path{$i}=$path{$this_tn}."\t"."$i";
                        $current{$i}=1;
                        $include{$i}=1;
                        if ($i=~/^N\d+/) {$go_on++;}  
                    }
		}
                }    
		$current{$this_tn}=0;             
        }
    }

    my $T2Tpath;

    foreach my $i(keys %path){
    my @temp=split(/\t/,$path{$i});

    if (($temp[0]=~/^T\d+/) && ($temp[-1]=~/^T\d+/)){
	if ($temp[0] eq $temp[-1]) {next;}
        my $tn_number=@temp;
        my $node_number_dist=$tn_number-2;
        my $taxon_dist;
        my $node_dist;
 
        my $first=shift(@temp);
        my $second;
        my @dist;
        my $from=$first;
        my $to;

	while(@temp){
	    $second=shift(@temp);
            $to=$second;
            my $this_dist=$class->{matrix}->{$first}->{$second};
             
            push(@dist, $this_dist);
            $taxon_dist+=$this_dist;
            $first=$second;
        } 
        
        $node_dist=$taxon_dist-$dist[0]-$dist[-1];
        
        $node_dist=sprintf("%.05f",$node_dist);
        $taxon_dist=sprintf("%.5f",$taxon_dist); 
        
             $T2Tpath->{$path{$i}}->{node_number_dist}=$node_number_dist; 
             $T2Tpath->{$path{$i}}->{taxon_dist}=$taxon_dist;
             $T2Tpath->{$path{$i}}->{node_dist}=$node_dist;
         
    }
}
    return $T2Tpath;

}



sub path_between_two_taxa{
   my ($class,$taxon_a,$taxon_b)=@_;
   my $matrix_ref=$class->{matrix};
    my %include;
    my %current;
    my %path;

    my $go_on=1;  

    $include{$taxon_a}=1;
    $current{$taxon_a}=1;    
    $path{$taxon_a}=$taxon_a;

    while($go_on){

        $go_on=0;

	foreach my $this_tn(keys %current){
	    if ($current{$this_tn}){

		foreach my $i(keys %{$matrix_ref->{$this_tn}}){
		    if ($matrix_ref->{$this_tn}->{$i} != -99999) {
			if ($include{$i}) {next;}
                        $path{$i}=$path{$this_tn}."\t"."$i";
                        $current{$i}=1;
                        $include{$i}=1;
                        if ($i=~/^N\d+/) {$go_on++;}
                        if ($i eq $taxon_b) {$go_on=0;}  
                    }
		}
                }    
		$current{$this_tn}=0;             
        }
    }

   my $return_path;
   foreach my $path(keys %path){
       if($path{$path}=~/^($taxon_a).+($taxon_b)$/){
       $return_path=$path{$path};
       last; 
   }
   }

   return $return_path; 

}




### node_number($tree)
### return: 0 if the tree is not in the right format
###         non 0 for the number of nodes in the tree

sub node_number
{
    my $tree=shift;

    my $left=0;
    my $right=0;
   
    if (!( ($tree=~/^\(/) &&($tree=~/\)$/) )) {return 0;}
    
    while($tree=~/\(/g){
	$left++;
     
    }
    while($tree=~/\)/g){
	$right++;
  
    }

    if ($left == $right) {return $left;}
    else {return 0;}
}


### all_taxons($tree)
### return ref to an array of taxons
#neva 11.03.11
sub all_taxons{
    my $tree=shift;
    $tree=~s/\)-?\d+\.?\d*:(-?\d+\.?\d*)/\):$1/g; # if support of clades was calculated
    $tree=~s/,/ /g;
    $tree=~s/\(/ /g;
    $tree=~s/\)/ /g;

    $tree=~s/^\s+//;
    $tree=~s/\s+$//;

    my @taxons;
    my @temp=split(/\s+/,$tree);
    foreach my $taxon(@temp){
	if ($taxon=~/^:/) {next;}
        else{$taxon=~s/:(-*\d*\.*\d*$|\s*NaN$)//;push(@taxons,$taxon);}
    }   
    return \@taxons;
}

### matrix_to_cluster
### input a hash reference 
### output an reference to an array of arrays devide the elements in the matrix into clusters
## based on single linkage clustering algorithm
##
sub matrix_to_cluster
{
    my $matrix_ref=shift;
    my $i;
    my $j;
    my %pos;
    my %include;
    my %all;

    
    my $id=0;

    foreach $i(keys %$matrix_ref){
        $all{$i}=1; 
	foreach $j(keys %{$matrix_ref->{$i}}){
	    my $dist=$matrix_ref->{$i}->{$j};

            $id++;
            if($dist=~/^N/) {next;}             
            if ($dist == -99999) {next;}              
            
	    my $old_i_id=$pos{$i};
            my $old_j_id=$pos{$j};
        
            if (!$old_i_id) {
		$pos{$i}=$id;
            }
            else {
		foreach my $key(keys %pos){
		    if ( $pos{$key} == $old_i_id ) {$pos{$key}=$id;} 
		}
	    }

            if (!$old_j_id) {
		$pos{$j}=$id;
            }
            else {
		foreach my $key(keys %pos){
		    if ( $pos{$key} == $old_j_id ) {$pos{$key}=$id;} 
		}
	    }
       }  
    }    

    my %cluster;

    foreach my $key(keys %pos) {
	$id=$pos{$key};
        $include{$key}=1; 
        push(@{$cluster{$id}},$key);
    }

    my @cluster;

    foreach my $key(keys %cluster){
	push(@cluster, $cluster{$key});
    }

    foreach my $key(keys %all){
    
	if (!$include{$key}){
        my @single=($key);
          push(@cluster,\@single); 
	}

}   
    


    return \@cluster;
}


## input reference to master matrix and reference to an array of taxon and nodes
## output a reference to the portion of matrix that only include the taxons and nodes on the array

## note: the master matrix is the unaltered matrix  
##       the array must came from matrix_to_cluster sub routine


sub cluster2matrix {
    my ($master_ref, $array_ref)=@_;
    my $sub_ref;

    my %tx_include;

    foreach my $i(@$array_ref){
	$tx_include{$i}=1;
    }

## get the value and the matrix
    foreach my $i(keys %$master_ref){
	foreach my $j(keys %{$master_ref->{$i}}){
            if ($tx_include{$i} && $tx_include{$j}){
	    $sub_ref->{$i}->{$j}=$master_ref->{$i}->{$j};
	}
	}
    }

   
## get ride of the hanging nodes
    
    my %false_taxon;
    my %real_taxon;

    foreach my $i(keys %$sub_ref){
    my $link=0;

    foreach my $j(keys %{$sub_ref->{$i}}){
	if ($sub_ref->{$i}->{$j} != -99999) {$link++;}
}   
    if ($link == 1) {$false_taxon{$i}=1;}

    }

    foreach my $i(keys %false_taxon){
   
    my $link=0;

    foreach my $j(keys %{$master_ref->{$i}}){
	if ($master_ref->{$i}->{$j} != -99999) {$link++;}
     }

   
    if ($link == 1) {
    $false_taxon{$i}=0;
    $real_taxon{$i}=1;
   }

    }

    my $temp_ref;

    foreach my $i(keys %$sub_ref){
	foreach my $j(keys %{$sub_ref->{$i}}){
	    if ($false_taxon{$i} || $false_taxon{$j}) {next;}
            else {
	     
            $temp_ref->{$i}->{$j}=$sub_ref->{$i}->{$j};   
            }
	}
    }
    $sub_ref=$temp_ref;
    undef $temp_ref;
          return $sub_ref;    
}


sub copy_matrix{
    my $old_max=shift;
    my $max;
    foreach my $i(keys %$old_max){
	foreach my $j(keys %{$old_max->{$i}}){
	    $max->{$i}->{$j}=$old_max->{$i}->{$j};
	}
    }
    return $max;

}


sub get_immediate_nodes_for_a_taxon{
    my ($class,$in_taxon)=@_;
    my $id=$class->{tn_id}->{$in_taxon};
    my $current=$id;;
    foreach my $key(keys %{$class->{matrix}->{$id}}){
    if((length($class->{matrix}->{$id}->{$key})>=1)&&($class->{matrix}->{$id}->{$key} != -99999)){
    $current=$key;
    last; 
     }    
    }
    return $current;
}




package DNATRIM;


sub trim_all_gap{
    my ($class,$gde0)=@_;
    my %seq;
    my $acc;

    my @t=split(/\n/,$gde0);
    foreach my $line(@t){
	if($line=~/^>|\#|\%/){
	    ($acc)=split(/\s+/,$line);
	}  
        else{
        $line=uc $line;
        $line=~s/\s+//g;
        $seq{$acc}.=$line; 
	}
    }
    
    my $len;
    foreach $acc(keys %seq){
	$len=length($seq{$acc});
       
    }
  
    my $mask="0"x$len;
   
    my @mask=split(//,$mask);

    foreach $acc(keys %seq){
	my @seq=split(//,$seq{$acc});
        for my $i(0..$len-1){
	    if ($seq[$i]=~/[A-Z]/) {$mask[$i]=1;}
	}
    }

    my $gde1;

       foreach $acc(keys %seq){
	   $gde1.=$acc."\n";
	my @seq=split(//,$seq{$acc});
           for my $i(0..$len-1){
	       if($mask[$i]){$gde1.=$seq[$i];
	}
       }
	   $gde1.="\n";
           delete $seq{$acc};
}

 
    undef @mask;
    undef $mask;
   

    return $gde1;

}

sub trim_gap_to_reference{
    my ($class,$gde0,$opt_t,$opt_G,$opt_e)=@_;
  
    my $gde1;
    my %seq;
    my @line=split(/\n/,$gde0);
    my $id;
    my @id;
   
    my $symbol;

    while(@line) {
	my $line=shift(@line);
        if ($line=~/^(\%|\#|>)/){
	    if(!$symbol){$symbol=substr($line,0,1);}
           ($id)=split(/\s+/,$line);
	   $id=substr($id,1);
	    if($id=~/\w+/){
           push(@id,$id); 
           }
	    next;
        }
        else {
            $line=~s/\s+//g;
            if($line){
	    $seq{$id}.=$line;
	}
        }
    }

   

    if (($opt_t) && ($seq{$opt_t}))  {       
       if ((!(length($opt_G) > 0) ) || (!($opt_G >=0))) {$opt_G=10;}     
       foreach my $key(keys %seq){
       $seq{$key}=~s/\s+//g; 
       }
     
       my $seq_t=$seq{$opt_t};
       my $length_t=length($seq_t);

       my $mask='1'x$length_t;
      



       while($seq_t=~/([^ATGCUatgcu]+)/g) {
	   my $pos2=pos($seq_t);
           my $gap=length($1);
           my $pos1=$pos2-$gap;
           
           if ($gap >= $opt_G) {
           substr($mask,$pos1,$gap)='0'x$gap;  
           }   
       }



       if($opt_e){
          
           my $start=0;
           my $end=0;
	   while($mask=~/^(0+)/g){
	       $start=length($1);
               last;
	   }
           while($mask=~/(0+)$/g){
               $end=length($1);
               last;
           } 

	   $mask=~tr/0/1/;
           if (($start >= 1) && ($start<=length($mask))){
           substr($mask,0,$start)='0'x$start; 
	   }

           if (($end >= 1) && ($end<=length($mask))){
           substr($mask,length($mask)-$end, $end)='0'x$start; 
	   }
       }
    


       my @keep=split(//,$mask);

      foreach my $key(keys %seq){
      my @this_seq=split(//,$seq{$key});
      my $this_seq='';
      my $i=0;

      foreach (@this_seq){
	  if ($keep[$i] >= 1) {
	      $this_seq.=$this_seq[$i];
          }
	  $i++;
      }
      $seq{$key}=$this_seq;
     } 

    
       foreach my $key (@id){
         
	   $gde1.=$symbol.$key."\n";
           $gde1.=$seq{$key}."\n";
       }
       
       return $gde1;
        
    }
    else{
	return $gde0;
    }
  }




sub trim_gde_by_alignment_score  {
    my ($class,$gde,$matrix_type,$cutoff,$gap_cutoff)=@_; ##opt_e trim from ends only (on)

 
   my $self={};
    $self->{ori_gde}=$gde;

    

    if (!$matrix_type) {$matrix_type='IUB';}
    if (!(length($cutoff) >= 1)) {$cutoff=0;}   ## column score cutoff 0
    if (!(length($gap_cutoff) >= 1)) {$gap_cutoff=100;} ## trim column with above 100% gaps
  
 

 
    $self->{cutoff}=$cutoff;
    $self->{gap_cutoff}=$gap_cutoff;
    $self=&get_matrix($self,$matrix_type);

     
    $self=&transform_seq($self);
    my $return_trim_ali=$self->{trim};

    foreach my $key(keys %$self){
	delete $self->{$key};
    }
    undef $self;

    return $return_trim_ali;
}



sub calscore{
    my ($class,$col)=@_;
    my $score;
    my %matrix=%{$class->{matrix}};
    my @symble=@{$class->{symble}};
    my %check;
    foreach my $key(@symble){
	$check{$key}=1;
    }
    $check{'-'}=1; 
    

    my @order=split(//,$col);
    foreach my $key(@order){
	if (!$check{$key}) {$key='N';}
    }
   
    my %freq;
    my $non_gap;
    my $total=length($col);

    foreach my $temp(@order){
	if ($temp ne '-'){$freq{$temp}++;
                          $non_gap++;}
    }

    if (!$non_gap) {return -1;}
    if (!$total) {return -1;}
    my $ratio=sprintf("%.2f",$non_gap*100/$total);
    if ($ratio < (100-$class->{gap_cutoff})) {return -1;} 


   
    my %con;   ## consensuse on r dimensions
    
    foreach my $d(@symble){
	$con{$d}=0;
        foreach my $i(keys %freq){
	    my $pair=$d.$i;
	 
        $con{$d}+= $matrix{$pair} * $freq{$i}; 
        }
       
	$con{$d}/=$total;
    }
    
    my @dist;

   foreach my $a(@order){
## a residue came in
       if ($a eq '-') {next;}
       
       my $dist;
 
       foreach my $d(@symble){
       my $pair=$a.$d;
       my $std=$matrix{$pair};
       my $diff=$con{$d}-$std;
            
          $dist+=$diff*$diff;            
       }
     
       $dist=sqrt($dist);
       push(@dist,$dist);
   } 
   
  
    my $num;
    foreach my $dist(@dist){
 	$score+=$dist;
      $num++; 
    }
    $score/=$num;

  

    $score=exp(-$score/1.5)*100*$non_gap/$total;
 
   $score=sprintf("%.0f",$score);

    return $score;
}



sub transform_seq{
    my $class=shift;
    my $gde=$class->{ori_gde};

 
    if($class->{cutoff}==0 && $class->{gap_cutoff}==100){
    $class->{trim}=$gde;
    return $class; 
    }
    
    
    $gde=~s/^\s+//;
    my $symbol=substr($gde,0,1);
    my @temp=split(/%|\#|>/,$gde);
    my @order;
    my %seq;
    my $len;

    foreach my $temp(@temp){
        if(!($temp=~/\w+/)){next;}
	my @ln=split(/\n/,$temp);
        my $seqname=shift(@ln);
       
        my $seq=join("",@ln);
	$seq=~s/\s+//g;
        $seq=uc($seq);
        push(@order,$seqname);
        if (!$len) {$len=length($seq);}
		    my @seq=split(//,$seq);
		    $seq{$seqname}=\@seq;
        
    }
    my @col;

    my $i;
    for($i=0;$i<$len;$i++){
	my $col='';
	foreach my $seqname(@order){
	    $col.=$seq{$seqname}->[$i];
	}
	push (@col,$col);

    }

    my @local_score;   
    foreach my $col(@col){
	my $score=&calscore($class,$col);
	if (!$score){$score=0;}
	push(@local_score,$score);
    }

   my $cutoff=$class->{cutoff};    
   my $mask=&mask(\@local_score,$cutoff);
     

   my @mask=split(//,$mask);

   my $trim_gde; 

  
   foreach my $id(@order) {
     
       $trim_gde.=$symbol.$id."\n";
       my $trim_seq;
   my @seq=@{$seq{$id}};
       my $i=0;
       while(@seq){
       my $seq=shift(@seq);
       if ($mask[$i]) {$trim_seq.=$seq;}
       $i++; 
       }

   for($i=0;$i<length($trim_seq);$i+=80){
       my $str=substr($trim_seq,$i,80);
       $trim_gde.=$str."\n";
   }
     
  } 
   
    $class->{trim}=$trim_gde;

    return $class;
}
   

sub mask{
    my ($ref,$cutoff)=@_;
    my @score=@$ref;
    my %local_score;
    my %column_score;
    my $mask;
    my $seqlength=@score;

    for my $i(0..$seqlength-1){
        if (!$score[$i]){$score[$i]=1;}
	$column_score{$i}=$score[$i];
    }

 

    for my $i ( 0 .. $seqlength - 1 ) {

	if ($column_score{$i} == -1) {$mask.=0;next;}

    if ($i<=2 or $i>= $seqlength-3) {
			$local_score{$i} = $column_score{$i};
		}
    else {
 $local_score{$i} = ($column_score{$i-3} + 2*$column_score{$i-2} + 3*$column_score{$i-1} + 4*$column_score{$i} + 3*$column_score{$i+1} + 2*$column_score{$i+2} + 1*$column_score{$i+3})/16;
		}
		if ($local_score{$i}/$column_score{$i} > 3 ) {
  my $score_l = $column_score{$i-3} + $column_score{$i-2} + $column_score{$i-1};
  my $score_r = $column_score{$i+1} + $column_score{$i+2} + $column_score{$i+3};

  if($score_r && $score_l){
		if ($score_l/$score_r > 20 or $score_l/$score_r < 0.05) {
				$local_score{$i} = $column_score{$i};
			}
	    }
 

		}
              
             
####	print "###########$cutoff \n";
		if ($local_score{$i} >= $cutoff) {
			$mask .= "1";
		}
		else {
			$mask .= "0";
		}
				
    
	}

   
       return $mask;



}


sub get_matrix{
    my ($class,$type)=@_;    
    my $matrix="nt_order = \"ABCDGHKMNRSTUVWXY\";

MATRIX IUB[]={
1.9,
0,  1.9,
0,  0,  1.9,
0,  0,  0,  1.9,
0,  0,  0,  0,  1.9,
0,  0,  0,  0,  0,  1.9,
0,  0,  0,  0,  0,  0,  1.9,
0,  0,  0,  0,  0,  0,  0,  1.9,
0,  1.9,0,  1.9,0,  1.9,1.9,1.9,1.9,
0,  0,  0,  0,  0,  0,  0,  0,  1.9,1.9,
0,  0,  0,  0,  0,  0,  0,  0,  1.9,0,  1.9,
0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1.9,
0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1.9,
0,  0,  0,  0,  0,  0,  0,  0,  1.9,0,  0,  0,  0,  1.9,
0,  0,  0,  0,  0,  0,  0,  0,  1.9,0,  0,  0,  0,  0,  1.9,
0,  1.9,0,  1.9,0,  1.9,1.9,1.9,1.9,1.9,1.9,0,  0,1.9,1.9,1.9,
0,  0,  0,  0,  0,  0,  0,  0,  1.9,0,  0,  0,  0,  0,  0,  0,  1.9,";


    if ($type && ($type eq 'CLUSTALW')) {
$matrix="nt_order = \"ABCDGHKMNRSTUVWXY\";

MATRIX clustalvdnamt[]={
  1,
  0,  0,
  0,  0,  1,
  0,  0,  0,  0,
  0,  0,  0,  0,  1,
  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,";
}



    my @nt;
    my %hash;
    my $i;
    my $j;
    my @temp=split(/\n/,$matrix);
    my @order;

    foreach my $temp(@temp){
	if ($temp=~/.*=.*\"([A-Z]+).*;/){
        @order=split(//,$1);
        next; 
        }
        if ($temp=~/\]=\{/) {$i=0;next;}
        if ($temp=~/\d*.*\d+,/){
        $temp=~s/\s+//g;
        my @score=split(/,/,$temp);
    	   for($j=0;$j<=$i;$j++){
           my $pair=$order[$i].$order[$j];
           $hash{$pair}=$score[$j];           
           $pair=$order[$j].$order[$i];
           $hash{$pair}=$score[$j];
	   }
	    $i++;
            next;
        }

        } 

       
	$class->{matrix}=\%hash;
    $class->{symble}=\@order;
        return $class; 

}






1;
