#!/usr/bin/perl
    use XML::TreePP;
    use Data::Dumper;
    use Bio::Phylo::IO;
    use Bio::Phylo::Forest::TreeRole;
    use Bio::Phylo::Taxa::Taxon;
    use Bio::Phylo::Factory;
    use DBI;
    use Bio::Phylo::IO 'parse';
    use XML::XSLT;
    use Path::Class qw(dir);
    use Bio::Phylo::Treedrawer;
    use XML::LibXSLT;
    use XML::LibXML;
    use Class::Struct;
	use Bio::Phylo::Util::CONSTANT 'looks_like_object';



    
  
  ## INPUT:   $dir - path to directory with STAP output folders
  #           (option) $chrom_list - path to file containing a list of chromatogram names (these must match the names of STAP output folders)
  #           by default this script will try to process all folders in $dir, assuming their names and contents comply with the aforementioned rules
  # STAP output folder X must contain X.ab1.cluster.fasta, X.mltree2.tab, X.mltree2.with_taxonomy, X.results;
  # this script produces X.mltree2.with_taxonomy_and_names (newick format), X.mltree2.svg (image) in $dir.X, 
  # and also results.xml in the X folder.
    
  my ($dir, $chrom_list) = @ARGV;
  my @chrom_array;
  my $dbh = DBI->connect('dbi:mysql:greengenes','weidewind','lotus82') or warn "Connection Error: $DBI::errstr\n";


   
  ## check if input directories exist 
    if (!-e $dir) {
    	print STDERR "No such directory: ".$dir."\n";
    	print STDERR "Exiting";
    	exit;
    }

   #print "_1_\n";
   
    if (!defined $chrom_list){
    	opendir DIR, $dir or die "Can't open directory $dir: $!";
    	
        while (my $folder = readdir(DIR)) {
        	if (-d "$dir/$folder"){
        		if (($folder ne ".") & ($folder ne "..")){	
        		push @chrom_array, $folder;
        		print $folder."\t";
        		}
       		}
        }
        closedir DIR;
    }
	
    else { 
		if(-e $chrom_list){
			open CHROMLIST, $chrom_list or die "Can't open file $chrom_list: $!";
			while (<DIRLIST>){
				my $folder = $_;
				$folder =~ s/^\s+//;
				$folder =~ s/\s+$//;
				push @chrom_array, $folder;
			}
    	
		}
		else {
			print STDERR "No such file: ".$chrom_list."\n";
			print STDERR "Exiting";
			exit;
		}
	}
	
    #print "_2_\n";
    ## 
 
    
	## assume that we are in standard BCV project directory, 
	## so we can get chromatogram length from  Input/X/ps_prj/poly_dir/X.ab1.bqs

	my $mydir = dir($dir);
	my $input_dir = $mydir->parent."/Input/"; 
	$dir = $dir."/";  
    print 	$input_dir."\n";

	
	## form xml structure: collect info for each chromatogram
	foreach my $chrom_name(@chrom_array){

    my $current_dir = $dir.$chrom_name."/";
    my $chrom_length;
    
    my $aggregate;
    my @chrom_stap_taxes;
    my @chrom_blast_taxes;
    
    ##
 
    my $clust_file =  $current_dir.$chrom_name.".ab1.cluster.fasta";
    open clusters, "< $clust_file" or die "Невозможно открыть файл ".$clust_file;

    while (my $row = <clusters>) {
	           	       my @splitter = split(/\s+/, $row);
            	       $row = <clusters>;
            	       $_ = $row;
					     my $exp = $splitter[1];
            	        $exp=~ s/\n$//;
                        my $name1 = (substr $splitter[0], 1)."_".(sprintf("%.3f",$exp));
           	      
				  
				  if ($exp >= $expectation_threshold && -e $current_dir.$name1.".results"){
				  
            	       my $tree_file = $current_dir.$name1.".mltree2.with_taxonomy";
					   print $tree_file."\n";
                    
 
 
 
                       ## print newick tree file with prokmsa names

 						my $tree_string = treeAsString($tree_file);
 						$tree_string =~ s/PROK(\d+)(\w+):/PROK.$1.$2."_".clean_string(get_prokmsa_name($1)).":"/ge;
						open NEW_TREE, "> $tree_file.with_names" or die "Невозможно открыть файл ".$tree_file.".with_names\n";
 						print NEW_TREE $tree_string;
						close NEW_TREE;
 						
 					  ##
 					
 					## draw tree with prokmsa names
 					my $treedrawer = Bio::Phylo::Treedrawer->new(
    					-width  => 2200,
    					-height => 2000,
    					-shape  => 'rect', 
   						-mode   => 'clado', 
    					-format => 'SVG'
 					);
					
					my $test_tree = parseStringTree($tree_string);
					print $test_tree;
#					if ( looks_like_object( $test_tree, _TREE_ ) ) {
#      print "ookay";
 #}

 					$treedrawer->set_tree($test_tree);
 					$treedrawer->set_text_width(1500);

                    my $svg_file =   $current_dir.$name1.".mltree2.svg";
					open FIG, "> $svg_file";
 					print FIG $treedrawer->draw;
 					close FIG;	
 				    ##	
            	}
				}      
	
	      

             close clusters;
} 


   
      
sub aggregate_tree{
	my ($stap_taxonomies, $blast_taxonomies)=@_;
     
     my $tax_tree2=TaxNode->new;
     $tax_tree2->stap_count(0);   
     $tax_tree2->blast_count(0);  
     $tax_tree2->stap_children({});
     $tax_tree2->blast_children({});
     
     my @taxonomy;
     
     for my $tax_string(@{$stap_taxonomies}){
     	@taxonomy = split (/\|/, $tax_string); 
     	build_taxtree_2($tax_tree2,\@taxonomy, 's');
     }
     
     for my $tax_string(@{$blast_taxonomies}){
     @taxonomy = split (/\|/, $tax_string); 
     build_taxtree_2($tax_tree2,\@taxonomy, 'b');
     }
     #
     
     my $I=taxtree2str_2 ($tax_tree2);
     my @tax_levels = split /\n/, $I;	
     
     return \@tax_levels;
}      
      
             
      
        
# from neva stap_summary.pl
sub build_taxtree{
	my ($T,$taxonomy)=@_;
	foreach my $taxa(@{$taxonomy}){
		my $t=$T->children()->{$taxa};
		if(!defined $t){
			$t=TaxNode->new;
			$t->count(1);
			$t->children({});
			$T->children->{$taxa}=$t;
		}else{
			$t->count($t->count+1);
		};
		$T=$t;
	};
};

sub build_taxtree_2{
	my ($T, $taxonomy, $type)=@_;
	if ($type eq 's'){
		foreach my $taxa(@{$taxonomy}){
			my $t=$T->stap_children()->{$taxa};
			if(!defined $t){
				$t=$T->blast_children()->{$taxa};
				if(!defined $t){
					$t=TaxNode->new;
					$t->blast_count(0);
					$t->blast_children({});
				}
				$t->stap_count(1);
				$t->stap_children({});
				$T->stap_children->{$taxa}=$t;
			}
			else{
				$t->stap_count($t->stap_count+1);
			};
			$T=$t;
		};
	}
	
	if ($type eq 'b'){
		foreach my $taxa(@{$taxonomy}){
			my $t=$T->blast_children()->{$taxa};
			if(!defined $t){
				$t=$T->stap_children()->{$taxa};
				if(!defined $t){
					$t=TaxNode->new;
					$t->stap_count(0);
					$t->stap_children({});					
				}
				$t->blast_count(1);
				$t->blast_children({});
				$T->blast_children->{$taxa}=$t;
			}
			else{
				$t->blast_count($t->blast_count+1);
			};
			$T=$t;
		};
	}	
};


sub taxtree2str_2{
	my $T=shift;
	my $outstr;
	my @stack;
	push @stack,$T;
	my %nodcnt;
	my $I=0;
	while(@stack>0){
		$T=$stack[-1];
		my $str="";
		for(my $i=0;$i<$I;$i++){
			$str.="----";
		};
		
		#merge two hashes
		my %hash_all = %{$T->stap_children};
		my %hash_blast = %{$T->blast_children};
		@hash_all{ keys %hash_blast } = values %hash_blast;
		
		foreach my $t(keys %hash_all){
			my $node=$T->stap_children->{$t};
			if (!defined $node) {
				$node=$T->blast_children->{$t};
			}
			if(!defined $nodcnt{$node}){
				$nodcnt{$node}=0;
				my $cnt_stap=$node->stap_count;
				my $cnt_blast=$node->blast_count;
				$str.="$t\t($cnt_stap;	$cnt_blast)\n";
				$outstr.=$str;
			}elsif($nodcnt{$node}==-1){
				next;
			};
			
			my %temp_hash_all = %{$node->stap_children};
		    my %temp_hash_blast = %{$node->blast_children};
		    @temp_hash_all{ keys %temp_hash_blast } = values %temp_hash_blast;
			
			my @tmp=keys %temp_hash_all;
			if($nodcnt{$node}==@tmp){
				#do not vizit this node next time
				$nodcnt{$node}=-1;
				#delete elements from stack
				$#stack--;
				$I--;
			}else{
				$I++;
				$nodcnt{$node}++;
				push @stack,$node;
			};
			last;
		};
	};
	return $outstr;
};


sub taxtree2str{
	my $T=shift;
	my $outstr;
	my @stack;
	push @stack,$T;
	my %nodcnt;
	my $I=0;
	while(@stack>0){
		$T=$stack[-1];
		my $str="";
		for(my $i=0;$i<$I;$i++){
			$str.="----";
		};
		
		
		foreach my $t(keys %{$T->children}){
			my $node=$T->children->{$t};
			if(!defined $nodcnt{$node}){
				$nodcnt{$node}=0;
				my $cnt=$node->count;
				$str.="$t\t($cnt)\n";
				$outstr.=$str;
			}elsif($nodcnt{$node}==-1){
				next;
			};
			my @tmp=keys %{$node->children()};
			if($nodcnt{$node}==@tmp){
				#do not vizit this node next time
				$nodcnt{$node}=-1;
				#delete elements from stack
				$#stack--;
				$I--;
			}else{
				$I++;
				$nodcnt{$node}++;
				push @stack,$node;
			};
			last;
		};
	};
	return $outstr;
};



sub clean_string{
	my $string = @_[0];
	$string =~ s/[\s+,:;'")(]+/_/g;
	return $string;
}

sub LCP{
                    	my ($seq1, $seq2)=@_;
                    	my $n =0;
                        my @array1 = split(//, $seq1);
						my @array2 = split(//, $seq2);
						while($array1[$n] eq $array2[$n] && ($n < scalar @array1 && $n < scalar @array2) ){
 						$n++; 
						}
						return $n;
                    }
                    
sub common_prefix {
	my ($seq1, $seq2)=@_;
    my $n =0;
    my @array1 = split(//, $seq1);
	my @array2 = split(//, $seq2);
	while($array1[$n] eq $array2[$n] && ($n < scalar @array1 && $n < scalar @array2) ){
 		$n++; 
	}
	return substr $seq1, 0, $n;
} 

sub multiple_common_prefix {
	my @array = @{$_[0]};
	my $size = scalar @array;
	if ($size == 2) {
		return common_prefix($array[0], $array[1]);
	}
	my $common = $array[0];
	for ($i = 1; $i < $size; $i++){
		$common = common_prefix($common, $array[$i]);
	}
	return $common;
}              
                    
                   
                    
sub difference{
                    	 my @array1 = @{$_[0]};
                    	 my @array2 = @{$_[1]};

                    	 my @difference; 
                    	 my %count = ();
                         foreach my $element (@array1, @array2) { $count{$element}++ }
                    
    					 foreach my $element (keys %count) {

    					 	if( $count{$element} == 1){
    					 		push @difference, $element;
    					 	}
    					 } 
    					 
                    	 return \@difference;
                    } 
                    
 sub hash_arr_difference{
 	my %hash = %{$_[0]};
 	my @array = @{$_[1]};
    
    
    my %difference = %hash;
    foreach(@array){
    	if(exists($difference{$_})){
    		delete $difference{$_};
    	}
    }
    return \%difference;
 }                  
                    
                    
sub treeAsString{
	                my $tree_file = $_[0];
					open TREE, "< $tree_file" or die "Невозможно открыть файл ".$tree_file."\n";
					# get a newick string from some source
 					my $tree_string;
 					my $t;
 					while($t = <TREE>){
 						$t =~ s/\n$//;
 						$tree_string = $tree_string.$t;
 					}
 					 close TREE;
 					 return $tree_string
}



sub parseTree{

 					# Call class method parse from Bio::Phylo::IO
 					my $tree = Bio::Phylo::IO->parse(
  					  -string => treeAsString($_[0]),
   					  -format => 'newick'
 					)->first;

 # note: newick parser returns 'Bio::Phylo::Forest'
 # Call ->first to retrieve the first tree of the forest.
 					return $tree;
	
}   

sub parseStringTree{

 					# Call class method parse from Bio::Phylo::IO
 					my $tree = Bio::Phylo::IO->parse(
  					  -string => $_[0],
   					  -format => 'newick'
 					)->first;

 # note: newick parser returns 'Bio::Phylo::Forest'
 # Call ->first to retrieve the first tree of the forest.
 					
					return $tree;
	
}  

sub find_distance0{
	

	   				my $tree = $_[0];
                    my @taxa = @{$_[1]};
                    	foreach (@taxa){
						#print $_;
					}
					my @nodes = @{ $tree->get_nodes_for_taxa(\@taxa) };
				
					my $mrca = $tree->get_mrca(\@nodes);
					my $dist1;
					my $t1 = $nodes[0];
					while ($t1 ne $mrca){
						$dist1 += $t1->get_branch_length();
						$t1 = $t1->get_parent();
					}
					my $dist2;
					my $t2 = $nodes[1];
					while ($t2 ne $mrca){
						$dist2 += $t2->get_branch_length();
						$t2 = $t2->get_parent();
					}
					return $dist1+$dist2;
}                 
 
 
sub find_distance{
	
	my $tree = $_[0];
	my $taxa_string1 = $_[1];
	my $taxa_string2 = $_[2];
 
    my $node1 = $tree->get_by_name($taxa_string1);
   # print "taxa1 ".$taxa_string1;
    my $node2 = $tree->get_by_name($taxa_string2);
   # print " taxa2 ".$taxa_string2;
  
    my $patristic_distance = $node1->calc_patristic_distance($node2);
    return $patristic_distance;
} 

sub min_distance_less_than{
	
	my $tree = $_[0];
	my @taxon1 = @{ $_[1] };
	my @taxon2 = @{ $_[2] };
	my $true = 0;
	#print scalar @taxon1;
	#print scalar @taxon2;
	OUTER: for my $i(0..(scalar @taxon1)-1){
		for my $j(0..(scalar @taxon2)-1){
		my $dist = find_distance($tree, $taxon1[$i], $taxon2[$j]);
		   if ($dist < $distance_threshold){
			  $true = 1;
			  last OUTER;
		   }
	    }
	}
	#print " that is true ".$true;
	return $true;

}

sub max_distance_less_than{
	
	my $tree = $_[0];
	my @taxon1 = @{ $_[1] };
	my @taxon2 = @{ $_[2] };
	my $true = 1;
	#print scalar @taxon1;
	#print scalar @taxon2;
	OUTER: for my $i(0..(scalar @taxon1)-1){
		for my $j(0..(scalar @taxon2)-1){
		my $dist = find_distance($tree, $taxon1[$i], $taxon2[$j]);
		   if ($dist > $distance_threshold){
			  $true = 0;
			  last OUTER;
		   }
	    }
	}
	#print " that is true ".$true;
	return $true;

}

 
## returns ProkMSA name for unique ProkMSA id, or empty string, if greengenes database doesn't contain such an id. 
 
sub get_prokmsa_name{
    my $sql = "select PROKMSANAME from id_to_name where PROKMSA_ID=".@_[0];
    my $sth = $dbh->prepare($sql);
    $sth->execute or die "SQL Error: $DBI::errstr\n";
    if (my @row = $sth->fetchrow_array){
    	return $row[0];
    }
    else {
    	return "";
    }

} 
 
