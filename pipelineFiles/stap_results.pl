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
	use vars qw($opt_i $opt_l $opt_x);
	use Getopt::Std;
	use XML::Simple;

$| = 1;
    
  
  ## INPUT:   -i $dir - path to directory with STAP output folders
  #           -l $chrom_list - path to file containing a list of chromatogram names (these must match the names of STAP output folders)
  #           by default this script will try to process all folders in $dir, assuming their names and contents comply with the aforementioned rules
  # STAP output folder X must contain X.ab1.cluster.fasta, X.mltree2.tab, X.mltree2.with_taxonomy, X.results;
  # this script produces X.mltree2.with_taxonomy_and_names (newick format), X.mltree2.svg (image) in $dir.X, 
  # and also results.xml in the X folder.
    
	my ($dir, $chrom_list);
  
	getopt('ilx');
  
	if (!$opt_i){
		print STDERR "Please, specify directory with STAP output folders";
		exit 4;
	}
	else {
		$dir = $opt_i;
	}
  
	if ($opt_l) {
		$chrom_list = $opt_l;
	}
	
	    ## magic constants
  my $xsl; ## path to xsl file which forms detailed output
  my $simple_xsl; ## path to xsl file which forms simple output
  my $genera_file; ## path to list of prok genera
  my $expectation_threshold;
  my $confidence_threshold;
  my $distance_threshold;
	
 my $email;
	if ($opt_x)	{
	print STDOUT " Using cfg\n";
		my $cfg = XML::Simple->new();
	
		my $data   = $cfg->XMLin($opt_x, 
		ForceArray => ['Flag','Read','Primer','PrimerDirection'],
		SuppressEmpty => 1,
		ContentKey => '-content' );

		$distance_threshold=$data->{DistanceThreshold};
		$xsl=$data->{XSLPath};
		$simple_xsl=$data->{SimpleXSLPath};
		$genera_file=$data->{GeneraFilePath};
		$confidence_threshold=$data->{ConfidenceThreshold};
		$expectation_threshold=$data->{ExpectThreshold};
		$email=$data->{Email};
	}
  else {
  print STDOUT " No cfg\n";
    ## magic constants
	$xsl = "/home/bcviss/pipelineFiles/results_with_aggregates.xsl"; ## path to xsl file which forms detailed output
	$simple_xsl = "/home/bcviss/pipelineFiles/simple.xsl"; ## path to xsl file which forms simple output
	$genera_file = "/home/bcviss/pipelineFiles/genera.txt"; ## path to list of prok genera
	$expectation_threshold = 0.05;
	$confidence_threshold = 0.80;
	$distance_threshold = 0.03;
  }
  
  
  
  my @chrom_array;
  
  my %valid_genera;
 # my %prokids_for_clstr;

  my $mode = "cluster"; # otu/gi/cluster

   
  ## check if input directories exist 
    if (!-e $dir) {
    	print STDERR "STAP results: No such directory: ".$dir."\n";
    	print STDERR "Exiting";
    	exit 2;
    }

   
    if (!defined $chrom_list){
    	opendir DIR, $dir or die "Can't open directory $dir: $!";
    	
        while (my $folder = readdir(DIR)) {
        	if (-d "$dir/$folder"){
        		if (($folder ne ".") & ($folder ne "..")){	
        		push @chrom_array, $folder;
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
			print STDERR "STAP results: No such file: ".$chrom_list."\n";
			print STDERR "Exiting";
			exit 2;
		}
	}
	
	if (-e $genera_file){
		open GENERA, "< $genera_file";
		while (<GENERA>){
			chomp;
			$valid_genera{$_} = 1;
		}
		close GENERA;
	}
	
    ## 
 
    
	## assume that we are in standard BCV project directory, 
	## so we can get chromatogram length from  Input/X/ps_prj/poly_dir/X.ab1.bqs

	my $mydir = dir($dir);
	my $input_dir = $mydir->parent."/Input/"; 
	my %xml_hash;
	my @all_chromatograms;
	my @all_taxonomies; #this array is used to form aggregated view
	my @all_blast_taxonomies; #this array is used to form aggregated view
	$dir = $dir."/";  
	
	 struct TaxNode => {
	 stap_count => '$',
	 blast_count => '$',
	 stap_children => '$',
	 blast_children => '$',
      };  
	
	

	
	
	## form xml structure: collect info for each chromatogram
	CHROM: foreach my $chrom_name(@chrom_array){

    my $current_dir = $dir.$chrom_name."/";
	
	
	 my @splitter =  grep { /\S/ } split(/[\/\\]/, $dir);
	 my $path_length = scalar @splitter;
	 my $application_current_dir = "";
	 
	 $application_current_dir = "./".$chrom_name."/";
	
	
	
    my $chrom_length;
    
    my $aggregate;
    my @chrom_stap_taxes;
    my @chrom_blast_taxes;
    
    print STDOUT $current_dir."\n";
    ## count lines in ab1.bqs to find chromatogram length
    my $count_lines_file = $input_dir.$chrom_name."/ps_prj/poly_dir/".$chrom_name.".ab1.bqs";
    if (-e $count_lines_file){
    open FILE, "< $count_lines_file";
    my $counter = 0;
    while (<FILE>){
    	$counter++;
    }
    close FILE;
    $chrom_length = $counter-2;
    }

    ##
 
    my $clust_file =  $current_dir.$chrom_name.".ab1.cluster.fasta";
    open clusters, "< $clust_file" or do {
		print STDERR "Can't open file ".$clust_file."\n";
		print STDERR "STAP results: skipping ".$chrom_name." because there is no $chrom_name.ab1.cluster.fasta\n";
		next;
	};

    my $total_expectation;

    my @candidates;
	my %hash_simple_view_1;
	
	
    while (my $row = <clusters>) {

           	       my @splitter = split(/\s+/, $row);

            	       $row = <clusters>;
            	       $_ = $row;
            	        my $exp = $splitter[1];
            	        $exp=~ s/\n$//;
						
						my $name1;
						if ($exp !~ /^\d\.\d*$/){
							$exp = "ND";
							$name1 = (substr $splitter[0], 1);
						}
						
                        else {
							$name1 = (substr $splitter[0], 1)."_".(sprintf("%.3f",$exp));
						}	
            	       
					   
            	        if (($exp >= $expectation_threshold || $exp eq "ND") && -e $current_dir.$name1.".mltree2.tab" && -e $current_dir.$name1.".results"){

							my $gaps += s/-//ig;
							if ($exp ne "ND"){
								$total_expectation += $exp;
							}
            	      
							my $cand_file = $current_dir.$name1.".mltree2.tab";
							my $tree_file = $current_dir.$name1.".mltree2.with_taxonomy";
							my $blast_file = $current_dir.$name1.".results";
            	       
							my %prokids_for_clstr;
							if ($mode eq "cluster"){
								my $clstr_file = $current_dir.$name1.".temp.clstr";
								open CLSTR, "<$clstr_file" or die "cannot open cluster file $clstr_file";
								my $str = <CLSTR>;
								while($str = <CLSTR>){
									my $centroid_id;
									my @ids;
									while(! ($str =~ m/^>/)){
										@t = split /[K_]/, $str;
										my $id = $t[1];
										if ($str =~ m/\*$/){
											$centroid_id = $id;
										}
										push @ids, $id;
										$str = <CLSTR>;
										if (!$str) {last;}
									}
									
									$prokids_for_clstr{$centroid_id} = \@ids;
							}
							close CLSTR;
							##unlink $clstr_file;
						}
					   
            	       my %id_to_name;
            	       
            	       ##all neighbours for this query sequence

            	       unless (open blast_file, "< $blast_file") { warn "cannot open ".$blast_file; next;}
            	       my $blast_line = <blast_file>;
            	       my @splitter =  split /[\t=]/, $blast_line;
            	       my $blast_id = "";
            	       my $blast_tax = "";
            	       my $blast_identity = "ND";
            	       if ($splitter[4] =~ m/^[\d\.]+$/) {$blast_identity = $splitter[4];}
            	       my $blast_covered = "ND";
            	       if ($splitter[6] =~ m/^[0-9]+$/) {$blast_covered = $splitter[6]}; ## be careful: this value must be divided by query length to get coverage
            	       my $blast_evalue = "ND";
            	       if ($splitter[8] =~ m/^[\d\.\-Ee\+]+$/) {$blast_evalue = $splitter[8]}
            	       my $blast_score = "ND";
            	       if ($splitter[10] =~ m/^[\d\.]+$/) {$blast_score = $splitter[10]}
            	       if ($splitter[2] =~ m/^PROK([0-9]+)_PROK[0-9\.]+\|([\w\|\]\[\-]+)/){
                       $blast_id = $1;
            	       $blast_tax = fineliner($2);
            	       }
            	       
            	       close blast_file;
            	       
            	       
            	       unless (open candidate, "< $cand_file") {warn "cannot open ".$cand_file; next;}
            	       my @neighbours;
					   my %proknames;
            	       my @closest_neighbours;
            	       my $go_on = 1;
            	       my $pivotal_node;
            	        
            	       while (my $hit = <candidate>) {
            	       my @hit_splitter = split(/\t/, $hit);
            	       my @prok_id = split(/[K_]/, $hit_splitter[1]);
            	       my $taxonomy = $hit_splitter[7];
            	       $taxonomy=~ s/\n$//;
            	       my $confidence = $hit_splitter[6];
					   if ($confidence =~/^PROK/) {
						$confidence = "undefined"; # it's a leaf (outgroup), not a node.  
					   }
            	       my $node = $hit_splitter[5];
            	       my $distance_in_nodes = $hit_splitter[4];
            	       my $distance_between_nodes = $hit_splitter[3];
            	       my $distance_between_leafs = $hit_splitter[2];
            	       my $chosen_one = "yes"; #boolean: if this hit is chosen for the final taxonomy analysis 
            	       my $prok_name = get_prokmsa_name($prok_id[1], \%prokids_for_clstr);
            	       
            	       $id_to_name{$prok_id[1]} = $prok_name;
            	      
            	       
            	       if (scalar @hit_splitter > 7){ ## else this sequence has no confidence, so I suppose user don't need it
            	        	if ($go_on == 1 || $node == $pivotal_node) {
            	        		$chosen_one = "yes";
								if ($mode eq "gi"){
									put_genus($prok_name, \%proknames);
								}
								elsif ($mode eq "otu" || $mode eq "cluster"){ 
									my @splitter = split /,/, $prok_name;
									foreach $org_string(@splitter){
										@name_count = split /[\)\(]/, $org_string;
										my $org_name = $name_count[0];
										$org_name =~ s/^\s+//g;
										$org_name =~ s/\s+$//g;
										my $org_count = $name_count[1];
										if (exists $proknames{$org_name}){
											my $count = $proknames{$org_name}+$org_count;
											$proknames{$org_name} = $count;
										}
										else {
											$proknames{$org_name} = $org_count;
										}
										
									}
								}
								

            	        	}
            	       		else {
            	        		$chosen_one = "no";
            	       		}

            	       	
            	      		 my $neighbour = {
            	       			'-taxonomy' => fineliner($taxonomy),
            	       			'-PROKid' => $prok_id[1],
            	       			'-confidence' => $confidence,
            	       			'-chosen_one' => $chosen_one,
            	       			'-PROKname' => $prok_name,
            	       			'-distance_between_leafs' => $distance_between_leafs,
            	       			'-distance_between_nodes' => $distance_between_nodes,
            	       			'-distance_in_nodes' => $distance_in_nodes,
            	       			'-branching_position' => $node
            	       		};
            	       		push @neighbours, $neighbour;
            	       	
            	    
            	       		if ($go_on == 1){ ## collect closest neighbours, which are not separated from the query on the given confidence level 
            	           		if ($confidence >= $confidence_threshold){
            	       	   			$go_on = 0;
            	       	   			$pivotal_node = $node;
            	           		}
            	           		push @closest_neighbours, $neighbour;
            	        	} 
            	       
            	        	elsif ($node == $pivotal_node){
            	       			push @closest_neighbours, $neighbour;
            	        	}

            	        }
            	        
            	        }
            	       
            	        close candidate;

            	     ##find neighbours with different taxonomy - this part is awful
            	    
            	     my %uniq;
                     foreach my $item (@closest_neighbours) {
                     	if (exists($uniq{$item->{'-taxonomy'}})) {
					 		push( @{ $uniq { $item->{'-taxonomy'} } }, $item->{'-PROKid'}); 
                     	}
						else{
					 		my @prokids;
					 		push @prokids, $item->{'-PROKid'};
					 		$uniq{$item->{'-taxonomy'}} = \@prokids;
					 	}
                     } 
            	     
                     ##
                     my @taxonomy;
                     my @taxonomies = keys %uniq;
                     my $number = scalar @taxonomies;
                    
                   
                     #taxonomy for candidate sequence
                     if ($number == 1) {
                     push @taxonomy, $taxonomies[0];
                     }
                     
                     elsif ($number > 1){
                     for  my $i (0 .. ($number-2))
                      {
                      	for my $j ($i+1 .. ($number-1)){
                      		my $diff = &trimmedLCP($taxonomies[$i], $taxonomies[$j]);  	
                      		my $sub1 = substr ($taxonomies[$i], $diff);
                      		my $sub2 = substr $taxonomies[$j], $diff;

                      	}
                      	
              
                      }
                      
                      ## собрали хэш %uniq_best, в котором хранятся все таксономии соседей, за исключением тех, для которых есть более глубокий аналог

                     my %uniq_best = %uniq; ##added: min to max

                     my @taxonomies_best = keys %uniq_best;
                     my $number_best = scalar @taxonomies_best;
                  
                     if ($number_best == 1){
                     	 push @taxonomy, $taxonomies_best[0];
                     }
                     
                     else {
                     	my @taxons;
                     	foreach my $tax (keys %uniq_best){
                     		    my @prok_ids = @{ $uniq_best{$tax} };
                     		    $tax =~ s/\|/_/g;   
                     		    my @ids_for_taxon; 
                     		    foreach my $prok_id(@prok_ids){
                     		    	push @ids_for_taxon, "PROK".$prok_id."_".$tax;
                     		    }	
                     			push @taxons, \@ids_for_taxon;
                     	}
                     	
                     	
                     	my $tree = parseTree($tree_file);
                     	my $are_close_relatives = 1;
                     	
                     	
                     	
                     	OUTER: for my $i (0 .. ($number_best-2)) {
                      	for my $j ($i+1 .. ($number_best-1)){

                      		my $are_close = max_distance_less_than($tree, \@{ $taxons[$i]} , \@{$taxons[$j]});
                      		if (!$are_close){
                      			$are_close_relatives = 0;
                            	last OUTER;
                      		}

                      	}

                      }
                    
                      if ($are_close_relatives == 1){
                      	#выдать всех
                      	foreach (@taxonomies_best){
                      		push @taxonomy, (substr $_, index($_, '_')+1);
                      	}
                      }
                      
                      else {
                      	#выдать общую подстроку
                      	my $t = $taxonomies_best[0];
                      	for my $i (1 .. ($number_best-1)){
                      		$t = substr ($t, 0, trimmedLCP($t, $taxonomies_best[$i]));
                      	}
                      	
                      	$t =~ s/\|$//;
                      	push @taxonomy, $t;
                      }
                     	

                     }
                     #
                     }
                     
                   my $coverage = "ND";
                   if ($blast_covered ne "ND") {$coverage = sprintf("%.2f", $blast_covered*100/(length($row)-$gaps));}
				   
				   my $proknames_string;
				   foreach my $name(keys %proknames){
					$proknames_string .= $name." (".$proknames{$name}."), ";
				   }
				   chop($proknames_string);
				   chop($proknames_string);
				   
				   foreach my $tx(@taxonomy){
				   push( @{ $hash_simple_view_1{$tx} }, $name1); 
				   }
				   
            	       my $candidate = {
                       '-name' => $name1,
                       '-expectation' =>  $exp,
                       '-length' => length($row)-$gaps,
					  '-newick_ref' => $application_current_dir.$name1.".mltree2.with_taxonomy.with_names",
                     ##  '-newick_ref' => $tree_file.".with_names",
					   '-img_ref' => $application_current_dir.$name1.".mltree2.html",
                      ## '-img_ref' => $current_dir.$name1.".mltree2.html",
					 
                       	'beast' => {
                       		'blast' => {
                       			'-prokID' => $blast_id,
								'-PROKname' => get_prokmsa_name($blast_id, 0),
            	                '-taxonomy' => $blast_tax,
            	                '-identity' => $blast_identity,
            	                '-coverage' => $coverage,
            	                '-e-value' => $blast_evalue,
            	                '-bit-score' => $blast_score
                       		},
                            'taxonomy' => \@taxonomy,
							'PROKnames' => $proknames_string,
							
                       		'-confidence' => $closest_neighbours[-1]->{'-confidence'},
                       		'neighbour' => \@neighbours
                       	}
            	       };                  
                     
                      push @candidates, $candidate;
                      push @all_blast_taxonomies, $blast_tax;
                      push @chrom_blast_taxes, $blast_tax;
                      
    ## one line taxonomy for STAP                  
	 	if ((scalar @taxonomy) == 1){
	 		push @all_taxonomies, $taxonomy[0];
	 		push @chrom_stap_taxes, $taxonomy[0];
	 	}
 	 else {
 	 		my $size = scalar @taxonomy;
 	 	    my $common_prefix = multiple_common_prefix(\@taxonomy);
 	 	    my $prefix_length = length($common_prefix);
 	 	    my $ambi_string;
 	 	    my $counter = 0;
 	 	    foreach my $str(@taxonomy){
 	 	    	$ambi_string .= (substr ($str, $prefix_length));
				$ambi_string =~ s/\|/_/g;
 	 	    	$counter++;
 	 	    	if ($counter < $size) {
 	 	    		$ambi_string .= " or ";
 	 	    	}
 	 	    }
	 		push @all_taxonomies, $common_prefix.$ambi_string;
	 		push @chrom_stap_taxes, $common_prefix.$ambi_string;
	 	}
                      
                      
                      
                      
                      ## print newick tree file with prokmsa names

 						my $tree_string = treeAsString($tree_file);
 						$tree_string =~ s/PROK(\d+)([\w\]\[\-]+):/PROK.$1.$2."_".clean_string($id_to_name{$1}).":"/ge;
 						
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


 					$treedrawer->set_tree(parseStringTree($tree_string));
 					$treedrawer->set_text_width(1500);

                    my $svg_file =   $current_dir.$name1.".mltree2.html";
					open FIG, "> $svg_file";
					print FIG "<!DOCTYPE html>\n";
					print FIG '<html lang="en">';
					print FIG "\n";
					print FIG "<head>\n";
					print FIG '<meta charset="UTF-8" />';
					print FIG "\n";
					print FIG "</head>\n";
					print FIG "<body>\n";
 					print FIG $treedrawer->draw;
					print FIG "</body>\n";
					print FIG "</html>\n";
 					close FIG;	
 				    ##						  

            	}      
            }
	
			if (scalar @candidates > 0) {
	
	
	       	my %candidate_for_name;
			foreach my $cand(@candidates){
			$candidate_for_name{$cand->{'-name'}} = $cand;
			}
			

			#create simple view: first we select only best candidates (we get rid of enterobacteriaceae, if e|shigella or e|unclassified is present)
			my @deepest;
			

			
			##update 1-08-2014
			
		#	foreach my $tx(keys %hash_simple_view_1){
		#		my $checker = 1;
		#	#	print STDOUT $tx." from hash1\n";
		#		foreach my $t(keys %hash_simple_view_1){
		#			if ($t=~ m/^\Q$tx\E\|.*/) {
		#			$checker = 0;
		#			last;
		#			}
		#		}
		#		if ($checker == 1){
		#			push @deepest, $tx;
		#		}
		#	}	
			
			my %tax_to_exp;
			
			foreach $tax(keys %hash_simple_view_1){
			        my $exp = 0;
					foreach my $cand_name(@{ $hash_simple_view_1{$tax} }){
					    my $t_exp = $candidate_for_name{$cand_name} -> {'-expectation'};
						my $taxons_count = scalar @{$candidate_for_name{$cand_name} -> {'beast'} -> {'taxonomy'}};
						if ($t_exp ne "ND"){
							$exp += $t_exp/$taxons_count;
						}
					}	
					$tax_to_exp{$tax} = $exp;
			}
			
			
			
			my %low_taxons_for_higher;
			my @sorted_taxons = keys %hash_simple_view_1;
			@sorted_taxons = sort {length($a) <=> length($b)} @sorted_taxons;
			

		 for (my $i = 0; $i < scalar @sorted_taxons; $i++){
			my $total_exp;
			
			my $deep = 1;
			my @lower_taxons;
			
				for (my $j = $i+1; $j < scalar @sorted_taxons; $j++){
				
						if ($sorted_taxons[$j]=~ m/^\Q$sorted_taxons[$i]\E\|.*/){
								$total_exp += $tax_to_exp{$sorted_taxons[$j]};
								push @lower_taxons, $sorted_taxons[$j];
								$deep = 0;
						}
				}
				
				if ( scalar (split /\|/, $sorted_taxons[$i]) > 4){
				
				foreach my $lower(@lower_taxons){
					my $t_exp = $tax_to_exp{$lower};
					my $add = $tax_to_exp{$sorted_taxons[$i]}*($t_exp/$total_exp);
					$tax_to_exp{$lower} = $t_exp + $add;
				}
			
				}
	
			if ($deep == 1){
				push @deepest, $sorted_taxons[$i];
			}	
	
			}

			my @array_simple_view;
			my %hash_simple_view;
			
			foreach my $tx(@deepest){
				my $best_cand;
				my $best_info = 0;
				
				my @cands_for_taxon;
				foreach my $cand_name(@{ $hash_simple_view_1{$tx} }){
					push @cands_for_taxon, $candidate_for_name{$cand_name};
				}	

				my @sorted_cands_for_taxon = sort {5*($b -> {'beast'} -> {'blast'} -> {'-identity'}) 
												   + ($b -> {'beast'} -> {'blast'} -> {'-coverage'}) <=>
												   5*($a -> {'beast'} -> {'blast'} -> {'-identity'})
												   + ($a -> {'beast'} -> {'blast'} -> {'-coverage'})} @cands_for_taxon;		
				
				$best_cand = $sorted_cands_for_taxon[0];
				my $stap_prokmsa_names = $best_cand ->{'beast'} -> {'PROKnames'};

				my $i = 1;
				while ($i < scalar @cands_for_taxon && 
					   ($best_cand -> {'beast'} -> {'blast'} -> {'-identity'} -
					   $sorted_cands_for_taxon[$i]-> {'beast'} -> {'blast'} -> {'-identity'} 
					   < 0.5) 
					   &&
					   ($best_cand -> {'beast'} -> {'blast'} -> {'-coverage'} -
					   $sorted_cands_for_taxon[$i]-> {'beast'} -> {'blast'} -> {'-coverage'} 
					   < 3) 
					   ){
					$stap_prokmsa_names = sum_proknames($stap_prokmsa_names, $sorted_cands_for_taxon[$i] ->{'beast'} -> {'PROKnames'});
					$i++;
					
				}
				
				
				
				my $simple_view = {
				    '-expectation' => $tax_to_exp{$tx},
					'-stap_taxonomy' => $tx,
					'-PROKnames' => $stap_prokmsa_names,
					'-blast_taxonomy' => $best_cand ->{'beast'} -> {'blast'} -> {'-taxonomy'},
					'-blast_coverage' => $best_cand -> {'beast'} -> {'blast'} -> {'-coverage'},
					'-blast_identity' => $best_cand -> {'beast'} -> {'blast'} -> {'-identity'},
					'-blast_PROKid' => $best_cand -> {'beast'} -> {'blast'} -> {'-prokID'},
					'-blast_PROKname' => $best_cand -> {'beast'} -> {'blast'} -> {'-PROKname'},
				};
				
				##upd 19-08
				my $identification_line = $simple_view -> {'-blast_PROKid'}." ".$simple_view -> {'-expectation'}." ".$simple_view -> {'-blast_identity'}." ".$simple_view -> {'-blast_coverage'};
				my @temp;
				if (exists $hash_simple_view{$identification_line}){
					@temp = @{$hash_simple_view{$identification_line}};
				}
				push @temp, $simple_view;
				$hash_simple_view{$identification_line} = \@temp;

			}
			
			##upd 19-08
			foreach my $id_line(keys %hash_simple_view){
				my $simple_view = ${$hash_simple_view{$id_line}}[0];
				if (scalar @{$hash_simple_view{$id_line}} > 1){
					my $stap_tax;
					my $exp = 0;
					foreach my $view(@{$hash_simple_view{$id_line}}){
						$stap_tax .= $view -> {'-stap_taxonomy'}." or ";
						if ($view -> {'-expectation'} ne "ND"){
							$exp += $view -> {'-expectation'};
						}
					}
					$stap_tax = substr $stap_tax, 0, length($stap_tax)-4;
					$simple_view -> {'-stap_taxonomy'} = $stap_tax;
					$simple_view -> {'-expectation'} = $exp;
				}
				push @array_simple_view, \%$simple_view;
			}
			#
			
			
			my @chrom_tax_levels =  @{aggregate_tree(\@chrom_stap_taxes, \@chrom_blast_taxes)};
		
            my $chrom = {
            	 '-name' => $chrom_name,
                 '-length' => $chrom_length,
				 '-fasta_ref' => $application_current_dir.$chrom_name.".ab1.cluster.fasta",
                # '-fasta_ref' => $clust_file, 
                 '-total_expect'  => $total_expectation,
            	 'candidate' => \@candidates,
                 'tree' => { 
                 	 'item' => \@chrom_tax_levels
                 },  
				 'simple' => \@array_simple_view
			};
			
            push @all_chromatograms, $chrom;
			
			}
			
			## 21 may 2015 upd
			else { ## if @candidates is empty
				my $chrom = {
            	 '-name' => $chrom_name,
                 '-length' => $chrom_length,
				 '-fasta_ref' => $application_current_dir.$chrom_name.".ab1.cluster.fasta",
                # '-fasta_ref' => $clust_file, 
                 '-total_expect'  => $total_expectation,
            	 'candidate' => "No valid sequences with acceptable BLAST results was produced",
                 'tree' => "-",                  
				 'simple' => "-"
				};
				push @all_chromatograms, $chrom;
			
			}

             close clusters;
} 


     my @all_chromatograms_sorted =  sort { $a->{'-name'} <=> $b->{'-name'} } @all_chromatograms;
	 $xml_hash ->{'case'} -> {'chromatogram'} = \@all_chromatograms_sorted;

	 if (@all_taxonomies && @all_blast_taxonomies){
		my @tax_levels = @{aggregate_tree(\@all_taxonomies, \@all_blast_taxonomies)};
		$xml_hash ->{'case'} -> {'tree'} -> {'item'} = \@tax_levels;
	 }
   
     ##   
   
     ## create xml file    
     my $tpp = XML::TreePP->new();
     # my $decl = "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>"."\n"."<?xml-stylesheet type=\"text/xsl\" href=\"".$xsl."\"?>";
     $tpp->set( xml_decl => $decl );
     $tpp->set( indent => 2 );
     $tpp->set( first_out => [ '-name', '-length', '-taxonomy', 'taxonomy'] );
     $tpp->set( last_out => [ '-fasta_ref'] );
     my $xml = $tpp->write( $xml_hash );
     my $xml_file =  $dir."results.xml";
     open xml, "> $xml_file" or die "Cannot open file ".$xml_file;
     print xml $xml;
     close xml;
     ##
     
     ##convert xml file into html
     
     # initialize the parser and XSLT processor
	 my $parser = XML::LibXML->new( );
	 my $xslt = XML::LibXSLT->new( );
	 my $stylesheet_detailed = $xslt->parse_stylesheet_file( $xsl );
	 my $stylesheet_simple = $xslt->parse_stylesheet_file( $simple_xsl );
	 

	 # parse, transform, print out result

  	 my $source_doc = $parser->parse_file( $xml_file );
  	 my $result_detailed = $stylesheet_detailed->transform( $source_doc );
	 my $result_simple = $stylesheet_simple->transform( $source_doc );
  	 my $html_detailed = $dir."detailed_results.html";
	 my $html_simple = $dir."simple_results.html";
	 
  	 open HTMLD, ">$html_detailed" or die "Cannot open ".$html_detailed;
  	 print HTMLD $stylesheet_detailed->output_string( $result_detailed );
  	 close HTMLD;
	 
	 open HTMLS, ">$html_simple" or die "Cannot open ".$html_simple;
  	 print HTMLS $stylesheet_simple->output_string( $result_simple );
  	 close HTMLS;
  	 ##

      
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
                    
 # for Bacteria|Planococcus and Bacteria|Planobacillus returns 9 (first mismatch index)         
 sub trimmedLCP   {
 	my ($seq1, $seq2) = @_;
 	if ($seq1 eq $seq2){
 		return length($seq1);
 	} 
 	else { ##else charAt(LCP) must be |, if no Planococcus interferes )
 		my $LCP = LCP($seq1, $seq2);
 		my $char = substr($seq1, $LCP-1, 1);
 		if ($char eq '|'){
 			return $LCP;
 		}
 		else {
 			return rindex((substr $seq1, 0, $LCP), '|')+1;
 		}
 	}
 }                
                  
                    
sub common_prefix {
	my ($seq1, $seq2) = @_;
    my $n = 0;
    my @array1 = split(//, $seq1);
	my @array2 = split(//, $seq2);
	while($array1[$n] eq $array2[$n] && ($n < scalar @array1 && $n < scalar @array2) ){
 		$n++; 
	}
	return substr $seq1, 0, $n;
} 

sub trimmed_common_prefix {
	my ($seq1, $seq2) = @_;
	return substr $seq1, 0, trimmedLCP($seq1, $seq2);
}

sub multiple_common_prefix {
	my @array = @{$_[0]};
	my $size = scalar @array;
	if ($size == 2) {
		return trimmed_common_prefix($array[0], $array[1]);
	}
	my $common = $array[0];
	for ($i = 1; $i < $size; $i++){
		$common = trimmed_common_prefix($common, $array[$i]);
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
my $string = treeAsString($_[0]);
$string =~ tr/][/}{/; ## bio phylo parser can't handle [TISSIERELLACEAE]

 					# Call class method parse from Bio::Phylo::IO
 					my $tree = Bio::Phylo::IO->parse(
  					  -string => $string,
   					  -format => 'newick'
 					)->first;

 # note: newick parser returns 'Bio::Phylo::Forest'
 # Call ->first to retrieve the first tree of the forest.
 					return $tree;
	
}   

sub parseStringTree{
my $string = $_[0];
$string =~ tr/][/}{/; ## bio phylo parser can't handle [TISSIERELLACEAE]
 					# Call class method parse from Bio::Phylo::IO
 					my $tree = Bio::Phylo::IO->parse(
  					  -string => $string,
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
	$taxa_string1  =~ tr/][/}{/;
	$taxa_string1  =~ s/\s//g;
	my $taxa_string2 = $_[2];
	$taxa_string2  =~ tr/][/}{/;
	$taxa_string2  =~ s/\s//g;
 
    my $node1 = $tree->get_by_name($taxa_string1);
    my $node2 = $tree->get_by_name($taxa_string2);

    my $patristic_distance = $node1->calc_patristic_distance($node2);
    return $patristic_distance;
} 

sub min_distance_less_than{
	
	my $tree = $_[0];
	my @taxon1 = @{ $_[1] };
	my @taxon2 = @{ $_[2] };
	my $true = 0;
	OUTER: for my $i(0..(scalar @taxon1)-1){
		for my $j(0..(scalar @taxon2)-1){
		my $dist = find_distance($tree, uc $taxon1[$i], uc $taxon2[$j]);
		   if ($dist < $distance_threshold){
			  $true = 1;
			  last OUTER;
		   }
	    }
	}

	return $true;

}

sub max_distance_less_than{
	
	my $tree = $_[0];
	my @taxon1 = @{ $_[1] };
	my @taxon2 = @{ $_[2] };
	my $true = 1;
	OUTER: for my $i(0..(scalar @taxon1)-1){
		for my $j(0..(scalar @taxon2)-1){

		my $dist = find_distance($tree, uc $taxon1[$i], uc $taxon2[$j]);
		   if ($dist > $distance_threshold){
			  $true = 0;
			  last OUTER;
		   }
	    }
	}
	
	return $true;

}

 
## returns ProkMSA name for unique ProkMSA id, or empty string, if greengenes database doesn't contain such an id. 
 
#sub get_prokmsa_name{
#	my $dbh = DBI->connect('dbi:mysql:greengenes','weidewind','lotus82') or die "Connection Error: $DBI::errstr\n";
#    my $sql = "select PROKMSANAME from id_to_name where PROKMSA_ID=".@_[0];
#    my $sth = $dbh->prepare($sql);
#    $sth->execute or die "SQL Error: $DBI::errstr\n";
#    if (my @row = $sth->fetchrow_array){
#    	return $row[0];
#    }
#    else {
#    	return "";
#    }
#} 

sub get_prokmsa_name{
	my $dbh = DBI->connect("DBI:mysql:database=bcviss;host=db",'bcviss','tyShmij2') or die "Connection Error: $DBI::errstr\n";
	if ($mode eq "gi" || (defined @_[1] && @_[1] == 0)){
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
	
	if ($mode eq "otu"){
    	my $sql = "select proknames from 99_otu_map where otu=".@_[0];
    	my $sth = $dbh->prepare($sql);
    	$sth->execute or die "SQL Error: $DBI::errstr\n";
    	if (my @row = $sth->fetchrow_array){
    		return $row[0];
    	}
    	else {
    		return "";
    	}
	}
	
	if ($mode eq "cluster"){
		my %prokids_for_clstr = %{$_[1]};
		my @ids = @{$prokids_for_clstr{$_[0]}};
		my %proknames;
		foreach my $id (@ids){
			my $sql = "select PROKMSANAME from id_to_name where PROKMSA_ID=".$id;
			my $sth = $dbh->prepare($sql);
			$sth->execute or die "SQL Error: $DBI::errstr\n";
			if (my @row = $sth->fetchrow_array){
				put_genus($row[0], \%proknames);
			}
		}
		my $proknames_string;
		foreach my $name(keys %proknames){
			$proknames_string .= $name." (".$proknames{$name}."), ";
		}
		chop($proknames_string);
		chop($proknames_string);
		
		return $proknames_string;
		
	}

   else {return ""};
} 

sub sum_proknames{
my ($str1, $str2) = @_;
my %proknames;
my $proknames_string;
my @arr1 = split /,/, $str1;
my @arr2 = split /,/, $str2;
foreach my $org(@arr1){
	my ($org_name, $org_count) = split /[\)\(]/, $org;
	$org_name =~ s/^\s+//;
	$org_name =~ s/\s+$//;
	$proknames{$org_name} = $org_count;
}
foreach my $org(@arr2){
	my ($org_name, $org_count) = split /[\)\(]/, $org;
	$org_name =~ s/^\s+//;
	$org_name =~ s/\s+$//;
	if (exists $proknames{$org_name}){
	my $t = $proknames{$org_name} + $org_count;
	$proknames{$org_name} = $t;
	}
	else {
	$proknames{$org_name} = $org_count;
	}
}

foreach my $name(keys %proknames){
	$proknames_string .= $name." (".$proknames{$name}."), ";
}
chop($proknames_string);
chop($proknames_string);

return $proknames_string;
}

sub fineliner{
my @full_taxonomy = split /\|/, $_[0];
my $fine_taxonomy;
for(my $i = 0; $i < scalar @full_taxonomy; $i++){
    my $tax_node = $full_taxonomy[$i];
	if ($tax_node eq "UNCLASSIFIED" || $i == 6){
		$fine_taxonomy .= lc $tax_node;
	}
	else {
		$fine_taxonomy .= substr $tax_node, 0, 1;
		my $temp = substr $tax_node, 1, length($tax_node)-1;
		$fine_taxonomy .= lc $temp;
		$fine_taxonomy .= "|";
	}
}
return $fine_taxonomy
}

sub put_genus {
				my ($prok_name, $proknames) = @_;
				if ($prok_name =~ m/(^([A-Z]([a-z]+|\.))\s([a-z\-_\.]+))/){
					if ($2 =~ m/^[A-Z]\.$/ || exists $valid_genera{$2}){
						if ($4 eq "subsp." || $4 eq "genomosp."){
							my @splitter = split /\s+/, $prok_name;
							my $tname = $splitter[0]." ".$splitter[1]." ".$splitter[2];
							if (exists $proknames -> {$tname}) {
								my $counter = $proknames -> {$tname};
								$counter++;
								$proknames -> {$tname} = $counter;
							}
							else {
								$proknames -> {$tname} = 1;
							}
						}
						else {
								if (exists $proknames -> {$1}) {
									my $counter = $proknames -> {$1};
									$counter++;
									$proknames -> {$1} = $counter;
								}
								else {
									$proknames -> {$1} = 1;
							}
						}	
					}
				
				}

}
 
