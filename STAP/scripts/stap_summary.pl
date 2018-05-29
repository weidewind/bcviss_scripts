#!/usr/bin/perl
#This script calculates summary statistics for STAP results
#OPTIONS
#dirlist [TXT] - a list of STAP output directories
#outdir - a path for output
#Output files:
#summary.results - union of STAP results
#summary.blast.txt - Taxonomy report
#summary.tree1.txt - Taxonomy report
#summary.tree2.txt - Taxonomy report
use strict;
use File::Find;
use File::Basename;
use Class::Struct;

my ($dirlist,$outdir)=@ARGV;
$outdir=~s/\s+\/?$//;
$outdir.="/";
my @results_files;
my @indirs;
my @rescnt;
sub get_results_file_name{
	my ($name,$dir,$ext) = fileparse($File::Find::name,'\.[^\.]*$');
	if($ext=~/\.results/){
		push @results_files,$File::Find::name;
		$rescnt[-1]++;
	};
};
open INF, "< $dirlist" or die "\nUnable to open input file $dirlist!";
while(<INF>){
	if(/\S+/){
		my $dir=$_;
		$dir=~s/\s+\/?$//;
		$dir.="/";
		push @rescnt,0;
		push @indirs,$dir;
		find (\&get_results_file_name,($dir));
	};
};
close INF;

my $ofile=$outdir."summary.results";
open SUMF, "> $ofile" or die "\nUnable to open file $ofile!";
my $I=0;

# initialize a structure
struct TaxNode => {
	count => '$',
	children => '$',
};
my $tax_blast=TaxNode->new;
$tax_blast->count(0);
$tax_blast->children({});

my $tax_tree1=TaxNode->new; 
$tax_tree1->count(0);    
$tax_tree1->children({});

my $tax_tree2=TaxNode->new;
$tax_tree2->count(0);    
$tax_tree2->children({});

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
			$str.="\t";
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

for( my $i=0;$i<@indirs;$i++){
	my $dir=$indirs[$i];
	print SUMF "indir=\"$dir\"\n";
	for(my $j=$I;$j-$I<$rescnt[$i];$j++){
		my $ifile=$results_files[$j];
		open RESF, "< $ifile";
		my $str=<RESF>;
		close RESF;
		print SUMF $str;
		my @cols=split /\s+/, $str;
		my @taxonomy;
		if(@cols>1){
			#blast
			$str=$cols[1];
			if($str=~/BLAST/){
				$str=~s/^\w+=//;
				$str=~s/[\w.]+\|//;
				@taxonomy=split /\|/, $str;
				build_taxtree($tax_blast,\@taxonomy);
			};
		};
		######
		if(@cols>2){
			#tree1
			$str=$cols[2];
			if($str=~/TREE1/){
				$str=~s/^\w+=//;
				$str=~s/[\w.]+\|//;
				@taxonomy=split /\|/, $str;
				build_taxtree($tax_tree1,\@taxonomy);
			};
		};
		######
		if(@cols>3){
			#tree2
			$str=$cols[3];
			if($str=~/TREE2/){
				$str=~s/^\w+=//;
				$str=~s/[\w.]+\|//;
				@taxonomy=split /\|/, $str;
				build_taxtree($tax_tree2,\@taxonomy);
			};
		};
		######
	};
	$I+=$rescnt[$i];
};
close SUMF;
	
$ofile=$outdir."summary.blast.txt";
open SUMF, "> $ofile" or die "\nUnable to open file $ofile!";
$I=taxtree2str $tax_blast;
print SUMF $I;
close SUMF;

$ofile=$outdir."summary.tree1.txt";
open SUMF, "> $ofile" or die "\nUnable to open file $ofile!";
$I=taxtree2str $tax_tree1;
print SUMF $I;
close SUMF;

$ofile=$outdir."summary.tree2.txt";
open SUMF, "> $ofile" or die "\nUnable to open file $ofile!";
$I=taxtree2str $tax_tree2;
print SUMF $I;
close SUMF;

