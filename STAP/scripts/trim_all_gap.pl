#!/usr/bin/perl


use strict;

my %opt=@ARGV;
my $opt_i=$opt{'-i'};
my $opt_o=$opt{'-o'};

my %mask;
my $length;

my $temp_seq;
my $temp_name;
my $seq;
my $name;
open(OUT, ">$opt_o") || die "cannot output to file $opt_o \n";
open(IN,$opt_i) || die "cannot open file $opt_i \n";
while(<IN>){
    if(/>/){
	$name=$temp_name;
        $seq=$temp_seq;
        $temp_name=$_;
        $temp_seq='';
        if($seq){&add_mask();} 
    }
    else{
        $_=~s/\s+//g; 
	$temp_seq.=$_;
    }
}
close IN;
	$name=$temp_name;
        $seq=$temp_seq;       
        if($seq){&add_mask();} 



$temp_seq='';
$temp_name='';
$seq='';
$name='';

open(IN,$opt_i) || die "cannot open file $opt_i \n";
while(<IN>){
    if(/>/){
	$name=$temp_name;
        $seq=$temp_seq;
        $temp_name=$_;
        $temp_seq='';
        if($seq){&output();} 
    }
    else{
        $_=~s/\s+//g; 
	$temp_seq.=$_;
    }
}
close IN;
	$name=$temp_name;
        $seq=$temp_seq;       
        if($seq){&output();} 

close OUT;


sub output{
    print OUT $name;
    my @t=split(//,$seq);
    my $i;
for $i(0..$length-1){
    if($mask{$i}){print OUT $t[$i];}
}
print OUT "\n";   
}



sub add_mask{
   
    my @t=split(//,$seq);

    
    my $count;

    if($length){$count=$length;}
    else{
    $length=@t;
    $count=$length; 
    }
    $count--;
    my $i;
    
    
    for $i(0..$count){

	if($t[$i]=~/[A-Za-z]/){
	    $mask{$i}=1;
	}



    }

   undef @t;
   undef $seq;  

}






