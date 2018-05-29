#!/usr/bin/perl


while(<>){
if(/>/){
($acc)=split(/\s+/);
$acc=~s/[^A-Za-z0-9-_]//g;
$acc.=".fa";
close OUT;
open(OUT,">$acc");
 }
print OUT $_;

}

