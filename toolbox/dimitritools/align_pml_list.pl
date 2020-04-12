#! /usr/bin/perl -w
##########################################################
# pdbchm.pl
# structure charmm minimization on input pdb.
##########################################################
use strict;
use File::Basename;
use Getopt::Long;

my $aln="/home/dbeglov/cryptosite_3/pml_align_sal.pl";
my $home = $ENV{'HOME'};

if (scalar @ARGV < 3)
{
        printf ("usage: %s pdb_to ch list\n", $0);
        exit;
}

my $pdb_to=shift @ARGV;
my $ch_to=shift @ARGV;
my $list=shift @ARGV;

open FI, "< $list" or die "couldn't open file $list: $!\n";
my @lines = <FI>;
chomp @lines;
close FI;

for my $l(@lines)
{
    my @ar=split('\t',$l);
    my $pdbto=lc($ar[1]);
    $pdbto=~s/ //g;
    $pdbto=$pdbto.".pdb";
    my $chains=$ar[3];
    $chains=~s/,//g;
    my @br=split(' ',$chains);
    for my $b(@br)
    {
       print $pdbto." ".$b."\n";
       system "$aln $pdb_to $ch_to $pdbto $b"; 
       
    }
}
