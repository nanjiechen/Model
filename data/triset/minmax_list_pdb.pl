#! /usr/bin/perl -w
##########################################################
# pdbchm.pl
# structure charmm minimization on input pdb.
##########################################################
use strict;
use File::Basename;
use Getopt::Long;

my $home = $ENV{'HOME'};

if (scalar @ARGV != 1)
{
        printf ("usage: %s pdb_list\n", $0);
        exit;
}

my $list=shift @ARGV;

open FI, "< $list" or die "couldn't open file $list: $!\n";
my @lines = <FI>;
chomp @lines;
close FI;

my $xmin=999999;
my $xmax=-999999;
my $ymin=$xmin;
my $ymax=$xmax;
my $zmin=$xmin;
my $zmax=$xmax;
for my $l1(@lines)
{
    open FI, "< $l1" or die "couldn't open file $l1: $!\n";
    my @lines2 = <FI>;
    chomp @lines2;
    close FI;
    for my $l(@lines2)
    {
       if(substr($l,0,4) eq "ATOM")
       {
	    my $x=substr($l,30,8);
	    if($x<$xmin){$xmin=$x;}
	    if($x>$xmax){$xmax=$x;}
	    my $y=substr($l,38,8);
	    if($y<$ymin){$ymin=$y;}
            if($y>$ymax){$ymax=$y;}
	    my $z=substr($l,46,8);
	    if($z<$zmin){$zmin=$z;}
            if($z>$zmax){$zmax=$z;}
       }
    }
}
$xmin=sprintf("%8.3f",$xmin);
print "ATOM      1  CA  ALA A   1    $xmin$ymin$zmin  1.00  0.00      A\n";
print "ATOM      2  CA  ALA A   2    $xmax$ymax$zmax  1.00  0.00      A\n";
print "END\n";
