#! /usr/bin/perl -w
##########################################################
# pdbchm.pl
# structure charmm minimization on input pdb.
##########################################################
use strict;
use File::Basename;
use Getopt::Long;

my $home = $ENV{'HOME'};

if (scalar @ARGV < 4)
{
        printf ("usage: %s pdb_to ch,ch/- pdb_from ch,ch/-\n", $0);
        exit;
}

my $pdb_to=shift @ARGV;
my ($name, $dir, $suffix) = fileparse ($pdb_to, qr/\.[^.]*/);
my $tnam=$name;
my $ch_to=shift @ARGV;
my $pdb_from=shift @ARGV;
($name, $dir, $suffix) = fileparse ($pdb_from, qr/\.[^.]*/);
my $fnam=$name;
my $ch_fr=shift @ARGV;
my @chf=split('\+',$ch_fr);
my @chs=split('\+',$ch_to);

open FI, ">temp.pml"or die "couldn't open file temp.pml: $!\n";
print FI "load $pdb_to\n";
if($ch_to ne "-")
{
   my $c0=shift @chs;
   print FI "select $tnam and chain $c0";
   for my $c(@chs){print FI "+$c";}
   print FI "\n";
}
else
{
   print FI "select $tnam\n";
}
print FI "cmd.create(None,\"sele\",zoom=0)\n";
print FI "delete sele\n";
print FI "load $pdb_from\n";
if($ch_fr ne "-")
{
   my $cf0=shift @chf;
   print FI "select $fnam and chain $cf0";
   for my $c(@chf){print FI "+$c";}
   print FI "\n";
}
else
{
   print FI "select $fnam\n";
}
print FI "cmd.enable('sele')\n";
print FI "cmd.create(None,\"sele\",zoom=0)\n";
print FI "cmd.align(\"polymer and name ca and (obj02)\",\"polymer and name ca and (obj01)\",quiet=0,object=\"aln_obj02_to_obj01\",reset=1)\n";
print FI "save $fnam\_$ch_fr\_$tnam\_$ch_to.aln, aln_obj02_to_obj01\n";
print FI "cmd.remove(\"(solvent and (all))\")\n";
print FI "save $fnam\_$ch_fr\_$tnam\_$ch_to.pdb,((obj02))\n";
print FI "save $tnam\_$ch_to.pdb,((obj01))\n";
print FI "quit\n";
close FI;
system "pymol -c temp.pml";
system "rm temp.pml";
