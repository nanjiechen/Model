#! /usr/bin/perl -w

##########################################################
# gen_moved_pdb.pl
# -------------------
# copyright : (C) 2006 by Ryan Brenke
# email : rbrenke@gmail.com
##########################################################

use strict;

if (scalar @ARGV < 10)
{
	print "usage: $0 tindx tindy tindz dimx dimy dimz csize orx ory orz \n";
	exit;
}

# command line args
my $tindx = shift @ARGV;
my $tindy = shift @ARGV;
my $tindz = shift @ARGV; # datfile index (0 based)
my $dimx = shift @ARGV;
my $dimy = shift @ARGV;
my $dimz = shift @ARGV;
my $csize =shift @ARGV;
my $orx = shift @ARGV;
my $ory = shift @ARGV;
my $orz = shift @ARGV;

my $ix =-($tindx-$orx)/$csize + $dimx/2;
my $iy =-($tindy-$ory)/$csize + $dimy/2;
my $iz =-($tindz-$orz)/$csize + $dimz/2;

print "$ix $iy $iz\n";
$ix= sprintf "%.0f", $ix;
$iy= sprintf "%.0f", $iy;
$iz= sprintf "%.0f", $iz;

print "$ix $iy $iz\n";
my $tind1=$iz+$dimz*($iy+$dimy*$ix);

print "$tind1\n";
