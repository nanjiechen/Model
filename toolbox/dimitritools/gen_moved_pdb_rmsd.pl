#! /usr/bin/perl -w

##########################################################
# gen_moved_pdb.pl
# -------------------
# copyright : (C) 2006 by Ryan Brenke
# email : rbrenke@gmail.com
##########################################################

use strict;

if (scalar @ARGV < 13)
{
	print "usage: $0 ROTPRM rotind tind dimx dimy dimz csize orx ory orz PDB_IFILE PDB_OFILE PDB_REF\n";
	exit;
}

# command line args
my $rotprm = shift @ARGV;
my $rotind = shift @ARGV;
my $tind = shift @ARGV; # datfile index (0 based)
my $dimx = shift @ARGV;
my $dimy = shift @ARGV;
my $dimz = shift @ARGV;
my $csize =shift @ARGV;
my $orx = shift @ARGV;
my $ory = shift @ARGV;
my $orz = shift @ARGV;
my $ligand = shift @ARGV;
my $ofile = shift @ARGV; # moved pdb ofile (- is STDOUT)
my $refile= shift @ARGV;

if ($ligand eq $ofile)
{
	print "PDB_IFILE is identical to PDB_OFILE, exiting\n";
	exit;
}

# get relevant values
my $roti = $rotind; # rotation index into rotprm
my $xc=$dimx/2;
my $yc=$dimy/2;
my $zc=$dimz/2;
my $iz=$tind % $dimz;
my $iy=(($tind-$iz)/$dimz) % $dimy;
my $ix=(($tind-$iz)/$dimz - $iy)/$dimy;
my $tind1=$iz+$dimz*($iy+$dimy*$ix);
#print "$tind $tind1\n";
my $tx = ($ix-$xc)*$csize;
my $ty = ($iy-$yc)*$csize;
my $tz = ($iz-$zc)*$csize;
$tx = -$tx +$orx;
$ty = -$ty +$ory;
$tz = -$tz +$orz;
#print "$tx, $ty, $tz\n";

# read in rotprm
open (ROTFILE, "< $rotprm") or die "could not open $rotprm: $!\n";
my @rotlines = <ROTFILE>;
close ROTFILE;
if (scalar @rotlines < 1)
{
	die "error: $rotprm contains no data\n";
}
my $rotline = $rotlines[$roti]; # get rotation roti
#print "$rotline\n";
$rotline =~ s/^\s+//; # remove leading whitespace
# get values
my $rotseti;
my $a11; my $a12; my $a13;
my $a21; my $a22; my $a23;
my $a31; my $a32; my $a33;
($rotseti, $a11,$a12,$a13,$a21,$a22,$a23,$a31,$a32,$a33) = split(/\s+/, $rotline);
#print "$a11 $a12 $a13 $a21 $a22 $a23 $a31 $a32 $a33\n";

my $tmppdb="$$.pdb";
my $vs_line = `center.pl $ligand $tmppdb`;
#print "$vs_line\n";
my @vs=split /\s+/, $vs_line;
# read in ligand pdb file
open (LIG, "< $tmppdb") or die "couldn't open $ligand: $!\n";
my @liglines = <LIG>;
close LIG;
chomp @liglines; # remove newlines
system "rm $tmppdb";

open(REF, "< $refile") or die "couldn't open $refile: $!\n";
my @relines = <REF>;
close REF;
chomp @relines; # remove newlines
my @xref;
my @yref;
my @zref;
my @aaref;
my $iaref=0;
for my $line (@relines)
{
	if ($line =~/^ATOM/)
        {
		$xref[$iaref]=substr ($line, 30, 8);
		$yref[$iaref]=substr ($line, 38, 8);
		$zref[$iaref]=substr ($line, 46, 8);
		my $a=substr ($line, 12,4);
		$a=~s/ //g;
		$aaref[$iaref]=$a;

		$iaref++;
	}
}


# print the ofile
open (OFILE, "> $ofile") or die "could not open $ofile: $!\n";
my $ind=0;
my $rmsd=0;
my $hind=0;
for my $line (@liglines)
{
	next if ($line !~ /^ATOM/);

	my $prefix = substr ($line, 0, 30);
	my $suffix = substr ($line, 54);

	my $x = substr ($line, 30, 8);
	my $y = substr ($line, 38, 8);
	my $z = substr ($line, 46, 8);
#print "$x, $y, $z\n";

	# rotate and translate
	my $newx = $a11*$x + $a12*$y + $a13*$z + $tx + $vs[0];
	my $newy = $a21*$x + $a22*$y + $a23*$z + $ty + $vs[1];
	my $newz = $a31*$x + $a32*$y + $a33*$z + $tz + $vs[2];

	printf OFILE ("%30s%8.3f%8.3f%8.3f%-50s\n", $prefix, $newx, $newy, $newz, $suffix);
	my $at=substr ($line,12,4);
	$at=~s/ //g;
	if(substr ($at,0,1) ne "H")
	{
		if($at ne $aaref[$ind])
		{
			print "atom mismatch with ref $at $aaref[$ind]\n";
			exit;
		}
		#		print "atom $at\n";
		my $d=$newx-$xref[$ind];
		$rmsd=$rmsd+$d*$d;
		$d=$newy-$yref[$ind];
		$rmsd=$rmsd+$d*$d;
		$d=$newz-$zref[$ind];
                $rmsd=$rmsd+$d*$d;
		$hind++;
	}
	$ind++;
}
$rmsd=sqrt($rmsd/$hind);
printf OFILE ("END\n");


close OFILE;

print "rmsd= $rmsd $hind\n";

