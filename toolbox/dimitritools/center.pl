#!/usr/bin/perl

if (scalar @ARGV < 2)
{
	print "usage: $0  PDB1_TO_CENTER  PDB1_OUT  (optional args: PDB2_TO_TRANSLATE_SIMILARLY  PDB2_OUT)\n";
	exit;
}

my $pdb1 = shift @ARGV;
my $pdb1out = shift @ARGV;
my $pdb2 = shift @ARGV;
my $pdb2out = shift @ARGV;


open (PDB1, "< $pdb1") or die "couldn't open $pdb1: $!\n";
my @lines = <PDB1>;
close PDB1;

my $sumx = 0;
my $sumy = 0;
my $sumz = 0;

my $i = 0;
for my $line (@lines)
{
	if ($line =~ /ATOM/)
	{
		$x1[$i] = substr($line, 30, 8);
		$y1[$i] = substr($line, 38, 8);
		$z1[$i] = substr($line, 46, 8);
		$x1[$i] =~ s/\s+//g;
		$y1[$i] =~ s/\s+//g;
		$z1[$i] =~ s/\s+//g;
		$begin[$i] = substr($line, 0, 30);
		$end[$i] = substr($line, 54, 50);
		$sumx += $x1[$i];
		$sumy += $y1[$i];
		$sumz += $z1[$i];
		$count++;
		$i++;
	}
}


$x0 = $sumx/$count;
$y0 = $sumy/$count;
$z0 = $sumz/$count;

print "$x0 $y0 $z0\n";

$i = 0;
open (PDB1OUT, "> $pdb1out") or die "couldn't open $pdb1out: $1\n";
for my $line (@lines)
{
	if ($line =~ /ATOM/)
	{
		$x1[$i] -= $x0;
		$y1[$i] -= $y0;
		$z1[$i] -= $z0;
		$end[$i] =~ s/\n//g;
		printf PDB1OUT ("%30s%8.3f%8.3f%8.3f%-50s\n", $begin[$i], $x1[$i], $y1[$i], $z1[$i], $end[$i]);
		$i++;
	}
}
close PDB1OUT;


exit if (! defined $pdb2);
if (! defined $pdb2out)
{
	print "not enough args for pdb2out\n";
	exit;
}

open (PDB2, "< $pdb2") or die "couldn't open $pdb2: $!\n";
@lines = <PDB2>;
close PDB2;

open (PDB2OUT, "> $pdb2out") or die "couldn't open $pdb2out: $1\n";
for my $line (@lines)
{
	if ($line =~ /ATOM/)
	{
		$x1 = substr($line, 30, 8);
		$y1 = substr($line, 38, 8);
		$z1 = substr($line, 46, 8);
		$x1 =~ s/\s+//g;
		$y1 =~ s/\s+//g;
		$z1 =~ s/\s+//g;
		$begin = substr($line, 0, 30);
		$end = substr($line, 54, 50);

		$x1 -= $x0;
		$y1 -= $y0;
		$z1 -= $z0;
		$end =~ s/\n//g;
		printf PDB2OUT ("%30s%8.3f%8.3f%8.3f%-50s\n", $begin, $x1, $y1, $z1, $end);
	}
}
close PDB2OUT;
