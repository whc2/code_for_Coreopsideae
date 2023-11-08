#!/usr/bin/perl -w
use strict;
use Data::Dumper;

# This script is used to calculate coverage of contigs. The resolution should be set before calculation.
# Then, contigs can be classified into four classes: 1x, 2x, 3x, 4x for a tetraploid species.

my $len_f = shift;
my $bed_f = shift;
#my $resolution = shift;	# bin size used to compute coverage for a bin.
#my $coverage = shift;	# whole genome coverage

my %Len;
open IN, $len_f || die "fail open $len_f\n";
while(<IN>) {
	chomp;
	my ($id, $len) = split /\t/,$_;
	$Len{$id} = $len;
}
close IN;

my %Bed;
open IN, $bed_f || die "fail open $bed_f\n";
while(<IN>) {
	chomp;
	my ($id, $beg, $end, $cov) = split /\t/,$_;	# begin, end, coverage
	$Bed{$id} += $cov * ($end - $beg);
}
close IN;

foreach my $ctg (sort keys %Len) {
	my $ctg_len = $Len{$ctg};
	my $cov_total = $Bed{$ctg};
	my $mean_cov = $cov_total / $ctg_len;
	print "$ctg\t$ctg_len\t$mean_cov\n";
}
