#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $gff_f = shift;

my %LTR;

open IN, $gff_f || die "fail open $gff_f\n";
while(<IN>) {
	chomp;
	next if($_ !~ /long_terminal_repeat/);

	# ID=lLTR_1;Parent=repeat_region_1;
	my @t = split /\t/,$_;
	my ($start, $end, $strand) = ($t[3]-1, $t[4]-1, $t[6]);

	my $id = $1 if($t[8] =~ /^ID=(\S+?);/);
	my $class = $1 if($t[8] =~ /Classification=(\S+?);/);

	my $proto = $id;
	$proto =~ s/^\w//g;

	$id =~ s/\_//g;
	push @{$LTR{$class}{$proto}}, join("_", $t[0], $start, $end, $id, 1, $strand);
}
close IN;
#print Dumper \%LTR;

foreach my $class (keys %LTR) {
	my $tmp = $LTR{$class};

	foreach my $ltr (keys %$tmp) {
		my $filename = $class . "_" . $ltr;
		$filename =~ s/\//\_/;

		open OUT, ">$filename.bed" || die "fail open $filename\n";
		foreach my $pair (@{$LTR{$class}{$ltr}}) {
			$pair =~ s/\_/\t/g;
			print OUT "$pair\n";
		}
		close OUT;
	}
}
