#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $copia_f = shift;
my $gypsy_f = shift;
my $unknown_f = shift;

my %Ltr;

open IN, $copia_f || die "fail open $copia_f\n";
while(<IN>) {
	chomp;
	$_ =~ s/"//g;
	my @t = split /\,/,$_;
	next if($t[0] eq "");
	push @{$Ltr{"Copia"}}, $t[3];
}
close IN;
#print Dumper \%Ltr;


open IN, $gypsy_f || die "fail open $gypsy_f\n";
while(<IN>) {
	chomp;
	$_ =~ s/"//g;
	my @t = split /\,/,$_;
	next if($t[0] eq "");
	push @{$Ltr{"Gypsy"}}, $t[3];
}
close IN;


open IN, $unknown_f || die "fail open $unknown_f\n";
while(<IN>) {
	chomp;
	$_ =~ s/"//g;
	my @t = split /\,/,$_;
	next if($t[0] eq "");
	push @{$Ltr{"unknown"}}, $t[3];
}
close IN;

print "LTR\tInsertion_time\n";
foreach my $type (sort keys %Ltr) {
	my $tp = $Ltr{$type};
	foreach my $tim (@$tp) {
		print "$type\t$tim\n";
	}
}
