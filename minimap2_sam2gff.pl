#!/usr/bin/perl

=head1 Name



=head1 Description



=head1 Version

  Author: Hengchao Wang, hengchaobioinf@gmail.com
  Version: 1.0,  Date: 2024-2-26
  Note:

=head1 Usage

  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple



=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my ($Verbose,$Help);
GetOptions(
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);

my $hisat2_bam = shift;

my $mark = 1;

my $total_bases;
my %Genes;

#open IN, "samtools view $hisat2_bam |" || die "fail ";
open IN, $hisat2_bam  || die "fail ";
while (<IN>) {
	next if($_ =~ /^@/);
	my ($gene_id, $gene_start, $map_cigar) = ($1,$2,$3) if(/^\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)/);
	next if($map_cigar =~ /\*/);
	next if($gene_id =~ /\*/);

	my @mat = $map_cigar =~ /(\d+)M/g;
	my $total_len = 0;

	foreach my $frag (@mat) {
		$total_len += $frag;
		$total_bases += $frag;
	}

	if(exists $Genes{$gene_id}) {
		$Genes{$gene_id} += $total_len;
	} else {
		$Genes{$gene_id} = $total_len;
	}
}
close IN;

print STDERR "Total_bases: $total_bases\n";

foreach my $g (sort keys %Genes) {
	my $bases = $Genes{$g};
	print "$g\t$bases\n";
}

####################################################
################### Sub Routines ###################
####################################################
