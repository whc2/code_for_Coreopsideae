#!/usr/bin/perl

=head1 Name



=head1 Description

This script uses tsv file as input to calculate the expression of genes in TPM.

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

my $bases_count_file = shift;
my $len_file = shift;
my $tsv_file = shift;

##get the Total_bases value
my $Total_bases;
open IN, $bases_count_file || die "fail";
while (<IN>) {
	$Total_bases = $1 if(/Total_bases:\s+(\d+)/);
}
close IN;

##get CDS length for each gene
my %GeneLen;
read_len($len_file,\%GeneLen);
#print Dumper \%CDSlen;

#get expression count for each gene
my %GeneExpress;
open IN, $tsv_file || die "fail";
while (<IN>) {
	next if(/^\#/);
	chomp;
    my ($gene, $bases) = split /\t/,$_;
    $GeneExpress{$gene} = $bases;
}
close IN;
#print Dumper \%GeneExpress;

my $Sum_tpk;
my %Tpk;
print "gene_id\tcds_len\tTotal_bases\tTPM\n";
foreach my $gene_id (sort keys %GeneExpress) {
	my $expression_count = $GeneExpress{$gene_id};
	my $cds_len = $GeneLen{$gene_id};

	my $tpk = $expression_count / ($cds_len / 1000);
	$Tpk{$gene_id} = $tpk;
	$Sum_tpk += $tpk / 1000000;
#	print "$gene_id\t$cds_len\t$Total_bases\t$expression_count\t$rpkm\n";
}

foreach my $gene_id (sort keys %GeneLen) {
	my $cds_len = $GeneLen{$gene_id};

    if(exists $Tpk{$gene_id}) {
        my $tpk = $Tpk{$gene_id};
        my $tpm = $tpk / $Sum_tpk;
        print "$gene_id\t$cds_len\t$Total_bases\t$tpm\n";
    } else {
        print "$gene_id\t$cds_len\t$Total_bases\t0\n";
    }
}

######################################

sub read_len{
	my $file=shift;
	my $hash=shift;

	open (IN,$file) || die ("fail open $file\n");
	while (<IN>) {
		chomp;
        my ($g, $len) = split /\t/,$_;
        $hash->{$g} = $len;
	}
	close(IN);	
}
