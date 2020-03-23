#!/usr/bin/env perl

=head1 NAME

  combine_tbl_cols.pl

=head1 DESCRIPTION

 combining columns of a table that have the same names.

 Rename the following taxa:
  s/^Paenibacillus_Cluster_1/Paenibacillus_polymyxa_Paenibacillus_peoriae_Paenibacillus_borealis_Paenibacillus_durus_Paenibacillus_barcinonensis_Paenibacillus_chibensis_Paenibacillus_glacialis_Paenibacillus_lautus_Paenibacillus_amylolyticus_Paenibacillus_graminis_Paenibacillus_agarexedens_Paenibacillus_lactis_Paenibacillus_nanensis_Paenibacillus_cineris_Paenibacillus_xylanexedens_Paenibacillus_cookii/g for @colNames;
  s/^Bacillus_Cluster_1/Bacillus_subtilis_Bacillus_vallismortis_Bacillus_mojavensis_Bacillus_amyloliquefaciens_Bacillus_tequilensis_Bacillus_licheniformis_Bacillus_sonorensis_Bacillus_velezensis/g for @colNames;
  s/^Pseudomonas_Cluster_1/Pseudomonas_reactans_Pseudomonas_fluorescens_Pseudomonas_salomonii_Pseudomonas_moraviensis_Pseudomonas_veronii_Pseudomonas_arsenicoxydans_Pseudomonas_koreensis_Pseudomonas_mandelii_Pseudomonas_marginalis_Pseudomonas_gessardii_Pseudomonas_jessenii/g for @colNames;
  s/^Lactobacillus_jensenii_1/Lactobacillus_jensenii/g
  s/^Lactobacillus_jensenii_3/Lactobacillus_jensenii/g
For vaginal samples, merge together the following taxa:
  s/^Lactobacillus_crispatus_Lactobacillus_helveticus/Lactobacillus_crispatus/g
  s/^Lactobacillus_gasseri_Lactobacillus_johnsonii/Lactobacillus_gasseri/g
  s/^Lactobacillus_gasseri_johnsonii/Lactobacillus_gasseri/g
  s/^Proteobacteria_bacterium/TM7_BVAB/g
  s/^Lactobacillus_kitasatonis/Lactobacillus_crispatus/g
  s/^Lactobacillus_acidophilus/Lactobacillus_crispatus/g
  s/^f_Beggiatoaceae/BVAB_TM7/g
  s/^Lactobacillus_kitasatonis_Lactobacillus_gallinarum_Lactobacillus_crispatus/Lactobacillus_crispatus/g
  s/^Atopobium_sp_Atopobium_parvulum_Atopobium_rimae_Atopobium_vaginae/Atopobium_vaginae/g
  s/^Lactobacillus_acidophilus_Lactobacillus_sp_Lactobacillus_taiwanensis_Lactobacillus_johnsonii/Lactobacillus_johnsonii/g
  s/^Lactobacillus_crispatus_helveticus/Lactobacillus_crispatus/g 
=head1 SYNOPSIS

  combine_tbl_cols.pl -i <input file> -o <output file> [Options]

=head1 OPTIONS

=over

=item B<--input-file, -i>
  Input file.

=item B<--output-file, -o>
  Output file.

=item B<--verbose, -v>
  Prints content of some output files.

=item B<--debug>
  Prints system commands

=item B<--dry-run>
  Print commands to be executed, but do not execute them.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE

  cd /Users/pgajer/projects/M_and_M/data

  combine_tbl_cols.pl -i MM_dada2_abundance_table_no_ASV.csv -o MM_dada2_abundance_table_spp.csv

=cut

use strict;
use warnings;
use diagnostics;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

# use Cwd qw(abs_path);
# use File::Temp qw/ tempfile /;
use List::Util qw( sum );

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################
my $vaginal = 0;
GetOptions(
    "input-file|i=s"  => \my $inFile,
    "output-file|o=s" => \my $outFile,
    "vaginal!"        => \$vaginal,
    "verbose|v"       => \my $verbose,
    "debug"           => \my $debug,
    "dry-run"         => \my $dryRun,
    "help|h!"         => \my $help,
) or pod2usage( verbose => 0, exitstatus => 1 );

if ($help) {
    pod2usage( verbose => 2, exitstatus => 0 );
    exit 1;
}

if ( !$inFile ) {
    print "\n\n\tERROR: Missing input file\n\n";
    pod2usage( verbose => 2, exitstatus => 0 );
    exit 1;
} elsif ( !$outFile ) {
    print "\n\n\tERROR: Missing output file\n\n";
    pod2usage( verbose => 2, exitstatus => 0 );
    exit 1;
}

if ( !-e $inFile ) {
    warn "\n\n\tERROR: $inFile does not exist";
    print "\n\n";
    exit;
}

####################################################################
##                               MAIN
####################################################################

open( IN, "<", "$inFile" ) or die "Cannot open $inFile for reading: $OS_ERROR";
my $header = <IN>;
chomp $header;
my ( $id, @colNames ) = split ",", $header;

## Substitution block --- alter known taxon names within the column header array to improve merging##
## The order of substitution is important --- longest names first and then short.
if ($vaginal) {
s/^Lactobacillus_crispatus_Lactobacillus_helveticus$/Lactobacillus_crispatus/g
      for @colNames;
    s/^Lactobacillus_gasseri_Lactobacillus_johnsonii$/Lactobacillus_gasseri/g
      for @colNames;
    s/^Lactobacillus_gasseri_johnsonii$/Lactobacillus_gasseri/g for @colNames;
    s/^Proteobacteria_bacterium$/BVAB_TM7/g                     for @colNames;
    s/^Lactobacillus_kitasatonis$/Lactobacillus_crispatus/g     for @colNames;
    s/^Lactobacillus_acidophilus$/Lactobacillus_crispatus/g     for @colNames;

    s/^f_Beggiatoaceae$/BVAB_TM7/g for @colNames;
s/^Lactobacillus_kitasatonis_Lactobacillus_gallinarum_Lactobacillus_crispatus$/Lactobacillus_crispatus/g
      for @colNames;
s/^Atopobium_sp_Atopobium_parvulum_Atopobium_rimae_Atopobium_vaginae$/Atopobium_vaginae/g
      for @colNames;
s/^Lactobacillus_acidophilus_Lactobacillus_sp_Lactobacillus_taiwanensis_Lactobacillus_johnsonii$/Lactobacillus_johnsonii/g
      for @colNames;
    s/^Lactobacillus_crispatus_helveticus$/Lactobacillus_crispatus/g
      for @colNames;
}
s/^Paenibacillus_polymyxa_Paenibacillus_peoriae_Paenibacillus_borealis_Paenibacillus_durus_Paenibacillus_barcinonensis_Paenibacillus_chibensis_Paenibacillus_glacialis_Paenibacillus_lautus_Paenibacillus_amylolyticus_Paenibacillus_graminis_Paenibacillus_agarexedens_Paenibacillus_lactis_Paenibacillus_nanensis_Paenibacillus_cineris_Paenibacillus_xylanexedens_Paenibacillus_cookii$/Paenibacillus_Cluster_1/g
  for @colNames;
s/^Bacillus_subtilis_Bacillus_vallismortis_Bacillus_mojavensis_Bacillus_amyloliquefaciens_Bacillus_tequilensis_Bacillus_licheniformis_Bacillus_sonorensis_Bacillus_velezensis$/Bacillus_Cluster_1/g
  for @colNames;
s/^Pseudomonas_reactans_Pseudomonas_fluorescens_Pseudomonas_salomonii_Pseudomonas_moraviensis_Pseudomonas_veronii_Pseudomonas_arsenicoxydans_Pseudomonas_koreensis_Pseudomonas_mandelii_Pseudomonas_marginalis_Pseudomonas_gessardii_Pseudomonas_jessenii$/Pseudomonas_Cluster_1/g
  for @colNames;
s/^Lactobacillus_jensenii_1$/Lactobacillus_jensenii/g for @colNames;
s/^Lactobacillus_jensenii_3$/Lactobacillus_jensenii/g for @colNames;

my @uqTxs = unique( \@colNames );
my %txIdx;
for my $tx (@uqTxs) {
    my @idx = grep { $tx eq $colNames[$_] } 0 .. $#colNames;
    push @{ $txIdx{$tx} }, @idx;
}

if ($debug) {
    print "txIdx\n";
    for my $tx ( keys %txIdx ) {
        print "$tx: ";
        my @idx = @{ $txIdx{$tx} };
        print_array( \@idx );
    }
}

open OUT, ">$outFile" or die "Cannot open $outFile for writing: $OS_ERROR\n";
print OUT "sampleID";
map { print OUT ",$_" } @uqTxs;
print OUT "\n";

while (<IN>) {
    chomp;
    my ( $id, @counts ) = split ",";
    print OUT $id;
    for my $tx (@uqTxs) {
        my @idx = @{ $txIdx{$tx} };
        print OUT "," . sum( @counts[@idx] );
    }
    print OUT "\n";
}
close IN;
close OUT;

print "Output written to $outFile\n";

####################################################################
##                               SUBS
####################################################################

# print array to stdout
sub print_array {
    my $a = shift;
    map { print "$_ " } @{$a};
    print "\n";
}

# extract unique elements from an array
sub unique {
    my $a = shift;
    my %saw;
    my @out = grep( !$saw{$_}++, @{$a} );

    return @out;
}

exit 0;
