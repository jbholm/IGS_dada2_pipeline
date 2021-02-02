#!/usr/bin/env perl

=head1 NAME

  combine_tx_for_ASV.pl

=head1 DESCRIPTION

  From the SILVA, PECAN, and other additional files, combine the taxonomic annotations of these multiple files giving
  preference to SILVA classifications, except when taxonomy = "g_Lactobacillus" or "Shuttleworthia", in which case
  the taxonomy from the PECAN file is used. The -x flag allows for an additional 2-column table to be considered after 
  SILVA and PECAN taxonomies are considered (BLAST output, for example). 

=head1 SYNOPSIS


=head1 OPTIONS

=over

=item B<--pecan-file, -p>
  A 2-column, PECAN taxonomy table 
  (ASV in column 1)

=item B<--silva-file, -s>
  The default, multi-level SILVA output table
  (ASV in column 1)

=item B<--homd-file>

=item B<--other-file, -o>
  Any other table with ASV in column 1 and preferred taxonomy in column 2.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE


=cut

use strict;
use warnings;
my $scriptsDir;
my $pipelineDir;

BEGIN {
    use File::Spec::Functions;
    use File::Basename;

    $scriptsDir = dirname(__FILE__);
}
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use List::Util qw( sum );
use Data::Dumper qw(Dumper);
use Cwd;
use File::Basename;

$OUTPUT_AUTOFLUSH = 1;
####################################################################
##                             OPTIONS
####################################################################
GetOptions(
    "PECAN-taxonomy|p=s"      => \my $pecanFile,
    "taxonomy|t=s"            => \my $taxonomy,
    "vaginal"                 => \my $vaginal,
    "ASV-count-table|c=s"     => \my $countTblFile,
    "quiet"                   => \my $quiet,
    "verbose|v"               => \my $verbose,
    "debug"                   => \my $debug,
    "dry-run"                 => \my $dryRun,
    "help|h!"                 => \my $help,
) or pod2usage( verbose => 0, exitstatus => 1 );

if ($help) {
    pod2usage( verbose => 2, exitstatus => 0 );
    exit 0;
}

if ( !$pecanFile ) {
    print "Please provide input taxonomy\n";
    pod2usage( verbose => 2, exitstatus => 0 );
    exit 0;
}
if ( !$taxonomy ) {
    print "Please provide name of the taxonomy\n";
    pod2usage( verbose => 2, exitstatus => 0 );
    exit 0;
}

if ($vaginal) {
    print "---Using vaginal-specific merge rules for taxa\n";
}

####################################################################
##                               MAIN
####################################################################

my %pecan = readTbl($pecanFile);
my %cmbTx;
my $count = 0;

foreach my $x ( keys %pecan ) {
    $count++;
    my $tx = $pecan{$x};
    $cmbTx{$x} = $tx;
}

# Sequence-specific name substitutions
my %subs = (
"GTACGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTACTGATTTGCTTAATTGCACCACATGTGTTTTTCTTTGAAACAAACTTGCTTTGGCGGTGGGCCCAGCCTGCCGCCAGAGGTCTAAACTTACAACCAATTTTTTATCAACTTGTCACACCAGATTATTACTAATAGTCAAAACTTTCAACAACGGATCTCTTGGTTCTCGCATCGATGAAGAACGCAGCGAAATGCGATACGTAATATGAATTGCAGATATTCGTGAATCATCGAGTCTTTGAACCATGAT"
      => "Candida_albicans", 
              "TGAGGAATATTCCACAATGGGCGAAAGCCTGATGGAGCAATGCCGCGTGCAGGATGAAGGCCCTCGGGTCGTAAACTGCTTTTATTAGAGAAGAATATGACGGTAACTAATGAATAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTCATACGTAGGTCCCAAGCGTTATCCGGAGTGACTGGGCGTAAAGAGTTGCGTAGGCGGCTAAGTAAGCGAGTAATGAAAACTATCGGCTCAACCGGTAGCCTGTTATTCGAACTGCTTGGCTCGAGATTATCAGAGGTCGCTGGAATTCCTAGTGTAGCAGTGAAATGCGTAGATATTAGGAAGAACACCAATGGCGTAGGCAGGCGACTGGGGTATTTCTGACGCTAAGGCACGAAAGCGTGGGGAGCGAACCGG" => "BVAB_TM7",
        "GTATTTAGCCGGAGCTTCTTAGTAAGGTACCGTCATTTTCTTCCCTTCTGATAGAGCTTTACATACCGAAATACTTCTTCGCTCACGCGGCGTCGCTGCATCAGGCTTTCGCCCATTGTGCAATATTCCCCACTGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCAGTCTCTCAACTCGGCTACTGATCTTCGCTTTGGTAGGCTTTTACCCCACCAACCGGCTAATCAGACGCGGGTCCATCCTATACCACCGGAGTTTTTCACACCATGTCATGCGACATTCGTGCGCTTATGCGGTATTATCAGCCGTTTCCGGCTGCTATCCCCCGGTACAGGGCAGGTTCCCCACGCGTTACTCACCCGTCCGCCACTAAGTAACTACATCTTCCGTCCGAAAACTTCCGTCGTAGCACTTCGTTCGACTTGCAT" => "BVAB1",
        "GTATTTAGCCGGAGCTTCTTAGTAAGGTACCGTCATTTTCTTCCCTTCTGATAGAGCTTTACATACCGAAATACTTCTTCGCTCACGCGGCGTCGCTGCATCAGGCTTTCGCCCATTGTGCAATATTCCCCACTGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCAGTCTCTCAACTCGGCTACTGATCTTCGCTTTGGTAGGCTTTTACCCCACCAACTGGCTAATCAGACGCGGGTCCATCCTATACCACCGGAGTTTTTCACACCATGTCATGCGACATTCGTGCGCTTATGCGGTATTATCAGCCGTTTCCGGCTGCTATCCCCCGGTACAGGGCAGGTTCCCCACGCGTTACTCACCCGTCCGCCACTAAGTAACTACATCTTCCGTCCGAAAACTTCCGTCGTAGCACTTCGTTCGACTTGCAT" =>"BVAB1",
        "GTATTTAGCCGGAGCTTCTTAGTAAGGTACCGTCATTTTCTTCCCTTCTGATAGAGCTTTACATACCGAAATACTTCTTCGCTCACGCGGCGTCGCTGCATCAGGCTTTCGCCCATTGTGCAATATTCCCCACTGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCAGTCTCTCAACTCGGCTACTGATCTTCGCTTTGGTAGGCTTTTACCCCACC" => "BVAB1",);

if (%subs) {

    # First search the ASV name-sequence reference:
    open my $refFH, "<all_runs_dada2_ASV.fasta";
    my $previous;    # contents of previous line
    while ( my $line = <$refFH> ) {
        chomp $line;
        ( my $match ) = grep ( /$line/, keys %subs );
        if ($match) {
            my $asvName = substr $previous, 1;
            chomp $asvName;
            $cmbTx{$asvName} = $subs{$match};
        }
        $previous = $line;
    }
}

print "---Combined taxonomy written for $count ASVs to $countTblFile\n";
my @suffixes = (".csv");
my $Prefix   = basename( $countTblFile, @suffixes );

print "---Adding taxonomy to $countTblFile\n";
my $cntWtx = "${Prefix}.${taxonomy}.asvs+taxa.csv.tmp";
open ALL, ">$cntWtx", or die "Cannot open $cntWtx for writing: $OS_ERROR\n";

my $cnttxon = "${Prefix}.${taxonomy}.taxa.csv.tmp";
open TXON, ">$cnttxon", or die "Cannot open $cnttxon for writing: $OS_ERROR\n";

my @list1 = get_file_data($countTblFile);
my $line  = 0;
my @asvs;
my @taxa;
my $i = 0;
foreach my $line1 (@list1) {
    chomp $line1;
    $line++;

    if ( $line == 1 ) {
        my @ASVs = split( ",", $line1 );
        shift @ASVs;
        foreach my $asv (@ASVs) {
            chomp $asv;
            $asv =~ s/[\r\n]+//;
            $i++;
            if ( !exists $cmbTx{$asv} ) {
                print "HD not present: #$i $asv\n";
            } else {
                push( @asvs, "$asv" );
                push( @taxa, "$cmbTx{$asv}" );
            }
        }
    }
}
print ALL "," . join( ",", @asvs ) . "\n";
my $taxaHead = "sampleID" . "," . join( ",", @taxa ) . "\n";
print ALL "$taxaHead";
print TXON "$taxaHead";
$line = 0;
foreach my $line1 (@list1) {
    chomp $line1;
    $line++;
    if ( $line > 1 ) {
        print ALL "$line1\n";
        print TXON "$line1\n";
    }
}

close TXON;
close ALL;
print "TOTAL ASVs: $i\n";

@suffixes = (".csv.tmp");
$Prefix   = basename( $cnttxon, @suffixes );
my $merged = $Prefix . "-merged.csv.tmp";

my $debugParam = $debug ? "--debug" : "";
my $vParam = $vaginal ? "--vaginal" : "";

my $cmd =
"$scriptsDir/combine_tbl_cols.pl -i $cnttxon -o $merged $debugParam $vParam";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0
    or die "system($cmd) failed with exit code: $?"
    if !$dryRun;


####################################################################
##                               SUBS
####################################################################
sub readTbl {

    my $file = shift;

    if ( !-f $file ) {
        warn "\n\n\tERROR in readTbl(): $file does not exist";
        print "\n\n";
        exit 1;
    }

    my %tbl;
    open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
    foreach (<IN>) {
        chomp;
        my ( $id, $t ) = split /\s+/, $_;
        $tbl{$id} = $t;
    }
    close IN;

    return %tbl;
}

sub read2colTbl {

    my $file = shift;
    my %tbl;

    open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
    foreach (<IN>) {
        chomp;
        next if $_ eq "";
        my ( $id, $t ) = split( /\s+/, $_, 2 );
        $tbl{$id} = $t;
    }
    close IN;

    return %tbl;
}

sub read3colTbl {

    my $file = shift;
    my %tbl;

    open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
    foreach (<IN>) {
        chomp;
        shift;
        next if $_ eq "";
        my ( $id, $t ) = split( /\s+/, $_, 2 );
        $tbl{$id} = $t;
    }
    close IN;

    return %tbl;
}

sub get_file_data {

    my ($filename) = @_;

    use strict;
    use warnings;

    # Initialize variables
    my @filedata = ();

    unless ( open( GET_FILE_DATA, $filename ) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}

exit;
