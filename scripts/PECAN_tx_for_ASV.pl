#!/usr/bin/env perl

=head1 NAME

  PECAN_tx_for_ASV.pl

=head1 DESCRIPTION

  Combine the ASV count table resulting from dada2 analyses with PECAN classifications for each ASV.

=head1 SYNOPSIS


=head1 OPTIONS

=over

=item B<--pecan-file, -p>
  A 2-column, PECAN taxonomy table 
  Usually looks like:  MC_order7_results.txt
  (ASV in column 1)

=item B<--count-table, -c>
  The .csv file resulting from the dada2 pipeline. 
  Usually looks like: (PROJECT)_all_runs_dada2_abundance_table.csv

=item B<--vaginal>
  Use this flag if your dataset is from vaginal samples.
  The flag will rename concatenated taxa to the single taxon known for the vaginal environment.
  Example: Lactobacillus_crispatus_Lactobacillus_amylovorus ==> Lactobacillus_crispatus

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE
  <local for Johanna> 
  cd ~/Documents/CONTRA
  PECAN_tx_for_ASV.pl -p MC_order7_results.txt -c CONTRA_all_runs_dada2_abundance_table.csv --vaginal

=cut

use strict;
use warnings;
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
  "PECAN-taxonomy|p=s"  => \my $pecanFile,
  "ASV-count-table|c=s" => \my $countTblFile,
  "vaginal"             => \my $vaginal,
  "quiet"               => \my $quiet,
  "verbose|v"           => \my $verbose,
  "debug"               => \my $debug,
  "dry-run"             => \my $dryRun,
  "help|h!"             => \my $help,
  )
  or pod2usage(verbose => 0,exitstatus => 1);

if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 0;
}

if (!$pecanFile) 
{
  print "Please provide input PECAN taxonomy\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 0;
}

if ($vaginal)
{
  print "---Using vaginal-specific merge rules for taxa\n";
}

####################################################################
##                               MAIN
####################################################################

my %pecan = readTbl($pecanFile);
my %cmbTx;
my $count=0;

foreach my $x (keys %pecan)
{ 
  $count++;
  my $tx = $pecan{$x};
  $cmbTx{$x} = $tx;
  if ($x =~ /GTACGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTACTGATTTGCTTAATTGCACCACATGTGTTTTTCTTTGAAACAAACTTGCTTTGGCGGTGGGCCCAGCCTGCCGCCAGAGGTCTAAACTTACAACCAATTTTTTATCAACTTGTCACACCAGATTATTACTAATAGTCAAAACTTTCAACAACGGATCTCTTGGTTCTCGCATCGATGAAGAACGCAGCGAAATGCGATACGTAATATGAATTGCAGATATTCGTGAATCATCGAGTCTTTGAACCATGAT/)
  {
    $cmbTx{$x} = "Candida_albicans";
  }
}

print "---Combined taxonomy written for $count ASVs to $countTblFile\n";
my @suffixes = (".csv");
my $Prefix = basename($countTblFile, @suffixes);

print "---Adding taxonomy to $countTblFile\n";
my $cntWtx = "$Prefix"."_w_taxa.csv";
open ALL, ">$cntWtx", or die "Cannot open $cntWtx for writing: $OS_ERROR\n";

my $cnttxon = "$Prefix"."_PECAN_taxa_only.csv";
open TXON, ">$cnttxon", or die "Cannot open $cnttxon for writing: $OS_ERROR\n";

my @list1 = get_file_data($countTblFile); 
my $line = 0;
my $head1;
my $head2;
my $i=0;
foreach my $line1 (@list1) 
{
  chomp $line1;
  $line++;

  if ($line == 1) 
  {  
    my @ASVs = split (",", $line1);
    shift @ASVs;
    foreach my $asv (@ASVs) 
    { 
      chomp $asv;
      $asv =~ s/[\r\n]+//;
      $i++;
      if (!exists $cmbTx{$asv}) 
      {
        print "HD not present: #$i $asv\n";
      }
      else
      {
        $head1 .= "$asv\,";
        $head2 .= "$cmbTx{$asv}\,";
      }
    }
  }
}
print ALL ",$head1\n";
print ALL ",$head2\n";
print TXON ",$head2\n";
$line = 0;
foreach my $line1 (@list1) 
{
  chomp $line1;
  $line++;
  if ($line > 1) 
  {
    print ALL "$line1\n";
    print TXON "$line1\n";
  } 
}

close TXON;
close ALL;
print "TOTAL ASVs: $i\n";

@suffixes = (".csv");
$Prefix = basename($cnttxon, @suffixes);
my $merged = $Prefix . "_merged.csv";

if ($debug)
{
  if ($vaginal)
  {
    my $cmd = "/home/jholm/bin/vaginal_combine_tbl_cols.pl -i $cnttxon -o $merged --debug";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }
  else
  {
    my $cmd = "/home/jholm/bin/combine_tbl_cols.pl -i $cnttxon -o $merged --debug";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }
}

if (!$debug)
{
  if ($vaginal)
  {
    my $cmd = "/home/jholm/bin/vaginal_combine_tbl_cols.pl -i $cnttxon -o $merged";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }
  else
  {
    my $cmd = "/home/jholm/bin/combine_tbl_cols.pl -i $cnttxon -o $merged";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }
}

####################################################################
##                               SUBS
####################################################################
sub readTbl{

  my $file = shift;

  if ( ! -f $file )
  {
    warn "\n\n\tERROR in readTbl(): $file does not exist";
    print "\n\n";
    exit 1;
  }

  my %tbl;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    my ($id, $t) = split /\s+/,$_;
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

sub read2colTbl{

  my $file = shift;
  my %tbl;

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    next if $_ eq "";
    my ($id, $t) = split (/\s+/,$_, 2);
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

sub read3colTbl{

  my $file = shift;
  my %tbl;

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    shift;
    next if $_ eq "";
    my ($id, $t) = split (/\s+/,$_, 2);
    $tbl{$id} = $t;
  }
  close IN;

  return %tbl;
}

sub get_file_data {

    my($filename) = @_;

    use strict;
    use warnings;

    # Initialize variables
    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}


exit;