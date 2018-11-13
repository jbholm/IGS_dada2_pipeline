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

=item B<--other-file, -o>
  Any other table with ASV in column 1 and preferred taxonomy in column 2.

=item B<-h|--help>
  Print help message and exit successfully.

=back


=head1 EXAMPLE


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
  "SILVA-taxonomy|s=s"  => \my $silvaFile,
  "UNITE-taxonomy|u=s"  => \my $uniteFile,
  "EZBio-taxonomy|e=s"  => \my $ezbioFile,
  "vaginal"             => \my $vaginal,
  "other-taxonomy|o=s"  => \my $otherFile,
  "ASV-count-table|c=s" => \my $countTblFile,
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

if ($pecanFile) 
{
  print "---Using PECAN taxonomy\n";
}

if ($uniteFile) 
{
  print "---Using UNITE taxonomy\n";
}

if ($silvaFile) 
{
  print "---Using SILVA taxonomy\n";
}

####################################################################
##                               MAIN
####################################################################
my %silva;
my %unite;
my %pecan;
my %other;
my %ezbio;

if ($silvaFile)
{
  %silva = readSILVATbl($silvaFile);
  my $silvaCond = "silva_condensed.txt";
  print "---Writing condensed SILVA taxonomy to $silvaCond\n";
  open OUT, ">$silvaCond", or die "Cannot open $silvaCond for writing: $OS_ERROR\n";
  foreach my $x (keys %silva)
  {
    print OUT "$x\t$silva{$x}\n";
  }
  close OUT;
}
if ($ezbioFile)
{
  %ezbio = readEzBioTbl($ezbioFile);
  my $ezbioCond = "ezbio_condensed.txt";
  print "---Writing condensed EzBio taxonomy to $ezbioCond\n";
  open OUT, ">$ezbioCond", or die "Cannot open $ezbioCond for writing: $OS_ERROR\n";
  foreach my $x (keys %ezbio)
  {
    print OUT "$x\t$ezbio{$x}\n";
  }
  close OUT;
}
if ($uniteFile)
{
  %unite = readUNITETbl($uniteFile);
  my $uniteCond = "unite_condensed.txt";
  print "---Writing condensed UNITE taxonomy to $uniteCond\n";
  open OUT, ">$uniteCond", or die "Cannot open $uniteCond for writing: $OS_ERROR\n";
  foreach my $x (keys %unite)
  {
    print OUT "$x\t$unite{$x}\n";
  }
  close OUT;
}
if ($pecanFile)
{
  %pecan = readTbl($pecanFile);
  print "---Using SILVA taxonomy except for species level classifications (PECAN)\n";
}
if ($otherFile)
{
  %other = readOtherTbl($otherFile);
  print "---Including $otherFile taxonomy\n";
  my $otherCond = "other_condensed.txt";
  print "---Writing condensed other taxonomy to $otherCond\n";
  open OUT, ">$otherCond", or die "Cannot open $otherCond for writing: $OS_ERROR\n";
  if ($otherFile)
  {
    foreach my $x (keys %other)
  {
    my $tx = "NA";
    if ($other{$x})
    {
      print OUT "$x\t$other{$x}\n";
    }
    else
    {
      print OUT "$x\t$tx\n";
    }
  }
    close OUT;
  }
}

my %cmbTx;
my $count=0;
my $countspp=0;

print "---Combining taxonomic assignments of all source files (includes " . scalar(keys %silva)." ASVs)\n";

if($silvaFile)
{
  foreach my $x (keys %silva)
  { 
    $count++;
    my $tx = $silva{$x};
    #print "The Silva taxonomy of $x is $tx\n";
    $cmbTx{$x} = $tx;

    if ($pecanFile)
    {
      if ($pecan{$x} && $pecan{$x} !~ /d_/ && $pecan{$x} !~ /p_/ && $pecan{$x} !~ /c_/ && $pecan{$x} !~ /o_/ && $pecan{$x} !~ /f_/ && $pecan{$x} !~ /g_/ )
      {
        ##print "The PECAN taxonomy of $x is $pecan{$x}\n";
        $countspp++;
        $tx = $pecan{$x};
        $cmbTx{$x} = $tx;
      }
    }
    if ($ezbioFile)
    {
      if ($ezbio{$x} && $ezbio{$x} !~ /d_/ && $ezbio{$x} !~ /p_/ && $ezbio{$x} !~ /c_/ && $ezbio{$x} !~ /o_/ && $ezbio{$x} !~ /f_/ && $ezbio{$x} !~ /g_/ )
      {
        ##print "The PECAN taxonomy of $x is $pecan{$x}\n";
        $countspp++;
        $tx = $ezbio{$x};
        $cmbTx{$x} = $tx;
      }
    }
    if ($otherFile)
    {
      $tx = $other{$x};
      if ($tx)
      {
        $cmbTx{$x} = $tx;
      }
    }
    if ($x =~ /GTACGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTACTGATTTGCTTAATTGCACCACATGTGTTTTTCTTTGAAACAAACTTGCTTTGGCGGTGGGCCCAGCCTGCCGCCAGAGGTCTAAACTTACAACCAATTTTTTATCAACTTGTCACACCAGATTATTACTAATAGTCAAAACTTTCAACAACGGATCTCTTGGTTCTCGCATCGATGAAGAACGCAGCGAAATGCGATACGTAATATGAATTGCAGATATTCGTGAATCATCGAGTCTTTGAACCATGAT/)
    {
      $cmbTx{$x} = "Candida_albicans";
    }
    if ($x =~ /TGAGGAATATTCCACAATGGGCGAAAGCCTGATGGAGCAATGCCGCGTGCAGGATGAAGGCCCTCGGGTCGTAAACTGCTTTTATTAGAGAAGAATATGACGGTAACTAATGAATAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTCATACGTAGGTCCCAAGCGTTATCCGGAGTGACTGGGCGTAAAGAGTTGCGTAGGCGGCTAAGTAAGCGAGTAATGAAAACTATCGGCTCAACCGGTAGCCTGTTATTCGAACTGCTTGGCTCGAGATTATCAGAGGTCGCTGGAATTCCTAGTGTAGCAGTGAAATGCGTAGATATTAGGAAGAACACCAATGGCGTAGGCAGGCGACTGGGGTATTTCTGACGCTAAGGCACGAAAGCGTGGGGAGCGAACCGG/)
    {
      $cmbTx{$x} = "BVAB_TM7";
    }
  }
}
print "\n$countspp ASVs assigned species via PECAN\n";
if ($uniteFile)
{
  foreach my $x (keys %unite)
  { 
    $count++;
    my $tx = $unite{$x};
    #print "The Silva taxonomy of $x is $tx\n";
    $cmbTx{$x} = $tx;
  }
}
my $cmbTblFile = "cmb_tx.txt";
open OUT, ">$cmbTblFile", or die "Cannot open $cmbTblFile for writing: $OS_ERROR\n";
foreach my $x (keys %cmbTx)
{
  print OUT "$x\t$cmbTx{$x}\n";
}
close OUT;

print "---Combined taxonomy written for $count ASVs to $countTblFile\n";
my @suffixes = (".csv");
my $Prefix = basename($countTblFile, @suffixes);

print "---Adding taxonomy to $countTblFile\n";
my $cntWtx = "$Prefix"."_w_taxa.csv";
open ALL, ">$cntWtx", or die "Cannot open $cntWtx for writing: $OS_ERROR\n";

my $cnttxon = "$Prefix"."_taxa_only.csv";
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
    my $cmd = " /home/jholm/bin/vaginal_combine_tbl_cols.pl -i $cnttxon -o $merged";
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
sub readEzBioTbl{

  my $file = shift;
  my %tbl;
  my $line=0;
  my $count=0;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    $line++;
    if ($line == 1)
    {
      next;
    }
    my ($seq, $d, $p, $c, $o, $f, $g, $s) = split (/[\s+,]/,$_);
    my $t;
    if ($s && $s ne "NA")
    {
      $count++;
      $t="$s";
    }
    elsif ($g ne "NA")
    {
      $t="g_$g";
    }
    elsif ($f ne "NA")
    {
      $t = "f_$f";
    }
    elsif ($o ne "NA")
    {
      $t = "o_$o";
    }
    elsif ($c ne "NA")
    {
      $t = "c_$c";
    }
    elsif ($p ne "NA")
    {
      $t = "p_$p";
    }
    else
    {
      $t = "d_$d";
    }
    $tbl{$seq} = $t;
  }
  close IN;
print "\n$count ASVs assigned species via EzBio\n";
  return %tbl;
}
sub readSILVATbl{

  my $file = shift;
  my %tbl;
  my $line=0;
  my $count=0;
  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    $line++;
    if ($line == 1)
    {
      next;
    }
    my ($seq, $d, $p, $c, $o, $f, $g, $s) = split (/[\s+,]/,$_);
    my $t;
    if ($s && $s ne "NA")
    {
      $count++;
      $t="$g"."_$s";
    }
    elsif ($g ne "NA")
    {
      $t="g_$g";
    }
    elsif ($f ne "NA")
    {
      $t = "f_$f";
    }
    elsif ($o ne "NA")
    {
      $t = "o_$o";
    }
    elsif ($c ne "NA")
    {
      $t = "c_$c";
    }
    elsif ($p ne "NA")
    {
      $t = "p_$p";
    }
    else
    {
      $t = "d_$d";
    }
    $tbl{$seq} = $t;
  }
  close IN;
print "\n$count ASVs assigned species via SILVA\n";
  return %tbl;
}

sub readUNITETbl{

  my $file = shift;
  my %tbl;
  my $line=0;

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    chomp;
    $line++;
    if ($line == 1)
    {
      next;
    }
    my ($seq, $k, $p, $c, $o, $f, $g, $s) = split (/[\s+,]/,$_);
    my $t;
    if ($s ne "NA")
    {
      $s =~ s/s__//g;
      $g =~ s/g__//g;
      $t="$g"."_$s";
    }
    elsif ($g ne "NA")
    {
      $g =~ s/g__//g;
      $t="g_$g";
    }
    elsif ($f ne "NA")
    {
      $f =~ s/f__//g;
      $t="f_$f";
    }
    elsif ($o ne "NA")
    {
      $o =~ s/o__//g;
      $t="o_$o";
    }
    elsif ($c ne "NA")
    {
      $c =~ s/c__//g;
      $t="c_$c";
    }
    elsif ($p ne "NA")
    {
      $p =~ s/p__//g;
      $t="p_$p";
    }
    elsif ($k ne "NA")
    {
      $k =~ s/k__//g;
      $t="k_$k";
    }
    else
    {
      $t="unclassified";
    }
    $tbl{$seq} = $t;
  }
  close IN;

  return %tbl;
}

sub readOtherTbl{

  my $file = shift;
  my %tbl;
  my $line=0;
  my %scores;

  open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
  foreach (<IN>)
  {
    my $final;
    my $t;
    chomp;
    my @line = split (/[\s+;]/,$_);
    my $seq = $line[0];
    shift @line;
    my $count=0;
    foreach my $piece (@line)
    {
      ## evaluate score, if > 90, make value in tbl
      my ($tax, $score) = split (/[\(\)]/,$piece);
      print "$tax has a score of $score\n";
      if ($score > 90)
      {
        $count++;
        $t = $tax;
      }
    }
    if ($count == 1)
    {
      $final= "d_$t";
    }
    if ($count == 2)
    {
      $final= "p_$t";
    }
    if ($count == 3)
    {
      $final= "c_$t";
    }
    if ($count == 4)
    {
      $final= "o_$t";
    }
    if ($count == 5)
    {
      $final= "f_$t";
    }
    if ($count == 6)
    {
      $final= "g_$t";
    }
    if ($count == 7)
    {
      $final= "$t";
    }
    else
    {
      print "No taxonomy for $seq\n";
    }
    $tbl{$seq} = $final;
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