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

my $scriptsDir;
my $pipelineDir;

BEGIN {
    use File::Spec::Functions;
    use File::Basename;

    $scriptsDir = dirname(__FILE__);

}
use strict;
use warnings;

use Pod::Usage;
use English      qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use List::Util   qw( sum );
use Data::Dumper qw(Dumper);
use Cwd;

$OUTPUT_AUTOFLUSH = 1;
####################################################################
##                             OPTIONS
####################################################################
my @fullTaxonomies = ();
GetOptions(
    "PECAN-taxonomy|p=s"      => \my $pecanFile,
    "SILVA-taxonomy|s=s"      => \my $silvaFile,
    "full-taxonomy|f=s"       => \@fullTaxonomies,
    "UNITE-taxonomy|u=s"      => \my $uniteFile,
    "homd-file=s"             => \my $homdFile,
    "ezBioCloud-taxonomy|e=s" => \my $ezBioFile,
    "vaginal"                 => \my $vaginal,
    "other-taxonomy|o=s"      => \my $otherFile,
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

my @taxNames;

if ($silvaFile) {
    print "---Using SILVA taxonomy\n";
    push @taxNames, "SILVA";
}
if ($pecanFile) {
    print "Using PECAN taxonomy\n";
    push @taxNames, "PECAN";
}
if ($ezBioFile) {
    print "---Using only ezBioCloud taxonomy\n";
    push @taxNames, "exBioCloud";
}
if ($homdFile) {
    print "---Using HOMD taxonomy\n";
    push @taxNames, "HOMD";
}
if ($uniteFile) {
    print "---Using UNITE taxonomy\n";
    push @taxNames, "UNITE";
}

####################################################################
##                               MAIN
####################################################################
my %silva;
my %homd;
my %unite;
my %pecan;
my %other;
my %ez;

foreach ( ( $silvaFile, $homdFile, $ezBioFile, $uniteFile, $otherFile ) ) {

}
if ($silvaFile) {
    %silva = readFullTaxonomyTbl($silvaFile);
    my $silvaCond = "silva_condensed.txt";
    print "---Writing condensed SILVA taxonomy to $silvaCond\n";
    open OUT, ">$silvaCond",
        or die "Cannot open $silvaCond for writing: $OS_ERROR\n";
    foreach my $x ( keys %silva ) {
        print OUT "$x\t$silva{$x}\n";
    }
    close OUT;
}

if ($homdFile) {
    %homd = readFullTaxonomyTbl($homdFile);
    my $homdCond = "homd_condensed.txt";
    print "---Writing condensed HOMD taxonomy to $homdCond\n";
    open OUT, ">$homdCond",
        or die "Cannot open $homdCond for writing: $OS_ERROR\n";
    foreach my $x ( keys %homd ) {
        print OUT "$x\t$homd{$x}\n";
    }
    close OUT;
}

if ($ezBioFile) {
    %ez = readezBio($ezBioFile);
    my $ezCond = "ezBio_condensed.txt";
    print "---Writing condensed ezBio taxonomy to $ezCond\n";
    open OUT, ">$ezCond",
        or die "Cannot open $ezCond for writing: $OS_ERROR\n";
    foreach my $x ( keys %ez ) {
        print OUT "$x\t$ez{$x}\n";
    }
    close OUT;
}

if ($uniteFile) {
    %unite = readUNITETbl($uniteFile);
    my $uniteCond = "unite_condensed.txt";
    print "---Writing condensed UNITE taxonomy to $uniteCond\n";
    open OUT, ">$uniteCond",
        or die "Cannot open $uniteCond for writing: $OS_ERROR\n";
    foreach my $x ( keys %unite ) {
        print OUT "$x\t$unite{$x}\n";
    }
    close OUT;
}
if ($pecanFile) {
    %pecan = read2ColTbl($pecanFile);
    print
        "---Using SILVA taxonomy except for Lactobacillus and BV species (PECAN)\n";
}
if ($otherFile) {
    %other = readOtherTbl($otherFile);
    print "---Including $otherFile taxonomy\n";
    my $otherCond = "other_condensed.txt";
    print "---Writing condensed other taxonomy to $otherCond\n";
    open OUT, ">$otherCond",
        or die "Cannot open $otherCond for writing: $OS_ERROR\n";
    if ($otherFile) {
        foreach my $x ( keys %other ) {
            my $tx = "NA";
            if ( $other{$x} ) {
                print OUT "$x\t$other{$x}\n";
            }
            else {
                print OUT "$x\t$tx\n";
            }
        }
        close OUT;
    }
}

my %cmbTx;
my $count = 0;

print "---Combining taxonomic assignments of all source files (includes "
    . scalar( keys %silva )
    . " ASVs)\n";

if ( $silvaFile || $homdFile ) {
    if (%homd) {
        %silva = %homd
            ; # process homd and silva identically. too lazy to rename variables
    }
    foreach my $x ( keys %silva ) {
        $count++;
        my $tx = $silva{$x};

        #print "The Silva taxonomy of $x is $tx\n";
        $cmbTx{$x} = $tx;

        if ($pecanFile) {

          # If the SILVA annotation is one of these genuses, use PECAN
          # NB The exact matches that yield species assignments, done during
          # assign_taxonomies.R, take precedence over PECAN. They are however,
          # excluded by these very strict regex patterns

           # Why can't we simply look for where genera match between SILVA and
           # PECAN assignments? Because Saccharibacteria, Saccharimonadia,
           # Saccharimonadales, Saccharimonadaceae, Candidatus_Saccharimonas,
           # and Mageeibacillus are not contained in our PECAN models

           # can we log when these are reassigned (uncomment that print statmt
           # below) so I can see what Saccharibacteria and Mageeibacillus ASVs
           # are assigned by speciateIT?
            if (   $tx && $tx =~ /^g_Lactobacillus$/
                || $tx && $tx =~ /^g_Shuttleworthia$/
                || $tx && $tx =~ /^p_Saccharibacteria$/
                || $tx && $tx =~ /^c_Saccharimonadia$/
                || $tx && $tx =~ /^o_Saccharimonadales$/
                || $tx && $tx =~ /^f_Saccharimonadaceae$/
                || $tx && $tx =~ /^g_Candidatus_Saccharimonas$/
                || $tx && $tx =~ /^g_Gardnerella$/
                || $tx && $tx =~ /^g_Prevotella$/
                || $tx && $tx =~ /^g_Sneathia$/
                || $tx && $tx =~ /^g_Atopobium$/
                || $tx && $tx =~ /^g_Mageeibacillus$/ )
            {
                #print "The PECAN taxonomy of $x is $pecan{$x}\n";
                $tx = $pecan{$x};
                $cmbTx{$x} = $tx;
            }
        }
        if ($otherFile) {
            $tx = $other{$x};
            if ($tx) {
                $cmbTx{$x} = $tx;
            }
        }
        if ($ezBioFile) {
            $tx = $ez{$x};
            if (   $tx && $tx =~ /Lactobacillus/
                || $tx && $tx =~ /Saccharibacteria/
                || $tx && $tx =~ /Gardnerella/
                || $tx && $tx =~ /Prevotella/
                || $tx && $tx =~ /Sneathia/
                || $tx && $tx =~ /Atopobium/ )
            {
                #print "The PECAN taxonomy of $x is $pecan{$x}\n";
                $tx = $ez{$x};
                $cmbTx{$x} = $tx;
            }
        }

    }

}
elsif ($ezBioFile) {
    foreach my $x ( keys %ez ) {
        $count++;
        my $tx = $ez{$x};

        $cmbTx{$x} = $tx;
    }
}
elsif ($uniteFile) {
    foreach my $x ( keys %unite ) {
        $count++;
        my $tx = $unite{$x};

        $cmbTx{$x} = $tx;
    }
}

my $cmbTblFile = "cmb_tx.txt";
open OUT, ">$cmbTblFile",
    or die "Cannot open $cmbTblFile for writing: $OS_ERROR\n";
foreach my $x ( keys %cmbTx ) {
    print OUT "$x\t$cmbTx{$x}\n";
}
close OUT;

print "---Combined taxonomy written for $count ASVs to $cmbTblFile\n";
my @suffixes = (".csv");
my $Prefix   = basename( $countTblFile, @suffixes );

print "---Adding taxonomy to $countTblFile\n";
my $taxString = join( "-", @taxNames );
my $cntWtx    = "${Prefix}.${taxString}.asvs+taxa.csv";
open ALL, ">$cntWtx", or die "Cannot open $cntWtx for writing: $OS_ERROR\n";

my $cnttxon = "${Prefix}.${taxString}.taxa.csv";
open TXON, ">$cnttxon",
    or die "Cannot open $cnttxon for writing: $OS_ERROR\n";

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
            }
            else {
                push( @asvs, "$asv" );
                push( @taxa, "$cmbTx{$asv}" );
            }
        }
    }
}
print ALL "," . join( ",", @asvs ) . "\n";
my $taxaHead = "," . join( ",", @taxa ) . "\n";
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

@suffixes = (".csv");
$Prefix   = basename( $cnttxon, @suffixes );
my $merged = $Prefix . "-merged.csv";

my $debugParam = $debug   ? "--debug"   : "";
my $vParam     = $vaginal ? "--vaginal" : "";

my $cmd
    = "$scriptsDir/combine_tbl_cols.pl -i $cnttxon -o $merged $debugParam $vParam";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

####################################################################
##                               SUBS
####################################################################
sub read2ColTbl {

    my $file = shift;

    if ( !-f $file ) {
        warn "\n\n\tERROR in read2ColTbl(): $file does not exist";
        print "\n\n";
        exit 1;
    }

    my %tbl;
    open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
    foreach (<IN>) {
        chomp;
        my ( $id, $t ) = split( /\s+/, $_ );
        $tbl{$id} = $t;
    }
    close IN;

    return %tbl;
}

sub readFullTaxonomyTbl {

    my $file = shift;
    my %tbl;
    my $line = 0;

    open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
    my @levelHints  = ();
    my @speciesCols = ();
    foreach (<IN>) {
        chomp;
        $line++;
        if ( $line == 1 ) {

            @levelHints = split( /,/, $_ );   # contains full taxonomic levels
            @levelHints = @levelHints[ 1 .. $#levelHints ];
            @speciesCols
                = grep { $levelHints[$_] =~ /^species$/i } 0 .. $#levelHints;

            # Create taxon level hints that consist of the first leter of the
            # level, lowercase, plus an underscore
            foreach (@levelHints) {
                $_ = lc( substr $_, 0, 1 );
                $_ = $_ . "_";
            }
            next;
        }

        my ( $asv, @taxa ) = split( /,\s*/, $_ );
        my $t;
        if ( @speciesCols && $taxa[ $speciesCols[0] ] ne "NA" ) {
            my $s = $taxa[ $speciesCols[0] ];
            my $g = $taxa[ $speciesCols[0] - 1 ];
            $t = $g . "_$s";
        }
        elsif (@taxa) {
            my $taxon = "NA";
            while ( $taxon eq "NA" && scalar @taxa > 0 ) {
                $taxon = pop @taxa;
            }
            if ( !@taxa ) {
                $t = "unclassified";
            }
            $t = $levelHints[ scalar @taxa ] . $taxon;
        }
        else {
            $t = "unclassified";
        }
        $tbl{$asv} = $t;
    }
    close IN;

    return %tbl;
}

sub readezBio {

    my $file = shift;
    my %tbl;
    my $line = 0;

    open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
    foreach (<IN>) {
        chomp;
        $line++;
        if ( $line == 1 ) {
            next;
        }
        my ( $seq, $d, $p, $c, $o, $f, $g, $s ) = split( /[\s+,]/, $_ );
        my $t;
        if ( $s && $s ne "NA" ) {
            $t = "$s";
        }
        elsif ( $g ne "NA" ) {
            $t = "g_$g";
        }
        elsif ( $f ne "NA" ) {
            $t = "f_$f";
        }
        elsif ( $o ne "NA" ) {
            $t = "o_$o";
        }
        elsif ( $c ne "NA" ) {
            $t = "c_$c";
        }
        elsif ( $p ne "NA" ) {
            $t = "p_$p";
        }
        else {
            $t = "d_$d";
        }
        $tbl{$seq} = $t;
    }
    close IN;

    return %tbl;
}

sub readUNITETbl {

    my $file = shift;
    my %tbl;
    my $line = 0;

    open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
    foreach (<IN>) {
        chomp;
        $line++;
        if ( $line == 1 ) {
            next;
        }
        my ( $seq, $k, $p, $c, $o, $f, $g, $s ) = split( /[\s+,]/, $_ );
        my $t;
        if ( $s ne "NA" ) {
            $s =~ s/s__//g;
            $g =~ s/g__//g;
            $t = "$g" . "_$s";
        }
        elsif ( $g ne "NA" ) {
            $g =~ s/g__//g;
            $t = "g_$g";
        }
        elsif ( $f ne "NA" ) {
            $f =~ s/f__//g;
            $t = "f_$f";
        }
        elsif ( $o ne "NA" ) {
            $o =~ s/o__//g;
            $t = "o_$o";
        }
        elsif ( $c ne "NA" ) {
            $c =~ s/c__//g;
            $t = "c_$c";
        }
        elsif ( $p ne "NA" ) {
            $p =~ s/p__//g;
            $t = "p_$p";
        }
        elsif ( $k ne "NA" ) {
            $k =~ s/k__//g;
            $t = "k_$k";
        }
        else {
            $t = "unclassified";
        }
        $tbl{$seq} = $t;
    }
    close IN;

    return %tbl;
}

sub readOtherTbl {

    my $file = shift;
    my %tbl;
    my $line = 0;
    my %scores;

    open IN, "$file" or die "Cannot open $file for reading: $OS_ERROR\n";
    foreach (<IN>) {
        my $final;
        my $t;
        chomp;
        my @line = split( /[\s+;]/, $_ );
        my $seq  = $line[0];
        shift @line;
        my $count = 0;
        foreach my $piece (@line) {
            ## evaluate score, if > 90, make value in tbl
            my ( $tax, $score ) = split( /[\(\)]/, $piece );
            print "$tax has a score of $score\n";
            if ( $score > 90 ) {
                $count++;
                $t = $tax;
            }
        }
        if ( $count == 1 ) {
            $final = "d_$t";
        }
        if ( $count == 2 ) {
            $final = "p_$t";
        }
        if ( $count == 3 ) {
            $final = "c_$t";
        }
        if ( $count == 4 ) {
            $final = "o_$t";
        }
        if ( $count == 5 ) {
            $final = "f_$t";
        }
        if ( $count == 6 ) {
            $final = "g_$t";
        }
        if ( $count == 7 ) {
            $final = "$t";
        }
        else {
            print "No taxonomy for $seq\n";
        }
        $tbl{$seq} = $final;
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
