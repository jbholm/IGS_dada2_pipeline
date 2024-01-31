#!/usr/local/packages/perl-5.30.2/bin/perl

=head1 NAME

  part2_illumina_dada2_devel.pl  

=head1 DESCRIPTION
  
  This is the second part of the MD Genomics 16S analysis pipeline. This script
  will create a project from one or more runs, remove bimeras, produce sequence
  (ASV) count tables, assign taxonomy to each ASV, create count tables using all
  taxonomic annotations, and create an HTML report.
  
  Most parameters are stored in config.json.
  
  Recommended to submit to SGE grid with request for 64G memory.

=head1 SYNOPSIS

  part2_illumina_dada2.pl -r <run> [-r <run2> … ] -p 16stest -m <map> [--pacbio] 
    [-v {V3V4, V4, ITS, FULL-LENGTH, OMPA}] [--oral]
    [--tax {SILVA, PECAN, PECAN-SILVA, UNITE, HOMD, SILVA138forPB}]

=head1 OPTIONS

=over

=item B<--run, -r> RUN
  Input run names (name of the run's directory)

=item B<--run-storage> PATH
  Full path to the parent directory where runs will be found. The default is in 
  config.json.

=item B<--project-ID, -p>
  Provide the project ID. A timestamp in the format "YYMMDD_" will be prepended.

=item B<--map, -m>
  A tab-delimited file listing the samples to add to this project. The file must
  have two columns named RUN.PLATEPOSITION and sampleID. RUN.PLATEPOSITION
  must be identical to the sample names in the "corrected" map created by Part 1.
  Default: project_map.txt.

=item B<--pacbio> 
  Use this flag for PacBio runs. Sets B<--variable-region> to "FULL-LENGTH" and 
  uses PacBio-specific chimera removal parameters.

=item B<--variable-region>|B<-v> {V3V4, V4, ITS, FULL-LENGTH}
  The targeted variable region, V3V4 by default. The variable region determines 
  the base taxonomic annotation scheme. --tax can be given to perform other 
  annotation schemes. See the mapping below:
  
  Variable region       Annotation scheme
  V3V4                  PECAN-SILVA
  V4                    SILVA
  ITS                   UNITE
  FULL-LENGTH           SILVA138forPB

=item B<--tax, -t> {PECAN-SILVA, SILVA, SILVA138forPB, PECAN, HOMD, UNITE}
  Perform additional taxonomic annotation schemes. For all references 
  besides PECAN, annotation is done by DADA2's implementation of the RDP 
  classifier.
  
  PECAN-SILVA - Initially assigns ASVs by the default SILVA database. Any ASVs assigned to the 
    genuses Lactobacillus, Gardnerella, Prevotella, Sneathia, Atopobium, 
    Shuttleworthia, or Saccharibacteria are classified by PECAN instead. 
    PECAN-only abundance tables are also produced.
  PECAN - PECAN is the only classifier run.
  SILVA - Assignments from the default SILVA database down to genus level.
  SILVA138forPB - Assignments from the SILVA138 database down to species level. 
    Use for PacBio (full-length) reads only.
  HOMD - Assigns from the Human Oral Microbiome Database, usually to species 
    level.
  UNITE - Assigns from the UNITE database, usually to species level. Use only 
    for the fungal ITS region.

=item B<--notVaginal>
  Use when project is NOT just vaginal sequences to prevent renaming of several
  SILVA-derived taxonomic assignments, including two that are renamed to 
  BVAB_TM7.

=item B<--oral>
  Changes the base taxonomic annotation scheme to "HOMD", regardless of
  --variable-region, and silently adds B<--notVaginal>. This can be used in 
  combination with B<--pacbio>, or other B<--tax> parameters. To get both 
  SILVA138forPB and HOMD annotations for a PacBio run, use the following 
  parameters in any order:

  B<--pacbio --oral --tax=SILVA138forPB>

=item B<--nocsts>
  Flag to skip CST assignment. By default, CSTs are assigned if the V3V4 variable region is being analyzed and samples are non-oral. Assignment data are written to *_PECAN_taxa-merged_StR_CST.csv.


=item B<--noreport>
  Flag to skip creating HTML report.

=item B<--verbose>
  Print extra messages to the log file.

=item B<--dry-run>
  Unsupported, probably doesn't work.

=item B<-h|--help>
  Print help message and exit successfully.

=back

=head1 EXAMPLE

  cd /local/scratch/MM
  part2_illumina_dada2.pl -i MM_01,MM_03,MM_05,MM_21 -v V3V4 -p MM
  
=cut

use strict;
use warnings;
my $scriptsDir;
my $pipelineDir;

BEGIN {
    use File::Spec::Functions;
    use File::Basename;

    $pipelineDir = dirname(__FILE__);
    $scriptsDir  = catdir( $pipelineDir, "scripts" );

}
use lib $scriptsDir;    # .pm files in ./scripts/ can be loaded
use lib "/usr/local/packages/perl-5.30.2/lib/site_perl/5.30.2/";

require Stats_gen;
require Version;

use Pod::Usage;
use English      qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
require Cwd;
use Cwd        qw(abs_path getcwd);
use File::Temp qw/ tempfile /;
require POSIX;
require IO::Tee;
use Data::Dumper;

#use Email::MIME;
#use Email::Sender::Simple qw(sendmail);
use File::Spec;
require IO::Tee;
use JSON;
use Time::Piece;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################
my $config_hashref = config();

my $region      = "V3V4";
my $csts        = 1;
my $pacbio      = 0;
my $report      = 1;
my $map_file    = "";
my $run_path    = $config_hashref->{"run_storage_path"};
my $multithread = 1;

# this is the way it is only to preserve the interface of --notVaginal. In the future, please change to --no-vaginal
my $vaginal    = 1;
my $notVaginal = 0;
my @schemes    = ();
GetOptions(
    "run|r=s"                                    => \my @runs,
    "run-storage=s"                              => \$run_path,
    "variable-region|v=s"                        => \$region,
    "project-ID|p=s"                             => \my $project,
    "map|m=s"                                    => \$map_file,
    "help|h!"                                    => \my $help,
    "verbose!"                                   => \my $verbose,
    "dry-run"                                    => \my $dryRun,
    "notVaginal|notvaginal|novaginal|no-vaginal" => \$notVaginal,
    "tax|t=s"                                    => \@schemes,
    "oral"                                       => \my $oral,
    "csts!"                                      => \$csts,
    "pacbio!"                                    => \$pacbio,
    "report!"                                    => \$report,
    "multithread=i"                              => \$multithread
    )

    or pod2usage( verbose => 0, exitstatus => 1 );

### REASONS FOR STOPPING IMMEDIATELY
if ($help) {
    pod2usage( verbose => 2, exitstatus => 0 );
    exit 1;
}
if ( !$project ) {
    die "Please provide a project name\n\n";
}

# make the year directory if needed, and add yymmdd to project name
my $t       = Time::Piece->new();
my $year    = $t->year;
my $yeardir = catdir( $config_hashref->{"project_storage_path"}, $year );
if ( !-d $yeardir ) {
    system("install -g igs -m 2770 -d $yeardir")
        or print "Didn't create $yeardir directory, $!\n";
}
$project =~ s/_/-/g;
$project = $t->strftime("\%y\%m\%d") . "_$project";
my $projDir = catdir( $yeardir, $project );
system( "install -g igs -m 2770 -d " . "\"$projDir\"" )
    or print "Didn't create $projDir directory, $!";

$run_path = abs_path($run_path) or die "Can't find run storage path.\n$!";

if ( !@runs ) {
    die "Please provide (a) run ID(s)\n\n";
}

##Compute the absolute path to each run
foreach (@runs) {
    my $pathsep = File::Spec->catfile( '', '' );
    if ( $_ =~ qr/$pathsep/ ) {
        die
            "$_ contains a path separator, and is not a valid run ID. Run IDs are the names"
            . " of directories found by default in:\n$run_path\n";
    }
    $_ = catdir( $run_path, $_ );
}

# Refine and validate the run paths (no duplicates, must exist on filesystem)
my %seen;
foreach ( @runs, $map_file ) {
    if ($_) {    # If the variable is a non-empty string...
        my $copy    = $_;
        my $relPath = $copy;
        $copy = abs_path($relPath);    # get abs path
        if ( !defined $copy || !-e $copy ) {
            die $_
                . " not found.\n"
                . "Current working directory: "
                . getcwd() . "\n";
        }
        $copy =~ s/\/$//;              # remove any trailing slash
        if ( $seen{ basename($copy) }++ ) {
            die "$_ is a duplicate run name and not allowed\n";
        }
        $_ = $copy
            ;    # external variable referenced by path has now been edited.
    }
}
chdir $projDir;

# PRINT TO LOG ASAP
my $log = "$project" . "_part2_16S_pipeline_log.txt";
open my $logFH, ">>$log" or die "Cannot open $log for writing: $OS_ERROR";
my $old = select $logFH;
$OUTPUT_AUTOFLUSH = 1;
select $old;

print $logFH "This file logs the progress of "
    . scalar(@runs)
    . " runs for $project 16S amplicon sequences through the illumina_dada2.pl pipeline.\n";

my $logTee = new IO::Tee( \*STDOUT, $logFH );
print $logTee "PIPELINE VERSION: " . Version::version() . "\n";
my $time = POSIX::strftime( "%Y-%m-%d %H:%M:%S", localtime(time) );
print $logTee "$time\n\n";

local $SIG{__WARN__} = sub {
    print $logTee "WARNING: $_[0]";
    print "WARNING: $_[0]";
    print STDERR "WARNING: $_[0]";
};    # Warnings go to log file, stdout, and stderr

$ENV{'LD_LIBRARY_PATH'} = "/usr/lib64/:/usr/local/packages/gcc/lib64";

my $python2 = $config_hashref->{"python2"};

$ENV{'PYTHONPATH'} = "";

if ($map_file) {
    my $orig_map_file = $map_file;

# rather roundabout...all I want is to rename this to <project>_sample_plate_map.txt
# Nextflow when
    if ( $map_file !~ /^${projDir}/ ) {
        $map_file = copy_to_project( $project, $map_file );
    }
    else {
        $map_file = move_to_project( $project, $map_file );
    }
    rename $map_file, catfile( dirname($map_file), "sample_plate_map.txt" );
    $map_file = move_to_project( $project,
        catfile( dirname($map_file), "sample_plate_map.txt" ) );
    $logTee->print("MAP: $orig_map_file => $map_file\n");
}
else {
    $logTee->print(
        "MAP: NONE (all samples from specified runs will be processed)\n");
}

$logTee->print("\n");
$logTee->print("RUN(s):\n");
map { $logTee->print("$_\n") } @runs;
$logTee->print("\n");

if ($oral) {
    $notVaginal = 1;
}
if ($notVaginal) {
    $logTee->print("VAGINAL: NO\n");
    $vaginal = 0;
}
else {
    $logTee->print("VAGINAL: YES\n");
}
if ( !$region ) {
    die
        "Please provide a variable region (-v), V3V4, V4, ITS, or FULL-LENGTH\n";
}
else {
    $logTee->print("REGION: $region\n");
}
$region = uc $region;
if (List::Util::none { $_ eq $region }
    ( 'V3V4', 'V4', 'ITS', 'FULL-LENGTH', 'OMPA' )
    )
{
    die
        "Illegal --variable-region (-v). V3V4, V4, ITS, ompA, and FULL-LENGTH are accepted.";
}

my %taxonomy_flags;

my $base_scheme;
my %region_to_scheme = (
    V3V4          => "PECAN-SILVA",
    V4            => "SILVA",
    ITS           => "UNITE",
    "FULL-LENGTH" => "SILVA138forPB",
    OMPA          => "(none)"
);

my $vaginal_adj = "vaginal";
if ( !$vaginal ) {
    $vaginal_adj = "non-vaginal";
    $region_to_scheme{"V3V4"} = "SILVA";
}

if ($oral) {
    $base_scheme = "HOMD";
    $logTee->print(
        "Using $base_scheme as the base annotation scheme because the samples were oral.\n"
    );
}
else {
    $base_scheme = $region_to_scheme{$region};
    my $message
        = "Using $base_scheme as the base annotation scheme because $region was targeted";
    if ( $region eq "V3V4" ) {
        $message .= " and the samples were $vaginal_adj";
    }
    $logTee->print("$message.\n");
}
unshift( @schemes, $base_scheme );

my $ps   = grep ( /^PECAN-SILVA$/, @schemes );
my $s    = grep ( /^SILVA$/,       @schemes );
my $s138 = grep( /^SILVA138forPB$/, @schemes );
my $p    = grep ( /^PECAN$/, @schemes );
my $h    = grep( /^HOMD$/,     @schemes );
my $u    = grep( /^UNITE$/,    @schemes );
my $na   = grep( /^\(none\)$/, @schemes );
if ( $ps + $s + $s138 + $p + $h + $u + $na == scalar @schemes ) { }
else {
    die
        "Illegal taxonomic annotation scheme given by --tax. Legal values for --tax|-t are "
        . "PECAN-SILVA, SILVA, SILVA138forPB, PECAN, HOMD, and UNITE.";
}

$logTee->print("Annotation schemes: ");
$logTee->print( join( " ", @schemes ) . "\n" );

# Which taxonomies need to be assigned from? (this is necessary bc PECAN-SILVA
# specifies two taxonomies)
%taxonomy_flags
    = ( SILVA138forPB => 0, SILVA => 0, HOMD => 0, UNITE => 0, PECAN => 0 );
if ($ps) {
    if ($csts) {
        print $logTee "Using Valencia to assign samples to CSTs\n";
    }
    else {
        print $logTee "Skipping CST assignment.\n";
    }
    $taxonomy_flags{PECAN} += $ps;
    $taxonomy_flags{SILVA} += $ps;
}
$taxonomy_flags{SILVA} += $s;
if ($p) {
    if ($csts) {
        print $logTee "Using Valencia to assign samples to CSTs\n";
    }
    else {
        print $logTee "Skipping CST assignment.\n";
    }
    $taxonomy_flags{PECAN} += $p;
}
$taxonomy_flags{SILVA138forPB} += $s138;
$taxonomy_flags{HOMD}          += $h;
$taxonomy_flags{UNITE}         += $u;

print $logTee "\n";
####################################################################
##                               MAIN
####################################################################
my @classifs = ();
my $projabund;

# the un-annotated abundance table signals dada2 already completed
# my $cmd = "rm -f *-dada2_abundance_table.rds";
# execute_and_log( $cmd, *STDOUT, $dryRun, "" );
my $cmd;

my $vopt = "";
if ($vaginal) {
    $vopt = "--vaginal";
}

##loop over array to copy the file to the main current working directory
## using the array string to also add a name
print $logTee "---Combining "
    . scalar(@runs)
    . " abundance table(s) and removing bimeras.\n";

my ( $abundRds, $abund, $stats, $fasta ) = (
    "all_runs_dada2_abundance_table.rds",
    "all_runs_dada2_abundance_table.csv",
    "DADA2_stats.txt",
    "all_runs_dada2_ASV.fasta"
);
( $abundRds, $abund, $stats, $fasta )
    = map { move_to_project( $project, $_, 1 ) }
    ( $abundRds, $abund, $stats, $fasta );

if ( List::Util::any { !-e $_ } ( $abundRds, $abund, $fasta, $stats ) ) {
    print $logTee "Combining runs and removing bimeras.\n";
    ( $abundRds, $abund, $fasta, $stats )
        = dada2_combine( $pacbio, \@runs, $map_file, $multithread );

    ( $abundRds, $abund, $fasta, $stats )
        = map { move_to_project( $project, $$_ ) }
        ( \$abundRds, \$abund, \$fasta, \$stats );

    $logTee->print("Outputs:\n");
    foreach ( ( $abundRds, $abund, $fasta, $stats ) ) {
        print $logTee "$_\n";
    }
    $logTee->print("\n");
}
else {
    print $logTee "Chimeras already removed. Skipping...\n";
}

# log original run paths
my $run_info_hash = combine_run_metadata(@runs);

# Give ASV's unique and easy-to-look-up IDs
$cmd = join(
    " ",
    (   "$python2",   "\"$scriptsDir/rename_asvs.py\"",
        "\"$fasta\"", "-p",
        "\"$project\""
    )
);
execute_and_log( $cmd, $logTee, $dryRun,
    "\n---Renaming ASVs in FASTA, abundance tables, and classification key(s)."
);

my $species_file;
my @full_classif_csvs;
my @two_col_classifs;

my @dada2_taxonomy_list;
foreach ( keys %taxonomy_flags ) {
    if ( $taxonomy_flags{$_} > 0 ) {
        if ( $_ eq "SILVA" ) {
            push @dada2_taxonomy_list, $config_hashref->{"SILVA_default"};
        }
        elsif ( $_ ~~ [ keys %{ $config_hashref->{"taxonomy_dbs"} } ] ) {
            push @dada2_taxonomy_list, $_;
        }
    }
}
if ( scalar @dada2_taxonomy_list > 0 ) {
    my @dada2Output
        = map { "$project" . "_$_.classification.csv" } @dada2_taxonomy_list;
    if ( List::Util::any { !-e $_ } @dada2Output ) {
        print $logTee
            "Using DADA2's RDP Classifier implementation to classify amplicon sequence variants (ASVs) with taxonomies:\n";
        print $logTee join( ",", @dada2_taxonomy_list ) . "\n";

        @dada2Output
            = dada2_classify( $fasta, \@dada2_taxonomy_list, $multithread )
            ;    # assignment
        ### @dada2Output now has every file that needs to be renamed
        print $logTee "---Renaming dada2 output files for project...\n";
        @dada2Output = map { move_to_project( $project, $_ ) } @dada2Output;
    }
    else {
        print $logTee
            "DADA2 has already been used to classify ASVs with taxonomies:\n";
        print $logTee join( ",", @dada2_taxonomy_list ) . "\n";
        print $logTee "Skipping.\n";
    }
    push( @full_classif_csvs, @dada2Output );

}

my $pecanFile = "MC_order7_results.clean.txt";
if ( $taxonomy_flags{PECAN} ) {
    if ( !-e $pecanFile ) {
        my $cmd = join(
            " ",
            (   "\"$config_hashref->{'python3.9'}\"",    "-s",
                "\"$scriptsDir/pecan_wrapper.py\"",      "-d",
                "\"$config_hashref->{'pecan_models'}\"", "-i",
                "\"$fasta\"",                            "-o",
                ".",                                     $vopt
            )
        );
        execute_and_log( $cmd, $logTee, 0,
                  "---Classifying ASVs with $region PECAN models (located in "
                . $config_hashref->{"pecan_models"}
                . ")" );

    }
    else {
        print $logTee
            "ASVs have already been classified with V3V4 PECAN models.\n";
        print $logTee "Located at:\n";
        print $logTee "\t"
            . File::Spec->catfile( $projDir, $pecanFile ) . "\n";
    }
    push @two_col_classifs, $pecanFile;
}

#################################################################
#########APPLY CLASSIFICATIONS TO COUNT TABLES ###################
#################################################################

foreach (@full_classif_csvs) {

# Apply classifications from the various classification files to the count
# tables. This is an error-prone task because we have different classification
# formats from the previous step. We reply on the filenames to tell us what the
# file format is, then take the appropriate action.
    my $db;
    if ( $_ =~ /SILVA138forPB/ ) {
        $db  = "SILVA138forPB";
        $cmd = join( " ",
            "\"$scriptsDir/combine_tx_for_ASV.pl\"",
            "-s", "\"$_\"", "-c", "\"$abund\"", $vopt );
    }
    elsif ( $_ =~ /(SILVA[^\.]*)\.classification/ ) {
        $db  = $1;
        $cmd = join(
            " ",
            (   "\"$scriptsDir/combine_tx_for_ASV.pl\"",
                "-s", "\"$_\"", "-c", "\"$abund\"", $vopt
            )
        );
    }
    elsif ( $_ =~ /HOMD/ ) {
        $db  = "HOMD";
        $cmd = join(
            " ",
            (   "\"$scriptsDir/combine_tx_for_ASV.pl\"",
                "--homd-file", "\"$_\"", "-c", "\"$abund\"", $vopt
            )
        );
    }
    elsif ( $_ =~ /UNITE/ ) {
        $db  = "UNITE";
        $cmd = join(
            " ",
            (   "\"$scriptsDir/combine_tx_for_ASV.pl\"",
                "-u", "\"$_\"", "-c", "\"$abund\"", $vopt
            )
        );
    }
    my $outfile = join( "_",
        ( $project, basename( $abund, ".csv" ), $db, "taxa-merged.csv" ) );
    if ( !-e $outfile ) {
        my $msg = "\n---Labeling count table with $db taxa.";
        execute_and_log( $cmd, $logTee, $dryRun, $msg );
    }

}
foreach (@two_col_classifs) {
    my $db;
    if ( $_ =~ /MC_order7_results/ ) {
        $db = "PECAN";
    }
    elsif ( $_ =~ /SPINGO/ ) {
        $db = "SPINGO";
    }
    my $outfile = join( "_",
        ( $project, basename( $abund, ".csv" ), $db, "taxa-merged.csv" ) );
    if ( !-e $outfile ) {
        my $cmd = join(
            " ",
            (   "\"$scriptsDir/PECAN_tx_for_ASV.pl\"",
                "-p", "\"$_\"", "-c", "\"$abund\"", "-t", "\"$db\"", $vopt
            )
        );

        my $msg = "---Labeling count table with $db taxa.";
        execute_and_log( $cmd, $logTee, $dryRun, $msg );
    }

}
if ( grep ( /^PECAN-SILVA$/, @schemes ) ) {
    ### APPLY PECAN+SILVA CLASSIFICATIONS TO COUNT TABLE (V3V4) ####
    #################################################################
    my $db      = "PECAN-SILVA";
    my $outfile = join( "_",
        ( $project, basename( $abund, ".csv" ), $db, "taxa-merged.csv" ) );
    my $msg       = "---Labeling count table with $db taxa.";
    my $silvaFile = "";
    foreach (@full_classif_csvs) {
        if ( $_ =~ /SILVA[^\.]*\.classification/ ) {
            $silvaFile = $_;
        }
    }
    if ( !$silvaFile ) {
        die "Can't find SILVA classification table from among:\n"
            . join( "\n", @full_classif_csvs ) . "\n";
    }
    $cmd = join(
        " ",
        (   "\"$scriptsDir/combine_tx_for_ASV.pl\"", "-p",
            "\"$pecanFile\"",                        "-s",
            "\"$silvaFile\"",                        "-c",
            "\"$abund\"",                            $vopt
        )
    );
    execute_and_log( $cmd, $logTee, $dryRun, $msg );
}
rename_temps("");

if ( $csts && $taxonomy_flags{PECAN} ) {

# $ENV{'LD_LIBRARY_PATH'} =
#   $ENV{'LD_LIBRARY_PATH'} . ':/usr/local/packages/python-3.5/lib'; #  This path doesn't exist; guess we don't need it

    my $pecanCountTbl
        = basename( $abund, (".csv") ) . ".PECAN.taxa-merged.csv";

    if ( !-e basename( $pecanCountTbl, (".csv") ) . "_StR_CST.csv" ) {
        print $logTee "---Assigning CSTs with Valencia\n";
        my $cmd
            = $config_hashref->{"valencia_python"}
            . " -s $scriptsDir/valencia_wrapper.py "
            . catdir( $pipelineDir,
            $config_hashref->{"valencia_coefficients"} )
            . " $pecanCountTbl";
        print $logTee "\tcmd=$cmd\n" if $dryRun || $verbose;

        system($cmd) == 0
            or die "system($cmd) failed with exit code: $?"
            if !$dryRun;
    }
    else {
        print $logTee "---Count table with CSTs already exists.\n";
    }
}

my $final_merge = glob("*.taxa-merged.csv");
$final_merge = defined $final_merge ? $final_merge : "(none)";
unlink glob "*.taxa.csv";
my $final_ASV_taxa = glob("*.asvs+taxa.csv");
$final_ASV_taxa = defined $final_ASV_taxa ? $final_ASV_taxa : "(none)";

$cmd = join(
    " ",
    (   "\"$pipelineDir/report/report16s.sh\"",
        "\"$projDir\"", "--map", "\"$map_file\"", "--runs", "@runs",
        "--project",    "\"$project\""
    )
);
execute_and_log( $cmd, $logTee, $dryRun, "Creating report..." );
print $logTee "---Final files succesfully produced!\n";
print $logTee
    "Final merged read count table: $final_merge\nFinal ASV table with taxa: $final_ASV_taxa\nFinal ASV count table: $abund\nASV sequences: all_runs_dada2_ASV.fasta\n"
    ;    #Read survival stats: $finalStats\n";
$logTee->close;

####################################################################
##                               SUBS
####################################################################
sub read_json {
    my $file = shift;
    my $mode = shift;

    my $json;
    {
        local $/;    #Enable 'slurp' mode
        if ( -e $file ) {
            open( my $FH, $mode, $file );
            seek $FH, 0, 0 or die;
            $json = <$FH>;
            close $FH;
        }
        else {
            return undef;
        }

    }

    my $hash = {};
    eval { $hash = decode_json($json); };

    return $hash;
}

sub config {
    my %arg             = @_;
    my $new_params      = delete $arg{new_params} // {};
    my $orig_param_file = "$pipelineDir/config.json";

    my $params = read_json( $orig_param_file, "<" );

    return $params;
}

sub combine_run_metadata {
    my $all_run_info = {};
    foreach my $rundir (@_) {
        my $run = basename($rundir);
        $all_run_info->{$run}
            = read_json( catfile( $rundir, ".meta.json" ), "<" );
        $all_run_info->{$run}->{"path"} = $rundir;
        if ( exists $all_run_info->{$run}->{"samples"}
            && ref $all_run_info->{$run}->{"samples"} eq ref {} )
        {
            foreach
                my $sample ( keys %{ $all_run_info->{$run}->{"samples"} } )
            {
                if (ref $all_run_info->{$run}->{"samples"}->{$sample} eq
                    ref {} )
                {
                    foreach my $filepath (
                        values
                        %{ $all_run_info->{$run}->{"samples"}->{$sample} } )
                    {
# an absent file will be decoded from JSON to either an empty hashref or undef
                        if ( !ref($filepath) && defined $filepath ) {
                            $filepath = catfile( $rundir, $filepath );
                        }
                    }
                }
                else {
                    warn "Malformed metadata: $run: samples: $sample.\n";
                }
            }
        }
        else {
            warn "Malformed metadata: $run: samples.\n";
        }
    }

    my $metadata = { runs => $all_run_info };
    open my $metadataFH, ">.meta.json";
    print $metadataFH encode_json($metadata);
    close $metadataFH;

    # symlink runs into project directory for part 3
    foreach my $rundir (@_) {

        # symlink if the run directory is not already a subdirectory of the
        # project directory (runs are external, by default)
        if ( $rundir !~ /^${projDir}/ ) {
            my $symlink = catdir( $projDir, basename($rundir) );
            if ( -e $symlink && -l $symlink ) {
                unlink $symlink;
                $logTee->print( "Relinking "
                        . basename($rundir)
                        . " into the project directory.\n" );
            }
            else {
                $logTee->print( "Linking "
                        . basename($rundir)
                        . " into the project directory.\n" );
            }
            symlink( $rundir, catdir( $projDir, basename($rundir) ) )
                or die "Problem creating "
                . catdir( $projDir, basename($rundir) )
                . " -> $rundir\n$!\n";
        }
    }
    return $metadata;
}

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
        my ( $id, $t ) = split( /\s+/, $_ );
        $tbl{$id} = $t;
    }
    close IN;

    return %tbl;
}

sub R {
    my $script = shift;
    my $args   = shift;

    my $R   = $config_hashref->{'R'};
    my $cmd = "${R}script --verbose $script $args";

    $logTee->print("\$ $cmd\n");

 # Get stdout and stderr; print stderr to the logfile/stdout because R sends
 # normal status messages through STDERR. Save STDOUT as file because all my R
 # scripts give the paths to their output files as STDOUT.
    if ($verbose) {
        $cmd .= " --verbose";
    }
    my $stdout_file = basename($script)
        . ".stdout"
        ;    # I tell R to print the names of its output files into STDOUT
    my $stderr_file = basename($script)
        . ".stderr"
        ;   # R's normal output, including error messages, comes out in STDERR
    $cmd .= " 2>$stderr_file 1>$stdout_file";

    system($cmd) == 0
        or die
        "system($cmd) failed with exit code $?. See this script's STDERR file for error description.";
    my $out = do {
        local $/ = undef;
        open my $fh, "<", $stdout_file
            or die "Could not open STDOUT captured from "
            . basename($script)
            . " to $stdout_file: $OS_ERROR\n";
        <$fh>;
    };
    unlink $stdout_file;

# the above solution would work, except that on our SGE, a slew of nonsense
# gets printed to STDERR whenever a new shell is started. So instead of
# redirecting R's STDERR to STDOUT, we give R our log file to print to.
# One disadvantage is, it doesn't go to perl's STDOUT.
# my $out = `$cmd`;
# if ($?)
# {
#     die
#       "system($cmd) failed with exit code $?. See the part2 LOG file for error description.";
# }

    print("\n");
    return ($out);
}

# Combines all part1 dada2 counts and stats into single files
# $_[0] Evaluated in a boolean context, gives whether PacBio full-length 16S
#   data is being used
# Returns the names of the three output file: abundance table (as RDS), abundance table (as
#   CSV), and stats
sub dada2_combine {
    my $pacbio      = shift;
    my $rundirs     = shift;
    my $map_file    = shift;
    my $multithread = shift;
    my @rundirs     = @$rundirs;

    my $sequencer = $pacbio ? "PACBIO" : "ILLUMINA";

    my $script = catfile( $pipelineDir, "scripts", "remove_bimeras.R" );

    my $args = "--seq=$sequencer --map $map_file --multithread $multithread";
    if ($verbose) {
        $args .= " --verbose";
    }
    $args .= " -- " . join( " ", @rundirs );

    return ( split( /\s/, R( $script, $args ) ) );
}

sub dada2_classify {
    my $fasta       = shift;
    my $taxonomies  = shift;
    my $multithread = shift;

    my @taxonomies = @$taxonomies;
    @taxonomies = map {"--tax=$_"} @taxonomies;
    $taxonomies = join( ' ', @taxonomies );

    my $script = catfile( $pipelineDir, "scripts", "assign_taxonomies.R" );

# The way this script is written, it writes to log instead of .stdout and .stderr . Can we tee to both?
    my $args = "$fasta $taxonomies --log $log --multithread $multithread";
    if ($verbose) {
        $args .= " --verbose";
    }
    my @taxonomic_refs = split( "\n", R( $script, $args ) );
    return @taxonomic_refs;
}

################################################################################
# Execute the given commands;
# The third to last argument is a user-readable message.
# the second to last argument is for dry run
# the last argument is the filehandle where messages should be logged (undef silences logging)
sub execute_and_log {
    my $cuteMsg   = pop @_;
    my $dryRun    = pop @_;
    my $lexicalFH = pop @_;

    my @ans = ();

    #$cuteMsg =~ s/^\n+|\n+$//g;    # remove leading and trailing newlines
    print $lexicalFH "$cuteMsg\n" if $lexicalFH;

    # CAREFUL, $cmd holds a reference to a variable in the caller!
    foreach my $cmd (@_) {

        # print each command
        print $lexicalFH "> $cmd\n" if $lexicalFH;

        if ( !$dryRun ) {

            # this gets the return code, STDOUT and STDERR separately
            my $res = system("$cmd 1>.perlstdout 2>.perlstderr");

            my $err;
            {
                open( my $fh, '<', '.perlstderr' );

                #my @stat = stat
                my $filesize = -s '.perlstderr';
                read $fh, $err, $filesize;
                close $fh;
            }

            if ( $res != 0 ) {

                die "system($cmd) failed with exit code: $res\n$err";
            }

            my $out;
            {
                open my $fh, '<', '.perlstdout';
                my $filesize = -s $fh;
                read $fh, $out, $filesize;
                close $fh;
            }

            if ($verbose) {
                $logTee->print("Output files:\n$out\n");
            }

            push @ans, { "stdout" => $out, "stderr" => $err };
        }

    }
    print $lexicalFH "\n" if $lexicalFH;

    unlink ".perlstdout", ".perlstderr";
    return (@ans);
}

# Moves the file  to the project directory and prefixes it with the project name
# separated by underscore
# Returns the new filename
sub move_to_project {
    my $proj   = shift;    # global
    my $file   = shift;
    my $dryRun = shift;

    my $base = fileparse($file);
    my $dest = File::Spec->catfile( $projDir, $proj . "_$base" );

    my $cmd = join( " ", ( "mv", "\"$file\"", "\"$dest\"" ) );
    $logTee->print($cmd) if $verbose;

    if ( !$dryRun ) {
        execute_and_log( $cmd, undef, 0, "" );
    }

    return $dest;
}

sub copy_to_project {
    my $proj = shift;    # global
    my $file = shift;

    my $base = fileparse($file);
    my $dest = File::Spec->catfile( $projDir, $proj . "_$base" );
    my $cmd  = join( " ", ( "cp", "\"$file\"", "\"$dest\"" ) );
    execute_and_log( $cmd, undef, 0, "" );
    return $dest;
}

# Renames each file ending in ".tmp" by removing ".tmp" and inserting $_[0]
# before the last remaining extension (found by ".") (or use $_[1] as 2 to insert
# before the last two extensions)
sub rename_temps {
    my $preExtension           = shift;
    my $insertBeforeExtensions = shift;
    if ( !defined $insertBeforeExtensions ) {
        $insertBeforeExtensions = 1;
    }

    my @temps = glob("*.tmp");
    foreach (@temps) {
        my $base = basename( $_, ".tmp" );

        my $newName;
        if ($preExtension) {
            my @parts = split( /\./, $base );

            my @firstPart = @parts[ 0 .. (
                    ( scalar @parts ) - $insertBeforeExtensions - 1 ) ];
            my @secondPart
                = @parts[ ( ( scalar @parts ) - $insertBeforeExtensions )
                .. ( ( scalar @parts ) - 1 ) ];

            $newName
                = join( ".", ( @firstPart, $preExtension, @secondPart ) );
        }
        else {
            $newName = $base;
        }
        system("mv $_ $newName") == 0
            or die "Failed to rename temporary file: $_. Exit code $?";
    }
}

exit 0;
