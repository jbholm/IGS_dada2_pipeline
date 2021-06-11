#!/usr/local/packages/perl-5.30.2/bin/perl

=head1 NAME

  part2_illumina_dada2_devel.pl  

=head1 DESCRIPTION

  Given a comma-separated list of Illumina amplicon runs and the variable region for the amplicons,
  and assuming these runs are also the names of the directories containing sequence abundance tables 
  written by R(dada2) saveRDS(seqtab, "dada2_abundance_table.rds") for those runs, this script will:

    -Read in each .rds file to R
    -Remove chimeras
    -Write a chimera-removed abundance table
    -Assign taxonomy via SILVA
    -Assign taxonomy via RDP
    -Classify sequences with PECAN

    -Default: 
      Classify ASVs with PECAN and SILVA
      For count table, write only PECAN classifications using vaginal merging rules.

  *** If providing multiple runs, this script cannot be q-submitted at this time.***
  *** Run from qlogin within the project directory produced in Part 1.***
  *** Recommended running w/multiple runs:

      From RHEL8:
      screen
      qlogin -P jravel-lab -l mem_free=16G -q interactive.q
      part2_illumina_dada2.pl -i <comma-separated-input-run-names> -v <variable-region> -p <project-ID>

      You can close the terminal at any time. 
      To return to this process the next time you login simply type:
      screen -r

=head1 SYNOPSIS

  part2_illumina_dada2.pl -i <input-run-names> -v <variable-region> -p <project-ID>

=head1 OPTIONS

=over

=item B<--input-run-names, -i>
  Comma-separated list of input run names (directories)

=item B<--variable-region>|B<-v> {V3V4, V4, ITS, FULL-LENGTH}
  The targeted variable region. V3V4 sets the default B<--tax> to PECAN-SILVA. 
  V4 sets the default B<--tax> to SILVA. ITS sequences are taxonomically assigned 
  from UNITE (overriding B<--tax> and B<--oral>). To reiterate,
  B<-v> ITS causes B<--tax> to be ignored.

=item B<--project-ID, -p>
  Provide the project ID

=item B<--map, -m>
  A tab-delimited file listing the samples to add to this project. The file must
  have two columns named RUN.PLATEPOSITION and sampleID. RUN.PLATEPOSITION
  must be identical to the sample names in the "corrected" map created by Part 1.
  Default: project_map.txt.

=item B<--pacbio> 
  Sets the B<--tax> default to SILVA138forPB and uses PacBio-specific chimera removal
  parameters.

=item B<--tax, -t> {PECAN-SILVA, SILVA, SILVA138forPB, PECAN, HOMD, UNITE}
  Override the default taxonomic reference(s) used for classifying ASVs. For all 
    references besides PECAN and SPECIES, assignment is done by DADA2's implementation of
    the RDP classifier.
  PECAN-SILVA - Initially assigns ASVs by SILVAv132. Any ASVs assigned to the 
    genuses Lactobacillus, Gardnerella,
    Prevotella, Sneathia, Atopobium, Shuttleworthia, or Saccharibacteria are
    classified by PECAN instead. PECAN-only abundance tables are also produced.
  PECAN - PECAN is the only classifier run.
  SILVA - Assignments from SILVA132 down to genus level.
  SILVA138forPB - Assignments from the SILVA138 database down to species level. Use for PacBio (full-length) reads only.
  HOMD - Assigns from the Human Oral Microbiome Database, usually to species level.
  UNITE - Assigns from the UNITE database, usually to species level. Appropriate only for the fungal ITS region.

=item B<--oral>
  Overrides B<-v>, B<--tax> and B<--pacbio>, using HOMD for all taxonomic assignment

=item B<--nocsts>
  Flag to skip CST assignment. By default, CSTs are assigned if the V3V4 variable region is being analyzed and samples are non-oral. Assignment data are written to *_PECAN_taxa-merged_StR_CST.csv.


=item B<--notVaginal>
  Use when project is NOT just vaginal sequences to prevent renaming of several
  SILVA-derived taxonomic assignments, including two that are renamed to 
  BVAB_TM7.

=item B<--noreport>
  Flag to skip creating HTML report.

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

BEGIN
{
    use File::Spec::Functions;
    use File::Basename;

    $pipelineDir = dirname(__FILE__);
    $scriptsDir  = catdir($pipelineDir, "scripts");

}
use lib $scriptsDir;    # .pm files in ./scripts/ can be loaded
use lib "/usr/local/packages/perl-5.30.2/lib/site_perl/5.30.2/";

require Stats_gen;
require Version;

use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
require Cwd;
use Cwd qw(abs_path getcwd);
use File::Temp qw/ tempfile /;
require POSIX;
require IO::Tee;
use Data::Dumper;

#use Email::MIME;
#use Email::Sender::Simple qw(sendmail);
use File::Spec;
require IO::Tee;
use JSON;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

my $csts     = 1;
my $pacbio   = 0;
my $report   = 1;
my $map_file = "";

# this is the way it is only to preserve the interface of --notVaginal. In the future, please change to --no-vaginal
my $vaginal    = 1;
my $notVaginal = 0;
my @strategies = ();
GetOptions(
           "input-runs|i=s"      => \my $inRuns,
           "variable-region|v=s" => \my $region,
           "project-ID|p=s"      => \my $project,
           "map|m=s"             => \$map_file,
           "overwrite|o=s"       => \my $overwrite,
           "help|h!"             => \my $help,
           "debug"               => \my $debug,
           "verbose!"            => \my $verbose,
           "dry-run"             => \my $dryRun,
           "skip-err-thld"       => \my $skipErrThldStr,
           "notVaginal"          => \$notVaginal,
           "tax|t=s"             => \@strategies,
           "oral"                => \my $oral,
           "csts!"               => \$csts,
           "pacbio!"             => \$pacbio,
           "report!"             => \$report
          )

  or pod2usage(verbose => 0, exitstatus => 1);

### REASONS FOR STOPPING IMMEDIATELY
if ($help)
{
    pod2usage(verbose => 2, exitstatus => 0);
    exit 1;
}
if (!$project)
{
    die "Please provide a project name\n\n";
}

my $projDir = Cwd::cwd;

if (!$inRuns)
{
    die "Please provide (a) run ID(s)\n\n";
}

##split the list of runs to an array
my @runs = split(",", $inRuns);

# Refine and validate the run paths (no duplicates, must exist on filesystem)
my %seen;
foreach (@runs, $map_file)
{
    if ($_)
    {    # If the variable is a non-empty string...
        my $copy = $_;
        $copy =~ s/\/$//;    # remove any trailing slash
        my $relPath = $copy;
        $copy = abs_path($relPath);    # get abs path
            # At the same time, abs_path returns undef if the path doesn't exist
            # so we can verify the existence of each file and directory
        if (!defined $copy)
        {
            print die $_
              . " not found.\n"
              . "Current working directory: "
              . getcwd() . "\n";
        }
        if ($seen{basename($copy)}++)
        {
            die "$_ is a duplicate run name and not allowed\n";
        }
        $_ = $copy;  # external variable referenced by path has now been edited.
    }
}

# PRINT TO LOG ASAP
my $log = "$project" . "_part2_16S_pipeline_log.txt";
open my $logFH, ">>$log" or die "Cannot open $log for writing: $OS_ERROR";
print $logFH "This file logs the progress of "
  . scalar(@runs)
  . " runs for $project 16S amplicon sequences through the illumina_dada2.pl pipeline.\n";

my $logTee = new IO::Tee(\*STDOUT, $logFH);
print $logTee "PIPELINE VERSION: " . Version::version() . "\n";
my $time = POSIX::strftime("%Y-%m-%d %H:%M:%S", localtime(time));
print $logTee "$time\n\n";

local $SIG{__WARN__} = sub {
    print $logFH "WARNING: $_[0]";
    print "WARNING: $_[0]";
    print STDERR "WARNING: $_[0]";
};    # Warnings go to log file, stdout, and stderr

my $params_hashref = config();

if (!exists $ENV{"LD_LIBRARY_PATH"})
{
    $ENV{"LD_LIBRARY_PATH"} = "";
}
$ENV{'LD_LIBRARY_PATH'} =
  $ENV{'LD_LIBRARY_PATH'} . ":/usr/local/packages/gcc/lib64";
my $output = `python2 --version 2>&1`;
if ($? == -1)
{
    $ENV{'PATH'} = $ENV{'PATH'} . ":" . dirname($params_hashref->{"python2"});
}
$ENV{'PYTHONPATH'} = "";

if ($notVaginal)
{
    $vaginal = 0;
}

if ($pacbio && !$region)
{
    $region = "FULL-LENGTH";
}
if (!$region)
{
    die
      "Please provide a variable region (-v), V3V4, V4, ITS, or FULL-LENGTH\n";
}
if (
    List::Util::none {$_ eq $region}
    ('V3V4', 'V4', 'ITS', 'FULL-LENGTH', 'full-length')
   )
{
    die
      "Illegal --variable-region (-v). V3V4, V4, ITS, and FULL-LENGTH are accepted.";
}

my $models;
my %taxonomy_flags;
my $pecan;
{

    if ($region eq "ITS")
    {
        @strategies = ("UNITE");
        print $logTee
          "Using the UNITE reference database because ITS region was targeted.\n";
    } elsif ($oral)
    {
        @strategies = ("HOMD");
        print $logTee
          "Using the HOMD reference database because samples are oral.\n";
    } elsif (scalar @strategies == 0)
    {
        if ($pacbio)
        {
            @strategies = ("SILVA138forPB");
        } elsif ($region eq "V3V4")
        {
            @strategies = ("PECAN-SILVA");
        } elsif ($region eq "V4")
        {
            @strategies = ("SILVA");
        } else
        {
            die
              "Unsure what taxonomic reference to use because known values of --tax, --region not given, and --oral and --pacbio absent.";
        }
    }

    my $ps = grep (/^PECAN-SILVA$/, @strategies);
    my $s  = grep (/^SILVA$/,       @strategies);
    my $s138 = grep(/^SILVA138forPB$/, @strategies);
    my $p    = grep (/^PECAN$/, @strategies);
    my $h    = grep(/^HOMD$/, @strategies);
    my $u    = grep(/^UNITE$/, @strategies);
    if ($ps + $s + $s138 + $p + $h + $u == scalar @strategies) { }
    else
    {
        die "Illegal taxonomic reference. Legal values for --tax|-t are "
          . "PECAN-SILVA, SILVA, SILVA138forPB, PECAN, HOMD, and UNITE.";
    }

# Which taxonomies need to be assigned from using DADA2? Which need to be assigned by SPINGO?
    %taxonomy_flags =
      (SILVA138 => 0, SILVA => 0, HOMD => 0, UNITE => 0, PECAN => 0);

    if ($ps)
    {
        print $logTee "Using the PECAN-SILVA assignment scheme.\n";
        $models = "/local/projects-t2/jholm/PECAN/v1.0/V3V4/merged_models/";
        if ($csts)
        {
            print $logTee "Using Valencia to assign samples to CSTs\n";
        } else
        {
            print $logTee "Skipping CST assignment.\n";
        }
        $taxonomy_flags{PECAN} = 1;
        $taxonomy_flags{SILVA} = 1;
    }
    if ($s)
    {
        print $logTee "Using SILVA-only assignment scheme.\n";
        $taxonomy_flags{SILVA} = 1;
    }
    if ($p)
    {
        print $logTee "Using the PECAN-only assignment scheme.\n";
        $models = "/local/projects-t2/jholm/PECAN/v1.0/V3V4/merged_models/";
        if ($csts)
        {
            print $logTee "Using Valencia to assign samples to CSTs\n";
        } else
        {
            print $logTee "Skipping CST assignment.\n";
        }
        $taxonomy_flags{PECAN} = 1;
    }
    if ($s138)
    {
        print $logTee "Using SILVA_v138-for-PacBio assignment scheme.\n";
        $taxonomy_flags{SILVA138forPB} = 1;
    }
    if (grep(/^HOMD$/, @strategies))
    {
        print $logTee "Using HOMD taxonomy\n";
        $taxonomy_flags{HOMD} = 1;
    }
    if (grep(/^UNITE$/, @strategies))
    {
        print $logTee "Using UNITE taxonomy\n";
        $taxonomy_flags{UNITE} = 1;
    }
}

print $logTee "\n";
####################################################################
##                               MAIN
####################################################################
$logTee->print("Run(s):\n");
map {$logTee->print("$_\n")} @runs;
$logTee->print("\n");

my $run_info_hash = combine_run_metadata(@runs);

my @classifs = ();
my $projabund;

# the un-annotated abundance table signals dada2 already completed
# my $cmd = "rm -f *-dada2_abundance_table.rds";
# execute_and_log( $cmd, *STDOUT, $dryRun, "" );
my $cmd;

##loop over array to copy the file to the main current working directory
## using the array string to also add a name
print $logTee "---Combining "
  . scalar(@runs)
  . " abundance table(s) and removing bimeras.\n";

my ($abundRds, $abund, $stats, $fasta) = (
                                          "all_runs_dada2_abundance_table.rds",
                                          "all_runs_dada2_abundance_table.csv",
                                          "DADA2_stats.txt",
                                          "all_runs_dada2_ASV.fasta"
                                         );
($abundRds, $abund, $stats, $fasta) =
  map {move_to_project($project, $_, 1)} ($abundRds, $abund, $stats, $fasta);

if (List::Util::any {!-e $_} ($abundRds, $abund, $fasta, $stats))
{
    print $logTee "Combining runs and removing bimeras.\n";
    ($abundRds, $abund, $fasta, $stats) = dada2_combine($pacbio, \@runs);

    ($abundRds, $abund, $fasta, $stats) = map {move_to_project($project, $$_)}
      (\$abundRds, \$abund, \$fasta, \$stats);

    $logTee->print("Outputs:\n");
    foreach (($abundRds, $abund, $fasta, $stats))
    {
        print $logTee "$_\n";
    }
} else
{
    print $logTee "Chimeras already removed. Skipping...\n";
}

# Give ASV's unique and easy-to-look-up IDs
$cmd = "$scriptsDir/rename_asvs.py $fasta -p $project";
execute_and_log($cmd, $logTee, $dryRun,
     "---Renaming ASVs in FASTA, abundance tables, and classification key(s).");

my $species_file;
my @full_classif_csvs;
my @two_col_classifs;

my @dada2_taxonomy_list;
foreach ("SILVA138forPB", "SILVA", "SILVA-PECAN", "UNITE", "HOMD")
{
    if ($taxonomy_flags{$_})
    {
        if ($_ eq "SILVA_PECAN")
        {
            push @dada2_taxonomy_list, "SILVA";
        } else
        {
            push @dada2_taxonomy_list, $_;
        }
    }
}
if (scalar @dada2_taxonomy_list > 0)
{
    my @output =
      map {"$project" . "_$_.classification.csv"} @dada2_taxonomy_list;
    my @dada2Output;
    if (List::Util::any {!-e $_} @output)
    {
        print $logTee
          "Using DADA2's RDP Classifier implementation to classify amplicon sequence variants (ASVs) with taxonomies:\n";
        print $logTee join(",", @dada2_taxonomy_list) . "\n";

        @dada2Output =
          dada2_classify($fasta, \@dada2_taxonomy_list);    # assignment

        ### @dada2Output now has every file that needs to be renamed
        print $logTee "---Renaming dada2 output files for project...\n";
        @dada2Output = move_to_project($project, @dada2Output);
    } else
    {
        print $logTee
          "DADA2 has already been used to classify ASVs with taxonomies:\n";
        print $logTee join(",", @dada2_taxonomy_list) . "\n";
        print $logTee "Skipping.\n";
    }
    push(@full_classif_csvs, @dada2Output);

}

my $pecanFile = "MC_order7_results.txt";
if ($taxonomy_flags{PECAN})
{
    if (!-e $pecanFile)
    {

        my $cmd =
          "/local/projects/pgajer/devel/MCclassifier/bin/classify -d $models -i $fasta -o .";
        execute_and_log(
            $cmd,
            $logTee,
            0,
            "---Classifying ASVs with $region PECAN models (located in $models)"
        );

    } else
    {
        print $logTee
          "ASVs have already been classified with V3V4 PECAN models.\n";
        print $logTee "Located at:\n";
        print $logTee "\t" . File::Spec->catfile($projDir, $pecanFile) . "\n";
    }
    push @two_col_classifs, $pecanFile;
}

#################################################################
#########APPLY CLASSIFICATIONS TO COUNT TABLES ###################
#################################################################

my $vopt = "";
if ($vaginal && !$oral)
{
    $vopt = "--vaginal";
}
foreach (@full_classif_csvs)
{

 # Apply classifications from the various classification files to the count
 # tables. This is an error-prone task because we have different classification
 # formats from the previous step. We reply on the filenames to tell us what the
 # file format is, then take the appropriate action.
    my $db;
    if ($_ =~ /SILVA138forPB/)
    {
        $db  = "SILVA138forPB";
        $cmd = "$scriptsDir/combine_tx_for_ASV.pl -s $_ -c $abund $vopt";
    } elsif ($_ =~ /SILVA\./)
    {
        $db  = "SILVA";
        $cmd = "$scriptsDir/combine_tx_for_ASV.pl --s $_ -c $abund $vopt";
    } elsif ($_ =~ /HOMD/)
    {
        $db = "HOMD";
        $cmd =
          "$scriptsDir/combine_tx_for_ASV.pl --homd-file $_ -c $abund $vopt";
    } elsif ($_ =~ /UNITE/)
    {
        $db  = "UNITE";
        $cmd = "$scriptsDir/combine_tx_for_ASV.pl -u $_ -c $abund $vopt";
    }
    my $outfile = join(
                       "_",
                       (
                        $project, basename($abund, ".csv"),
                        $db,      "taxa-merged.csv"
                       )
                      );
    if (!-e $outfile)
    {
        my $msg = "---Labeling count table with $db taxa.";
        execute_and_log($cmd, $logTee, $dryRun, $msg);
    }

}
foreach (@two_col_classifs)
{
    my $db;
    if ($_ =~ /MC_order7_results/)
    {
        $db = "PECAN";
    } elsif ($_ =~ /SPINGO/)
    {
        $db = "SPINGO";
    }
    my $outfile = join(
                       "_",
                       (
                        $project, basename($abund, ".csv"),
                        $db,      "taxa-merged.csv"
                       )
                      );
    if (!-e $outfile)
    {
        my $cmd =
          "$scriptsDir/PECAN_tx_for_ASV.pl -p $_ -c $abund -t $db $vopt";

        my $msg = "---Labeling count table with $db taxa.";
        execute_and_log($cmd, $logTee, $dryRun, $msg);
    }

}
if (grep (/^PECAN-SILVA$/, @strategies))
{
    ### APPLY PECAN+SILVA CLASSIFICATIONS TO COUNT TABLE (V3V4) ####
    #################################################################
    my $db = "PECAN-SILVA";
    my $outfile = join(
                       "_",
                       (
                        $project, basename($abund, ".csv"),
                        $db,      "taxa-merged.csv"
                       )
                      );
    my $msg       = "---Labeling count table with $db taxa.";
    my $silvaFile = "";
    foreach (@full_classif_csvs)
    {
        if ($_ =~ /SILVA\./)
        {
            $silvaFile = $_;
        }
    }
    $cmd =
      "$scriptsDir/combine_tx_for_ASV.pl -p $pecanFile -s $silvaFile -c $abund $vopt";
    execute_and_log($cmd, $logTee, $dryRun, $msg);
}
rename_temps("");

if ($csts && $pecan)
{
    $ENV{'LD_LIBRARY_PATH'} =
      $ENV{'LD_LIBRARY_PATH'} . ':/usr/local/packages/python-3.5/lib';
    my $pecanCountTbl =
      basename($projabund, (".csv")) . "_PECAN_taxa-merged.csv";

    if (!-e basename($pecanCountTbl, (".csv")) . "_StR_CST.csv")
    {
        print $logTee "---Assigning CSTs with Valencia\n";
        my $cmd =
          "python3 -s $scriptsDir/valencia_wrapper.py "
          . catdir($pipelineDir, "ext", "valencia",
                   "CST_profiles_jan28_mean.csv")
          . " $pecanCountTbl";
        print $logTee "\tcmd=$cmd\n" if $dryRun || $debug;

        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;
    } else
    {
        print $logTee "---Count table with CSTs already exists.\n";
    }
}

my $final_merge = glob("*.taxa-merged.csv");
unlink glob "*.taxa.csv";
my $final_ASV_taxa = glob("*.asvs+taxa.csv");

$cmd =
  "$pipelineDir/report/report16s.sh '$projDir' --runs @runs --project $project";
execute_and_log($cmd, $logTee, $dryRun, "Creating report...");
print $logTee "---Final files succesfully produced!\n";
print $logTee
  "Final merged read count table: $final_merge\nFinal ASV table with taxa: $final_ASV_taxa\nFinal ASV count table: $abund\nASV sequences: all_runs_dada2_ASV.fasta\n"
  ;    #Read survival stats: $finalStats\n";
$logTee->close;

####################################################################
##                               SUBS
####################################################################
sub read_json
{
    my $file = shift;
    my $mode = shift;

    my $json;
    {
        local $/;    #Enable 'slurp' mode
        if (-e $file)
        {
            open(my $FH, $mode, $file);
            seek $FH, 0, 0 or die;
            $json = <$FH>;
            close $FH;
        } else
        {
            return undef;
        }

    }

    my $hash = {};
    eval {$hash = decode_json($json);};

    return $hash;
}

sub config
{
    my %arg             = @_;
    my $new_params      = delete $arg{new_params} // {};
    my $orig_param_file = "$pipelineDir/config.json";

    my $params = read_json($orig_param_file, "<");

    return $params;
}

sub combine_run_metadata
{
    my $all_run_info = {};
    foreach my $rundir (@_)
    {
        my $run = basename($rundir);
        $all_run_info->{$run} = read_json(catfile($rundir, ".meta.json"), "<");
        $all_run_info->{$run}->{"path"} = $rundir;
        if (exists $all_run_info->{$run}->{"samples"}
            && ref $all_run_info->{$run}->{"samples"} eq ref {})
        {
            foreach my $sample (keys %{$all_run_info->{$run}->{"samples"}})
            {
                if (ref $all_run_info->{$run}->{"samples"}->{$sample} eq ref {})
                {
                    foreach my $filepath (
                        values %{$all_run_info->{$run}->{"samples"}->{$sample}})
                    {
  # an absent file will be decoded from JSON to either an empty hashref or undef
                        if (!ref($filepath) && defined $filepath)
                        {
                            $filepath = catfile($rundir, $filepath);
                        }
                    }
                } else
                {
                    warn "Malformed metadata: $run: samples: $sample.\n";
                }
            }
        } else
        {
            warn "Malformed metadata: $run: samples.\n";
        }
    }

    my $metadata = {runs => $all_run_info};
    open my $metadataFH, ">.meta.json";
    print $metadataFH encode_json($metadata);
    close $metadataFH;

    return $metadata;
}

sub readTbl
{

    my $file = shift;

    if (!-f $file)
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
        my ($id, $t) = split(/\s+/, $_);
        $tbl{$id} = $t;
    }
    close IN;

    return %tbl;
}

sub R
{
    my $script = shift;
    my $args   = shift;

    open my $scriptFH, "<$script", or die "cannot read header of $script: $!\n";
    my $shebang = <$scriptFH>;
    chomp $shebang;
    $shebang =~ s/#!//;    # Remove shebang itself
    my $pathsep = catfile('', '');

    # my $outR    = catfile( $projDir, basename($script) . "out" );

    my $cmd = "$shebang --verbose $script $args";
    my ($stdout) = execute_and_log($cmd, undef, 0, "");

    # check_R_for_error(\$stdout); execute_and_log should die on error
    return ($stdout);
}

# Combines all part1 dada2 counts and stats into single files
# $_[0] Evaluated in a boolean context, gives whether PacBio full-length 16S
#   data is being used
# Returns the names of the three output file: abundance table (as RDS), abundance table (as
#   CSV), and stats
sub dada2_combine
{
    my $pacbio  = shift;
    my $rundirs = shift;
    my @rundirs = @$rundirs;

    my $sequencer = $pacbio ? "PACBIO" : "ILLUMINA";

    my $script = catfile($pipelineDir, "scripts", "remove_bimeras.R");
    my $args   = "--seq=$sequencer --map $map_file -- " . join(" ", @rundirs);

    return (split(/\s/, R($script, $args)));
}

sub dada2_classify
{
    my $fasta      = shift;
    my $taxonomies = shift;

    my @taxonomies = @$taxonomies;
    @taxonomies = map {"--tax=$_"} @taxonomies;
    $taxonomies = join(' ', @taxonomies);

    my $script = catfile($pipelineDir, "scripts", "assign_taxonomies.R");
    my $args   = "$fasta $taxonomies";
    my @taxonomic_refs = split("\n", R($script, $args));
    return @taxonomic_refs;
}

sub check_R_for_error
{
    my $outR = shift;
    open IN, "$outR" or die "Cannot open $outR for reading: $OS_ERROR\n";
    my $exitStatus = 1;

    foreach my $line (<IN>)
    {
        if ($line =~ /Error/)
        {
            print "R script crashed. Check $outR for details\n";
            print "";
            $exitStatus = 0;
            exit;
        }
    }
    close IN;
}
################################################################################
# Execute the given commands;
# the second to last argument must evaluate false to prevent loggin to log file.
# the last argument is 1 if this is a dry run
sub execute_and_log
{
    my $cuteMsg   = pop @_;
    my $dryRun    = pop @_;
    my $lexicalFH = pop @_;

    # Print to STDOUT if filehandle was not provided
    if (!$lexicalFH)
    {
        $lexicalFH = *STDOUT;
    }

    my @ans = ();
    $cuteMsg =~ s/^\n+|\n+$//g;    # remove leading and trailing newlines
    print $lexicalFH "$cuteMsg\n";

    # CAREFUL, $cmd holds a reference to a variable in the caller!
    foreach my $cmd (@_)
    {

        # print each command
        print $lexicalFH "> $cmd\n";

        if (!$dryRun)
        {
            my $res = system("$cmd 1>.perlstdout 2>.perlstderr");
            if ($res != 0)
            {
                open(my $fh, '<', '.perlstderr');

                #my @stat = stat
                my $filesize = -s '.perlstderr';
                my $err;
                read $fh, $err, $filesize;
                close $fh;
                die "system($cmd) failed with exit code: $?\n$err";
            }
            open my $fh, '<', '.perlstdout';
            my $filesize = -s $fh;
            my $out;
            read $fh, $out, $filesize;
            if ($verbose)
            {
                print $out;
                open(my $fh, '<', '.perlstderr');
                my $filesize = -s '.perlstderr';
                my $err;
                read $fh, $err, $filesize;
                close $fh;
            }
            close $fh;
            push @ans, $out;
        }

    }
    print $lexicalFH "\n";

    # unlink ".perlstdout", ".perlstderr";
    return (@ans);
}

# Moves the file  to the project directory and prefixes it with the project name
# separated by underscore
# Returns the new filename
sub move_to_project
{
    my $proj   = shift;    # global
    my $file   = shift;
    my $dryRun = shift;

    my $base = fileparse($file);
    my $dest = File::Spec->catfile($projDir, $proj . "_$base");
    if (!$dryRun)
    {
        my $cmd = "mv $file $dest";
        execute_and_log($cmd, undef, 0, "");
    }

    return $dest;
}

sub copy_to_project
{
    my $proj = shift;    # global
    my $file = shift;

    my $base = fileparse($file);
    my $dest = File::Spec->catfile($projDir, $proj . "_$base");
    my $cmd  = "cp $file $dest";
    execute_and_log($cmd, undef, 0, "");
    return $dest;
}

# Renames each file ending in ".tmp" by removing ".tmp" and inserting $_[0]
# before the last remaining extension (found by ".") (or use $_[1] as 2 to insert
# before the last two extensions)
sub rename_temps
{
    my $preExtension           = shift;
    my $insertBeforeExtensions = shift;
    if (!defined $insertBeforeExtensions)
    {
        $insertBeforeExtensions = 1;
    }

    my @temps = glob("*.tmp");
    foreach (@temps)
    {
        my $base = basename($_, ".tmp");

        my $newName;
        if ($preExtension)
        {
            my @parts = split(/\./, $base);

            my @firstPart =
              @parts[0 .. ((scalar @parts) - $insertBeforeExtensions - 1)];
            my @secondPart = @parts[((scalar @parts) - $insertBeforeExtensions)
              .. ((scalar @parts) - 1)];

            $newName = join(".", (@firstPart, $preExtension, @secondPart));
        } else
        {
            $newName = $base;
        }
        system("mv $_ $newName") == 0
          or die "Failed to rename temporary file: $_. Exit code $?";
    }
}

exit 0;
