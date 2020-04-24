#!/usr/bin/env perl

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

      From RHEL5:
      screen
      qlogin -P jravel-lab -l mem_free=500M -q interactive.q
      export LD_LIBRARY_PATH=/usr/local/packages/gcc/lib64
      source /usr/local/packages/usepackage/share/usepackage/use.bsh
      use python-2.7
      use qiime
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

=item B<--variable-region>|B<-v> {V3V4, V4, ITS}
  The targeted variable region. V3V4 sequences are taxonomically assigned from 
  SILVA and PECAN; abundance tables for the individual assignments as well as
  a table combining the two assignments, are produced. V4 sequences are 
  taxonomically assigned from SILVa. ITS sequences are taxonomically assigned 
  from UNITE.

=item B<--desire-pecan-models, -v>
  V4 or V3V4 or ITS

=item B<--project-ID, -p>
  Provide the project ID

=item B<--notVaginal>
  optional flag. Use when project is NOT
  just vaginal sequences.

=item B<--pecan-silva, --pecan+silva>
  Ignored. If B<-v V3V4> is given, both references are always used.

=item B<--nocsts>
  Flag to skip CST assignment. By default, CSTs are assigned if the V3V4 variable region is being analyzed and samples are non-oral. Assignment data are written to *_PECAN_taxa-merged_StR_CST.csv.

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

BEGIN {
    use File::Spec::Functions;
    use File::Basename;

    $pipelineDir = dirname(__FILE__);
    $scriptsDir  = catdir( $pipelineDir, "scripts" );

}
use lib $scriptsDir;    # .pm files in ./scripts/ can be loaded

require Stats_gen;
require Version;

use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
require Cwd;
use File::Temp qw/ tempfile /;

#use Email::MIME;
#use Email::Sender::Simple qw(sendmail);
use File::Spec;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

my $csts   = 1;
my $report = 1;

# this is the way it is only to preserve the interface of --notVaginal. In the future, please change to --no-vaginal
my $vaginal    = 1;
my $notVaginal = 0;
GetOptions(
    "input-runs|i=s"      => \my $inRuns,
    "variable-region|v=s" => \my $region,
    "project-ID|p=s"      => \my $project,
    "overwrite|o=s"       => \my $overwrite,
    "help|h!"             => \my $help,
    "debug"               => \my $debug,
    "dry-run"             => \my $dryRun,
    "skip-err-thld"       => \my $skipErrThldStr,
    "notVaginal"          => \$notVaginal,
    "pecan-silva"         => \my $pecanSilva,
    "oral"                => \my $oral,
    "csts!"               => \$csts,
    "report!"             => \$report
  )

  or pod2usage( verbose => 0, exitstatus => 1 );

if ($help) {
    pod2usage( verbose => 2, exitstatus => 0 );
    exit 1;
}

$ENV{'LD_LIBRARY_PATH'} =
  $ENV{'LD_LIBRARY_PATH'} . ":/usr/local/packages/gcc/lib64";
my $R = "/usr/local/packages/r-3.4.0/bin/R";

if ($notVaginal) {
    $vaginal = 0;
}
if ( !$region ) {
    print "Please provide a variable region (-v), V3V4 or V4\n";
    pod2usage( verbose => 2, exitstatus => 0 );
    exit 1;
}
if ( List::Util::none { $_ eq $region } ( 'V3V4', 'V4', 'ITS' ) ) {
    die "Illegal --variable-region (-v). V3V4, V4, and ITS are accepted.";
}
my $models;

my @taxonomies;
if ( $region eq 'V3V4' && !$oral ) {
    print "Using PECAN and SILVA taxonomies\n";
    $models = "/local/projects-t2/jholm/PECAN/v1.0/V3V4/merged_models/";
    if ($csts) {
        print "Using Valencia to assign samples to CSTs\n";
    } else {
        print "Skipping CST assignment.\n";
    }
} else {
    if ( $region eq 'V4' && !$oral ) {
        print "Using SILVA taxonomy\n";
    }
    if ( $region eq 'V4' && $oral ) {
        print "Using HOMD taxonomy only\n";
    }
    if ( $region eq 'ITS' ) {
        print "Using UNITE taxonomy only\n";
    }
}

if ( $region eq 'V3V4' && !$models ) {
    print "Please provide a valid variable region\n\n";
    pod2usage( verbose => 2, exitstatus => 0 );
    exit 1;
}

if ( !$project ) {
    print "Please provide a project name\n\n";
    pod2usage( verbose => 2, exitstatus => 0 );
    exit 1;
}

if ( !$inRuns ) {
    print "Please provide (a) run ID(s)\n\n";
    pod2usage( verbose => 2, exitstatus => 0 );
    exit 1;
}

####################################################################
##                               MAIN
####################################################################

my $projDir = Cwd::cwd;

##split the list of runs to an array
my @runs = split( ",", $inRuns );

my $log = "$project" . "_part2_16S_pipeline_log.txt";

open my $logFH, ">>$log" or die "Cannot open $log for writing: $OS_ERROR";

print $logFH "PIPELINE VERSION: " . Version::version() . "\n";
print $logFH "This file logs the progress of "
  . scalar(@runs)
  . " runs for $project 16S amplicon sequences through the illumina_dada2.pl pipeline.\n";

my $abundance = "all_runs_dada2_abundance_table.csv";
my $projabund = $project . "_" . $abundance;
my $silva     = "silva_classification.csv";
my $unite     = "unite_classification.csv";
my $homd      = "homd_classification.csv";    ##########################for HOMD
my $projSilva = $project . "_" . $silva;
my $projUNITE = $project . "_" . $unite;
my $projHOMD  = $project . "_" . $homd;       ##########################for HOMD
my @classifs  = ();
my $cmd;

# the un-annotated abundance table signals dada2 already completed
if ( !-e $projabund ) {
    my $cmd = "rm -f *-dada2_abundance_table.rds";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0
      or die "system($cmd) failed with exit code: $?"
      if !$dryRun;

    ##loop over array to copy the file to the main current working directory
    ## using the array string to also add a name
    if ( scalar(@runs) > 1 ) {
        print "---Copying "
          . scalar(@runs)
          . " abundance tables to this directory & combining\n";
        print $logFH "---Copying "
          . scalar(@runs)
          . " abundance tables to this directory & combining\n";
        print $logFH "Runs:\n";
    } else {
        print "---Copying 1 abundance table to this directory\n";
        print "---Proceeding to chimera removal for 1 run\n";
        print $logFH "---Copying 1 abundance table to this directory\n";
        print $logFH "---Proceeding to chimera removal for 1 run\n";
        print $logFH "Run:\n";
    }

    foreach my $i (@runs) {
        print $logFH "$i\n";
        my $currTbl = $i . "/dada2_abundance_table.rds";
        my $newTbl  = $project . "_" . $i . "-dada2_abundance_table.rds";
        my $cmd     = "cp $currTbl $newTbl";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;
        print $logFH "$cmd\n";

        my $currStats = $i . "/dada2_part1_stats.txt";
        my $newStats  = $project . "_" . $i . "-dada2_part1_stats.txt";
        $cmd = "cp $currStats $newStats";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;
        print $logFH "$cmd\n";
    }

    print
"---Performing chimera removal on merged tables and classifying amplicon sequence variants (ASVs)\n";
    print $logFH
"---Performing chimera removal on merged tables and classifying amplicon sequence variants (ASVs)\n";
    if ( $region eq 'ITS' ) {
        dada2_combine_and_classifyITS($inRuns);
        print
"---Merged, chimera-removed abundance tables written to all_runs_dada2_abundance_table.csv\n";
        print
          "---ASVs classified via UNITE written to unite_classification.csv\n";
        print "---dada2 completed successfully\n";
        print $logFH
"---Merged, chimera-removed abundance tables written to all_runs_dada2_abundance_table.csv\n";
        print $logFH
          "---ASVs classified via UNITE written to unite_classification.csv\n";
        print $logFH "---dada2 completed successfully\n";

        print $logFH "---Renaming dada2 files for project\n";
        $cmd = "mv $abundance $projabund";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;
        print $logFH "$cmd\n";

        print "---Renaming UNITE classification file for project\n";
        print $logFH "---Renaming UNITE classification file for project\n";
        $cmd = "mv $unite $projUNITE";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;
        print $logFH "$cmd\n";

    } elsif ($oral) {
        dada2_combine_and_classifyHOMD($inRuns);

        print
"---Merged, chimera-removed abundance tables written to all_runs_dada2_abundance_table.csv\n";
        print
          "---ASVs classified via HOMD written to homd_classification.csv\n";
        print
"---Final ASVs written to all_runs_dada2_ASV.fasta for classification via PECAN\n";
        print "---dada2 completed successfully\n";

        print $logFH
"---Merged, chimera-removed abundance tables written to all_runs_dada2_abundance_table.csv\n";
        print $logFH
          "---ASVs classified via HOMD written to homd_classification.csv\n";
        print $logFH
          "---ASVs classified via RDP written to rdp_classification.csv\n";
        print $logFH
"---Final ASVs written to all_runs_dada2_ASV.fasta for classification via PECAN\n";
        print $logFH "---dada2 completed successfully\n";

        print "---Renaming dada2 files for project\n";
        print $logFH "---Renaming dada2 files for project\n";
        $cmd = "mv $abundance $projabund";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;
        print $logFH "$cmd\n";

        print "---Renaming HOMD classification file for project\n";
        print $logFH "---Renaming HOMD classification file for project\n";
        $cmd = "mv $homd $projHOMD";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;
        print $logFH "$cmd\n";
    } else {
        dada2_combine_and_classify($inRuns);

        print
"---Merged, chimera-removed abundance tables written to all_runs_dada2_abundance_table.csv\n";
        print
          "---ASVs classified via silva written to silva_classification.csv\n";
        print "---ASVs classified via RDP written to rdp_classification.csv\n";
        print
"---Final ASVs written to all_runs_dada2_ASV.fasta for classification via PECAN\n";
        print "---dada2 completed successfully\n";

        print $logFH
"---Merged, chimera-removed abundance tables written to all_runs_dada2_abundance_table.csv\n";
        print $logFH
          "---ASVs classified via silva written to silva_classification.csv\n";
        print $logFH
          "---ASVs classified via RDP written to rdp_classification.csv\n";
        print $logFH
"---Final ASVs written to all_runs_dada2_ASV.fasta for classification via PECAN\n";
        print $logFH "---dada2 completed successfully\n";

        print "---Renaming dada2 files for project\n";
        print $logFH "---Renaming dada2 files for project\n";
        $cmd = "mv $abundance $projabund";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;
        print $logFH "$cmd\n";

        print "---Renaming SILVA classification file for project\n";
        print $logFH "---Renaming SILVA classification file for project\n";
        $cmd = "mv $silva $projSilva";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;
        print $logFH "$cmd\n";
    }
} else {
    print "DADA2 chimera removal already finished. Skipping...\n";
}

if ( $region eq 'ITS' ) {

    # SILVA file not renamed after DADA2 combine and classify R script
    push @classifs, ( $projUNITE, "silva_classification.csv" );
} elsif ($oral) {
    push @classifs, $projHOMD;
} else {
    push @classifs, $projSilva;
}

my @combinedStats = ( "stats_file_cmp.txt", "dada2_part2_stats.txt" );
if ( !List::Util::all { return -e $_ } @combinedStats ) {

    # Combine dada2 stats from all runs, and the overall project, into one file
    Stats_gen::combine_dada2_stats( $projDir, @runs );
} else {
    print "Combined DADA2 stats files already exist at:\n";
    print "\t$projDir/stats_file_cmp.txt\n";
    print "\t$projDir/dada2_part2_stats.txt\n";
}

my $projpecan = "";
if ( $region eq 'V3V4' && !$oral ) {
    my $pecan = "MC_order7_results.txt";
    $projpecan = $project . "_" . "MC_order7_results.txt";
    if ( !-e $projpecan ) {
        print
"---Classifying ASVs with $region PECAN models (located in $models)\n";
        print $logFH
"---Classifying ASVs with $region PECAN models (located in $models)\n";
        my $fasta;
        $fasta = "all_runs_dada2_ASV.fasta";

        $cmd =
"/local/projects/pgajer/devel/MCclassifier/bin/classify -d $models -i $fasta -o .";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;

        $cmd = "mv $pecan $projpecan";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;

    } else {
        print "ASVs have already been classified with V3V4 PECAN models.\n";
        print "Located at:\n";
        print "\t$projDir/$projpecan\n";
    }
}

# Give ASV's unique and easy-to-look-up IDs
# Figure out what to do when DADA2_combine_and_classify_ITS is run (both silva
# and unite classification csvs are created)
print
"---Renaming ASVs in FASTA, abundance tables, and SILVA classification key.\n";
print $logFH
"---Renaming ASVs in FASTA, abundance tables, and SILVA classification key.\n";

# BAD! Current version of the pipeline has exception case where both SILVA-PECAN
# and SILVA-only count tables are produced, but only the first is renamed,
# according to the order of @taxonomies above. The next major release will
# rename ASVs BEFORE taxa are applied to the count table.
$cmd =
"python2 $scriptsDir/rename_asvs.py -p $project -c @classifs --pecan $projpecan";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0
  or print "$cmd failed with exit code: $?. Continuing...\n"
  if !$dryRun;
print $logFH "$cmd\n";

#################################################################
#########APPLY CLASSIFICATIONS TO COUNT TABLE ###################
#################################################################

#### APPLY NON-PECAN CLASSIFICATIONS TO COUNT TABLE (V4) ####
##############################################################
if ( $region eq 'V4' || $region eq 'V3V4' ) {
    if ( !$oral ) {
        print "---Creating count table with SILVA classifications only\n";
        print $logFH
          "---Creating count table with SILVA classifications only\n";
        $cmd = "$scriptsDir/combine_tx_for_ASV.pl -s $projSilva -c $projabund";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;
        push @taxonomies, "SILVA";
    }
}

if ( $region eq 'V3V4' ) {
    if ( !$oral ) {

        my $vopt = "";
        if ($vaginal) {
            $vopt = "--vaginal";
        }

#### APPLY PECAN-ONLY CLASSIFICATIONS TO COUNT TABLE (V3V4) ####
################################################################
        print "---Creating count table with PECAN classifications\n";
        print $logFH "---Creating count table with PECAN classifications\n";

        $cmd =
          "$scriptsDir/PECAN_tx_for_ASV.pl -p $projpecan -c $projabund $vopt";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;

        push @taxonomies, "PECAN";

        if ($csts) {
            $ENV{'LD_LIBRARY_PATH'} =
              $ENV{'LD_LIBRARY_PATH'} . ':/usr/local/packages/python-3.5/lib';
            my $pecanCountTbl =
              basename( $projabund, (".csv") ) . "_PECAN_taxa-merged.csv";

            if ( !-e basename( $pecanCountTbl, (".csv") ) . "_StR_CST.csv" ) {
                print "---Assigning CSTs with Valencia\n";
                print $logFH "---Assigning CSTs with Valencia\n";
                $cmd = "$scriptsDir/valencia_wrapper.py "
                  . catdir( $pipelineDir, "ext", "valencia",
                    "CST_profiles_jan28_mean.csv" )
                  . " $pecanCountTbl";
                print "\tcmd=$cmd\n" if $dryRun || $debug;
                system($cmd) == 0
                  or die "system($cmd) failed with exit code: $?"
                  if !$dryRun;
            } else {
                print "---Count table with CSTs already exists.\n";
            }

        }

        #### APPLY PECAN+SILVA CLASSIFICATIONS TO COUNT TABLE (V3V4) ####
        #################################################################
        print "---Creating count table with SILVA+PECAN classifications\n";
        print $logFH
          "---Creating count table with SILVA+PECAN classifications\n";
        $cmd =
"$scriptsDir/combine_tx_for_ASV.pl -p $projpecan -s $projSilva -c $projabund $vopt";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;
        push @taxonomies, "SILVA-PECAN";

    }
}

if ( $region eq 'ITS' ) {
    print "---Creating count table with UNITE classifications\n";
    print $logFH "---Creating count table with UNITE classifications\n";
    $cmd = "$scriptsDir/combine_tx_for_ASV.pl -u $projUNITE -c $projabund";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0
      or die "system($cmd) failed with exit code: $?"
      if !$dryRun;
    push @taxonomies, "UNITE";

} elsif ($oral) {    # note: same logic as above: ITS prioritizes over oral

    print "---Creating count table with HOMD classifications\n";
    print $logFH "---Creating count table with HOMD classifications\n";
    $cmd =
      "$scriptsDir/combine_tx_for_ASV.pl --homd-file $projHOMD -c $projabund";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0
      or die "system($cmd) failed with exit code: $?"
      if !$dryRun;
    push @taxonomies, "HOMD";
}

my $final_merge = glob("*_taxa-merged.csv");
unlink glob "*_taxa.csv";
my $final_ASV_taxa = glob("*_asvs+taxa.csv");

if ($report) {
    print "\nCreating report...\n";
    $cmd = "$pipelineDir/report/report16s.sh '$projDir' --runs @runs";
    execute_and_log( $cmd, $logFH, $dryRun );
}

print $logFH "---Final files succesfully produced!\n";
print $logFH
"Final merged read count table: $final_merge\nFinal ASV table with taxa: $final_ASV_taxa\nFinal ASV count table: $projabund\nASV sequences: all_runs_dada2_ASV.fasta\n"
  ;    #Read survival stats: $finalStats\n";
close $logFH;

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
        my ( $id, $t ) = split( /\s+/, $_ );
        $tbl{$id} = $t;
    }
    close IN;

    return %tbl;
}

sub dada2_combine_and_classify {
    my ($inRuns) = shift;

    my $Rscript = qq~

library("dada2")
packageVersion("dada2")
path<-getwd()

## list all of the files matching the pattern
tables<-list.files(path, pattern="-dada2_abundance_table.rds", full.names=TRUE)
stats<-list.files(path, pattern="-dada2_part1_stats.txt", full.names=TRUE)

## get the run names using splitstring on the tables where - exists
sample.names <- sapply(strsplit(basename(tables), "-"), `[`, 1)
##sample.names
##names(tables) <- sample.names

runs <- vector("list", length(sample.names))
names(runs) <- sample.names
for(run in tables) {
  cat("Reading in:", run, "\n")
  runs[[run]] <- readRDS(run)
}

runstats <- vector("list", length(sample.names))
names(runstats) <- sample.names
for(run in stats) {
  cat("Reading in:", run, "\n")
  runstats[[run]] <- read.delim(run, )
}

unqs <- unique(c(sapply(runs, colnames), recursive=TRUE))
n<-sum(unlist(lapply(X=runs, FUN = nrow)))
st <- matrix(0L, nrow=n, ncol=length(unqs))
rownames(st) <- c(sapply(runs, rownames), recursive=TRUE)
colnames(st) <- unqs
for(sti in runs) {
  st[rownames(sti), colnames(sti)] <- sti
}
st <- st[,order(colSums(st), decreasing=TRUE)]

##st.all<-mergeSequenceTables(runs)
# Remove chimeras
seqtab <- removeBimeraDenovo(st, method="consensus", multithread=TRUE)
# Assign taxonomy
silva <- assignTaxonomy(seqtab, "$pipelineDir/taxonomy/silva_nr_v138_train_set.fa.gz", multithread=TRUE)
# Write to disk
saveRDS(seqtab, "all_runs_dada2_abundance_table.rds") # CHANGE ME to where you want sequence table saved
write.csv(seqtab, "all_runs_dada2_abundance_table.csv", quote=FALSE)
write.csv(silva, "silva_classification.csv", quote=FALSE)

fc = file("all_runs_dada2_ASV.fasta")
fltp = character()
for( i in 1:ncol(seqtab))
{
  fltp <- append(fltp, paste0(">", colnames(seqtab)[i]))
  fltp <- append(fltp, colnames(seqtab)[i])
}
writeLines(fltp, fc)
rm(fltp)
close(fc)

track<-as.matrix(rowSums(seqtab))
colnames(track) <- c("nonchimeric")
write.table(track, "dada2_part2_stats.txt", quote=FALSE, append=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
 ~;
    run_R_script($Rscript);
}

sub dada2_combine_and_classifyHOMD {
    my ($inRuns) = shift;

    my $Rscript = qq~

library("dada2")
packageVersion("dada2")
path<-getwd()

## list all of the files matching the pattern
tables<-list.files(path, pattern="-dada2_abundance_table.rds", full.names=TRUE)
stats<-list.files(path, pattern="-dada2_part1_stats.txt", full.names=TRUE)

## get the run names using splitstring on the tables where - exists
sample.names <- sapply(strsplit(basename(tables), "-"), `[`, 1)
##sample.names
##names(tables) <- sample.names

runs <- vector("list", length(sample.names))
names(runs) <- sample.names
for(run in tables) {
  cat("Reading in:", run, "\n")
  runs[[run]] <- readRDS(run)
}

runstats <- vector("list", length(sample.names))
names(runstats) <- sample.names
for(run in stats) {
  cat("Reading in:", run, "\n")
  runstats[[run]] <- read.delim(run, )
}

unqs <- unique(c(sapply(runs, colnames), recursive=TRUE))
n<-sum(unlist(lapply(X=runs, FUN = nrow)))
st <- matrix(0L, nrow=n, ncol=length(unqs))
rownames(st) <- c(sapply(runs, rownames), recursive=TRUE)
colnames(st) <- unqs
for(sti in runs) {
  st[rownames(sti), colnames(sti)] <- sti
}
st <- st[,order(colSums(st), decreasing=TRUE)]

##st.all<-mergeSequenceTables(runs)
# Remove chimeras
seqtab <- removeBimeraDenovo(st, method="consensus", multithread=TRUE)
# Assign taxonomy
homd <- assignTaxonomy(seqtab, "/home/jholm/bin/HOMD_v15.1_DADA2_taxonomy_final.txt", multithread=TRUE)
# Write to disk
saveRDS(seqtab, "all_runs_dada2_abundance_table.rds") # CHANGE ME to where you want sequence table saved
write.csv(seqtab, "all_runs_dada2_abundance_table.csv", quote=FALSE)
write.csv(homd, "homd_classification.csv", quote=FALSE)

fc = file("all_runs_dada2_ASV.fasta")
fltp = character()
for( i in 1:ncol(seqtab))
{
  fltp <- append(fltp, paste0(">", colnames(seqtab)[i]))
  fltp <- append(fltp, colnames(seqtab)[i])
}
writeLines(fltp, fc)
rm(fltp)
close(fc)

track<-as.matrix(rowSums(seqtab))
colnames(track) <- c("nonchimeric")
write.table(track, "dada2_part2_stats.txt", quote=FALSE, append=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
 ~;
    run_R_script($Rscript);
}

sub dada2_combine_and_classifyITS {
    my ($inRuns) = shift;

    my $Rscript = qq~

library("dada2")
packageVersion("dada2")
path<-getwd()

## list all of the files matching the pattern
tables<-list.files(path, pattern="-dada2_abundance_table.rds", full.names=TRUE)
stats<-list.files(path, pattern="-dada2_part1_stats.txt", full.names=TRUE)

## get the run names using splitstring on the tables where - exists
sample.names <- sapply(strsplit(basename(tables), "-"), `[`, 1)
##sample.names
##names(tables) <- sample.names

runs <- vector("list", length(sample.names))
names(runs) <- sample.names
for(run in tables) {
  cat("Reading in:", run, "\n")
  runs[[run]] <- readRDS(run)
}

runstats <- vector("list", length(sample.names))
names(runstats) <- sample.names
for(run in stats) {
  cat("Reading in:", run, "\n")
  runstats[[run]] <- read.delim(run, )
}

unqs <- unique(c(sapply(runs, colnames), recursive=TRUE))
n<-sum(unlist(lapply(X=runs, FUN = nrow)))
st <- matrix(0L, nrow=n, ncol=length(unqs))
rownames(st) <- c(sapply(runs, rownames), recursive=TRUE)
colnames(st) <- unqs
for(sti in runs) {
  st[rownames(sti), colnames(sti)] <- sti
}
st <- st[,order(colSums(st), decreasing=TRUE)]

##st.all<-mergeSequenceTables(runs)
# Remove chimeras
seqtab <- removeBimeraDenovo(st, method="consensus", multithread=TRUE)# THIS is the source of ASV names and sequences
# Assign taxonomy (requires colnames of seqtab to be ASV sequences)
unite <- assignTaxonomy(seqtab, "/home/jholm/bin/sh_general_release_dynamic_01.12.2017.fasta", multithread=TRUE)
silva <- assignTaxonomy(seqtab, "$pipelineDir/taxonomy/silva_nr_v138_train_set.fa.gz", multithread=TRUE)

# Name ASVs (Due to the constraints of assignTaxonomy, this step must occur after it)

# Write to disk
saveRDS(seqtab, "all_runs_dada2_abundance_table.rds") # CHANGE ME to where you want sequence table saved
write.csv(seqtab, "all_runs_dada2_abundance_table.csv", quote=FALSE)
write.csv(unite, "unite_classification.csv", quote=FALSE)
write.csv(silva, "silva_classification.csv", quote=FALSE)

fc = file("all_runs_dada2_ASV.fasta")
fltp = character()
for( i in 1:ncol(seqtab))
{
  fltp <- append(fltp, paste0(">", colnames(seqtab)[i])) # This reference FASTA is actually created from colnames(seqtab)
  fltp <- append(fltp, colnames(seqtab)[i])
}
writeLines(fltp, fc)
rm(fltp)
close(fc)

track<-as.matrix(rowSums(seqtab))
colnames(track) <- c("nonchimeric")
write.table(track, "dada2_part2_stats.txt", quote=FALSE, append=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
 ~;
    run_R_script($Rscript);
}

sub run_R_script {

    my $Rscript = shift;

    my $outFile = "rTmp.R";
    open OUT, ">$outFile", or die "cannot write to $outFile: $!\n";
    print OUT "$Rscript";
    close OUT;

    #local $ENV{LD_LIBRARY_PATH} = "/usr/local/packages/gcc/lib64";
    #system($cmd) == 0 or die "system($cmd) failed:$?\n";
    my $cmd = "$R CMD BATCH $outFile";
    system($cmd) == 0 or die "system($cmd) failed:$?\n";

    my $outR = $outFile . "out";
    open IN, "$outR" or die "Cannot open $outR for reading: $OS_ERROR\n";
    my $exitStatus = 1;

    foreach my $line (<IN>) {
        if ( $line =~ /learnErrors/ ) {
            next;
        } elsif ( $line =~ /Error/ ) {
            print "R script crashed at\n$line";
            print "check $outR for details\n";
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
sub execute_and_log {
    my $dryRun = pop @_;
    my $logFH  = pop @_;

    # CAREFUL, $cmd holds a reference to a variable in the caller!
    foreach my $cmd (@_) {
        print "\t$cmd\n" if $debug;
        if ($logFH) {
            print $logFH "\t$cmd\n";
        }
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;
    }
}
exit 0;
