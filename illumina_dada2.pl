#!/usr/local/packages/perl-5.30.2/bin/perl
#. /usr/local/packages/usepackage/share/usepackage/use.bsh

=head1 NAME

illumina_dada2.pl

=head1 SYNOPSIS

illumina_dada2.pl (-i <input directory> | -r1 <fwd reads> -r2 <rev reads> [-i1 <index 1> -i2 <index 2>]) --wd <directory> -m <map> -v <variable-region> [<options>]

=head1 DESCRIPTION

This is the MSL 16S pipeline for ILLUMINA runs. It can be launched from any 
location on the IGS server. When running the pipeline in full, the
essential inputs are a set of raw Illumina sequencing files (usually R1, R2, I1, 
and I2), a QIIME-formatted mapping file, the targeted variable region, and the
path to a working directory where the outputs will be placed.

The steps of the pipeline are as follows:

=over 

=item 1. Extracts barcodes from the raw files. Required inputs:

=over

=item Z<>* raw reads and indexes

=item Z<>* map

=back

=item 2. Demultiplexes the raw reads into fastq files containing reads specific to 
this project. Required inputs:

=over

=item Z<>* ./barcodes.fastq

=item Z<>* map

=back

=item 3. Produces individual .fastq files for each sample listed in the mapping file.
Required inputs:

=over

=item Z<>* ./fwdSplit/seqs.fastq

=item Z<>* ./revSplit/seqs.fastq

=back

=item 4. Performs tag-cleaning of each sample-specific file. Required inputs:

=over

=item Z<>* ./fwdSplit/split_by_sample_out/<sample_id>_*.fastq

=item Z<>* ./revSplit/split_by_sample_out/<sample_id>_*.fastq

=back

=item 5. Runs the forward and reverse reads through the dada2 pipeline for the 
V3V4 16S rRNA gene region. Alternatively, analysis of the V4 or ITS region may
be specified.

=over

=item Z<>* ./<sample_id>_*R1_tc.fastq

=item Z<>* ./<sample_id>_*R2_tc.fastq

=back

=back

Most pipeline products are stored in a directory named after the run. By default,
run directories are stored in /local/projects-t3/MSL/runs/. A log file is 
written at ./<RUN>_16S_pipeline_log.txt. Barcodes, raw libraries, and QC'd
reads are deleted if the DADA2 count table and stats exist when the pipeline
terminates.

=head1 OPTIONS

=head2 GENERAL

=over

=item B<--working-dir>=PATH, B<-wd> PATH

Indicate an existing directory in which to place output, and from which the
run name will be parsed. The last directory on PATH must be named after the run.
(Many of the result files will be named after the run.)

=item B<--1Step>

Use this flag if the data are prepared by 1-Step PCR (only forward & reverse read files
available). This processes the input files correctly and activates appropriate
parameters during adapter trimming, quality trimming/filtering, and denoising.

=item B<-h>, B<--help>

Print help message and exit successfully.

=item B<--qsub-project>=space, B<-qp> space

Indicate which qsub-project space should be used for all qsubmissions. The
default is jravel-lab.

=item B<--debug> {barcodes, demux, splitsamples, tagclean, dada2}

Runs the specified section of the pipeline. Multiple --debug options can be given
to run multiple consecutive parts of the pipeline, provided that the input to
the earliest requested step is present. Any non-consecutive steps will be 
ignored.

Usually helpful to add B<--noskip>.

=item B<--noskip>

Will not check for any possible existing output files when deciding whether to 
run a section of the pipeline.

=item B<--verbose>

Prints each command to STDOUT and to the log file.

=item B<--nodelete|--no-delete>

Don't delete intermediate files

=back

=head2 INPUT

=over

=item B<--raw-path>=path, B<-i> path

Single full path to directory containing raw files. File names must have the 
pattern *AA.fastq[.gz], where AA is either R1, R2, R3, R4, I1, or I2, and gzip 
compression is optional. B<--raw-path> is incompatible with B<-r1>, B<-r2>, B<-i1>, and B<-i2>.

Three combinations of files are allowed:

    DEFAULT (TWO-STEP PCR):
    R1, R2, I1, I2

    OLD TWO-STEP NAMING STYLE (automatically detected):
    R1, R2, R3, R4 
    (R1 and R4 are fwd and rev reads, respectively. R2 and R3 are index 1 
    and index 2, respectively.)

    ONE-STEP PCR (B<--1step> MUST BE GIVEN):
    R1, R2

=item B<-r1> file

Full path to raw foward read file (R1). Gzip compression optional.

=item B<-r2> file

Full path to raw reverse read file (R2, or R4 in old naming scheme). Gzip 
compression optional.

=item B<-i1> file

Full path to raw index 1 file (I1, or R2 in old naming scheme). Incompatible 
with B<--1step>. Gzip compression optional.

=item B<-i2> file

Full path to raw index 2 file (I2, or R3 in old naming scheme). Incompatible 
with B<--1step>. Gzip compression optional.

=item B<--map>=file, B<-m> file

The full path to the Qiime-formatted mapping file.

=back

=head2 BARCODE EXTRACTION AND DEMULTIPLEXING

=over

=item B<--bclen> LENGTH

Manually specify the length of forward and reverse barcodes. This many bases is 
removed from the index 1 and index 2 of each read, concatenated, and used to 
demultiplex the reads according to the provided map.

In our current sequencing configuration, the indexes ARE exactly this length.
By default, the pipeline sets the barcode length equal to the length of the
first index.

=item B<--troubleshoot_barcodes>

Try all four permutations of reverse-complementing and switching the 
concatenation order of the indexes. Whichever of the four permutations yields a
successful demux is used. Note: Performing these transformations on the indexes 
may coincidentally yield a barcode that seems to be correct, even though the 
overall demux is incorrect. 

=back

=head2 TRIMMING, FILTERING, AND DENOISING

=over

=item B<--var-reg>={V3V4, V4, ITS}, B<-v> {V3V4, V4, ITS}

The targeted variable region.

=item B<--dada2>="options"

Overrides the default DADA2 parameters used at the MSL. The following options
are allowed:

 --dada2-truncLen-f, -for (defaults: V3V4: 225 | V4: 200 | ITS: 0)
 --dada2-truncLen-r, -rev (defaults: V3V4: 225 | V4: 200 | ITS: 0)
 --dada2-maxN (default: 0)
 --dada2-maxEE (defaults: V3V4: 2 | V4: 2 | ITS: 0)
 --dada2-truncQ (default: 2)
 --dada2-rmPhix (default: TRUE)
 --dada2-maxLen (default: Inf)
 --dada2-minLen (default: V3V4: 20 | V4: 20 | ITS: 50)
 --dada2-minQ (default: 0)

Please see https://rdrr.io/bioc/dada2/man/filterAndTrim.html for descriptions
of the parameters. The parameters should be given within double quotes as shown 
below:

part1.sh --dada2="--dada2-maxEE 5 --dada2-minQ 10" ...

=item B<--dada2-mem> memory

The amount of memory to request for the DADA2 qsub job. Use Sun Grid Engine
syntax: 10G, 500M, etc.

=back 

=cut

use 5.010;    # for //
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
use lib catdir($pipelineDir, "lib", "perl5");
use lib catdir($pipelineDir, "lib");

require Version;

use File::Find;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Cwd qw(abs_path getcwd);
use File::Temp qw/ tempfile /;
use POSIX;
require File::Copy;
require List::Util;
use JSON;
require IO::Tee;
use File::Path qw( remove_tree );
use File::Copy::Recursive qw(dircopy );
use Data::Dumper;
use Hash::Merge qw( merge );
use Path::Tiny qw(path);
use Schedule::SGE;
use Capture::Tiny qw(capture_stderr);
use PathFinder;
use File::Copy::Recursive qw(rcopy_glob);
use Data::Dumper;
$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################
my @dbg;
my $oneStep       = 0;
my $dada2mem      = "1G";
my $tbshtBarcodes = 0;
my $delete        = 1;
my $trySkip       = 1;
GetOptions(
           "raw-path|i=s"           => \my $raw_dir,
           "r1=s"                   => \my $r1,
           "r2=s"                   => \my $r2,
           "i1=s"                   => \my $i1,
           "i2=s"                   => \my $i2,
           "map|m=s"                => \my $map,
           "var-reg|v=s"            => \my $var,
           "help|h!"                => \my $help,
           "d|debug=s"              => \@dbg,
           "verbose!"               => \my $verbose,
           "dry-run!"               => \my $dryRun,
           "skip!"                  => \$trySkip,
           "bclen=i"                => \my $bcLen,
           "dada2-truncLen-f|for=i" => \my $truncLenL,
           "dada2-truncLen-r|rev=i" => \my $truncLenR,
           "dada2-maxN=s"           => \my $maxN,
           "dada2-maxEE=s"          => \my $maxEE,
           "dada2-truncQ=s"         => \my $truncQ,
           "dada2-rmPhix!"          => \my $phix,
           "dada2-maxLen=s"         => \my $maxLen,
           "dada2-minLen=s"         => \my $minLen,
           "dada2-minQ=s"           => \my $minQ,
           "dada2-mem:s"            => \$dada2mem,
           "1Step!"                 => \$oneStep,
           "working-dir|wd=s"       => \my $wd,
           "qsub-project|qp=s"      => \my $qproj,
           "troubleshoot_barcodes!" => \$tbshtBarcodes,
           "delete!"                => \$delete,
          )
##add option for final resting place of important data

  or pod2usage(verbose => 0, exitstatus => 1);

if ($help)
{
    pod2usage(verbose => 2, exitstatus => 0);
    exit 1;
}

####################################################################
## PARAMETER CHECKING PT 1
####################################################################

# Check existence/compatibility of required parameters
if (   !$raw_dir
    && !((!$oneStep && $r1 && $r2 && $i1 && $i2) || ($oneStep && $r1 && $r2)))
{
    my $parameters = $oneStep ? "-r1, -r2" : "-r1, -r2, -i1, -i2";
    die "\n\tPlease provide the location of the raw sequencing files "
      . "(single directory => -i)\n\t\tOR \n\tFull paths to each raw file => "
      . "$parameters)\n\n";
}
if ($raw_dir && ($r1 || $r2 || $i1 || $i2))
{
    die
      "\n\tInput directory (-i) and raw files (-r*|-i*) were both specified. Please "
      . "provide one or the other.\n\n";
}
if (!$map)
{
    die "\n\tPlease provide the path to the run mapping file (-m)\n\n";
}
if ($oneStep && ($i1 || $i2))
{
    die "\n\t--1Step is incompatible with -i1 and -i2.\n\n";
}
if (!$var)
{
    die "\n\tPlease indicate the targeted variable region (-v {V3V4, V4, ITS,"
      . "ompA})\n\n";
}
if (!$wd)
{
    die "\n***Please choose a working directory (-wd).";
}

if ($truncLenL && !$truncLenR)
{
    die "***\nPlease provide truncation lengths for forward and reverse "
      . "reads\n";
}

# Refine and validate all variables that refer to the filesystem
my (@paths) = (
               {path => \$wd,      name => "Working directory"},
               {path => \$raw_dir, name => "Raw file directory"},
               {path => \$r1,      name => "Raw forward reads file"},
               {path => \$r2,      name => "Raw reverse reads file"},
               {path => \$i1,      name => "Raw forward index file"},
               {path => \$i2,      name => "Raw reverse index file"},
               {path => \$map,     name => "Mapping file"}
              );
foreach (@paths)
{
    if (${$$_{path}})
    {    # If the variable is a non-empty string...
        my $copy = ${$$_{path}};
        $copy =~ s/\/$//;    # remove any trailing slash
        my $relPath = $copy;
        $copy = abs_path($relPath);    # get abs path
            # At the same time, abs_path returns undef if the path doesn't exist
            # so we can verify the existence of each file and directory
        if (!defined $copy)
        {
            die $$_{name}
              . " not found. Looked for "
              . $relPath . ".\n"
              . "Current working directory: "
              . getcwd() . "\n";
        }
        ${$$_{path}} =
          $copy;    # external variable referenced by path has now been edited.
    }
}

# After shaving off potential trailing slash,
# Split on the last path separator, then check if the last string is a word
$wd =~ s/\/$//;
if ((File::Spec->splitpath("$wd"))[2] !~ /[^\/]/)
{
    die "Working directory (-wd) must be a path to a directory";
}

# now we can start logging to the log file
my $GLOBAL_PATHS = PathFinder->new(wd => $wd);

my ($readsForInput, $readsRevInput, $index1Input, $index2Input);
my %localNames;
if (!@dbg || grep(/^barcodes$/, @dbg) || grep(/^demux$/, @dbg))
{
    # Find the input files, and find their putative locations if the index files
    # were to be decompressed to our working directory.
    # In one-step runs, $index1Input and $index2Input are the SAME files
    # pointed to by $readsForInput and $readsRevInput
    ($readsForInput, $readsRevInput, $index1Input, $index2Input) = $raw_dir
      ? find_raw_files(
                       dir          => $raw_dir,
                       pf           => $GLOBAL_PATHS,
                       pcr_one_step => $oneStep
                      )
      : $oneStep ? ($r1, $r2, $r1, $r2)
      :            ($r1, $r2, $i1, $i2);
    @localNames{("readsFor", "readsRev", "index1", "index2")} =
      convert_to_local_if_gz($wd, $readsForInput, $readsRevInput, $index1Input,
                             $index2Input);
}

####################################################################
## INITIALIZATION
###################################################################
$GLOBAL_PATHS->{'raw_dir'} = $raw_dir;
$GLOBAL_PATHS->{'r1'}      = $r1;
$GLOBAL_PATHS->{'r2'}      = $r2;
$GLOBAL_PATHS->{'i1'}      = $i1;
$GLOBAL_PATHS->{'i2'}      = $i2;
$GLOBAL_PATHS->{'map'}     = $map;

my $config_hashref = config();

chdir $GLOBAL_PATHS->{'wd'};

# if we're in the global run directory, share all files with group
my $global_run_storage = abs_path($config_hashref->{"run_storage_path"});
if ($GLOBAL_PATHS->{'wd'} =~ /^${global_run_storage}/)
{
    umask(0002);
}

my $run  = basename($GLOBAL_PATHS->{'wd'});
my $time = strftime("%Y-%m-%d %H:%M:%S", localtime(time));

local $SIG{__WARN__} = sub {
    $GLOBAL_PATHS->print("WARNING: $_[0]");
    print STDERR "WARNING: $_[0]";
};    # Warnings go to log file, stdout, and stderr

# Now that we're done checking parameters, I want the more serious error
# messages to be printed to the log file.
$SIG{__DIE__} = sub {

    # if ($^S) {
    #     die @_;    # If die was called inside an eval, simply pass it along
    # } else {
    $GLOBAL_PATHS->print(@_ . "\n");
    return 1;

    # }
};

print "Logging to: " . $GLOBAL_PATHS->part1_log() . "\n";
$GLOBAL_PATHS->print("PIPELINE VERSION: " . Version::version() . "\n");
$GLOBAL_PATHS->print("$time\n");

my $map_log = catfile($GLOBAL_PATHS->part1_error_log(),
                      basename($GLOBAL_PATHS->{'map'}, ".txt") . ".log");

my @split;

my @errors;

my %number2step = (
                   0 => "map",
                   1 => "meta",
                   2 => "raw",
                   3 => "barcodes",
                   4 => "library",
                   5 => "samples",
                   6 => "tagcleaned",
                   7 => "dada2part1out"
                  );
my %step2number = (
                   "map"           => 0,
                   "meta"          => 1,
                   "raw"           => 2,
                   "barcodes"      => 3,
                   "library"       => 4,
                   "samples"       => 5,
                   "tagcleaned"    => 6,
                   "dada2part1out" => 7
                  );

my $models;
if ($var eq 'V3V4')
{
    $models = "/local/projects-t2/jholm/PECAN/v1.0/V3V4/merged_models/";
} else
{
    $models = "PECAN not used.\n";
}

if (!$qproj)
{
    print
      "\nqsub-project ID (--qsub-project, -qp) not provided. Using jravel-lab "
      . " as default\n";
    $qproj = "jravel-lab";
}

if (!-e $GLOBAL_PATHS->part1_error_log())
{
    mkdir $GLOBAL_PATHS->part1_error_log();
    print $GLOBAL_PATHS->part1_error_log()
      . " did not exist -> making "
      . $GLOBAL_PATHS->part1_error_log() . "\n";
} else
{
    # Remove old error logs
    my @oldLogs = glob($GLOBAL_PATHS->part1_error_log() . "/*.*");
    my @cmds;
    foreach my $log (@oldLogs)
    {
        my @dontDelete = (
                 $GLOBAL_PATHS->part1_error_log() . "/illumina_dada2.pl.stderr",
                 $GLOBAL_PATHS->part1_local_map(), $map_log
        );
        if (
            List::Util::none {$_ eq $log}
            @dontDelete
           )
        {
            push @cmds, "rm $log";
        }
    }
    if (@cmds)
    {
        execute_and_log(
                        cmds    => \@cmds,
                        dry_run => $dryRun,
                        msg     => "Removing stale error logs."
                       );
    }
}

if (!-e $GLOBAL_PATHS->part1_stdout_log())
{
    mkdir $GLOBAL_PATHS->part1_stdout_log();
    print $GLOBAL_PATHS->part1_stdout_log()
      . " did not exist -> making "
      . $GLOBAL_PATHS->part1_stdout_log() . "\n";
}

my $metadata = project_metadata();
$metadata =
  project_metadata(metadata => {"params" => {"platform" => "ILLUMINA"}});

my $qiime = "$GLOBAL_PATHS->{'wd'}/${run}_qiime_config.txt";

if ($tbshtBarcodes && @dbg)
{
    $GLOBAL_PATHS->print(
        "Disabling --debug sections because --troubleshoot_barcodes was given.\n"
    );
    @dbg = ();
}

if (@dbg)
{
    $GLOBAL_PATHS->print("DBG FLAGS: ");
    for (@dbg)
    {
        $GLOBAL_PATHS->print("$_ ");
    }

    my $e  = grep(/^barcodes$/,     @dbg);
    my $de = grep(/^demux$/,        @dbg);
    my $s  = grep(/^splitsamples$/, @dbg);
    my $t  = grep(/^tagclean$/,     @dbg);
    my $da = grep(/^dada2$/,        @dbg);
    if ($e + $de + $s + $t + $da != scalar @dbg)
    {
        die "Illegal debug option. Legal debug options are "
          . "barcodes, demux, splitsamples, tagclean, and dada2.";
    }
    if ($trySkip)
    {
        $trySkip = 0;
        $GLOBAL_PATHS->print(
                "Disabling skipping because --debug sections were selected.\n");
    }

    $GLOBAL_PATHS->print("\n");
}

$GLOBAL_PATHS->print(  "RUN: $run\nVARIABLE REGION: $var\n"
                     . "R VERSION: "
                     . $config_hashref->{'R'}
                     . "\nPECAN MODELS: $models\n");
if ($oneStep)
{
    $GLOBAL_PATHS->print("PCR PREPARATION METHOD: 1-Step\n");
} else
{
    $GLOBAL_PATHS->print("PCR PREPARATION METHOD: 2-Step\n");
}

my $skipMe;

if (   !@dbg
    || grep(/^barcodes$/,     @dbg)
    || grep(/^demux$/,        @dbg)
    || grep(/^splitsamples$/, @dbg))
{
    ###### BEGIN CHECK OF QIIME CONFIGURATION ###########
    #####################################################
    $GLOBAL_PATHS->print("\n");
    my $cmd = qiime_cmd('print_qiime_config.py', $config_hashref) . " > $qiime";
    execute_and_log(
                    cmds    => [$cmd],
                    logger  => $GLOBAL_PATHS,
                    dry_run => $dryRun,
                    msg     => "Printing QIIME configuration details\n"
                   );
    $GLOBAL_PATHS->print("See: $qiime\n");

    if ($trySkip)
    {
        $skipMe = skippable(
                        inputs  => [$GLOBAL_PATHS->{'map'}],
                        outputs => [$GLOBAL_PATHS->part1_local_map(), $map_log],
                        checksums_in  => $metadata->{"checkpoints"}{"map"},
                        checksums_out => $metadata->{"checkpoints"}{"meta"},
                        step_name     => "map validation",
                        previous_step_skipped => 1,
                        pf                    => $GLOBAL_PATHS
        );
    }

    if (!$skipMe)
    {
        File::Copy::copy($GLOBAL_PATHS->{'map'},
                         $GLOBAL_PATHS->part1_local_map());
        preprocess_map($GLOBAL_PATHS->part1_local_map(), $run);

        ###### BEGIN VALIDATION OF MAPPING FILE ###########
        ################################################
        $GLOBAL_PATHS->print(
                    "MAPPING FILE: " . $GLOBAL_PATHS->part1_local_map() . "\n");

        my $cmd =
            qiime_cmd('validate_mapping_file.py', $config_hashref) . " -m "
          . $GLOBAL_PATHS->part1_local_map()
          . " -s -o "
          . $GLOBAL_PATHS->part1_error_log();
        execute_and_log(
                        cmds    => [$cmd],
                        logger  => $GLOBAL_PATHS,
                        dry_run => $dryRun,
                        msg => "Validating map from $GLOBAL_PATHS->{'map'}\n"
                       );
        my $mappingError = glob($GLOBAL_PATHS->part1_error_log() . "/*.log");
        if ($mappingError)
        {
            open MAPERROR, "<$mappingError"
              or die "Cannot open $mappingError for " . "reading: $OS_ERROR";
            $_ = <MAPERROR>;    # gets first line

            chomp;
            if ($_ =~ /No errors or warnings found in mapping file./)
            {
                $GLOBAL_PATHS->print("---Map passed validation.\n");

            } else
            {
                while (<MAPERROR>)
                {
                    if ($_ =~ m/^Errors -----------------------------$/)
                    {
                        my $nextLine = <MAPERROR>;
                        if ($nextLine =~
                            m/^Warnings ---------------------------$/)
                        {
                            $GLOBAL_PATHS->print(
                                       "---Warnings during map validation. See "
                                         . "$mappingError for details.\n");
                        } else
                        {
                            close MAPERROR;
                            die
                              "***Error in mapping file. See $mappingError for "
                              . "details. Exiting.\n";
                        }
                    }
                }
                $GLOBAL_PATHS->print("---Map passed validation.\n");
            }
            close MAPERROR;
        } else
        {
            die
              "validate_mapping_file.py terminated but did not signal success. Normally a success message is printed in its error file.";
        }
        File::Copy::move($GLOBAL_PATHS->part1_local_map() . "_corrected.txt",
                         $GLOBAL_PATHS->part1_local_map());
        project_metadata(
             metadata => {"map" => {"file" => $GLOBAL_PATHS->part1_local_map()}}
        );

        ###### BEGIN EVALUATION OF SAMPLES VIA MAPPING FILE ###########
        ###############################################################
        open MAP, "<" . $GLOBAL_PATHS->part1_local_map()
          or die "Cannot open "
          . $GLOBAL_PATHS->part1_local_map()
          . " for reading: $OS_ERROR";
        my $extctrl     = 0;
        my $pcrpos      = 0;
        my $pcrneg      = 0;
        my $projSamples = 0;
        my $linecount   = 0;
        my $null        = 0;
        while (<MAP>)
        {
            chomp;
            ## don't count header as sample; don't count any line if it doesn't
            if ($. > 1)
            {
                ## start with four tab-separated fields; the first three must contain non-whitespace chars
                if ($_ =~ /^(\S+\t){3}/)
                {
                    if (   $_ =~ "EXTNTC"
                        || $_ =~ "NTC.EXT"
                        || $_ =~ "NTCEXT"
                        || $_ =~ "EXT.NTC"
                        || $_ =~ "NTC"
                        || $_ =~ "EXTNEG")
                    {
                        $extctrl++;
                    } elsif (   $_ =~ "PCRPOS"
                             || $_ =~ "PCR.pos"
                             || $_ =~ "pCR.pos"
                             || $_ =~ "POSCTRL"
                             || $_ =~ "POS.CTRL"
                             || $_ =~ "POSCON"
                             || $_ =~ "posctr")
                    {
                        $pcrpos++;
                    } elsif (   $_ =~ "PCRNTC"
                             || $_ =~ "PCR.NEG"
                             || $_ =~ "PCR.NTC"
                             || $_ =~ "PCRNEG"
                             || $_ =~ "PCRNEGCTRL"
                             || $_ =~ "ntcctr")
                    {
                        $pcrneg++;
                    } elsif ($_ =~ /NULL/)
                    {
                        $null++;
                    } else
                    {
                        $projSamples++;
                    }
                } elsif ($_ =~ /\S/)
                {

                    # QIIME's validate_mapping_file seems to already check this:
                    die
                      "In mapping file the line $. does not have four tab-separated fields.";
                }
            }
        }
        my $nSamples = $extctrl + $pcrpos + $pcrneg + $null + $projSamples;
        close MAP;

        $GLOBAL_PATHS->print("NO. SAMPLES: $nSamples\nNO. NULLS: $null\n"
            . "NO. EXTRACTION NEGATIVE CONTROLS: $extctrl\n"
            . "NO. PCR POSITIVE CONTROLS: $pcrpos\nNO. PCR NEGATIVE CONTROLS: $pcrneg\n"
        );

        if (!@dbg)
        {
            cacheChecksums(
                           files     => [$GLOBAL_PATHS->{'map'}],
                           step_name => "map",
                           pf        => $GLOBAL_PATHS
                          );
            cacheChecksums(
                 files     => [$mappingError, $GLOBAL_PATHS->part1_local_map()],
                 step_name => "meta",
                 pf        => $GLOBAL_PATHS
                          );
        }
    } else
    {

     # It is possible the first run had a map error, but the user failed to fix.
     # skippable would return true, but the log file would still inform us of
     # the error.
        my $mappingError =
          glob(catfile($GLOBAL_PATHS->part1_error_log(), "*.log"));

        open MAPERROR, "<", $mappingError
          or die "Cannot open $mappingError for " . "reading: $OS_ERROR";
        $_ = <MAPERROR>;
        close MAPERROR;
        chomp;
        if (!$_ =~ /No errors or warnings found in mapping file./)
        {
            die "***Unresolved error in mapping file. See $mappingError for "
              . "details. Exiting.\n";
        }
        $GLOBAL_PATHS->print(
                       "Mapping file has been validated already. Moving on.\n");
    }

} else
{
    $GLOBAL_PATHS->print("NOT USING QIIME\n");
}

###### BEGIN BARCODES ##########
#######################################
my $newSamNo;
if ($tbshtBarcodes)
{
    # make three other directories
    my %dirs_params = (
          "$GLOBAL_PATHS->{'wd'}_rc_fwd"     => "--rev_comp_bc1",
          "$GLOBAL_PATHS->{'wd'}_rc_rev"     => "--rev_comp_bc2",
          "$GLOBAL_PATHS->{'wd'}_rc_fwd_rev" => "--rev_comp_bc1 --rev_comp_bc2",
          "$GLOBAL_PATHS->{'wd'}"            => ""
    );

    # initialize directories by copying all contents from cwd
    for my $dir (keys %dirs_params)
    {

        sub init
        {
            my $dir = shift;
            if ($dir ne $GLOBAL_PATHS->{'wd'})
            {
                remove_tree($dir) if (-e $dir);
                mkdir $dir or die "Unable to create $dir\n";
                rcopy_glob("$GLOBAL_PATHS->{'wd'}/*", $dir)
                  or die
                  "Could not perform dircopy of $GLOBAL_PATHS->{'wd'} to $dir: $!";
            }
        }
        init($dir);
        $dir = "${dir}_switched";
        init($dir);
    }

    my $success = 0;
    $GLOBAL_PATHS->print(
        "Attempting all permutations of reverse-complemented indexes, with and without switching of indexes."
          . "\n");
    my $summary = "";

    for my $dir (keys %dirs_params)
    {

        sub barcodes_demux
        {
            my %arg             = @_;
            my $dir             = delete $arg{wd};
            my $barcode_ori     = delete $arg{barcode_ori};
            my $switch_barcodes = delete $arg{switch_barcodes} // 0;

       # clone this and replace wd because demux can take place in any directory
            my $paths = $GLOBAL_PATHS->clone();
            $paths->set_wd($dir);
            $GLOBAL_PATHS->print("Now working in: $dir.\n");
            $GLOBAL_PATHS->print("Please refer to " . $paths->part1_log());

            # Extract barcodes
            barcodes(
                     oriParams  => $barcode_ori,
                     pathfinder => $paths,
                     switch     => $switch_barcodes
                    );

            # demux, and check for success
            $success = demux(pathfinder => $paths, die_on_fail => 0);
            if ($success)
            {
                my $msg =
                  "Demux may have succeeded when extract_barcodes.py called with "
                  . "command-line options \"${barcode_ori}\"";
                $msg .= " and switched indexes" if $switch_barcodes;
                $msg .= "\n";
                $GLOBAL_PATHS->print("\n$msg\n");
                $summary .= $msg;
            } else
            {
                my $msg = "Demux failed when extract_barcodes.py called with "
                  . "command-line options \"${barcode_ori}\"";
                $msg .= " and switched indexes" if $switch_barcodes;
                $msg .= "\n";
                $summary .= $msg;
            }
            $paths->close();
        }
        barcodes_demux(wd => $dir, barcode_ori => $dirs_params{$dir});
        barcodes_demux(
                       wd              => "${dir}_switched",
                       barcode_ori     => $dirs_params{$dir},
                       switch_barcodes => 1
                      );

    }

    my $msg =
        "\nFinished troubleshooting barcodes.\n"
      . "Please inspect the split library logs in each trial directory to determine which demux was correct.\n"
      . "Correct files can be moved back to this run directory, and the pipeline can be continued using '--debug splitsamples --debug tagclean --debug dada2'.\n\n";
    $GLOBAL_PATHS->print($msg);
    $GLOBAL_PATHS->print("SUMMARY:\n$summary\n");
    exit 0;

} else
{
    if ((!@dbg) || grep(/^barcodes$/, @dbg))
    {
        barcodes(oriParams => "", pathfinder => $GLOBAL_PATHS);
    }
    if (!@dbg || grep(/^demux$/, @dbg))
    {
        demux(pathfinder => $GLOBAL_PATHS, die_on_file => 1);
    }

    # Then continue with code below the two subroutines
}

sub barcodes
{
    my %arg       = @_;
    my $oriParams = delete $arg{oriParams} // '';
    my $switch    = delete $arg{switch} // 0;
    my $paths     = delete $arg{pathfinder} // $GLOBAL_PATHS;
    my $barcodes  = "$paths->{'wd'}/barcodes.fastq";

    if ($switch)
    {    # Consider the i1's and i2's switched
        my $oldi1 = $index1Input;
        $index1Input = $index2Input;
        $index2Input = $oldi1;
    }

    my $nSamples = count_samples($paths->part1_local_map());

    ## change to full path to full barcodes (flag)
    my $count = 0;

    my @files;
    my @names;
    if ($oneStep)
    {
        $paths->print(  "Forward and reverse reads will be obtained from:\n"
                      . "\t$readsForInput\n\t$readsRevInput\n");
        @files = ($readsForInput, $readsRevInput);
        @names = ("RF", "RR");
    } else
    {
        $paths->print(  "Forward and reverse reads will be obtained from:\n"
                      . "\t$readsForInput\n\t$readsRevInput\n");
        $paths->print(  "Index 1 and Index 2 will be obtained from:\n"
                      . "\t$index1Input\n\t$index2Input\n");
        @files = ($readsForInput, $readsRevInput, $index1Input, $index2Input);
        @names = ("RF", "RR", "I1", "I2");
    }

    my $input  = \@files;
    my $output = [$barcodes];
    if ($trySkip)
    {
        $skipMe = skippable(
                        inputs        => $input,
                        outputs       => $output,
                        checksums_in  => $metadata->{"checkpoints"}{"raw"},
                        checksums_out => $metadata->{"checkpoints"}{"barcodes"},
                        step_name     => "barcode extraction",
                        in_names      => \@names,
                        previous_step_skipped => 1,
                        pf                    => $GLOBAL_PATHS
        );
    }

    my $step1;

    if (!$skipMe)
    {

        my $start = time;
        $step1 = "extract_barcodes.py";

        # Get index files as .fastq
        my @cmds = ();
        if ($index1Input !~ /.gz$/)
        {

            # if the index files aren't .gz, just read in place.
        } else
        {    # Otherwise...
                # decompress to our run directory
            push(@cmds,
                 "gzip --decompress --force < $index1Input > $localNames{\"index1\"}"
                );
            $paths->print(
                "---Decompressing $run barcode and index file from $index1Input to $localNames{\"index1\"}\n"
            );
        }
        if ($index2Input !~ /.gz$/)
        {       # same for reverse reads
        } else
        {
            push(@cmds,
                 "gzip --decompress --force < $index2Input > $localNames{\"index2\"}"
                );
            $paths->print(
                "---Decompressing $run barcode and index file from $index2Input to $localNames{\"index2\"}\n"
            );
        }

        # Execute the commands queued above
        if (scalar @cmds > 0)
        {
            execute_and_log(
                cmds    => \@cmds,
                dry_run => $dryRun,
                logger  => $paths,
                msg =>
                  "Decompressing raw files containing indexes to run directory...\n",
                verbose => 1
            );
            $paths->print("---Raw index file(s) decompressed.\n");
            @cmds = ();
        }

        # Sanity check
        my $rawLnsFwd = count_lines($localNames{"index1"});
        my $rawLnsRev = count_lines($localNames{"index2"});
        if ($rawLnsFwd == $rawLnsRev)
        {
            $paths->print(($rawLnsFwd / 4 . " barcoded reads to process.\n"));
        } else
        {
            $paths->print(
                 "Warning: raw index files have differing numbers of lines.\n");
        }

        # Was in original code, but maybe unnecessary?
        my $mapOpt = $oneStep ? "-m " . $paths->part1_local_map() . "" : "";

        if (!$bcLen)
        {
            if ($oneStep)
            {
                # Idk if auto-detecting the barcode length works under one-step
                # PCR. We've never done any one-step PCR runs while I've been
                # here
                $bcLen = 12;
            } else
            {
                $bcLen = find_index_length($localNames{"index1"});
                $paths->print("Detected barcode length as $bcLen \n");
            }
        }

        my $cmd = qiime_cmd('extract_barcodes.py', $config_hashref)
          . " -f $localNames{\"index1\"} -r $localNames{\"index2\"} -c barcode_paired_end --bc1_len $bcLen --bc2_len $bcLen $mapOpt -o $paths->{'wd'} $oriParams";

        execute_and_log(
                       cmds    => [$cmd],
                       logger  => $paths,
                       dry_run => $dryRun,
                       msg => "Waiting for barcode extraction to complete...\n",
                       verbose => 1
        );

        # Sanity check
        my $barcodeLns = count_lines("$barcodes");
        if ($rawLnsFwd != $barcodeLns)
        {
            die
              "Forward index file had $rawLnsFwd lines but barcodes.fastq had $barcodeLns lines.";
        }

        my $duration = time - $start;
        $paths->print(  "---Barcode extraction complete... "
                      . ($barcodeLns / 4)
                      . " barcodes extracted\n");
        $paths->print("---Duration of barcode extraction: $duration s\n");

        # Record hash of barcodes.fastq, and inputs
        if (!@dbg)
        {
            cacheChecksums(
                           files      => \@files,
                           step_name  => "raw",
                           file_names => \@names,
                           pf         => $paths
                          );
            cacheChecksums(
                           files     => [$barcodes],
                           step_name => "barcodes",
                           pf        => $paths
                          );
        }

        # Delete unneeded output of barcodes.py
        @cmds = ();
        push @cmds, "rm -rf $paths->{'wd'}/reads1.fastq";
        push @cmds, "rm -rf $paths->{'wd'}/reads2.fastq";

        # If two-step, delete the temporarily decompressed index files now
        if (!$oneStep)
        {
            if ($index1Input =~ /.gz$/)
            {
                push @cmds, "rm -rf $localNames{\"index1\"}";
            }
            if ($index2Input =~ /.gz$/)
            {
                push @cmds, "rm -rf $localNames{\"index2\"}";
            }
        }
        execute_and_log(
                        cmds    => \@cmds,
                        dry_run => $dryRun,
                        msg     => "Cleaning up after extract_barcodes.py...\n",
                        logger  => $paths
                       );

        if (@dbg && !grep(/^demux$/, @dbg))
        {
            $paths->print(  "Finished extracting barcodes. Terminated "
                          . "because --debug demux was not specified.\n");
            exit 0;
        }
    } else
    {
        $paths->print(  "Barcodes already extracted as "
                      . basename($barcodes)
                      . ". Moving on.\n");
    }
}

###### BEGIN SPLIT LIBRARIES ##########
#######################################

sub demux
{
    my %arg   = @_;
    my $paths = delete $arg{pathfinder} // $GLOBAL_PATHS;
    $paths->print("Demuxing in: $paths->{'wd'}\n");

    my $die_on_fail = delete $arg{die_on_fail};

    my $barcodes = "$paths->{'wd'}/barcodes.fastq";
    my $nSamples = count_samples($paths->part1_local_map());
    my $step2;
    my $step3;

    my $split_log = catfile($paths->fwd_demux_dir(), "split_library_log.txt");

    my @cmds;

    # split_libraries_fastq.py accepts fastq.gz, so no need to
    # convert_to_local_if_gz... Unless it's one-step, in which case we can save
    # time by using the already-decompressed raw files.
    if ($oneStep)
    {
        ($readsForInput, $readsRevInput) =
          convert_to_local_if_gz($paths->{'wd'}, $readsForInput,
                                 $readsRevInput);
    }

    if ($trySkip)
    {
        $skipMe = skippable(
            inputs => [
                $readsForInput, $readsRevInput,

                # FIXME need to copy local map to other directory?
                $paths->part1_local_map(), $barcodes
                      ],
            outputs      => [$paths->fwd_library(), $paths->rev_library()],
            checksums_in => {
                            "RF"  => $metadata->{"checkpoints"}{"raw"}->{"RF"},
                            "RR"  => $metadata->{"checkpoints"}{"raw"}->{"RR"},
                            "map" => $metadata->{"checkpoints"}{"map"},
                            "barcodes" => $metadata->{"checkpoints"}{"barcodes"}
            },
            checksums_out         => $metadata->{"checkpoints"}{"library"},
            step_name             => "demultiplexing",
            in_names              => ["RF", "RR", "map", "barcodes"],
            previous_step_skipped => 1,
            pf                    => $paths
                           );
    }

    if (!$skipMe)
    {
        my $start  = time;
        my $step2  = "split_libraries_fastq.py";
        my $script = qiime_cmd($step2, $config_hashref);

   # qiime's split_libraries_fastq accepts some keywords too, such as "golay_12"
        my $barcodeType = 2 * $bcLen;

        @cmds = ();
        push @cmds,
            "$config_hashref->{'executor'} -N $step2 -P $qproj -e "
          . $paths->part1_error_log() . " -o "
          . $paths->part1_stdout_log()
          . " $script -i $readsForInput -o "
          . $paths->fwd_demux_dir()
          . " -b $barcodes -m "
          . $paths->part1_local_map()
          . " --max_barcode_errors 1 --store_demultiplexed_fastq --barcode_type $barcodeType -r 999 -n 999 -q 0 -p 0.0001";
        push @cmds,
            "$config_hashref->{'executor'} -N $step2 -P $qproj -e "
          . $paths->part1_error_log() . " -o "
          . $paths->part1_stdout_log()
          . " $script -i $readsRevInput -o "
          . $paths->rev_demux_dir()
          . " -b $barcodes -m "
          . $paths->part1_local_map()
          . " --max_barcode_errors 1 --store_demultiplexed_fastq --barcode_type $barcodeType -r 999 -n 999 -q 0 -p 0.0001";

        execute_and_log(
                        cmds    => \@cmds,
                        logger  => $paths,
                        dry_run => $dryRun,
                        msg     => "Demultiplexing to get the run library...\n"
                       );

        $paths->print(
                     "---Waiting for fwd and rev seqs.fastq to complete....\n");
        $paths->print("---Monitoring $step2 error logs....\n");
        $paths->check_error_log(prefix => $step2);
        while (!(-e $paths->fwd_library()) || !(-e $paths->rev_library()))
        {
            sleep 1;

            # check_error_log allows pipeline to terminate if the qsubbed
            # split_libraries_fastq.py prints error
            $paths->check_error_log(prefix => $step2);
        }
        my $duration = time - $start;
        $paths->print(
             "---Duration of fwd and rev seqs.fastq production: $duration s\n");

        ##Check split_library log for 0's
        my @split =
          readSplitLog(file => $split_log, verbose => 1, logger => $paths);

        $newSamNo = $nSamples - scalar @split;
        if (scalar @split > 0 && scalar @split ne $nSamples)
        {
            $paths->print("---The following "
                . scalar @split
                . " samples returned 0 "
                . "reads, due to either a lack of barcodes or a filtering of those "
                . "reads in split_libraries_fastq.py:\n");
            foreach my $x (@split)
            {
                $paths->print("   $x\n");
            }
            $paths->print(
                         "---Number of samples after split_libraries_fastq.py: "
                           . "$newSamNo\n");
        } elsif (scalar @split == $nSamples)
        {
            if ($die_on_fail)
            {
                $paths->print(
                        "---No samples were successfully demultiplexed. Is the "
                          . "mapping file correct? Exiting.\n");
                die;
            } else
            {
                return 0;
            }
        } elsif (scalar @split == 0)
        {
            $paths->print("---Reads from all samples were demultiplexed.\n");
        }

        @cmds = ();
        push @cmds,
            catfile($scriptsDir, "get_split_library_stats.sh") . " < "
          . catfile($paths->fwd_demux_dir(), "split_library_log.txt") . " > "
          . catfile($paths->fwd_demux_dir(), "split_library_stats.txt");
        push @cmds,
            catfile($scriptsDir, "get_split_library_stats.sh") . " < "
          . $paths->rev_demux_dir()
          . "/split_library_log.txt > "
          . $paths->rev_demux_dir()
          . "/split_library_stats.txt";
        execute_and_log(
                        cmds    => \@cmds,
                        logger  => $paths,
                        dry_run => $dryRun,
                        msg     => "Extracting demultiplex stats...\n"
                       );

        if ($dada2mem eq "1G")
        {    # if it's the default
            open(my $sll, "<",
                 ${paths}->fwd_demux_dir() . "/split_library_log.txt")
              or die "Can't open fwd split library log: $!";
            while (my $line = <$sll>)
            {
                if ($line =~ m/^Total number seqs written\s+(\d+)/)
                {
                    $dada2mem =
                      List::Util::max(POSIX::ceil(4.2 * log($1) - 4.2) / 10, 1);
                    $dada2mem = "${dada2mem}G";
                }
            }
        }

        ###### BEGIN FASTQC ON SEQS.FASTQ #####
        #######################################
        # Replace this with calls to execute_and_log after merging in master
        my $binary = "/local/projects-t3/MSL/pipelines/bin/fastqc";
        @cmds = ();
        push @cmds,
            "$config_hashref->{'executor'} -l mem_free=300M -P $qproj -e "
          . $paths->part1_error_log() . " -o "
          . $paths->part1_stdout_log()
          . " $binary --limits $pipelineDir/ext/fastqc/limits.txt --outdir "
          . $paths->fwd_demux_dir() . " "
          . $paths->fwd_library();
        push @cmds,
            "$config_hashref->{'executor'} -l mem_free=300M -P $qproj -e "
          . $paths->part1_error_log() . " -o "
          . $paths->part1_stdout_log()
          . " $binary --limits $pipelineDir/ext/fastqc/limits.txt --outdir "
          . $paths->rev_demux_dir() . " "
          . $paths->rev_library();
        execute_and_log(
                        cmds    => \@cmds,
                        logger  => $paths,
                        dry_run => $dryRun,
                        msg     => "Assessing library quality...\n"
                       );

        # Remove temporarily decompressed files
        if ($oneStep)
        {
            my @cmds = ();
            push(@cmds, "rm -rf $readsForInput");
            push(@cmds, "rm -rf $readsRevInput");
            execute_and_log(
                cmds    => \@cmds,
                dry_run => $dryRun,
                logger  => $paths,
                msg => "---Removing decompressed raw files from run directory\n"
            );
        }

        if (!@dbg)
        {
            cacheChecksums(
                    files     => [$paths->fwd_library(), $paths->rev_library()],
                    step_name => "library",
                    pf        => $paths
                          );
        }
        project_metadata(
                  metadata => {"map" => {"file" => $paths->part1_local_map()}});

    } else
    {
        $paths->print(  "Library already produced as "
                      . basename($paths->fwd_library()) . " and "
                      . basename($paths->rev_library())
                      . ". Moving on.\n");
        my @split =
          readSplitLog(file => $split_log, verbose => 0, logger => $paths);

# Still need number of demuxed samples for later, even though we didn't actually demux this time
        $newSamNo = $nSamples - scalar @split;
    }

    if (@dbg && !grep(/^splitsamples$/, @dbg))
    {
        die "Finished demultiplexing libaries. Terminated "
          . "because -d splitsamples was not specified.";
    }
    return 1;
}

if (!@dbg || grep(/^splitsamples$/, @dbg))
{
    ###### BEGIN SPLIT BY SAMPLE ##########
    #######################################
    my $n_fq1 = 0;

    # need to look for both fastqs and fastq.gz's because after tagcleaning
    # all raw demuxed files are gzipped. Look for either
    my @fwdOutput =
      $GLOBAL_PATHS->splitsamples_output_fwd(exts => ["fastq", "fastq.gz"]);
    my @revOutput =
      $GLOBAL_PATHS->splitsamples_output_rev(exts => ["fastq", "fastq.gz"]);

    if ($trySkip)
    {
        $skipMe = skippable(
                         inputs => [
                                    File::Spec->abs2rel(
                                                   $GLOBAL_PATHS->fwd_library(),
                                                   $GLOBAL_PATHS->{'wd'}
                                    ),
                                    File::Spec->abs2rel(
                                                   $GLOBAL_PATHS->rev_library(),
                                                   $GLOBAL_PATHS->{'wd'}
                                    )
                                   ],
                         outputs       => [@fwdOutput, @revOutput],
                         checksums_in  => $metadata->{"checkpoints"}{"library"},
                         checksums_out => $metadata->{"checkpoints"}{"samples"},
                         step_name     => "sample splitting",
                         previous_step_skipped => 1,
                         pf                    => $GLOBAL_PATHS
        );
    }

    if (!$skipMe)
    {
        my $step3  = "split_sequence_file_on_sample_ids.py";
        my $script = qiime_cmd($step3, $config_hashref);
        my @cmds   = ();
        while (   !(-e $GLOBAL_PATHS->fwd_library())
               || !(-e $GLOBAL_PATHS->rev_library()))
        {
            sleep 1;
        }

        push @cmds,
          "$config_hashref->{'executor'} -N $step3 -l mem_free=5G -P $qproj -e "
          . $GLOBAL_PATHS->part1_error_log() . " -o "
          . $GLOBAL_PATHS->part1_stdout_log()
          . " $script -i "
          . $GLOBAL_PATHS->fwd_library()
          . " --file_type fastq -o "
          . $GLOBAL_PATHS->fwd_sample_dir();
        push @cmds,
          "$config_hashref->{'executor'} -N $step3 -l mem_free=5G -P $qproj -e "
          . $GLOBAL_PATHS->part1_error_log() . " -o "
          . $GLOBAL_PATHS->part1_stdout_log()
          . " $script -i "
          . $GLOBAL_PATHS->rev_library()
          . " --file_type fastq -o "
          . $GLOBAL_PATHS->rev_sample_dir();
        execute_and_log(
                    cmds    => \@cmds,
                    logger  => $GLOBAL_PATHS,
                    dry_run => $dryRun,
                    msg => "Splitting $run seqs.fastq " . "files by sample ID\n"
        );

        ## the $nSamples needs to be altered if a sample has 0 reads, because the sample-specific fastq won't be produced
        my $n_fq   = 0;
        my $nLines = 0;
        my @fwdOutput;
        my @revOutput;

        # Count the number of lines in fwdSplit/seqs.fastq
        open my $FWD, "<" . $GLOBAL_PATHS->fwd_library();
        while (<$FWD>) { }
        my $fwdLines = $.;
        close $FWD;

# if we know how many new files there are supposed to be, wait for them all to appear
        if (defined $newSamNo)
        {

            # Count the number of files in the directory
            while ($n_fq != $newSamNo)
            {
                $n_fq = 0;
                @fwdOutput =
                  $GLOBAL_PATHS->splitsamples_output_fwd(exts => ["fastq"]);
                $n_fq = scalar @fwdOutput;
                $GLOBAL_PATHS->check_error_log(prefix => $step3);
            }
        }
        print "---All samples ($n_fq) accounted for in "
          . $GLOBAL_PATHS->fwd_sample_dir() . "\n";

# Now check that all the reads are still contained in the sample-specific FASTQ's
        while ($nLines != $fwdLines)
        {
            $nLines = 0;
            @fwdOutput =
              $GLOBAL_PATHS->splitsamples_output_fwd(exts => ["fastq"]);
            if (@fwdOutput)
            {
                foreach my $file (@fwdOutput)
                {

                    # Also keep a running total of lines in the
                    # split_by_sample_out files
                    open SAMPLE, "<$file";
                    while (<SAMPLE>) { }
                    $nLines += $.;
                    close SAMPLE;
                }
            }
            $GLOBAL_PATHS->check_error_log(prefix => $step3);
        }

        print "---All reads (@{[$nLines / 4]}) accounted for in "
          . $GLOBAL_PATHS->fwd_sample_dir() . "\n";

        $n_fq   = 0;
        $nLines = 0;

        # Count the number of lines in R4split/seqs.fastq
        open my $REV, "<" . $GLOBAL_PATHS->rev_library();
        while (<$REV>) { }
        my $revLines = $.;
        close $REV;

# if we know how many new files there are supposed to be, wait for them all to appear
        if (defined $newSamNo)
        {
            $n_fq++;    # Count the number of files in the directory
            while ($n_fq != $newSamNo)
            {
                $n_fq = 0;
                @revOutput =
                  $GLOBAL_PATHS->splitsamples_output_rev(exts => ["fastq"]);
                $n_fq = scalar @revOutput;
                $GLOBAL_PATHS->check_error_log(prefix => $step3);
            }
        }

# Now check that all the reads are still contained in the sample-specific FASTQ's
        while ($nLines != $revLines)
        {
            $nLines = 0;
            @revOutput =
              $GLOBAL_PATHS->splitsamples_output_rev(exts => ["fastq"]);
            if (@revOutput)
            {
                foreach my $file (@revOutput)
                {

                    # Also keep a running total of lines in the
                    # split_by_sample_out files
                    open SAMPLE, "<$file";
                    while (<SAMPLE>) { }
                    $nLines += $.;
                    close SAMPLE;
                }
            }
            $GLOBAL_PATHS->check_error_log(prefix => $step3);
        }

        # Append _R1 or _R2 to each filename
        {
            $GLOBAL_PATHS->print(
                      "Appending _R1 or _R2 to each demuxed FASTQ filename.\n");
            foreach my $oldname (@fwdOutput)
            {
                my ($name, $path, $suffix) = fileparse($oldname, ".fastq");
                my $newname = catfile($path, "${name}_R1.fastq");
                File::Copy::move($oldname, $newname);
                $oldname = $newname;
            }
            foreach my $oldname (@revOutput)
            {
                my ($name, $path, $suffix) = fileparse($oldname, ".fastq");
                my $newname = catfile($path, "${name}_R2.fastq");
                File::Copy::move($oldname, $newname);
                $oldname = $newname;
            }
        }

        if (!@dbg)
        {
            cacheChecksums(
                           files     => [@fwdOutput, @revOutput],
                           step_name => "samples",
                           pf        => $GLOBAL_PATHS
                          );
        }
        print "--All samples ($n_fq) and reads (@{[$nLines / 4]}) accounted for"
          . " in "
          . $GLOBAL_PATHS->rev_sample_dir() . "\n";
    } else
    {
        $GLOBAL_PATHS->print(
            "Library already split to $newSamNo sample-specific files. Moving on.\n"
        );
    }

    # Remove temporarily decompressed files
    if ($oneStep)
    {
        my ($readsForInput, $readsRevInput, $index1Input, $index2Input) =
          $GLOBAL_PATHS->{'raw_dir'}
          ? find_raw_files(
                           pf           => $GLOBAL_PATHS,
                           dir          => $GLOBAL_PATHS->{'raw_dir'},
                           pcr_one_step => $oneStep
                          )
          : $oneStep ? (
                        $GLOBAL_PATHS->{'r1'}, $GLOBAL_PATHS->{'r2'},
                        $GLOBAL_PATHS->{'r1'}, $GLOBAL_PATHS->{'r2'}
                       )
          : (
             $GLOBAL_PATHS->{'r1'}, $GLOBAL_PATHS->{'r2'},
             $GLOBAL_PATHS->{'i1'}, $GLOBAL_PATHS->{'i2'}
            );

        my @cmds = ();
        push(@cmds, "rm -rf $readsForInput");
        push(@cmds, "rm -rf $readsRevInput");
        execute_and_log(
             cmds    => \@cmds,
             dry_run => $dryRun,
             msg =>
               "---Removing decompressed raw files from $GLOBAL_PATHS->{'wd'}\n"
        );
    }

    if (@dbg && !grep(/^tagclean$/, @dbg))
    {
        $GLOBAL_PATHS->print(
             "Finished splitting library to sample-specific FASTQs. Terminated "
               . "because --debug tagclean was not specified.\n");
        exit 0;
    }
}

###### BEGIN TAGCLEANING ##########
###################################

my $start = time;
if (!@dbg || grep(/^tagclean$/, @dbg))
{
    my $do = 1;

    my @inputsF =
      $GLOBAL_PATHS->splitsamples_output_fwd(exts => ["fastq", "fastq.gz"]);
    my @inputsR =
      $GLOBAL_PATHS->splitsamples_output_rev(exts => ["fastq", "fastq.gz"]);

    if (scalar @inputsF != scalar @inputsR)
    {
        die "Unequal numbers of forward and reverse sample FASTQ's.\n";
    }
    my $newSamNo = scalar @inputsF;

    my @fwdTcFiles = glob("$GLOBAL_PATHS->{'wd'}/*R1_tc.fastq.gz");
    my @revTcFiles = glob("$GLOBAL_PATHS->{'wd'}/*R2_tc.fastq.gz");

    if ($trySkip)
    {
        $skipMe = skippable(
              inputs                => [(@inputsF, @inputsR)],
              outputs               => [@fwdTcFiles, @revTcFiles],
              checksums_in          => $metadata->{"checkpoints"}{"samples"},
              checksums_out         => $metadata->{"checkpoints"}{"tagcleaned"},
              step_name             => "primer removal",
              previous_step_skipped => 1,
              pf                    => $GLOBAL_PATHS
                           );
    }

    if (!$skipMe)
    {
        my @cmds;

        $var = uc $var;
        if (grep(/^${var}$/, ("V3V4", "V4", "ITS")))
        {
            my $fwd_adapt;
            my $rev_adapt;
            my $base_cmd =
                "$config_hashref->{'executor'} -l mem_free=400M -P $qproj -e "
              . $GLOBAL_PATHS->part1_error_log() . " -o "
              . $GLOBAL_PATHS->part1_stdout_log()
              . " -N tagcleaner.pl 'perl "
              . $config_hashref->{'part1'}->{"tagcleaner"};

            if ($oneStep)
            {
                $base_cmd .= " -trim_within 50";
                if ($var eq "V3V4")
                {
                    $fwd_adapt = "GGACTACHVGGGTWTCTAAT";
                    $rev_adapt = "ACTCCTACGGGAGGCAGCAG";
                }
                if ($var eq "V4")
                {
                    $fwd_adapt = "GTGCCAGCMGCCGCGGTAA";
                    $rev_adapt = "ACTCCTACGGGAGGCAGCAG";
                }
            } else
            {
                if ($var eq "V3V4")
                {
                    $fwd_adapt = "ACTCCTACGGGAGGCAGCAG";
                    $rev_adapt = "GGACTACHVGGGTWTCTAAT";
                }
                if ($var eq "V4")
                {
                    $fwd_adapt = "GTGCCAGCMGCCGCGGTAA";
                    $rev_adapt = "GGACTACHVGGGTWTCTAAT";
                }
                if ($var eq "ITS")
                {
                    $fwd_adapt = "CTGCCCTTTGTACACACCGC";
                    $rev_adapt = "TTTCGCTGCGTTCTTCATCG";
                }
            }

            my $filename;
            foreach my $filename (@inputsF)
            {
                my $cmd = $base_cmd;
                my ($name, $dir, $ext) = fileparse($filename, qr{\.gz});
                if ($ext)
                {
                    $cmd =~ s/perl/gunzip $filename; perl/;
                }
                my @suffixes = (".fastq");
                my $Prefix   = basename($name, @suffixes);
                my $tc       = "$GLOBAL_PATHS->{'wd'}/$Prefix" . "_tc";
                $cmd =
                    $cmd
                  . " -fastq "
                  . catfile($dir, $name)
                  . " -out $tc -line_width 0 "
                  . "-verbose -tag5 $fwd_adapt -mm5 2 -nomatch 1'";
                push @cmds, $cmd;
            }

            foreach my $filename (@inputsR)
            {
                my $cmd = $base_cmd;
                my ($name, $dir, $ext) = fileparse($filename, qr{\.gz});
                if ($ext)
                {
                    $cmd =~ s/perl/gunzip $filename; perl/;
                }
                my @suffixes = (".fastq");
                my $Prefix   = basename($name, @suffixes);
                my $tc       = "$GLOBAL_PATHS->{'wd'}/$Prefix" . "_tc";
                $cmd =
                    $cmd
                  . " -fastq "
                  . catfile($dir, $name)
                  . " -out $tc -line_width 0 "
                  . "-verbose -tag5 $rev_adapt -mm5 2 -nomatch 1'";
                push @cmds, $cmd;
            }
        } elsif ($var eq "OMPA")
        {
            my $bbduk = $config_hashref->{"part1"}->{"bbduk.sh"};
            foreach my $in_F (@inputsF)
            {
                my ($name, $dir, $ext) = fileparse($in_F, qr{\.gz});
                if ($ext)
                {
                    my $cmd = "gunzip $in_F";
                    push @cmds, $cmd;

                    # $name would now contain <sample_name>.fastq
                }
                my @suffixes = ("_R1.fastq");
                my $prefix   = basename($name, @suffixes);
                my $in_R     = catfile($GLOBAL_PATHS->rev_sample_dir(),
                                   "${prefix}_R2.fastq");
                my $out_F =
                  catfile($GLOBAL_PATHS->{'wd'}, "${prefix}_R1_tc.fastq");
                my $out_R =
                  catfile($GLOBAL_PATHS->{'wd'}, "${prefix}_R2_tc.fastq");

                # trim fwd reads first (skipr2) then trim rev reads (skipr1)
                # max length 301 - (20 - 2) = 283
                my $cmd =
                  "$bbduk -Xmx4915m -Xms4327m in=$in_F in2=$in_R out=stdout.fastq  minlen=60 literal=TGGGATCGTTTTGATGTATT,TGGGATCGCTTTGATGTATT copyundefined rcomp=f skipr2=t restrictleft=30 ktrim=l k=20 hdist=2  mink=18 hdist2=0";
                $cmd .= " 2>>"
                  . catfile($GLOBAL_PATHS->part1_error_log(),
                            "bbduk.sh.stderr");
                $cmd .=
                  " | $bbduk -Xmx4915m -Xms4327m in=stdin.fastq int=t out=$out_F out2=$out_R  minlen=0 maxlen=283 literal=TAAACTTGCTTGCCACTCATG,ACTTGCTTGCCATTCATGGTA copyundefined rcomp=f skipr1=t restrictleft=30 ktrim=l k=20 hdist=2 mink=18 hdist2=0";
                $cmd .= " 2>>"
                  . catfile($GLOBAL_PATHS->part1_error_log(),
                            "bbduk.sh.stderr");
                push @cmds, $cmd;
            }
        }

        my $qsub = grep(/^${var}$/, ("V3V4", "V4", "ITS"));
        execute_and_log(
                        cmds    => \@cmds,
                        logger  => $GLOBAL_PATHS,
                        dry_run => $dryRun,
                        msg  => "Removing $var primers from all sequences.\n",
                        qsub => $qsub
                       );

        my @fwdFiles = glob("$GLOBAL_PATHS->{'wd'}/*R1_tc.fastq");
        for (my $i = 0; $i < scalar @fwdFiles; $i++)
        {
            my $revFile = $fwdFiles[$i] =~ s/R1_tc.fastq$/R2_tc.fastq/r;
            if (-z $fwdFiles[$i] || -z $revFile)
            {
                if (-z $fwdFiles[$i])
                {
                    $GLOBAL_PATHS->print(
                        "Primer trimming filtered all reads from $fwdFiles[$i]. Deleting both sample FASTQs.\n"
                    );
                } else
                {
                    $GLOBAL_PATHS->print(
                        "Primer trimming filtered all reads from $revFile. Deleting both sample FASTQs.\n"
                    );
                }
                unlink $fwdFiles[$i];
                unlink $revFile;
            }
        }

        my $to_gzip = join ",", "*_tc.fastq",
          catfile($GLOBAL_PATHS->fwd_sample_dir(), "*.fastq"),
          catfile($GLOBAL_PATHS->rev_sample_dir(), "*.fastq");
        my $cmd = "gzip -f9 {" . $to_gzip . "}";
        execute_and_log(
                        cmds    => [$cmd],
                        logger  => $GLOBAL_PATHS,
                        dry_run => $dryRun,
                        msg     => "Compressing tagcleaned FASTQ's...\n"
                       );
        $GLOBAL_PATHS->print("---Raw and tagcleaned FASTQ's compressed.\n");

        if (!@dbg)
        {
            cacheChecksums(
                  files => [
                            glob(catfile($GLOBAL_PATHS->fwd_sample_dir(), "*")),
                            glob(catfile($GLOBAL_PATHS->rev_sample_dir(), "*"))
                           ],
                  step_name => "samples",
                  pf        => $GLOBAL_PATHS
            );

            my @gzipped = tagcleaned_files();
            cacheChecksums(
                           files     => \@gzipped,
                           step_name => "tagcleaned",
                           pf        => $GLOBAL_PATHS
                          );
        }

        my $duration = time - $start;
        $GLOBAL_PATHS->print(
                       "---Primer sequences removed from $newSamNo samples.\n");
        $GLOBAL_PATHS->print("---Duration of tagcleaning: $duration s\n");
    } else
    {
        $GLOBAL_PATHS->print(
            "---Primers already removed from forward and reverse reads of $newSamNo samples. Moving on.\n"
        );
    }

    if (@dbg && !grep(/^dada2$/, @dbg))
    {
        $GLOBAL_PATHS->print(
            "Finished removing primers. Terminated because --debug dada2 was not "
              . "specified.\n");
        exit 0;
    }
}

my $stats_file  = "dada2_part1_stats.txt";
my $counts_file = "dada2_abundance_table.rds";

###### BEGIN DADA2 ##########
#############################
if ((!@dbg) || grep(/^dada2$/, @dbg))
{
    my @outputs = ($stats_file, $counts_file);
    my @inputs  = glob("$GLOBAL_PATHS->{'wd'}/*R[1|2]_tc.fastq.gz");

    if ($trySkip)
    {
        $skipMe = skippable(
                   inputs        => \@inputs,
                   outputs       => \@outputs,
                   checksums_in  => $metadata->{"checkpoints"}{"tagcleaned"},
                   checksums_out => $metadata->{"checkpoints"}{"dada2part1out"},
                   step_name     => "DADA2",
                   previous_step_skipped => 1,
                   pf                    => $GLOBAL_PATHS
        );
    }

    if (!$skipMe)
    {
        my $truncLen;

        if ($oneStep)
        {
            if ($var eq "V3V4")
            {

                $truncLenL = 255 if (!$truncLenL);
                $truncLenR = 255 if (!$truncLenR);

                if (!$maxN)
                {
                    $maxN = 0;
                }
                if (!$maxEE)
                {
                    $maxEE = "2";
                }
                if (!$truncQ)
                {
                    $truncQ = 2;
                }
                if (!$phix)
                {
                    $phix = "1";
                }
                if (!$maxLen)
                {
                    $maxLen = "Inf";
                }
                if (!$minLen)
                {
                    $minLen = 20;
                }
                if (!$minQ)
                {
                    $minQ = 0;
                }
            }

            if ($var eq "V4")
            {
                $truncLenL = 200 if (!$truncLenL);
                $truncLenR = 200 if (!$truncLenR);

                if (!$maxN)
                {
                    $maxN = 0;
                }
                if (!$maxEE)
                {
                    $maxEE = "2";
                }
                if (!$truncQ)
                {
                    $truncQ = 2;
                }
                if (!$phix)
                {
                    $phix = "1";
                }
                if (!$maxLen)
                {
                    $maxLen = "Inf";
                }
                if (!$minLen)
                {
                    $minLen = 20;
                }
                if (!$minQ)
                {
                    $minQ = 0;
                }
            }
        } else
        {
            if ($var eq "V3V4")
            {
                $truncLenL = 255 if (!$truncLenL);
                $truncLenR = 225 if (!$truncLenR);

                if (!$maxN)
                {
                    $maxN = 0;
                }
                if (!$maxEE)
                {
                    $maxEE = "2";
                }
                if (!$truncQ)
                {
                    $truncQ = 2;
                }
                if (!$phix)
                {
                    $phix = "1";
                }
                if (!$maxLen)
                {
                    $maxLen = "Inf";
                }
                if (!$minLen)
                {
                    $minLen = 20;
                }
                if (!$minQ)
                {
                    $minQ = 0;
                }
            }

            if ($var eq "V4")
            {
                $truncLenL = 200 if (!$truncLenL);
                $truncLenR = 200 if (!$truncLenR);

                if (!$maxN)
                {
                    $maxN = 0;
                }
                if (!$maxEE)
                {
                    $maxEE = "2";
                }
                if (!$truncQ)
                {
                    $truncQ = 2;
                }
                if (!$phix)
                {
                    $phix = "1";
                }
                if (!$maxLen)
                {
                    $maxLen = "Inf";
                }
                if (!$minLen)
                {
                    $minLen = 20;
                }
                if (!$minQ)
                {
                    $minQ = 0;
                }
            }

            if ($var eq "ITS")
            {
                $truncLenL = 0 if (!$truncLenL);
                $truncLenR = 0 if (!$truncLenR);

                if (!$maxN)
                {
                    $maxN = 0;
                }
                if (!$maxEE)
                {
                    $maxEE = "1";
                }
                if (!$truncQ)
                {
                    $truncQ = 2;
                }
                if (!$phix)
                {
                    $phix = "1";
                }
                if (!$maxLen)
                {
                    $maxLen = "Inf";
                }
                if (!$minLen)
                {
                    $minLen = 50;
                }
                if (!$minQ)
                {
                    $minQ = 0;
                }
            }
            if ($var eq "OMPA")
            {
                $truncLenL = 259 if (!$truncLenL);
                $truncLenR = 150 if (!$truncLenR);

                if (!$maxN)
                {
                    $maxN = 0;
                }
                if (!$maxEE)
                {
                    $maxEE = "2";
                }
                if (!$truncQ)
                {
                    $truncQ = 2;
                }
                if (!$phix)
                {
                    $phix = "1";
                }
                if (!$maxLen)
                {
                    $maxLen = "Inf";
                }
                if (!$minLen)
                {
                    $minLen = 20;
                }
                if (!$minQ)
                {
                    $minQ = 0;
                }
            }
        }
        my $r_out = dada2(
                          $run,    $truncLenL, $truncLenR, $maxN,
                          $maxEE,  $truncQ,    $phix,      $maxLen,
                          $minLen, $minQ,      $dada2mem
                         );

###### EVALUATING DADA2 OUTPUT ##########
#########################################
        my @removed;
        ## do a scan of this file for those samples not passing filtering (0 sequences surviving - and log this)

        open ROUT, "<${r_out}"
          or die "Cannot open ${r_out} for reading: $OS_ERROR";
        while (<ROUT>)
        {
            chomp;
            if ($_ =~ /The filter removed all reads/)
            {
                push @removed, $_;
            }
        }
        close ROUT;

        if (scalar @removed > 0)
        {
            $GLOBAL_PATHS->print(
                scalar @removed . " samples removed during dada2 filtering:\n");
            for my $x (@removed)
            {
                $GLOBAL_PATHS->print("$x\n");
            }
        }
        $GLOBAL_PATHS->print("\n");

        if (List::Util::all {-e $_} @outputs)
        {
            $GLOBAL_PATHS->print(
                        "\nFor $var region, dada2 used the following filtering "
                          . "requirements:\ntruncLen = $truncLenL, $truncLenR\n"
                          . "maxN = $maxN\nmaxEE = $maxEE\ntruncQ = $truncQ\n"
                          . "rm.phix = $phix\n");

            # "dada2 completed successfully" is a key phrase that causes
            # an appropriate message printed to STDOUT
            $GLOBAL_PATHS->print(
                    "dada2 completed successfully!\nAbundance table for "
                  . "run $run located at $GLOBAL_PATHS->{'wd'}/dada2_abundance_table.rds\n"
            );
            $GLOBAL_PATHS->print(
                        "See $outputs[0] for dada2 table of reads surviving by "
                          . "step\n\n");

            if (!@dbg)
            {
                cacheChecksums(
                               files     => \@outputs,
                               step_name => "dada2part1out",
                               pf        => $GLOBAL_PATHS
                              );
            }
        } else
        {
            $GLOBAL_PATHS->print(
                "---dada2 did not complete successfully, something went wrong!\n"
                  . "---Check R output file (usually R.stderr in SGE stderr directory).\n"
            );
            die;
        }

    } else
    {
        $GLOBAL_PATHS->print(
            "DADA2 already processed the same tagcleaned files during a previous run. Doing nothing.\n"
        );
    }
}

if ($delete && -e $stats_file && -e $counts_file)
{
    my @delete_files = (
                        "$GLOBAL_PATHS->{'wd'}/barcodes.fastq",
                        catfile($GLOBAL_PATHS->fwd_demux_dir(), "seqs.fna"),
                        catfile($GLOBAL_PATHS->rev_demux_dir(), "seqs.fna"),
                        $GLOBAL_PATHS->fwd_library(),
                        $GLOBAL_PATHS->rev_library(),
                        tagcleaned_files(),
                        glob(
                             catfile(
                                     $GLOBAL_PATHS->filtered_dir(),
                                     "*_filt\\.fastq\\.gz"
                                    )
                            )
                       );

    for my $file (@delete_files)
    {
        unlink $file;
    }
    rmdir $GLOBAL_PATHS->filtered_dir();
}

###### COMPLETING $logFH FILE ##############
#########################################

END
{
    if (defined ${global_run_storage})
    {
        if ($GLOBAL_PATHS->{'wd'} =~ /^${global_run_storage}/)
        {
            find(
                sub {
                    chmod(0664, $File::Find::name) unless -d;
                    chmod(0775, $File::Find::name) if -d;
                },
                $GLOBAL_PATHS->{'wd'}
                );
            chmod(0775, $GLOBAL_PATHS->{'wd'});
        }
    }

    if (defined $GLOBAL_PATHS)
    {
        $GLOBAL_PATHS->print("\n\n\n");
        $GLOBAL_PATHS->close();
    }

}

####################################################################
##                               SUBS
####################################################################
sub read_json
{
    my $file = shift;
    my $mode = shift;

    my $hash = {};

    if (-e $file)
    {
        # get existing metadata on filesystem
        my $json;
        {
            local $/;    #Enable 'slurp' mode
            open(my $FH, $mode, $file);
            seek $FH, 0, 0 or die;
            $json = <$FH>;
            close $FH;
        }

        eval {$hash = decode_json($json);};
    }

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

sub project_metadata
{
    my %arg          = @_;
    my $any_args     = scalar %arg;
    my $new_metadata = delete %arg{"metadata"} // {};
    my $replace      = delete $arg{"replace"} // 0;
    my $pf           = delete %arg{"pf"} // $GLOBAL_PATHS;

    my $metadataFile = "$pf->{'wd'}/.meta.json";
    my $old_metadata = read_json($metadataFile, "+<");

    # If no arguments, user wants to GET existing metadata
    if (!$any_args)
    {
        return $old_metadata;
    } else
    {
        # If arguments, user wants to UPDATE AND GET existing metadata
        my $combined_metadata;

        if (!$replace)
        {
            # recursively merge params and checkpoints with $project_metadata
            $combined_metadata = merge($old_metadata, $new_metadata);
        } else
        {
            # replace each present element of new_metadata
            $combined_metadata = $old_metadata;
            foreach my $step (keys %$new_metadata)
            {
                $pf->print("Replacing $step\n");
                $combined_metadata->{$step} = $new_metadata->{$step};
            }
        }

        open my $metadataFH, ">$metadataFile";
        print $metadataFH encode_json($combined_metadata);
        close $metadataFH;
        return $combined_metadata;
    }
}

sub preprocess_map
{
    my $map_filepath = shift;
    my $run_name     = shift;

    $run_name =~ s/[^a-zA-Z0-9\.]/\./g;

    rename($map_filepath, $map_filepath . '.bak');
    open(IN,  '<' . $map_filepath . '.bak') or die $!;
    open(OUT, '>' . $map_filepath)          or die $!;
    while (<IN>)
    {
        s/^([a-zA-Z0-9\.]+)\t/$run_name\.$1\t/g;
        print OUT $_;
    }
    close(IN);
    close(OUT);
    unlink($map_filepath . '.bak');

    return;
}

sub qiime_cmd
{
    my $script = shift;
    my $config = shift;

    my $python =
      catfile(abs_path($config->{'part1'}->{"qiime1/bin"}), "python");
    $script = catfile(abs_path($config->{'part1'}->{"qiime1/bin"}), $script);
    return "$python -s $script";
}

#' @_[0] Filenames of inputs to checksum
#' @_[1] Filenames of outputs to checksum (not guaranteed to be existing files)
#' @_[2] A hash pointer containing the checksums for inputs
#' @_[3] A hash pointer containing the checksums for outputs
#' @_[4] The name of the current step, for log messages.
#' @_[5] The names of the inputs (should match keys of @_[2])
#' @_[6] The names of the outputs (should match keys of @_[3])
sub skippable
{
    my %arg          = @_;
    my $inputsR      = delete %arg{"inputs"};
    my $outputsR     = delete %arg{"outputs"};
    my $ci           = delete %arg{"checksums_in"} // {};
    my %checksumsIn  = %{$ci};
    my $co           = delete %arg{"checksums_out"} // {};
    my %checksumsOut = %{$co};
    my $step         = delete %arg{"step_name"};
    my $inNames      = delete %arg{"in_names"};
    my $outNames     = delete %arg{"out_names"};
    my $prevSkip     = delete %arg{"previous_step_skipped"};
    my $pf           = delete %arg{"pf"};

#consider changing keys of %checksumsIn and %checksumsOut to absolute paths. Same
# with anything incoming $inputsR and $outputsR. whoever uses this function
# shouldn't be tasked with knowing the cache structure

    # Remove any key-value pairs with undef value. This is the necessary to
    # test if we have reference checksums for all the files present.
    delete $checksumsIn{$_}
      for grep {not defined $checksumsIn{$_}} keys %checksumsIn;
    delete $checksumsOut{$_}
      for grep {not defined $checksumsOut{$_}} keys %checksumsOut;

    if ($verbose)
    {
        $pf->print("Input files: " . Dumper($inputsR) . "\n");
        $pf->print("Cached files: " . Dumper(%checksumsIn) . "\n");
        $pf->print("Output files: " . Dumper($outputsR) . "\n");
        $pf->print("Cached files: " . Dumper(%checksumsOut) . "\n");
    }

    if (   scalar @$outputsR == 0
        || List::Util::none {defined $_} @$outputsR
        || List::Util::none {-e $_} @$outputsR)
    {
        $pf->print(ucfirst("$step has not been run yet. Running now...\n"));
        return 0;
    }    # no output files found by caller

    if (!List::Util::all {-e $_} @$outputsR)
    {

        # Caller named output files but not all existed
        $pf->print(
            ucfirst(
                "$step output is incomplete. Removing outputs files and executing $step...\n"
            )
        );
        foreach (@$outputsR)
        {
            if (-e $_)
            {
                unlink();
            }
        }
        return 0;
    }

    # When previously stored, the checksums should have been named consistently,
    # or named after the files themselves.
    if (!defined $inNames)
    {
        $inNames =
          [map {File::Spec->abs2rel($_, $GLOBAL_PATHS->{'wd'})} @$inputsR];
    }
    if (!defined $outNames)
    {
        $outNames =
          [map {File::Spec->abs2rel($_, $GLOBAL_PATHS->{'wd'})} @$outputsR];
    }

    if (!setIdent($inNames, [keys %checksumsIn]))
    {
        # The input files are named differently from last time (or pipeline
        # provided inconsistent file nicknames)
        $pf->print(
            "Input to $step does not match previous records. Discarding $step ",
            "outputs and repeating $step.\n"
        );
        unlink @$outputsR;
        return 0;
    }
    if (!setIdent($outNames, [keys %checksumsOut]))
    {

# The outputs files are named differently from last time (or pipeline provided inconsistent file nicknames)
        $pf->print(
            ucfirst($step)
              . " output does not match previous records. Discarding $step outputs and repeating $step...\n"
        );
        unlink @$outputsR;
        return 0;
    }

    $pf->print(ucfirst("$step") . " output exists. Trying to validate...\n");

    if (!$prevSkip)
    {

        # If the previous step was not skipped, the inputs must be checksummed
        if ($verbose)
        {
            $pf->print("Checksumming input files:\n" . Dumper($inputsR));
        }
        my $i          = -1;
        my $inputValid = List::Util::all
        {
            if ($verbose)
            {
                $pf->print("Checksumming input file: $_\n");
            }
            $i++;
            my @stdout   = split(/\s/, `sha512sum $_`);
            my $checksum = $stdout[0];

            $checksum eq $checksumsIn{$inNames->[$i]};
        }
        @$inputsR;

        # If all inputs are/nt valid give appropriate messages
        if ($inputValid)
        {
            $pf->print("Input to $step matches previous records.\n");
        } else
        {
            $pf->print(
                "Input to $step has changed since last execution. Discarding $step outputs and\n",
                "repeating $step.\n"
            );
            unlink @$outputsR;
            return 0;
        }
    }    # Otherwise, we know the inputs are valid
         # then do the same for all the outputs

    my $i = -1;
    if ($verbose)
    {
        $pf->print("Checksumming output files:\n" . Dumper($outputsR));
    }
    my $outputValid = List::Util::all
    {
        $i++;
        if ($verbose)
        {

            $pf->print("Checksumming output file:$_\n");
        }
        my @stdout   = split(/\s/, `sha512sum $_`);
        my $checksum = $stdout[0];

        $checksum eq $checksumsOut{$outNames->[$i]};
    }
    @$outputsR;

    # If all outputs are/nt valid give appropriate messages
    if ($outputValid)
    {
        $pf->print("Output to $step matches previous records.\n");
        return 1;
    } else
    {
        $pf->print(
            "Output to $step does not match previous records. Discarding $step outputs and\n",
            "repeating $step.\n"
        );
        unlink @$outputsR;
        return 0;
    }

    return 0;
}

#' @_[0] A pointer to an array containing filenames to checksum
#' @_[1] The name of the step whose products are being cached
#' @_[2] Optional pointer to array of names (include to avoid IDing files by filenames)
sub cacheChecksums
{
    my %arg   = @_;
    my $files = delete %arg{"files"};
    my $step  = delete %arg{"step_name"};
    my $pf    = delete %arg{"pf"};
    my $names = delete %arg{"file_names"}
      // [map {File::Spec->abs2rel($_, $pf->{'wd'})} @$files];

    $pf->print("Saving progress in run folder...\n");

    my $destroyNextStep = sub {
        if ($step2number{$step} + 1 < scalar keys %number2step)
        {

            # Reset the next set of checksums so the next step is forced to run
            $metadata->{"checkpoints"}{$number2step{$step2number{$step} + 1}} =
              {};
        }
    };
    my $update = 0;
    my @keys   = keys %{$metadata->{"checkpoints"}{$step}};
    if (!setIdent(\@keys, $names))
    {
        $update = 1;
        $destroyNextStep->()
          ;    # MUST DO THIS BEFORE ANY NEW CHECKSUMS ARE WRITTEN
         # PREVENTS "CHIMERIC" CHECKSUM FILE IN CASE OF PIPELINE INTERRUPTION HERE
    }
    my %newChecksums;

    for my $i (0 .. scalar @$files - 1)
    {
        my $file = $files->[$i];
        if (!-e -f $file)
        {
            die "Non-existant file passed to cacheChecksums(): $file\n";
        }
        my $name     = $names->[$i];
        my @output   = split(/\s/, `sha512sum $file`);
        my $checksum = $output[0];
        if (!$update)
        { # !$update means that $checksums{$step}->{ $file } is defined for all files
            if ($metadata->{"checkpoints"}{$step}->{$name} ne $checksum)
            {
                $destroyNextStep->();
            }
        }

        # Store the checksum of the files just produced.
        $newChecksums{$name} = $checksum;
    }

    $metadata->{"checkpoints"}{$step} = \%newChecksums;
    if (List::Util::all {defined $_} $metadata->{"checkpoints"})
    {
        project_metadata(
                      metadata => {"checkpoints" => $metadata->{"checkpoints"}},
                      pf       => $pf,
                      replace  => 1
        );
    }
}

sub setIdent
{
    my ($left, $right) = @_;
    return 0 if scalar @$left != scalar @$right;
    my %hash;
    @hash{@$left, @$right} = ();
    return scalar keys %hash == scalar @$left;
}

# @param 0 mapping file path
sub count_samples
{
    my $file = shift;
    my $ans  = count_lines($file) - 1;
    return $ans;
}

sub count_lines
{
    my $file = shift;
    open my $FH, "<$file" or die "Cannot open $file for reading: $OS_ERROR";
    while (<$FH>) { }
    my $ans = $.;
    close $FH;
    return $ans;
}

# @param 0 original input directory
# @param 1 whether one-step or not
# @param 2 a scalar file handle to write logs to
sub find_raw_files
{
    my %arg     = @_;
    my $pf      = delete %arg{"pf"};
    my $dir     = delete %arg{"dir"};
    my $oneStep = delete %arg{"pcr_one_step"};

    my $index1Input;
    my $index2Input;
    my $readsForInput;
    my $readsRevInput;

    # These are the only input files needed for 1step
    my @r1s = glob("$dir/*R1.fastq $dir/*R1.fastq.gz");
    my @r2s = glob("$dir/*R2.fastq $dir/*R2.fastq.gz");

    if ($oneStep)
    {
        if (scalar @r1s == 1 && scalar @r2s == 1)
        {
            $index1Input = $readsForInput = $r1s[0];
            $index2Input = $readsRevInput = $r2s[0];
        } else
        {
            my $files_found = "";
            my @file_arrays = (\@r1s, \@r2s);
            foreach my $file_array_ref (@file_arrays)
            {
                my @file_array = @$file_array_ref;
                if (scalar @file_array)
                {
                    $files_found .= (join "\n", @file_array);
                }
            }

            my $message =
                "Could not find a complete and exclusive set of raw files."
              . " Since --1step given, input directory must have exactly one"
              . " R1 and R2 file. Files found:\n"
              . $files_found;

            $pf->print($message);

            die $message;

        }
    } else
    {
        my @i1s = glob("$dir/*I1.fastq $dir/*I1.fastq.gz");
        my @i2s = glob("$dir/*I2.fastq $dir/*I2.fastq.gz");

        # also globbed for r1 and r2 files above
        my @r3s = glob("$dir/*R3.fastq $dir/*R3.fastq.gz");
        my @r4s = glob("$dir/*R4.fastq $dir/*R4.fastq.gz");

        if (   scalar @i1s == 1
            && scalar @i2s == 1
            && scalar @r1s == 1
            && scalar @r2s == 1
            && scalar @r3s == 0
            && scalar @r4s == 0)
        {
            $index1Input   = $i1s[0];
            $index2Input   = $i2s[0];
            $readsForInput = $r1s[0];
            $readsRevInput = $r2s[0];
        } elsif (   scalar @i1s == 0
                 && scalar @i2s == 0
                 && scalar @r1s == 1
                 && scalar @r2s == 1
                 && scalar @r3s == 1
                 && scalar @r4s == 1)
        {
            $index1Input   = $r2s[0];
            $index2Input   = $r3s[0];
            $readsForInput = $r1s[0];
            $readsRevInput = $r4s[0];
        } else
        {
            my $files_found = "";
            my @file_arrays = (\@i1s, \@i2s, \@r1s, \@r2s, \@r3s, \@r4s);
            foreach my $file_array_ref (@file_arrays)
            {
                my @file_array = @$file_array_ref;
                if (scalar @file_array)
                {
                    $files_found .= (join "\n", @file_array);
                }
            }
            my $message =
              "Could not find a complete and exclusive set of raw files. Input directory must"
              . " contain exactly one set of I1, I2, R1, and R2, OR R1, R2, R3, and R4.\n"
              . "File found in $dir:\n$files_found";
            $pf->print($message);
            die $message;
        }
    }
    return ($readsForInput, $readsRevInput, $index1Input, $index2Input);
}

# Given *R1.fastq(.gz), if it's .gz, then removes the extension and gives a
# filepath in the local directory
# @param 0 The run directory
# @params 1... The full path to the original file
sub convert_to_local_if_gz
{
    my $wd    = shift;
    my @files = @_;
    my @ans;
    foreach my $file (@files)
    {
        if ($file =~ /.gz$/)
        {

            # Rename *[R|I][1|2].fastq.gz to <WD>/<RUN>_[R|I][1|2].fastq
            my $dest = basename($file);

            my $suffix =
              substr($dest, (length($dest) - 11), 8);    # get suffix (sans .gz)

            $wd =~ s/\/$//;    # remove any trailing slash
            my @dirs = File::Spec->splitdir($wd);
            $file = "$wd/$dirs[scalar(@dirs) - 1]" . "_$suffix";
        }
        push(@ans, $file);
    }
    return @ans;
}

sub find_index_length
{
    my $file = shift;
    open my $reads, $file or die "Could not open $file: $!";
    my $first_index;
    while (<$reads>)
    {
        chomp;
        $first_index = $_ if $. == 2;
        last if $. == 2;
    }
    close $reads;

    return (length $first_index);
}

sub readSplitLog
{
    my %arg       = @_;
    my $pf        = delete %arg{"logger"};
    my $split_log = delete %arg{"file"};
    my $verbose   = delete %arg{"verbose"} // undef;

    ##Check split_library log for 0's
    open SPLIT, "<$split_log"
      or die "Cannot open $split_log for writing: " . "$OS_ERROR";
    my @split;
    while (<SPLIT>)
    {
        if ($_ =~ /\t/)
        {
            chomp;
            my ($sample, $nReads) = split, /\t/;
            chomp $nReads;
            if ($nReads eq "0")
            {
                push @split, $sample;
            } else
            {
                if ($verbose)
                {
                    $pf->print("$_\n");
                }
            }
        }
    }
    close SPLIT;
    return @split;
}

sub dada2
{
    my $run       = shift;
    my $truncLenL = shift;
    my $truncLenR = shift;
    my $maxN      = shift;
    my $maxEE     = shift;
    my $truncQ    = shift;
    my $phix      = shift;
    my $maxLen    = shift;
    my $minLen    = shift;
    my $minQ      = shift;
    my $dada2mem  = shift;

    my $Rscript =
      catfile($pipelineDir, "scripts", "filter_and_denoise_illumina.R");
    my $args =
      "--truncLenL=$truncLenL --truncLenR=$truncLenR --maxN=$maxN --maxEE=$maxEE --truncQ=$truncQ --rm.phix=$phix";

    chdir $GLOBAL_PATHS->{'wd'};

    my $exitStatus = 1;
    $dada2mem =~ s/[a-zA-Z]//g;
    my $previous_mem = 0;
    my $try          = 1;
    my $R_out        = catfile($GLOBAL_PATHS->part1_error_log(), "R.stderr");

    while ($exitStatus == 1 && ($dada2mem != $previous_mem || $try < 3))
    {
        my $cmd =
            "rm -rf "
          . $GLOBAL_PATHS->filtered_dir() . " "
          . catfile($GLOBAL_PATHS->{'wd'}, "dada2_part1_stats.txt ")
          . catfile($GLOBAL_PATHS->{'wd'}, "dada2_abundance_table.rds ")
          . catfile($GLOBAL_PATHS->{'wd'}, ".RData");
        execute_and_log(
            cmds    => [$cmd],
            dry_run => $dryRun,
            msg =>
              "Removing old filtered fastq files, stats, and Rout files from previous runs\n"
        );

        $cmd =
          "$config_hashref->{'executor'} -l mem_free=${dada2mem}G -P $qproj -e ${R_out} -o "
          . $GLOBAL_PATHS->part1_stdout_log()
          . " -N Rscript \"$config_hashref->{R}script $Rscript $args\"";
        execute_and_log(
            cmds    => [$cmd],
            logger  => $GLOBAL_PATHS,
            dry_run => $dryRun,
            msg =>
              "Running DADA2 with fastq files in $GLOBAL_PATHS->{'wd'} for $var region...\n",
            qsub => 1
        );

        if (-e "dada2_part1_stats.txt")
        {    # sign of success
            $GLOBAL_PATHS->print("---R script completed without errors.\n");
            $GLOBAL_PATHS->print(
                  "---DADA2-specific commands can be found in " . "$Rscript\n");
            $exitStatus = 0;    # Move on from DADA2
        } else
        {
            open IN, "<${R_out}"
              or die "Could not open ${R_out} to read DADA2 R log: $OS_ERROR\n";
            my $line = <IN>;
            while ($line)
            {
                if (            # signs of bad termination that might be fixed
                                # with more memory
                    (
                        $line =~ /Error in/
                     || $line =~ /Execution halted/
                     || $line =~ /encountered errors/
                     || $line =~ /Error:/
                     || $line =~ /Traceback:/
                    )
                    && $line !~ /learnErrors/    # False positives; don't match
                    && $line !~ /error rates/
                    && $line !~ /errors/
                   )
                {
                    # Don't change $exitStatus, so DADA2 is restarted

                    $GLOBAL_PATHS->print("---R script crashed at: $line\n");

                    # Preserve the last R log file that errored. Get rid of the
                    # old R output file, then run R again.
                    File::Copy::move("${R_out}", "${R_out}.old");
                    $GLOBAL_PATHS->print("---See ${R_out}.old for details.\n");
                    $GLOBAL_PATHS->print("---Attempting to restart R...\n");

                    $previous_mem = $dada2mem;
                    $dada2mem     = List::Util::min($dada2mem * 2, 1);

                    last;
                } elsif ($line =~ /(E|e)rror/ || $line =~ /ERROR/)
                {
                    my $error_msg = "Error in DADA2. See ${R_out} for details.";
                    $GLOBAL_PATHS->print($error_msg);
                    die $error_msg;
                }
                $line = <IN>;
            }
            close IN;

            $try += 1;
        }
    }
    if ($try > 2)
    {
        my $msg =
          "Won't attempt to restart R, as the same command has been attempted twice.";
        $GLOBAL_PATHS->print($msg);
        die $msg;
    }
    if ($exitStatus == 1)
    {
        my $msg =
          "R output doesn't indicate error, but analysis files were not produced.";
        $GLOBAL_PATHS->print($msg);
        die $msg;
    }

    return ${R_out};
}

sub source
{
    my $name = shift;

    open my $fh, "<$name" or die "could not open $name: $!";

    while (<$fh>)
    {
        chomp;
        my ($k, $v) = split /=/, $_, 2;
        $v =~ s/^(['"])(.*)\1/$2/;            #' fix highlighter
        $v =~ s/\$([a-zA-Z]\w*)/$ENV{$1}/g;
        $v =~ s/`(.*?)`/`$1`/ge;              #dangerous
        $ENV{$k} = $v;
    }
}

################################################################################
# Execute the given commands;
#  par 1: array ref containing commands
# par 2: Lexical filehandle that overrides printing to STDOUT. Pass 0 for no override.
# par 3: $dryRun parameter
# par 4: Human-friendly message to be printed to STDOUT.
# par 5: true if commands are to be qsubbed (causes return value to be array of SGE job IDs)
# Remaining args: cmds to be run through system()
sub execute_and_log
{
    my %arg      = @_;
    my $cmds     = delete %arg{"cmds"};
    my @commands = @$cmds;
    my $fh       = delete %arg{"logger"} // $GLOBAL_PATHS;
    my $dryRun   = delete %arg{"dry_run"} // 0;
    my $cuteMsg  = delete %arg{"msg"} // 0;
    my $qsub     = delete %arg{"qsub"} // 0;

    if (!@commands)
    {
        warn "execute_and_log() called with no commands\n";
    }

    $fh->print("$cuteMsg\n");

    my @pids;
    my $sge;
    if ($qsub)
    {
        $sge = Schedule::SGE->new(
                -executable => {
                    qsub  => '/usr/local/packages/sge-root/bin/lx24-amd64/qsub',
                    qstat => '/usr/local/packages/sge-root/bin/lx24-amd64/qstat'
                }
        );
    }

    foreach my $cmd (@commands)
    {
        if ($verbose)
        {    # print each command
            $fh->print("\$ $cmd\n");
        } else
        {
            STDOUT->print("\$ $cmd\n");
        }

        if ($qsub)
        {
            $sge->command($cmd);
            push @pids, $sge->execute();
        } else
        {

            system($cmd) == 0
              or die "system($cmd) failed with exit code: $?"
              if !$dryRun;
        }

    }

    if ($qsub)
    {
        my $done = 0;
        while (!$done)
        {
            $sge->status();
            $done = 1;
            foreach (@pids)
            {
                my @job_stats = @{$sge->brief_job_stats($_)};
                if (@job_stats)
                {
                    $done = 0;
                    last;
                }
            }
            sleep 1;
        }
    }
    $fh->print("Finished " . lcfirst($cuteMsg) . "\n");
}

sub tagcleaned_files
{
    return glob(catdir($GLOBAL_PATHS->{'wd'}, "*R[1|2]_tc.fastq.gz"));
}

exit 0;
