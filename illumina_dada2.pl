#!/usr/bin/env perl
#. /usr/local/packages/usepackage/share/usepackage/use.bsh

=head1 NAME

illumina_dada2.pl

=head1 DESCRIPTION

The script can be launched from any location on the IGS server, it automatically 
produces a directory in /local/groupshare/ravel named after the project and run
ID provided. 

Beginning with a set of raw Illumina sequencing files (usually R1, R2, I1, and I2),
a mapping file, a project ID, a run ID, and specifying the targeted variable
region, this script:
  1. Extracts barcodes from the raw files.
  2. Demultiplexes the raw reads into fastq files containing reads specific to 
  this project.
  3. Produces individual .fastq files for each sample listed in the mapping file
  4. Performs tag-cleaning of each file
  5. Runs the forward and reverse reads through the dada2 pipeline for the 
  specified 16S rRNA gene region.

A log file is written at <PROJECT>/<RUN>/<PROJECT>_<RUN>_16S_pipeline_log.txt

=head1 SYNOPSIS

THE FOLLOWING STEPS ARE REQUIRED BEFORE THIS SCRIPT WILL RUN:

1. Qlogin to RHEL7 if you ssh'd onto a non-RHEL7 host: 

  qlogin -P jravel-lab -l mem_free=500M -q interactive.q -V

2. Execute from bash: 

  export LD_LIBRARY_PATH=/usr/local/packages/gcc/lib64
  source /usr/local/packages/usepackage/share/usepackage/use.bsh
  use python-2.7
  use qiime

3. then run script as below: 

FOR 2-STEP:
qsub -cwd -b y -l mem_free=1G -P jravel-lab -q threaded.q -pe thread 4 -V
  -e <path_to_logs> -o <path_to_logs> /home/jholm/bin/illumina_dada2.pl 
  -i <path to raw files> -p <project name> -r <run id> -v <variable region> 
  -m <full path to mapping file> -sd <storage directory>

OR:
qsub -cwd -b y -l mem_free=1G -P jravel-lab -q threaded.q -pe thread 4 -V
  -e <path_to_logs> -o <path_to_logs> /home/jholm/bin/illumina_dada2.pl 
  -r1 <path_to_fwd_reads_file> -r2 <path_to_rev_reads_file> -i1 <path_to_index_1_file>
  -i2 <path_to_index_2_file> -p <project name> -r <run id> -v <variable region> 
  -m <full path to mapping file> -sd <storage directory>

FOR 1-STEP
qsub -cwd -b y -l mem_free=1G -P jravel-lab -q threaded.q -pe thread 4 -V
  -e <path_to_logs> -o <path_to_logs> illumina_dada2.pl -i <input directory>
  -p <project name> -r <run ID> -m <mapping file> -v <variable region> 
  -sd <storage directory> --1Step

qsub -cwd -b y -l mem_free=1G -P jravel-lab -q threaded.q -pe thread 4 -V
  -e <path_to_logs> -o <path_to_logs> illumina_dada2.pl -r1 <path_to_fwd_reads_file>
  -r2 <path_to_rev_reads_file> -p <project name> -r <run ID> -m <mapping file>
  -v <variable region> -sd <storage directory> --1Step

=head1 OPTIONS

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

=item B<--var-reg>={V3V4, V4, ITS}, B<-v> {V3V4, V4, ITS}

The targeted variable region.

=item B<--1Step>

Use this flag if the data are prepared by 1-Step PCR (only forward & reverse read files
available)

=item B<--working-dir>=PATH, B<-wd> PATH

Indicate an existing directory in which to place output, and from which the
project and run names will be parsed. The last directory on PATH should be named
after the run, and the second-to-last directory on PATH should be named after
the project.

=item B<-h>, B<--help>

Print help message and exit successfully.

=item B<--qsub-project>=space, B<-qp> space

Indicate which qsub-project space should be used for all qsubmissions. The
default is jravel-lab.

=item B<--debug> {validate, barcodes, demux, splitsamples, tagclean, dada2}

Runs the specified section of the pipeline. Multiple --debug options can be given
to run multiple consecutive parts of the pipeline, provided that the input to
the earliest requested step is present. Any non-consecutive steps will be 
ignored.

=item B<--verbose>

Prints each command to STDOUT.

=back 

=cut

use strict;
use warnings;
my $scriptsDir;
my $pipelineDir;

BEGIN {
    use File::Spec::Functions;
    require File::Basename;
    $pipelineDir = File::Basename::dirname(__FILE__);
    $scriptsDir  = catdir( $pipelineDir, "scripts" );

}
use lib $scriptsDir;    # .pm files in ./scripts/ can be loaded

require Version;

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

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

my @dbg;
my $oneStep = 0;
my $dada2mem = "1G";
GetOptions(
    "raw-path|i=s"           => \my $inDir,
    "r1=s"                   => \my $r1file,
    "r2=s"                   => \my $r2file,
    "i1=s"                   => \my $i1file,
    "i2=s"                   => \my $i2file,
    "map|m=s"                => \my $map,
    "var-reg|v=s"            => \my $var,
    "help|h!"                => \my $help,
    "d|debug=s"              => \@dbg,
    "verbose!"               => \my $verbose,
    "dry-run!"               => \my $dryRun,
    "noskip!"                => \my $noSkip,
    "dada2-truncLen-f|for=i" => \my $f,
    "dada2-truncLen-r|rev=i" => \my $r,
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
  )
##add option for final resting place of important data

  or pod2usage( verbose => 0, exitstatus => 1 );

$verbose = 1;

if ($help) {
    pod2usage( verbose => 2, exitstatus => 0 );
    exit 1;
}

####################################################################
## PARAMETER CHECKING PT 1
####################################################################

# Check existence/compatibility of required parameters
if (
    !$inDir
    && !(
           ( !$oneStep && $r1file && $r2file && $i1file && $i2file )
        || ( $oneStep && $r1file && $r2file )
    )
  )
{
    my $parameters = $oneStep ? "-r1, -r2" : "-r1, -r2, -i1, -i2";
    die "\n\tPlease provide the location of the raw sequencing files "
      . "(single directory => -i)\n\t\tOR \n\tFull paths to each raw file => "
      . "$parameters)\n\n";
}
if ( $inDir && ( $r1file || $r2file || $i1file || $i2file ) ) {
    die
"\n\tInput directory (-i) and raw files (-r*|-i*) were both specified. Please "
      . "provide one or the other.\n\n";
}
if ( !$map ) {
    print "\n\tPlease provide the path to the project mapping file (-m)\n\n";
    exit 1;
}
if ( $oneStep && ( $i1file || $i2file ) ) {
    die "\n\t--1Step is incompatible with -i1 and -i2.\n\n";
}
if ( !$var ) {
    print "\n\tPlease indicate the targeted variable region (-v V3V4 or -v V4"
      . " or -v ITS)\n\n";
    exit 1;
}
if ( !$wd ) {
    die "\n***Please choose a working directory (-wd).";
}

if (@dbg) {
    my $q  = grep( /^validate$/,     @dbg );
    my $e  = grep( /^barcodes$/,     @dbg );
    my $de = grep( /^demux$/,        @dbg );
    my $s  = grep( /^splitsamples$/, @dbg );
    my $t  = grep( /^tagclean$/,     @dbg );
    my $da = grep( /^dada2$/,        @dbg );
    if ( $q + $e + $de + $s + $t + $da == scalar @dbg ) { }
    else {
        die "Illegal debug option. Legal debug options are "
          . "validate, barcodes, demux, splitsamples, tagclean, and dada2.";
    }
    $noSkip = 1;
}


if ( $f && !$r ) {
    print "***\nPlease provide truncation lengths for forward and reverse "
      . "reads\n";
    die;
}

# Refine and validate all variables that refer to the filesystem
my (@paths) = (
    {
        path => \$wd,
        name => "Working directory"
    },
    {
        path => \$inDir,
        name => "Raw file directory"
    },
    {
        path => \$r1file,
        name => "Raw forward reads file"
    },
    {
        path => \$r2file,
        name => "Raw reverse reads file"
    },
    {
        path => \$i1file,
        name => "Raw forward index file"
    },
    {
        path => \$i2file,
        name => "Raw reverse index file"
    },
    {
        path => \$map,
        name => "Mapping file"
    }
);
foreach (@paths) {
    if ( ${$$_{path}} ) {    # If the variable is a non-empty string...
        my $copy = ${$$_{path}};
        $copy =~ s/\/$//;    # remove any trailing slash
        my $relPath = $copy;
        $copy = abs_path( $relPath );    # get abs path
            # At the same time, abs_path returns undef if the path doesn't exist
            # so we can verify the existence of each file and directory
        if ( !defined $copy ) {
            die $$_{name}
              . " not found. Looked for "
              . $relPath . ".\n"
              . "Current working directory: "
              . getcwd() . "\n";
        }
        ${$$_{path}} = $copy; # external variable referenced by path has now been edited.
    }
}


# Find the input files, and find their putative locations if the index files
# were to be decompressed to our working directory.
# In one-step runs, $index1Input and $index2Input are the SAME files
# pointed to by $readsForInput and $readsRevInput
my $index1Input;
my $index2Input;
my $readsForInput;
my $readsRevInput;
( $readsForInput, $readsRevInput, $index1Input, $index2Input ) =
    $inDir ? find_raw_files( $inDir, $oneStep )
  : $oneStep ? ( $r1file, $r2file, $r1file, $r2file )
  :            ( $r1file, $r2file, $i1file, $i2file );
my %localNames;
@localNames{ ( "readsFor", "readsRev", "index1", "index2" ) } =
  convert_to_local_if_gz( $wd, $readsForInput, $readsRevInput,
    $index1Input, $index2Input );

if (
    # Check if there are at least two directories in $wd
    scalar( File::Spec->splitdir( ( File::Spec->splitpath("$wd/") )[1] ) ) -
    2 < 2
  )
{
    die "Working directory (-wd) must have the pattern */PROJECT/RUN";
}

####################################################################
## INITIALIZATION
###################################################################

# Initialize log file
my $run     = File::Basename::basename($wd);
my $pd      = ( File::Basename::fileparse($wd) )[1];
my $project = File::Basename::basename($pd);

my $time = strftime( "%Y-%m-%d %H:%M:%S", localtime(time) );
my $log  = "$wd/$project" . "_" . $run . "_16S_pipeline_log.txt";
open my $logFH, ">>$log" or die "Cannot open $log for writing: $OS_ERROR";
my $logTee = new IO::Tee( \*STDOUT, $logFH );

local $SIG{__WARN__} = sub {
    print $logFH "WARNING: $_[0]";
    print "WARNING: $_[0]";
    print STDERR "WARNING: $_[0]";
};    # Warnings go to log file, stdout, and stderr

# Now that we're done checking parameters, I want the more serious error
# messages to be printed to the log file.
$SIG{__DIE__} = sub {
    if ($^S) {
        die @_;    # If die was called inside an eval, simply pass it along
    } else {
        print $logFH @_ . "\n";
        return 1;
    }
};

print "Logging to: $log\n";
print $logTee "PIPELINE VERSION: " . Version::version() . "\n";
print $logTee "$time\n";

my $error_log  = "$wd/qsub_error_logs";
my $stdout_log = "$wd/qsub_stdout_logs";
## Instead of working directory, give option for provided directory (or current
#  working directory)

my $localMap = "$error_log/" . File::Basename::basename( $map, ".txt" ) . "_corrected.txt";
my $mapLog   = "$error_log/" . File::Basename::basename( $map, ".txt" ) . ".log";

my $fwdProjDir = "$wd/fwdSplit";
my $revProjDir = "$wd/revSplit";

my $fwdSampleDir = "$fwdProjDir/split_by_sample_out";
my $revSampleDir = "$revProjDir/split_by_sample_out";

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
my $skip;

my $models;
if ( $var eq 'V3V4' ) {
    $models = "/local/projects-t2/jholm/PECAN/v1.0/V3V4/merged_models/";
}
if ( $var eq 'V4' || $var eq 'ITS' ) {
    $models = "PECAN not used.\n";
}

if ( !$qproj ) {
    print
      "\nqsub-project ID (--qsub-project, -qp) not provided. Using jravel-lab "
      . " as default\n";
    $qproj = "jravel-lab";
}

my $R = "/usr/local/packages/r-3.4.0/bin/R";

if ( !-e $error_log ) {
    mkdir $error_log;
    print "$error_log did not exist -> making $error_log\n";
} else {

    # Remove old error logs
    my @oldLogs = glob("$error_log/*.*");
    if (@oldLogs) {
        my @cmds;
        foreach my $log (@oldLogs) {
            my @dontDelete = ( "$error_log/illumina_dada2.pl.stderr", $localMap, $mapLog );
            if ( List::Util::none { $_ eq $log }
                @dontDelete )
            {
                push @cmds, "rm $log";
            }
        }
        execute_and_log( @cmds, 0, $dryRun, "Removing stale error logs." );
    }
}

if ( !-e $stdout_log ) {
    mkdir $stdout_log;
    print "$stdout_log did not exist -> making $stdout_log\n";
}

my $checksumFile = "$wd/.checkpoints.json";
my %checksums    = read_checkpoints($checksumFile);

my $qiime = "$wd/$project" . "_" . $run . "_" . "qiime_config.txt";

if (@dbg) {
    print $logTee "DBG FLAGS: ";
    for (@dbg) {
        print $logTee "$_ ";
    }
    print $logTee "\n";
}
print $logTee "PROJECT: $project\nVARIABLE REGION: $var\n"
  . "R VERSION: $R\nPECAN MODELS: $models\n";
if ($oneStep) {
    print $logTee "PCR PREPARATION METHOD: 1-Step\n";
} else {
    print $logTee "PCR PREPARATION METHOD: 2-Step\n";
}

if (   !@dbg
    || grep( /^barcodes$/,     @dbg )
    || grep( /^demux$/,        @dbg )
    || grep( /^splitsamples$/, @dbg ) )
{
    ###### BEGIN CHECK OF QIIME CONFIGURATION ###########
    #####################################################
    my $cmd = "print_qiime_config.py > $qiime";
    execute_and_log( $cmd, $logTee, $dryRun, "QIIME CONFIGURATION DETAILS:\nsee $qiime\n" );
} else {
    print $logTee "NOT USING QIIME\n";
}

if ( !$noSkip ) {
    $skip = skippable( [$map], [ $localMap, $mapLog ],
        $checksums{"map"}, $checksums{"meta"}, "map validation", ["map"] );
}

if ( !$skip ) {
    ###### BEGIN VALIDATION OF MAPPING FILE ###########
    ################################################
    print $logTee "MAPPING FILE: $map\n";

    my $cmd = "validate_mapping_file.py -m $map -s -o $error_log";
    execute_and_log( $cmd, $logTee, $dryRun, "Validating map from $map\n" );
    my $mappingError = glob("$error_log/*.log");
    if ($mappingError) {
        open MAPERROR, "<$mappingError"
          or die "Cannot open $mappingError for " . "reading: $OS_ERROR";
        $_ = <MAPERROR>;
        close MAPERROR;
        chomp;
        if ( $_ =~ /No errors or warnings found in mapping file./ ) {
            print $logTee "---Map passed validation.\n";

        } else {
            die "***Error in mapping file. See $mappingError for "
              . "details. Exiting.\n";
        }
    } else {
        die
"validate_mapping_file.py terminated but did not signal success. Normally a success message is printed in its error file.";
    }

    ###### BEGIN EVALUATION OF SAMPLES VIA MAPPING FILE ###########
    ###############################################################
    open MAP, "<$map" or die "Cannot open $map for reading: $OS_ERROR";
    my $extctrl     = 0;
    my $pcrpos      = 0;
    my $pcrneg      = 0;
    my $projSamples = 0;
    my $linecount   = 0;
    my $null        = 0;
    while (<MAP>) {
        chomp;
        ## don't count header as sample; don't count any line if it doesn't
        if ( $. > 1 ) {
            ## start with four tab-separated fields containing non-whitespace chars
            if ( $_ =~ /^\S+(\t\S+){3}/ ) {
                if (   $_ =~ "EXTNTC"
                    || $_ =~ "NTC.EXT"
                    || $_ =~ "NTCEXT"
                    || $_ =~ "EXT.NTC"
                    || $_ =~ "NTC"
                    || $_ =~ "EXTNEG" )
                {
                    $extctrl++;
                } elsif ( $_ =~ "PCRPOS"
                    || $_ =~ "PCR.pos"
                    || $_ =~ "pCR.pos"
                    || $_ =~ "POSCTRL"
                    || $_ =~ "POS.CTRL"
                    || $_ =~ "POSCON"
                    || $_ =~ "posctr" )
                {
                    $pcrpos++;
                } elsif ( $_ =~ "PCRNTC"
                    || $_ =~ "PCR.NEG"
                    || $_ =~ "PCR.NTC"
                    || $_ =~ "PCRNEG"
                    || $_ =~ "PCRNEGCTRL"
                    || $_ =~ "ntcctr" )
                {
                    $pcrneg++;
                } elsif ( $_ =~ /NULL/ ) {
                    $null++;
                } else {
                    $projSamples++;
                }
            } elsif ( $_ =~ /\S/ ) {

                # QIIME's validate_mapping_file seems to already check this:
                die
"In mapping file the line $. does not have four tab-separated fields.";
            }
        }
    }
    my $nSamples = $extctrl + $pcrpos + $pcrneg + $null + $projSamples;
    close MAP;

    print $logTee "NO. SAMPLES: $nSamples\nNO. NULLS: $null\n"
      . "NO. EXTRACTION NEGATIVE CONTROLS: $extctrl\n"
      . "NO. PCR POSITIVE CONTROLS: $pcrpos\nNO. PCR NEGATIVE CONTROLS: $pcrneg\n";

    if ( !@dbg ) {
        cacheChecksums( [$map], "map", ["map"] );
        cacheChecksums(
            [
                $mappingError,
                "$error_log/"
                  . File::Basename::basename( $map, ".txt" )
                  . "_corrected.txt"
            ],
            "meta"
        );
    }
} else {

    # It is possible the first run had a map error, but the user failed to fix.
    # skippable would return true, but the log file would still inform us of
    # the error.
    my $mappingError = glob("$error_log/*.log");

    open MAPERROR, "<", $mappingError
      or die "Cannot open $mappingError for " . "reading: $OS_ERROR";
    $_ = <MAPERROR>;
    close MAPERROR;
    chomp;
    if ( !$_ =~ /No errors or warnings found in mapping file./ ) {
        die "***Unresolved error in mapping file. See $mappingError for "
          . "details. Exiting.\n";
    }
    print $logTee "Mapping file has been validated already. Moving on.\n";
}

###### BEGIN BARCODES ##########
#######################################

if ( ( !@dbg ) || grep( /^barcodes$/, @dbg ) ) {
    my $barcodes = "$wd/barcodes.fastq";
    my $nSamples = count_samples($map);

    ## change to full path to full barcodes (flag)
    my $count = 0;

    my @files;
    my $ids;
    if ($oneStep) {
        print $logTee "Forward and reverse reads will be obtained from:\n"
          . "\t$readsForInput\n\t$readsRevInput\n";
        @files = ( $readsForInput, $readsRevInput );
        $ids   = [ "RF", "RR" ];
    } else {
        print $logTee "Forward and reverse reads will be obtained from:\n"
          . "\t$readsForInput\n\t$readsRevInput\n";
        print $logTee "Index 1 and Index 2 will be obtained from:\n"
          . "\t$index1Input\n\t$index2Input\n";
        @files = ( $readsForInput, $readsRevInput, $index1Input, $index2Input );
        $ids   = [ "RF", "RR", "I1", "I2" ];
    }

    my $input = [ $readsForInput, $readsRevInput, $index1Input, $index2Input ];
    my $inputNames = [ "RF", "RR", "I1", "I2" ];
    my $output     = [$barcodes];
    if ( !$noSkip ) {
        $skip =
          skippable( $input, $output, $checksums{"raw"}, $checksums{"barcodes"},
            "barcode extraction", $inputNames );
    }

    my $step1;

    if ( !$skip ) {

        my $start = time;
        $step1 = "extract_barcodes.py";

        # Get index files as .fastq
        my @cmds = ();
        if ( $index1Input !~ /.gz$/ ) {

            # if the index files aren't .gz, just read in place.
        } else {    # Otherwise...
                    # decompress to our project/run directory
            push( @cmds,
"gzip --decompress --force < $index1Input > $localNames{\"index1\"}"
            );
            print $logTee
"---Decompressing $project barcode and index file from $index1Input to $localNames{\"index1\"}\n";
        }
        if ( $index2Input !~ /.gz$/ ) {    # same for reverse reads
        } else {
            push( @cmds,
"gzip --decompress --force < $index2Input > $localNames{\"index2\"}"
            );
            print $logTee
"---Decompressing $project barcode and index file from $index2Input to $localNames{\"index2\"}\n";
        }

        # Execute the commands queued above
        if ( scalar @cmds > 0 ) {
            execute_and_log( @cmds, 0, $dryRun,
"Decompressing raw files containing indexes to run directory...\n"
            );
            print "---Raw index file(s) decompressed.\n";
            @cmds = ();
        }

        # Sanity check
        my $rawLnsFwd = count_lines( $localNames{"index1"} );
        my $rawLnsRev = count_lines( $localNames{"index2"} );
        if ( $rawLnsFwd == $rawLnsRev ) {
            print $logTee ( $rawLnsFwd / 4 . " barcoded reads to process.\n" );
        } else {
            print $logTee
              "Warning: raw index files have differing numbers of lines.\n";
        }

        my $mapOpt = "";
        my $bcLen  = 8;    # 2-step pcr
        if ($oneStep) {

            # Was in original code, but maybe unnecessary?
            $mapOpt = "-m $map";

            $bcLen = 12;    # 1-step pcr
        }
        my $cmd =
"extract_barcodes.py -f $localNames{\"index1\"} -r $localNames{\"index2\"} -c barcode_paired_end --bc1_len $bcLen --bc2_len $bcLen $mapOpt -o $wd";

        execute_and_log( $cmd, $logTee, $dryRun,
            "Waiting for barcode extraction to complete...\n" );

        # Sanity check
        my $barcodeLns = count_lines("$barcodes");
        if ( $rawLnsFwd != $barcodeLns ) {
            die
"Forward index file had $rawLnsFwd lines but barcodes.fastq had $barcodeLns lines.";
        }

        my $duration = time - $start;
        print $logTee "---Barcode extraction complete... "
          . ( $barcodeLns / 4 )
          . " barcodes extracted\n";
        print $logTee "---Duration of barcode extraction: $duration s\n";

        # Record hash of barcodes.fastq, and inputs
        if ( !@dbg ) {
            cacheChecksums( \@files, "raw", $ids );
            cacheChecksums( [$barcodes], "barcodes" );
        }

        # Delete unneeded output of barcodes.py
        @cmds = ();
        push @cmds, "rm -rf $wd/reads1.fastq";
        push @cmds, "rm -rf $wd/reads2.fastq";

        # If two-step, delete the temporarily decompressed index files now
        if ( !$oneStep ) {
            if ( $index1Input =~ /.gz$/ ) {
                push @cmds, "rm -rf $localNames{\"index1\"}";
            }
            if ( $index2Input =~ /.gz$/ ) {
                push @cmds, "rm -rf $localNames{\"index2\"}";
            }
        }
        execute_and_log( @cmds, 0, $dryRun,
            "Cleaning up after extract_barcodes.py...\n" );

        if ( @dbg && !grep( /^demux$/, @dbg ) ) {
            print $logTee "Finished extracting barcodes. Terminated "
              . "because --debug demux was not specified.\n";
            exit 0;
        }
    } else {
        print $logTee "Barcodes already extracted as "
          . File::Basename::basename($barcodes)
          . ". Moving on.\n";
    }
}

###### BEGIN SPLIT LIBRARIES ##########
#######################################
my $newSamNo;
if ( !@dbg || grep( /^demux$/, @dbg ) ) {
    my $barcodes = "$wd/barcodes.fastq";
    my $nSamples = count_samples($map);
    my $step2;
    my $step3;
    my $split_log = "$fwdProjDir/split_library_log.txt";

    my @cmds;
    my $rForSeqsFq = "$fwdProjDir/seqs.fastq";
    my $rRevSeqsFq = "$revProjDir/seqs.fastq";

    my ( $readsForInput, $readsRevInput, $index1Input, $index2Input ) =
        $inDir ? find_raw_files( $inDir, $oneStep, $logFH )
      : $oneStep ? ( $r1file, $r2file, $r1file, $r2file )
      :            ( $r1file, $r2file, $i1file, $i2file );

    # split_libraries_fastq.py accepts fastq.gz, so no need to
    # convert_to_local_if_gz... Unless it's one-step, in which case we can save
    # time by using the already-decompressed raw files.
    if ($oneStep) {
        ( $readsForInput, $readsRevInput ) =
          convert_to_local_if_gz( $wd, $readsForInput, $readsRevInput );
    }

    if ( !$noSkip ) {
        $skip = skippable(
            [ $readsForInput, $readsRevInput, $map, $barcodes ],
            [ $rForSeqsFq,    $rRevSeqsFq ],
            {
                "RF" => $checksums{"raw"}->{"RF"},
                "RR" => $checksums{"raw"}->{"RR"},
                %{ $checksums{"map"} },
                %{ $checksums{"barcodes"} }
            },
            $checksums{"library"},
            "demultiplexing",
            [ "RF", "RR", "map", File::Spec->abs2rel( $barcodes, $wd ) ]
        );
    }

    if ( !$skip ) {
        my $start = time;
        my $step2 = "split_libraries_fastq.py";

        ## including the stitch script before so that demultiplex happens by barcode NOT order
        my $barcodeType;
        if ($oneStep) {
            $barcodeType = "24";
        } else {
            $barcodeType = "16";
        }

        @cmds = ();
        push @cmds,
"qsub -cwd -b y -l mem_free=1G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log $step2 -i $readsForInput -o $fwdProjDir -b $barcodes -m $map --max_barcode_errors 1 --store_demultiplexed_fastq --barcode_type $barcodeType -r 999 -n 999 -q 0 -p 0.0001";
        push @cmds,
"qsub -cwd -b y -l mem_free=1G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log $step2 -i $readsRevInput -o $revProjDir -b $barcodes -m $map --max_barcode_errors 1 --store_demultiplexed_fastq --barcode_type $barcodeType -r 999 -n 999 -q 0 -p 0.0001";

        execute_and_log( @cmds, $logTee, $dryRun,
            "Demultiplexing to get the project library...\n" );

        print "---Waiting for fwd and rev seqs.fastq to complete....\n";
        print "---Monitoring $step2 error logs....\n";
        check_error_log( $error_log, $step2 );
        while ( !( -e $rForSeqsFq ) || !( -e $rRevSeqsFq ) ) {
            sleep 1;

            # check_error_log allows pipeline to terminate if the qsubbed
            # split_libraries_fastq.py prints error
            check_error_log( $error_log, $step2 );
        }
        my $duration = time - $start;
        print $logTee
          "---Duration of fwd and rev seqs.fastq production: $duration s\n";

        ##Check split_library log for 0's
        my @split = readSplitLog( $split_log, 1 );

        $newSamNo = $nSamples - scalar @split;
        if ( scalar @split > 0 && scalar @split ne $nSamples ) {
            print $logTee "---The following "
              . scalar @split
              . " samples returned 0 "
              . "reads, due to either a lack of barcodes or a filtering of those "
              . "reads in split_libraries_fastq.py:\n";
            foreach my $x (@split) {
                print $logTee "   $x\n";
            }
            print $logTee
              "---Number of samples after split_libraries_fastq.py: "
              . "$newSamNo\n";
        } elsif ( scalar @split == $nSamples ) {
            print $logTee
              "---No samples were successfully demultiplexed. Is the "
              . "mapping file correct? Exiting.\n";
            die;
        } elsif ( scalar @split == 0 ) {
            print $logTee "---Reads from all samples were demultiplexed.\n";
        }

        ###### BEGIN FASTQC ON SEQS.FASTQ #####
        #######################################
        # Replace this with calls to execute_and_log after merging in master
        my $cmd =
"qsub -cwd -b y -l mem_free=300M -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log fastqc --limits $pipelineDir/ext/fastqc/limits.txt --outdir $fwdProjDir $fwdProjDir/seqs.fastq";
        print "\tcmd=$cmd\n" if $verbose;
        system($cmd) == 0 
          or die "system($cmd) failed with exit code: $?" 
          if !$dryRun;
        $cmd =
"qsub -cwd -b y -l mem_free=300M -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log fastqc --limits $pipelineDir/ext/fastqc/limits.txt --outdir $revProjDir $revProjDir/seqs.fastq";
        print "\tcmd=$cmd\n" if $verbose;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;
        
        # Remove temporarily decompressed files
        if ($oneStep) {
            print;
            my @cmds = ();
            push( @cmds, "rm -rf $readsForInput" );
            push( @cmds, "rm -rf $readsRevInput" );
            execute_and_log( @cmds, 0, $dryRun,
                "---Removing decompressed raw files from run directory\n" );
        }

        if ( !@dbg ) {
            cacheChecksums( [ $rForSeqsFq, $rRevSeqsFq ], "library" );
        }

    } else {
        print $logTee "Library already produced as "
          . File::Basename::basename($rForSeqsFq) . " and "
          . File::Basename::basename($rRevSeqsFq)
          . ". Moving on.\n";
        my @split = readSplitLog( $split_log, 0 );

        # Still need number of demuxed samples for later, even though we didn't actually demux this time
        $newSamNo = $nSamples - scalar @split;
    }

    if ( @dbg && !grep( /^splitsamples$/, @dbg ) ) {
        die "Finished demultiplexing libaries. Terminated "
          . "because -d splitsamples was not specified.";
    }
}

if ( !@dbg || grep( /^splitsamples$/, @dbg ) ) {

    ###### BEGIN SPLIT BY SAMPLE ##########
    #######################################
    my $rForSeqsFq = "$fwdProjDir/seqs.fastq";
    my $rRevSeqsFq = "$revProjDir/seqs.fastq";

    my $n_fq1 = 0;

    my @fwdOutput = glob("$fwdSampleDir/*.fastq");
    my @revOutput = glob("$revSampleDir/*.fastq");

    if ( !$noSkip ) {
        $skip = skippable(
            [ $rForSeqsFq, $rRevSeqsFq ],
            [ @fwdOutput,  @revOutput ],
            $checksums{"library"}, $checksums{"samples"}, "sample splitting"
        );
    }

    if ( !$skip )
    {
        my $step3 = "split_sequence_file_on_sample_ids.py";
        my @cmds  = ();
        while ( !( -e $rForSeqsFq ) || !( -e $rRevSeqsFq ) ) { sleep 1; }

        push @cmds,
"qsub -cwd -b y -l mem_free=4G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log -V $step3 -i $rForSeqsFq --file_type fastq -o $fwdSampleDir";
        push @cmds,
"qsub -cwd -b y -l mem_free=4G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log -V $step3 -i $rRevSeqsFq --file_type fastq -o $revSampleDir";
        execute_and_log( @cmds, $logTee, $dryRun,
            "Splitting $project seqs.fastq " . "files by sample ID\n" );

        ## the $nSamples needs to be altered if a sample has 0 reads, because the sample-specific fastq won't be produced
        my $n_fq   = 0;
        my $nLines = 0;
        my @fwdOutput;
        my @revOutput;

        # Count the number of lines in fwdSplit/seqs.fastq
        open my $FWD, "<$fwdProjDir/seqs.fastq";
        while (<$FWD>) { }
        my $fwdLines = $.;
        close $FWD;

        # if we know how many new files there are supposed to be, wait for them all to appear
        if ( defined $newSamNo ) {
            # Count the number of files in the directory
            while ( $n_fq != $newSamNo ) {
                $n_fq      = 0;
                @fwdOutput = glob("$fwdSampleDir/*.fastq");
                $n_fq = scalar @fwdOutput;
                check_error_log( $error_log, $step3 );
            }
        }
        print "---All samples ($n_fq) accounted for in $fwdSampleDir\n";

        # Now check that all the reads are still contained in the sample-specific FASTQ's
        while ( $nLines != $fwdLines ) {
            $nLines    = 0;
            @fwdOutput = glob("$fwdSampleDir/*.fastq");
            if (@fwdOutput) {
                foreach my $file (@fwdOutput) {
                    # Also keep a running total of lines in the
                    # split_by_sample_out files
                    open SAMPLE, "<$file";
                    while (<SAMPLE>) { }
                    $nLines += $.;
                    close SAMPLE;
                }
            }
            check_error_log( $error_log, $step3 );
        }

        print "---All reads (@{[$nLines / 4]}) accounted for in $fwdSampleDir\n";

        $n_fq   = 0;
        $nLines = 0;

        # Count the number of lines in R4split/seqs.fastq
        open my $REV, "<$revProjDir/seqs.fastq";
        while (<$REV>) { }
        my $revLines = $.;
        close $REV;

        # if we know how many new files there are supposed to be, wait for them all to appear
        if ( defined $newSamNo ) {
            $n_fq++;    # Count the number of files in the directory
            while ( $n_fq != $newSamNo ) {
                $n_fq      = 0;
                @revOutput = glob("$revSampleDir/*.fastq");
                $n_fq = scalar @revOutput;
                check_error_log( $error_log, $step3 );
            }
        }

        # Now check that all the reads are still contained in the sample-specific FASTQ's
        while ( $nLines != $revLines ) {
            $nLines    = 0;
            @revOutput = glob("$revSampleDir/*.fastq");
            if (@revOutput) {
                foreach my $file (@revOutput) {
                    # Also keep a running total of lines in the
                    # split_by_sample_out files
                    open SAMPLE, "<$file";
                    while (<SAMPLE>) { }
                    $nLines += $.;
                    close SAMPLE;
                }
            }
            check_error_log( $error_log, $step3 );
        }

        if ( !@dbg ) {
            cacheChecksums( [ @fwdOutput, @revOutput ], "samples" );
        }
        print
          "--All samples ($n_fq) and reads (@{[$nLines / 4]}) accounted for"
          . " in $revSampleDir\n";
    } else {
        print $logTee
"Library already split to $newSamNo sample-specific files. Moving on.\n";
    }

    # Remove temporarily decompressed files
    if ($oneStep) {
        my @cmds = ();
        push( @cmds, "rm -rf $readsForInput" );
        push( @cmds, "rm -rf $readsRevInput" );
        execute_and_log( @cmds, 0, $dryRun,
            "---Removing decompressed raw files from $wd\n" );
    }

    if ( @dbg && !grep( /^tagclean$/, @dbg ) ) {
        print $logTee
"Finished splitting library to sample-specific FASTQs. Terminated "
          . "because --debug tagclean was not specified.\n";
        exit 0;
    }
}

###### BEGIN TAGCLEANING ##########
###################################

my $start = time;
if ( !@dbg || grep( /^tagclean$/, @dbg ) ) {
    my $do       = 1;
    my @inputsF  = glob("$fwdSampleDir/*.fastq");
    my @inputsR  = glob("$revSampleDir/*.fastq");
    my @outPatts = ( "_R1_tc", "_R2_tc" );

    if ( scalar @inputsF != scalar @inputsR ) {
        die "Unequal numbers of forward and reverse sample-specific FASTQ's.\n";
    }
    my $newSamNo = scalar @inputsF;
    
    my @fwdTcFiles = glob("$wd/*R1_tc.fastq");
    my @revTcFiles = glob("$wd/*R2_tc.fastq");

    # my @forFilenames = glob("$fwdSampleDir/*.fastq");
    # my @revFilenames = glob("$revSampleDir/*.fastq");

    if ( !$noSkip ) {
        $skip = skippable(
            [
                (
                    glob("$fwdSampleDir/*.fastq"), glob("$revSampleDir/*.fastq")
                )
            ],
            [ @fwdTcFiles, @revTcFiles ],
            $checksums{"samples"},
            $checksums{"tagcleaned"},
            "primer removal"
        );
    }

    my @cmds;
    my $cmd;
    if ( !$skip )
    {
        if ($oneStep) {
            if ( $var eq "V3V4" ) {
                my $filename;
                opendir R1, $fwdSampleDir
                  or die "Cannot open directory $fwdSampleDir\n";
                while ( $filename = readdir R1 ) {
                    if ( $filename =~ /.fastq$/ ) {
                        my @suffixes = ( ".fastq", ".fq" );
                        my $Prefix =
                          File::Basename::basename( $filename, @suffixes );
                        my $tc = "$wd/$Prefix" . "_R1_tc";
                        $cmd =
"qsub -cwd -b y -l mem_free=400M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $fwdSampleDir/$filename -out $tc -line_width 0 -verbose -tag5 GGACTACHVGGGTWTCTAAT -mm5 2 -trim_within 50";
                        push @cmds, $cmd;
                    }
                }
                close R1;

                opendir R4, $revSampleDir
                  or die "Cannot open directory $revSampleDir\n";
                while ( $filename = readdir R4 ) {
                    if ( $filename =~ /.fastq$/ ) {

                        my @suffixes = ( ".fastq", ".fq" );
                        my $Prefix =
                          File::Basename::basename( $filename, @suffixes );
                        my $tc = "$wd/$Prefix" . "_R2_tc";
                        $cmd =
"qsub -cwd -b y -l mem_free=400M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $revSampleDir/$filename -out $tc -line_width 0 -verbose -tag5 ACTCCTACGGGAGGCAGCAG -mm5 2 -trim_within 50";
                        push @cmds, $cmd;
                    }
                }
                close R4;
            }

            if ( $var eq "V4" ) {
                my $filename;
                opendir R1, $fwdSampleDir
                  or die "Cannot open directory $fwdSampleDir\n";
                while ( $filename = readdir R1 ) {
                    if ( $filename =~ /.fastq$/ ) {
                        my @suffixes = ( ".fastq", ".fq" );
                        my $Prefix =
                          File::Basename::basename( $filename, @suffixes );
                        my $tc = "$wd/$Prefix" . "_R1_tc";
                        $cmd =
"qsub -cwd -b y -l mem_free=400M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $fwdSampleDir/$filename -out $tc -line_width 0 -verbose -tag5 GTGCCAGCMGCCGCGGTAA -mm5 2";
                        push @cmds, $cmd;
                    }
                }
                close R1;

                opendir R4, $revSampleDir
                  or die "Cannot open directory $revSampleDir\n";
                while ( $filename = readdir R4 ) {
                    if ( $filename =~ /.fastq$/ ) {

                        my @suffixes = ( ".fastq", ".fq" );
                        my $Prefix =
                          File::Basename::basename( $filename, @suffixes );
                        my $tc = "$wd/$Prefix" . "_R2_tc";
                        $cmd =
"qsub -cwd -b y -l mem_free=400M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $revSampleDir/$filename -out $tc -line_width 0 -verbose -tag5 ACTCCTACGGGAGGCAGCAG -mm5 2";
                        push @cmds, $cmd;
                    }
                }
                close R4;
            }
        } else {
            if ( $var eq "V3V4" ) {
                my $filename;
                opendir R1, $fwdSampleDir
                  or die "Cannot open directory $fwdSampleDir\n";
                while ( $filename = readdir R1 ) {
                    if ( $filename =~ /.fastq$/ ) {
                        my @suffixes = ( ".fastq", ".fq" );
                        my $Prefix =
                          File::Basename::basename( $filename, @suffixes );
                        my $tc = "$wd/$Prefix" . "_R1_tc";

                        $cmd =
"qsub -cwd -b y -l mem_free=400M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $fwdSampleDir/$filename -out $tc -line_width 0 -verbose -tag5 ACTCCTACGGGAGGCAGCAG -mm5 2";
                        push @cmds, $cmd;
                    }
                }
                close R1;

                opendir R4, $revSampleDir
                  or die "Cannot open directory $revSampleDir\n";
                while ( $filename = readdir R4 ) {
                    if ( $filename =~ /.fastq$/ ) {
                        my @suffixes = ( ".fastq", ".fq" );
                        my $Prefix =
                          File::Basename::basename( $filename, @suffixes );
                        my $tc = "$wd/$Prefix" . "_R2_tc";

                        $cmd =
"qsub -cwd -b y -l mem_free=400M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $revSampleDir/$filename -out $tc -line_width 0 -verbose -tag5 GGACTACHVGGGTWTCTAAT -mm5 2";
                        push @cmds, $cmd;
                    }
                }
                close R4;
            }
            if ( $var eq "V4" ) {
                my $filename;
                opendir R1, $fwdSampleDir
                  or die "Cannot open directory $fwdSampleDir\n";
                while ( $filename = readdir R1 ) {
                    if ( $filename =~ /.fastq$/ ) {
                        my @suffixes = ( ".fastq", ".fq" );
                        my $Prefix =
                          File::Basename::basename( $filename, @suffixes );
                        my $tc = "$wd/$Prefix" . "_R1_tc";
                        $cmd =
"qsub -cwd -b y -l mem_free=400M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $fwdSampleDir/$filename -out $tc -line_width 0 -verbose -tag5 GTGCCAGCMGCCGCGGTAA -mm5 2";
                        push @cmds, $cmd;

                    }
                }
                close R1;

                opendir R4, $revSampleDir
                  or die "Cannot open directory $revSampleDir\n";
                while ( $filename = readdir R4 ) {
                    if ( $filename =~ /.fastq$/ ) {

                        my @suffixes = ( ".fastq", ".fq" );
                        my $Prefix =
                          File::Basename::basename( $filename, @suffixes );
                        my $tc = "$wd/$Prefix" . "_R2_tc";
                        $cmd =
"qsub -cwd -b y -l mem_free=400M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $revSampleDir/$filename -out $tc -line_width 0 -verbose -tag5 ACTCCTACGGGAGGCAGCAG -mm5 2";
                        push @cmds, $cmd;
                    }
                }
                close R4;
            }
            if ( $var eq "ITS" ) {
                my $filename;
                opendir R1, $fwdSampleDir
                  or die "Cannot open directory $fwdSampleDir\n";
                while ( $filename = readdir R1 ) {
                    if ( $filename =~ /.fastq$/ ) {
                        my @suffixes = ( ".fastq", ".fq" );
                        my $Prefix =
                          File::Basename::basename( $filename, @suffixes );
                        my $tc = "$wd/$Prefix" . "_R1_tc";
                        $cmd =
"qsub -cwd -b y -l mem_free=400M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $fwdSampleDir/$filename -out $tc -line_width 0 -verbose -tag5 CTGCCCTTTGTACACACCGC -mm5 2";
                        push @cmds, $cmd;
                    }
                }
                close R1;

                opendir R4, $revSampleDir
                  or die "Cannot open directory $revSampleDir\n";
                while ( $filename = readdir R4 ) {
                    if ( $filename =~ /.fastq$/ ) {

                        my @suffixes = ( ".fastq", ".fq" );
                        my $Prefix =
                          File::Basename::basename( $filename, @suffixes );
                        my $tc = "$wd/$Prefix" . "_R2_tc";
                        $cmd =
"qsub -cwd -b y -l mem_free=400M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $revSampleDir/$filename -out $tc -line_width 0 -verbose -tag5 TTTCGCTGCGTTCTTCATCG -mm5 2";
                        push @cmds, $cmd;
                    }
                }
                close R4;
            }
        }
        execute_and_log( @cmds, $logTee, $dryRun,
            "Removing $var primers from all sequences.\n" );

        my @fwdFiles = glob("$wd/*R1_tc.fastq");
        my $nFiles   = @fwdFiles;
        my $equalLns = 0;
        while ( $nFiles != $newSamNo || !$equalLns ) {
            @fwdFiles = glob("$wd/*R1_tc.fastq");
            $nFiles   = @fwdFiles;

# Ensure that each tagcleaned file has the same number of lines as its input file
            $equalLns = 1;
            foreach my $file (@fwdFiles) {
                my $basename =
                  File::Basename::basename( $file, ("_R1_tc.fastq") );
                if ( count_lines($file) !=
                    count_lines("$fwdSampleDir/$basename.fastq") )
                {
                    $equalLns = 0;
                }
            }
            check_error_log( $error_log, "perl" );
        }
        print $logTee "---All tagcleaned R1 samples accounted for in $wd\n";

        $equalLns = 0;
        my @revFiles = glob("$wd/*R2_tc.fastq");
        $nFiles = @revFiles;
        while ( $nFiles != $newSamNo || !$equalLns ) {
            @revFiles = glob("$wd/*R2_tc.fastq");
            $nFiles   = @revFiles;

# Ensure that each tagcleaned file has the same number of lines as its input file
            $equalLns = 1;
            foreach my $file (@revFiles) {
                my $basename =
                  File::Basename::basename( $file, ("_R2_tc.fastq") );
                if ( count_lines($file) !=
                    count_lines("$fwdSampleDir/$basename.fastq") )
                {
                    $equalLns = 0;
                }
            }
            check_error_log( $error_log, "perl" );
        }
        print $logTee
          "---All tagcleaned R4 (R2) samples accounted for in $wd\n";

        if ( !@dbg ) {
            cacheChecksums( [ @fwdFiles, @revFiles ], "tagcleaned" );
        }

        my $duration = time - $start;
        print $logTee "---Primer sequences removed from $newSamNo samples.\n";
        print $logTee "---Duration of tagcleaning: $duration s\n";
    } else {
        print $logTee
"---Primers already removed from forward and reverse reads of $newSamNo samples. Moving on.\n";
    }

    if ( @dbg && !grep( /^dada2$/, @dbg ) ) {
        print $logTee
"Finished removing primers. Terminated because --debug dada2 was not " .
"specified.\n";
        exit 0;
    }
}

###### BEGIN DADA2 ##########
#############################
if ( ( !@dbg ) || grep( /^dada2$/, @dbg ) ) {
    my $projrt = "$wd/$project" . "_" . $run . "_dada2_part1_rTmp.R";
    my $projrtout = "$wd/$project" . "_" . $run . "_dada2_part1_rTmp.Rout";

    my @outputs = ( "$projrt", "$projrtout" );
    my @inputs = glob("$wd/*R[1|2]_tc.fastq");

    if ( !$noSkip ) {
        $skip = skippable(
            \@inputs, \@outputs,
            $checksums{"tagcleaned"},
            $checksums{"dada2part1out"}, "DADA2"
        );
    }

    if ( !$skip ) {
        my $truncLen;
        if ( $f && $r ) {
            $truncLen = "c($f,$r)";
        }

        chdir $pd;
        if ($oneStep) {
            if ( $var eq "V3V4" ) {

                if ( !$truncLen ) {
                    $truncLen = "c(255,255)";
                }
                if ( !$maxN ) {
                    $maxN = 0;
                }
                if ( !$maxEE ) {
                    $maxEE = "c(2,2)";
                }
                if ( !$truncQ ) {
                    $truncQ = 2;
                }
                if ( !$phix ) {
                    $phix = "TRUE";
                }
                if ( !$maxLen ) {
                    $maxLen = "Inf";
                }
                if ( !$minLen ) {
                    $minLen = 20;
                }
                if ( !$minQ ) {
                    $minQ = 0;
                }
                dada2(
                    $run,  $truncLen, $maxN,   $maxEE, $truncQ,
                    $phix, $maxLen,   $minLen, $minQ,  $logFH
                );
            }

            if ( $var eq "V4" ) {
                if ( !$truncLen ) {
                    $truncLen = 200;
                }
                if ( !$maxN ) {
                    $maxN = 0;
                }
                if ( !$maxEE ) {
                    $maxEE = "c(2,2)";
                }
                if ( !$truncQ ) {
                    $truncQ = 2;
                }
                if ( !$phix ) {
                    $phix = "TRUE";
                }
                if ( !$maxLen ) {
                    $maxLen = "Inf";
                }
                if ( !$minLen ) {
                    $minLen = 20;
                }
                if ( !$minQ ) {
                    $minQ = 0;
                }
                dada2(
                    $run,  $truncLen, $maxN,   $maxEE, $truncQ,
                    $phix, $maxLen,   $minLen, $minQ,  $logFH
                );
            }
        } else {
            if ( $var eq "V3V4" ) {
                if ( !$truncLen ) {
                    $truncLen = "c(255,225)";
                }
                if ( !$maxN ) {
                    $maxN = 0;
                }
                if ( !$maxEE ) {
                    $maxEE = "c(2,2)";
                }
                if ( !$truncQ ) {
                    $truncQ = 2;
                }
                if ( !$phix ) {
                    $phix = "TRUE";
                }
                if ( !$maxLen ) {
                    $maxLen = "Inf";
                }
                if ( !$minLen ) {
                    $minLen = 20;
                }
                if ( !$minQ ) {
                    $minQ = 0;
                }
                dada2(
                    $run,  $truncLen, $maxN,   $maxEE, $truncQ,
                    $phix, $maxLen,   $minLen, $minQ,  $logFH
                );
            }

            if ( $var eq "V4" ) {
                if ( !$truncLen ) {
                    $truncLen = 200;
                }
                if ( !$maxN ) {
                    $maxN = 0;
                }
                if ( !$maxEE ) {
                    $maxEE = "c(2,2)";
                }
                if ( !$truncQ ) {
                    $truncQ = 2;
                }
                if ( !$phix ) {
                    $phix = "TRUE";
                }
                if ( !$maxLen ) {
                    $maxLen = "Inf";
                }
                if ( !$minLen ) {
                    $minLen = 20;
                }
                if ( !$minQ ) {
                    $minQ = 0;
                }
                dada2(
                    $run,  $truncLen, $maxN,   $maxEE, $truncQ,
                    $phix, $maxLen,   $minLen, $minQ,  $logFH
                );
            }

            if ( $var eq "ITS" ) {
                if ( !$truncLen ) {
                    $truncLen = 0;
                }
                if ( !$maxN ) {
                    $maxN = 0;
                }
                if ( !$maxEE ) {
                    $maxEE = "1";
                }
                if ( !$truncQ ) {
                    $truncQ = 2;
                }
                if ( !$phix ) {
                    $phix = "TRUE";
                }
                if ( !$maxLen ) {
                    $maxLen = "Inf";
                }
                if ( !$minLen ) {
                    $minLen = 50;
                }
                if ( !$minQ ) {
                    $minQ = 0;
                }
                dada2(
                    $run,  $truncLen, $maxN,   $maxEE, $truncQ,
                    $phix, $maxLen,   $minLen, $minQ,  $logFH
                );
            }
        }

        # Rename DADA2 R files
        my $rt     = "$wd/dada2_part1_rTmp.R";
        my @cmds   = ();
        push @cmds, "mv $rt $projrt";
        my $rtout     = "$wd/dada2_part1_rTmp.Rout";
        push @cmds, "mv $rtout $projrtout";
        eval {
            execute_and_log( @cmds, 0, $dryRun, "Renaming DADA2 R files.\n" );
        };
        warn $@ if $@;
        print $logTee "---DADA2-specific commands with output can be found in "
          . "$projrtout\n";

###### EVALUATING DADA2 OUTPUT ##########
#########################################
        my @removed;
## do a scan of this file for those samples not passing filtering (0 sequences surviving - and log this)
        open RTOUT, "<$projrtout"
          or die "Cannot open $rtout for reading: $OS_ERROR";
        while (<RTOUT>) {
            chomp;
            if ( $_ =~ /filterAndTrim/ ) {

                #print "Splitting $_\n";
                my @filtTrim = split( /\s/, $_ );
                for my $x (@filtTrim) {

                    #print "Searching $x\n";
                    if ( $x =~ /truncLen/ ) {
                        chomp;
                        $truncLen = $x;
                    }
                    if ( $x =~ /maxN/ ) {
                        chomp;
                        $maxN = $x;
                    }
                    if ( $x =~ /maxEE/ ) {
                        chomp;
                        $maxEE = $x;
                    }
                    if ( $x =~ /truncQ/ ) {
                        chomp;
                        $truncQ = $x;
                    }
                    if ( $x =~ /phix/ ) {
                        chomp;
                        $phix = $x;
                    }
                }
            }
            if ( $_ =~ /The filter removed all reads/ ) {
                push @removed, $_;
            }
        }
        close RTOUT;

        if ( scalar @removed > 0 ) {
            print $logTee scalar @removed
              . " samples removed during dada2 filtering:\n";
            for my $x (@removed) {
                print $logTee "$x\n";
            }
        }
        print $logTee "\n";

        if ( List::Util::all { -e $_ } @outputs ) {
            print $logTee
              "\nFor $var region, dada2 used the following filtering "
              . "requirements:\n$truncLen\n$maxN\n$maxEE\n$truncQ\n$phix\n";

            # "dada2 completed successfully" is a key phrase that causes
            # an appropriate message printed to STDOUT
            print $logTee "dada2 completed successfully!\nAbundance table for "
              . "$project run $run located at $wd/dada2_abundance_table.rds\n";
            print $logTee
              "See $outputs[0] for dada2 table of reads surviving by "
              . "step\n\n";

            if ( !@dbg ) {
                cacheChecksums( \@outputs, "dada2part1out" );
            }
        } else {
            print $logTee
              "---dada2 did not complete successfully, something went wrong!\n"
              . "---Check $projrtout.\n";
            die;
        }
        
    } else {
        print $logTee
"DADA2 already processed the same tagcleaned files during a previous run. Doing nothing.\n";
    }
}

###### COMPLETING $logFH FILE ##############
#########################################

END {

    if ( defined $logTee ) {
        print $logTee "\n\n\n";
        close $logFH;
    }

}

####################################################################
##                               SUBS
####################################################################

sub read_checkpoints {
    my $filepath = shift;

    my $json;
    {
        local $/;    #Enable 'slurp' mode
        if ( -e $filepath ) {
            open my $FH, "+<$filepath";
            seek $FH, 0, 0 or die;
            $json = <$FH>;
            close $FH;
        }

    }

    my @checkpoints;
    my $data = {};
    eval {
        $data = decode_json($json);

    };

    return %{$data};
}

#' @_[0] Filenames of inputs to checksum
#' @_[1] Filenames of outputs to checksum (not guaranteed to be existing files)
#' @_[2] A hash pointer containing the checksums for inputs
#' @_[3] A hash pointer containing the checksums for outputs
#' @_[4] The name of the current step, for log messages.
#' @_[5] The names of the inputs (should match keys of @_[2])
#' @_[6] The names of the outputs (should match keys of @_[3])
sub skippable {
    my $inputsR  = shift;
    my $outputsR = shift;
    my $par3     = shift;
    my %checksumsIn;
    if ( defined $par3 ) {

        # If the checksums json didn't contain any info about this step, the
        # value of the parameter will be undef. To maintain the same interface,
        # the default value is a reference an empty hash.
        %checksumsIn = %{$par3};
    }

    # Remove any key-value pairs with undef value. This is the necessary to
    # test if we have reference checksums for all the files present.
    delete $checksumsIn{$_}
      for grep { not defined $checksumsIn{$_} } keys %checksumsIn;

    my $par4 = shift;
    my %checksumsOut;
    if ( defined $par4 ) {
        %checksumsOut = %{$par4};
    }
    delete $checksumsOut{$_}
      for grep { not defined $checksumsOut{$_} } keys %checksumsOut;

    my $prevSkip = $skip;    # file-level my-variable
    my $step     = shift;
    my $inNames  = shift;
    my $outNames = shift;
    my $wd       = $wd;      # This is a file-level my-variable

    if ( scalar @$outputsR == 0 ) {
        return 0;
    }                        # no output files found by caller
    if ( List::Util::none { defined $_ } @$outputsR ) {
        return 0;
    }
    if ( List::Util::none { -e $_ } @$outputsR ) {

        # Caller named output files but none existed. No need to print message
        return 0;
    }
    if ( !List::Util::all { -e $_ } @$outputsR ) {

        # Caller named output files but not all existed
        print $logTee
"$step output is incomplete. Removing outputs files and executing $step...\n";
        foreach (@$outputsR) {
            if ( -e $_ ) {
                unlink();
            }
        }
        return 0;
    }

    # When previously stored, the checksums should have been named consistently,
    # or named after the files themselves.
    if ( !defined $inNames ) {
        $inNames = [ map { File::Spec->abs2rel( $_, $wd ) } @$inputsR ];
    }
    if ( !defined $outNames ) {
        $outNames = [ map { File::Spec->abs2rel( $_, $wd ) } @$outputsR ];
    }

    if ( !setIdent( $inNames, [ keys %checksumsIn ] ) ) {

        # The input files are named differently from last time (or pipeline
        # provided inconsistent file nicknames)
        print $logTee
          "Input to $step does not match previous records. Discarding $step ",
          "outputs and repeating $step.\n";
        unlink @$outputsR;
        return 0;
    }
    if ( !setIdent( $outNames, [ keys %checksumsOut ] ) ) {

# The outputs files are named differently from last time (or pipeline provided inconsistent file nicknames)
        print $logTee ucfirst($step)
          . " output does not match previous records. Discarding $step outputs and repeating $step...\n";
        unlink @$outputsR;
        return 0;
    }

    print ucfirst("$step") . " output exists. Trying to validate...\n";

    if ( !$prevSkip ) {

        # If the previous step was not skipped, the inputs must be checksummed

        my $i          = -1;
        my $inputValid = List::Util::all {
            $i++;
            my @stdout   = split( /\s/, `sha512sum $_` );
            my $checksum = $stdout[0];
            if ( $checksum eq $checksumsIn{ $inNames->[$i] } ) {
                return 1;
            } else {
                return 0;
            }
        }
        @$inputsR;

        # If all inputs are/nt valid give appropriate messages
        if ($inputValid) {
            print "Input to $step matches previous records.\n";
        } else {
            print $logTee
"Input to $step has changed since last execution. Discarding $step outputs and\n",
              "repeating $step.\n";
            unlink @$outputsR;
            return 0;
        }
    }    # Otherwise, we know the inputs are valid
         # then do the same for all the outputs

    my $i           = -1;
    my $outputValid = List::Util::all {
        $i++;
        my @stdout   = split( /\s/, `sha512sum $_` );
        my $checksum = $stdout[0];
        if ( $checksum eq $checksumsOut{ $outNames->[$i] } ) {
            return 1;
        } else {
            return 0;
        }
    }
    @$outputsR;

    # If all outputs are/nt valid give appropriate messages
    if ($outputValid) {
        print "Output to $step matches previous records.\n";
        return 1;
    } else {
        print $logTee
"Output to $step does not match previous records. Discarding $step outputs and\n",
          "repeating $step.\n";
        unlink @$outputsR;
        return 0;
    }

    return 0;
}

#' @_[0] A pointer to an array containing filenames to checksum
#' @_[1] The name of the step whose products are being cached
#' @_[2] Optional pointer to array of names (include to avoid IDing files by filenames)
sub cacheChecksums {
    my $files = shift;
    my $step  = shift;
    my $names = shift;
    my $wd    = $wd;     # This is a file-level $my variable

    if ( !defined $names ) {
        $names = [ map { File::Spec->abs2rel( $_, $wd ) } @$files ];
    }
    if ( defined $checksumFile ) {

        print $logTee "Saving progress in run folder...\n";

        my $destroyNextStep = sub {
            if ( $step2number{$step} + 1 < scalar keys %number2step ) {

             # Reset the next set of checksums so the next step is forced to run
                $checksums{ $number2step{ $step2number{$step} + 1 } } = {};
            }
        };
        my $update = 0;
        my @keys   = keys %{ $checksums{$step} };
        if ( !setIdent( \@keys, $names ) ) {
            $update = 1;
            $destroyNextStep->()
              ;    # MUST DO THIS BEFORE ANY NEW CHECKSUMS ARE WRITTEN
             # PREVENTS "CHIMERIC" CHECKSUM FILE IN CASE OF PIPELINE INTERRUPTION HERE
        }
        my %newChecksums;

        for my $i ( 0 .. scalar @$files - 1 ) {
            my $file     = $files->[$i];
            my $name     = $names->[$i];
            my @output   = split( /\s/, `sha512sum $file` );
            my $checksum = $output[0];
            if ( !$update )
            { # !$update means that $checksums{$step}->{ $file } is defined for all files
                if ( $checksums{$step}->{$name} ne $checksum ) {
                    $destroyNextStep->();
                }
            }

            # Store the checksum of the files just produced.
            $newChecksums{$name} = $checksum;
        }
        $checksums{$step} = \%newChecksums;

        open my $checksumFH, ">$checksumFile";
        if ( List::Util::all { defined $_ } %checksums ) {
            print $checksumFH encode_json( \%checksums );
        }
        close $checksumFH;

    }
}

sub setIdent {
    my ( $left, $right ) = @_;
    return 0 if scalar @$left != scalar @$right;
    my %hash;
    @hash{ @$left, @$right } = ();
    return scalar keys %hash == scalar @$left;
}

# @param 0 mapping file path
sub count_samples {
    my $map = shift;
    my $ans = count_lines($map) - 1;
    return $ans;
}

sub count_lines {
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
sub find_raw_files {
    my $index1Input;
    my $index2Input;
    my $readsForInput;
    my $readsRevInput;
    my ( $wd, $oneStep, $log ) = @_;

    # These are the only input files needed for 1step
    my @r1s = glob("$wd/*R1.fastq $wd/*R1.fastq.gz");
    my @r2s = glob("$wd/*R2.fastq $wd/*R2.fastq.gz");

    if ($oneStep) {
        if (   scalar @r1s == 1
            && scalar @r2s == 1 )
        {
            $index1Input = $readsForInput = $r1s[0];
            $index2Input = $readsRevInput = $r2s[0];
        } else {
            print $log "Couldn't find input files in $wd.\n";
            print $log "Files found:\n";

            my $printme = join "\n", @r1s;
            print $log "$printme\n" if $printme;
            $printme = join "\n", @r2s;
            print $log "$printme\n" if $printme;

            die
"Could not find a complete and exclusive set of raw files. Since --1step given, input directory must"
              . " have exactly one R1 and R2 file. See pipeline log for a list of the files"
              . " found.";
        }
    } else {
        my @i1s = glob("$wd/*I1.fastq $wd/*I1.fastq.gz");
        my @i2s = glob("$wd/*I2.fastq $wd/*I2.fastq.gz");

        # also globbed for r1 and r2 files above
        my @r3s = glob("$wd/*R3.fastq $wd/*R3.fastq.gz");
        my @r4s = glob("$wd/*R4.fastq $wd/*R4.fastq.gz");

        if (   scalar @i1s == 1
            && scalar @i2s == 1
            && scalar @r1s == 1
            && scalar @r2s == 1
            && scalar @r3s == 0
            && scalar @r4s == 0 )
        {
            $index1Input   = $i1s[0];
            $index2Input   = $i2s[0];
            $readsForInput = $r1s[0];
            $readsRevInput = $r2s[0];
        } elsif ( scalar @i1s == 0
            && scalar @i2s == 0
            && scalar @r1s == 1
            && scalar @r2s == 1
            && scalar @r3s == 1
            && scalar @r4s == 1 )
        {
            $index1Input   = $r2s[0];
            $index2Input   = $r3s[0];
            $readsForInput = $r1s[0];
            $readsRevInput = $r4s[0];
        } else {
            print $log "Couldn't find input files.\n";
            print $log "Files found in $wd:\n";
            my $printme = join "\n", @i1s;
            print $log "$printme\n" if $printme;
            $printme = join "\n", @i2s;
            print $log "$printme\n" if $printme;
            $printme = join "\n", @r1s;
            print $log "$printme\n" if $printme;
            $printme = join "\n", @r2s;
            print $log "$printme\n" if $printme;
            $printme = join "\n", @r3s;
            print $log "$printme\n" if $printme;
            $printme = join "\n", @r4s;
            print $log "$printme\n" if $printme;
            die
"Could not find a complete and exclusive set of raw files. Input directory must"
              . " contain exactly one set of I1, I2, R1, and R2, OR R1, R2, R3, and R4.\n"
              . "See pipeline log for a list of the files found.";
        }
    }
    return ( $readsForInput, $readsRevInput, $index1Input, $index2Input );
}

# Given *R1.fastq(.gz), if it's .gz, then removes the extension and gives a
# filepath in the local directory
# @param 0 The project/run directory
# @params 1... The full path to the original file
sub convert_to_local_if_gz {
    my $wd    = shift;
    my @files = @_;
    my @ans;
    foreach my $file (@files) {
        if ( $file =~ /.gz$/ ) {

            # Rename *[R|I][1|2].fastq.gz to <WD>/<PROJ>_<RUN>_[R|I][1|2].fastq
            my $dest = File::Basename::basename($file);

            my $suffix = substr( $dest, ( length($dest) - 11 ), 8 )
              ;    # get suffix (sans .gz)

            my @dirs = File::Spec->splitdir($wd);
            $file =
                "$wd/$dirs[scalar(@dirs) - 2]" . "_"
              . "$dirs[scalar(@dirs) - 1]"
              . "_$suffix";
        }
        push( @ans, $file );
    }
    return @ans;
}

sub readSplitLog {
    my $split_log = shift;
    my $verbose   = shift;
    ##Check split_library log for 0's
    open SPLIT, "<$split_log"
      or die "Cannot open $split_log for writing: " . "$OS_ERROR";
    my @split;
    while (<SPLIT>) {
        if ( $_ =~ /\t/ ) {
            chomp;
            my ( $sample, $nReads ) = split, /\t/;
            chomp $nReads;
            if ( $nReads eq "0" ) {
                push @split, $sample;
            } else {
                if ($verbose) {
                    print $logTee "$_\n";
                }
            }
        }
    }
    close SPLIT;
    return @split;
}

sub dada2 {
    my $run      = shift;
    my $truncLen = shift;
    my $maxN     = shift;
    my $maxEE    = shift;
    my $truncQ   = shift;
    my $phix     = shift;
    my $maxLen   = shift;
    my $minLen   = shift;
    my $minQ     = shift;

    my $Rscript = qq~

  library("dada2")
  packageVersion("dada2")
  cwd<-getwd()


  ## perform filtering and trimming
  path <- cwd
  filtpath <- file.path(path, "filtered")
  fastqFs <- sort(list.files(path, pattern="R1_tc.fastq"))
  fastqRs <- sort(list.files(path, pattern="R2_tc.fastq"))
  sample.names <- sapply(strsplit(basename(fastqFs), "_"), `[`, 1)
  filtFs<-(paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs<-(paste0(sample.names, "_R_filt.fastq.gz"))
  filtFs_files<-file.path(filtpath,filtFs)
  filtRs_files<-file.path(filtpath,filtRs)
  if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
  out<-filterAndTrim(fwd=file.path(path,fastqFs), filt=filtFs_files, rev=file.path(path,fastqRs), filt.rev=filtRs_files, truncLen=$truncLen, maxN=$maxN, maxEE=$maxEE, truncQ=$truncQ, rm.phix=$phix, compress=TRUE, multithread=TRUE, verbose=TRUE, matchIDs=TRUE)

## use this to produce the first 5 files to check the quality of the samples.
  ##make post-trimmed quality figures
  ##postscript("filt_forward_reads_quality.eps")
  ##plotQualityProfile(filtFs) ## can take some time if being produced for all samples
  ##dev.off()
  ##if you specify a range it makes only that number of plots, otherwise all plots are produced in grid fashion
  ##postscript("filt_reverse_reads_quality.eps")
  ##plotQualityProfile(filtRs)

  ## Learn errors
  filtFs <- list.files(filtpath, pattern="_F_filt.fastq.gz", full.names = TRUE)
  filtRs <- list.files(filtpath, pattern="_R_filt.fastq.gz", full.names = TRUE)
  sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
  sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)
  if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
  names(filtFs) <- sample.names
  names(filtRs) <- sample.namesR
  set.seed(100)
  # Learn forward error rates
  errF <- learnErrors(filtFs, nread=1e6, multithread=TRUE)
  # Learn reverse error rates
  errR <- learnErrors(filtRs, nread=1e6, multithread=TRUE)
  # Sample inference and merger of paired-end reads
  mergers <- vector("list", length(sample.names))
  names(mergers) <- sample.names
  for(sam in sample.names) {
    cat("Processing:", sam, "\n")
      derepF <- derepFastq(filtFs[[sam]])
      ddF <- dada(derepF, err=errF, multithread=TRUE)
      derepR <- derepFastq(filtRs[[sam]])
      ddR <- dada(derepR, err=errR, multithread=TRUE)
      merger <- mergePairs(ddF, derepF, ddR, derepR)
      mergers[[sam]] <- merger
  }

  rm(derepF); rm(derepR)

  ## Make sequence abundance table 
  seqtab <- makeSequenceTable(mergers)
  saveRDS(seqtab, "dada2_abundance_table.rds")

  getN <- function(x) sum(getUniques(x))
  ## track <- cbind(out, rowSums(seqtab))
  v<-rowSums(seqtab)
  v0<-numeric(nrow(out))
  track<-cbind(out, v0)
  rownames(track)<-Map(function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][1], rownames(track))
  track[names(v),3]<-v
  colnames(track) <- c("input", "filtered", "merged")
  write.table(track, "dada2_part1_stats.txt", quote=FALSE, append=FALSE, sep=\t, row.names=TRUE, col.names=TRUE)
  ~;
    run_R_script( $Rscript, $logFH );
}

sub run_R_script {
    my $Rscript = shift;
    my $logFH   = shift;
    chdir $wd;

    my $outFile = "dada2_part1_rTmp.R";
    open OUT, ">$outFile", or die "cannot write to $outFile: $!\n";
    print OUT "$Rscript";
    close OUT;

    my $exitStatus = 1;

    while ( $exitStatus == 1 ) {
        my $outR = $outFile . "out";
        my $cmd =
"rm -rf $wd/filtered $outR $wd/dada2_part1_stats.txt $wd/dada2_abundance_table.rds $wd/.RData";
        execute_and_log( $cmd, 0, $dryRun,
"---Removing old filtered fastq files, stats, and Rout files from previous runs\n"
        );

        $cmd =
"qsub -cwd -b y -l mem_free=1G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log -V $R CMD BATCH $outFile";
        execute_and_log( $cmd, $logTee, $dryRun,
            "Running DADA2 with fastq files in $wd for $var region...\n" );

        while ( !-e $outR ) {
            check_error_log( $error_log, "R" );
        }

        # Until DADA2 succeeds, look for the R output file and monitor for
        # signs of termination
        if ( -e $outR ) {
            my $decided = 0;    # flag to stop re-reading .Rout
            while ( !$decided ) {
                open IN, "<$outR"
                  or die "Cannot open $outR for reading: $OS_ERROR\n";
                my $line = <IN>;
                while ( $line && !$decided ) {
                    if (        # signs of bad termination
                        (
                               $line =~ /Error in/
                            || $line =~ /Execution halted/
                            || $line =~ /encountered errors/
                            || $line =~ /Error :/
                        )
                        && $line !~ /learnErrors/ # False positives; don't match
                        && $line !~ /error rates/
                        && $line !~ /errors/
                      )
                    {
                        $decided =
                          1;    # Leave $exitStatus = 1 so DADA2 is restarted
                        print $logTee "---R script crashed at: $line\n";

                     # Preserve the last R log file that errored. Get rid of the
                     # old R output file, then run R again.
                        File::Copy::move( "$outR", "$outR.old" );
                        print $logTee "---See $outR.old for details.\n";
                        print $logTee "---Attempting to restart R...\n";

                    } elsif ( $line =~ /proc.time()/ ) {    # sign of success

                        print $logTee "---R script completed without errors."
                          if $verbose;
                        $exitStatus = 0;                    # Move on from DADA2
                        $decided    = 1;                    # Stop reading .Rout

                    }
                    $line = <IN>;
                }
                close IN;
            }
        }
    }
}

sub source {
    my $name = shift;

    open my $fh, "<$name" or die "could not open $name: $!";

    while (<$fh>) {
        chomp;
        my ( $k, $v ) = split /=/, $_, 2;
        $v =~ s/^(['"])(.*)\1/$2/;            #' fix highlighter
        $v =~ s/\$([a-zA-Z]\w*)/$ENV{$1}/g;
        $v =~ s/`(.*?)`/`$1`/ge;              #dangerous
        $ENV{$k} = $v;
    }
}

# Check for qsub errors.
# @_[0] The error log directory
# @_[1] Prefix of error files to search.
sub check_error_log {
    my $error_log = shift;
    my $step      = shift;
    my $error;
    ## make $error array so that goes through all errors for that step.
  RECHECK:
    $error = glob("$error_log/$step.e*");
    if ( defined $error ) {
        open ERROR, "<$error" or die "Can't open $error.\n";
        while (<ERROR>) {
            if ( $_ =~ /error/ || $_ =~ /Error/ || $_ =~ /ERROR/ ) {
                print "Error in $step. See $error.\n";
                print $logTee "Error in $step. See $error.\n";
                seek( ERROR, 0, 0 );
                my $errorMessage = do { local $/; <ERROR> };
                die $errorMessage;
            }
        }
    } else {
        ##print "---Re-checking for $step error file\n";
        goto RECHECK;
    }

    close ERROR;
}

################################################################################
# Execute the given commands;
# 1st pop: Human-friendly message to be printed to STDOUT.
# 2nd pop: $dryRun parameter
# 3rd pop: Lexical filehandle that overrides printing to STDOUT. Pass 0 for no override.
# Remaining args: cmds to be run through system()
sub execute_and_log {
    my $cuteMsg   = pop @_;
    my $dryRun    = pop @_;
    my $lexicalFH = pop @_;

    if ( $lexicalFH ) { # print the human-friendly message once
        print $lexicalFH "$cuteMsg\n";    # Print to whatever handle was provided
    } else {
        print "$cuteMsg\n"; 
    }

    # CAREFUL, $cmd holds a reference to a variable in the caller!
    foreach my $cmd (@_) {
        if ( $lexicalFH && $verbose ) { # print each command
            print $lexicalFH "\$ $cmd\n"; 
        } elsif ( !$lexicalFH && $verbose ) {
            print "\$ $cmd\n";
        }

        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;
    }
}

exit 0;

