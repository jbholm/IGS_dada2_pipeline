#!/usr/bin/env perl
#. /usr/local/packages/usepackage/share/usepackage/use.bsh
use File::Basename;
my $pipelineDir = dirname(__FILE__);

=head1 NAME

illumina_dada2.pl

=head1 DESCRIPTION

The script can be launched from any location on the IGS server, it automatically 
produces a directory in /local/groupshare/ravel named after the project and run
ID provided. 

Beginning with the path to four raw Illumina sequencing files (R1, R2, R3, R4),
a mapping file, a project ID, a run ID, and specifying the targeted variable
region, this script:
  1. Produces individual .fastq files for each sample listed in the mapping file
  2. Performs tag-cleaning of each file
  3. Runs the R1 (forward) and R4 (later called R2, reverse) files through the dada2 pipeline for either the 16S rRNA gene V3V4 or V4 regions.

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
  -r1 <path_to_R1_file> -r2 <path_to_I1_file> -r3 <path_to_I2_file>
  -r4 <path_to_R2_file> -p <project name> -r <run id> -v <variable region> 
  -m <full path to mapping file> -sd <storage directory>

FOR 1-STEP
qsub -cwd -b y -l mem_free=1G -P jravel-lab -q threaded.q -pe thread 4 -V
  -e <path_to_logs> -o <path_to_logs> illumina_dada2.pl -i <input directory>
  -p <project name> -r <run ID> -m <mapping file> -v <variable region> 
  -sd <storage directory> --1Step

qsub -cwd -b y -l mem_free=1G -P jravel-lab -q threaded.q -pe thread 4 -V
  -e <path_to_logs> -o <path_to_logs> illumina_dada2.pl -r1 <path_to_R1_file>
  -r2 <path_to_R2_file> -p <project name> -r <run ID> -m <mapping file>
  -v <variable region> -sd <storage directory> --1Step

=head1 OPTIONS

=over

=item B<--raw-path>=path, B<-i> path

Single full path to directory containing raw R1, R2, R3, and R4 files, or R1, 
R2, I1, and I2 files. Incompatible with -r1, -r2, -r3, and -r4.

=item B<--r1-path>=file, B<-r1> file

Full path to raw R1 read file (forward read file, or r1) (.fastq.gz). 
Incompatible with b<--raw-path, -i>.

=item B<--r2-path>=file, B<-r2> file

Full path to raw R2 read file (barcode file, or i1) (.fastq.gz). Incompatible 
with b<--raw-path, -i>.

=item B<--r3-path>=file, B<-r3> file

Full path to raw R3 read file (barcode file, or i2) (.fastq.gz). Incompatible 
with b<--raw-path, -i>. Ignored when --1Step given.

=item B<--r4-path>=file, B<-r4> file

Full path to raw R4 read file (reverse read file, r2 or r4) (.fastq.gz). 
Incompatible with b<--raw-path, -i>. Ignored when --1Step given.

=item B<--project-name>=name, B<-p> name

Create the project folder with this name.

=item B<--run-ID>=name, B<-r> name

Create the run folder with this name.

=item B<--map>=file, B<-m> file

The full path to the Qiime-formatted mapping file.

=item B<--var-reg>={V3V4, V4, ITS}, B<-v> {V3V4, V4, ITS}

The targeted variable region.

=item B<--1Step>

Use this flag if the data are prepared by 1-Step PCR (only r1 & r2 raw files
available)

=item B<--storage-dir>=path, B<-sd> path

Indicate an existing directory in which to place the project directory.
"scratch" evaluates to "/local/scratch/" and "groupshare" evaluates to
"/local/groupshare/ravel"

=item B<-h>, B<--help>

Print help message and exit successfully.

=item B<--qsub-project>=space, B<-qp> space

Indicate which qsub-project space should be used for all qsubmissions. The
default is jravel-lab.

=item B<-dbg> {qiime_and_validation, extract_barcodes, demultiplex, tagclean, dada2}

Runs only one section of the pipeline and writes every qsub command to the log file.

=back 

=cut

use strict;
use warnings;
use Pod::Usage;
use English qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Cwd qw(abs_path);
use File::Temp qw/ tempfile /;
use File::Basename;
use POSIX;
use File::Spec::Functions;

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
    "raw-path|i=s"           => \my $inDir,
    "r1-path|r1=s"           => \my $r1file,
    "r2-path|r2=s"           => \my $r2file,
    "r3-path|r3=s"           => \my $r3file,
    "r4-path|r4=s"           => \my $r4file,
    "project-name|p=s"       => \my $project,
    "run-ID|r=s"             => \my $run,
    "map|m=s"                => \my $map,
    "var-reg|v=s"            => \my $var,
    "help|h!"                => \my $help,
    "dbg:s"                  => \my $dbg,
    "dry-run"                => \my $dryRun,
    "skip-err-thld"          => \my $skipErrThldStr,
    "dada2-truncLen-f|for:i" => \my $f,
    "dada2-truncLen-r|rev:i" => \my $r,
    "dada2-maxN:s"           => \my $maxN,
    "dada2-maxEE:s"          => \my $maxEE,
    "dada2-truncQ:s"         => \my $truncQ,
    "dada2-rmPhix"           => \my $phix,
    "dada2-maxLen:s"         => \my $maxLen,
    "dada2-minLen:s"         => \my $minLen,
    "dada2-minQ:s"           => \my $minQ,
    "1Step"                  => \my $oneStep,
    "storage-dir|sd=s"       => \my $sd,
    "qsub-project|qp:s"      => \my $qproj,
  )
##add option for final resting place of important data

  or pod2usage( verbose => 0, exitstatus => 1 );

if ($dbg) {
    if (   $dbg eq "qiime_and_validation"
        || $dbg eq "extract_barcodes"
        || $dbg eq "demultiplex"
        || $dbg eq "tagclean"
        || $dbg eq "dada2" )
    {
        print "DBG = $dbg";
    } else {
        die "Illegal debug option: -dbg $dbg";
    }
}

if ($help) {
    pod2usage( verbose => 2, exitstatus => 0 );
    exit 1;
}

if ( !$inDir && !( $r1file && $r2file && $r3file && $r4file ) ) {
    die "\n\tPlease provide the location of the raw sequencing files "
      . "(single directory => -i)\n\t\tOR \n\tFull paths to each R1, R2, R3, "
      . "and R4 file => -r1, -r2, -r3, -r4)\n\n";
}

if ( $inDir && ( $r1file || $r2file || $r3file || $r4file ) ) {
    die
"\n\tInput directory (-i) and raw files (-r*) were both specified. Please "
      . "provide one or the other.\n\n";
}

if ( !$map ) {
    print "\n\tPlease provide a full path to the project mapping file (-m)\n\n";
    exit 1;
}
if ( !$var ) {
    print "\n\tPlease indicate the targeted variable region (-v V3V4 or -v V4"
      . " or -v ITS)\n\n";
    exit 1;
}
my $models;
if ( $var eq 'V3V4' ) {
    $models = "/local/projects-t2/jholm/PECAN/v1.0/V3V4/merged_models/";
}
if ( $var eq 'V4' || $var eq 'ITS' ) {
    $models = "PECAN not used.\n";
}
if ( !$sd ) {
    die "\n***Please choose a storage directory (-sd), either 'scratch' or "
      . "'groupshare'.\nOR provide a full path to an existing directory.";
}

if ( $f && !$r ) {
    print "***\nPlease provide truncation lengths for forward and reverse "
      . "reads\n";
    die;
}

if ( !$qproj ) {
    print
      "\nqsub-project ID (--qsub-project, -qp) not provided. Using jravel-lab "
      . " as default\n";
    $qproj = "jravel-lab";
}

my $R = "/usr/local/packages/r-3.4.0/bin/R";

####################################################################
##                               DIRECTORY NAMES, VARIABLES
####################################################################
my $pd;
my $wd;
if ($inDir) {
    $inDir =~ s/\/$//;    # remove trailing slash
}

if ( $sd eq "scratch" ) {
    $pd = "/local/scratch/$project";
    $wd = "/local/scratch/$project/$run";
} elsif ( $sd eq "groupshare" ) {
    $pd = "/local/groupshare/ravel/$project";
    $wd = "/local/groupshare/ravel/$project/$run";
} else {
    $pd = "$sd/$project";
    $wd = "$sd/$project/$run";
}

my $error_log  = "$wd/qsub_error_logs";
my $stdout_log = "$wd/qsub_stdout_logs";
## Instead of working directory, give option for provided directory (or current
#  working directory)

my $rForSplit;
my $rRevSplit;

# my $r1 = "$wd/$project" . "_" . "$run" . "_"
#   . "R1.fastq";    ## change to full path (flag)
# my $r2 = "$wd/$project" . "_" . "$run" . "_"
#   . "R2.fastq";    ## change to full path (flag)
# my $r3 = "$wd/$project" . "_" . "$run" . "_"
#   . "R3.fastq";    ## change to full path - usr provides full path (flag)
# my $r4;

if ($oneStep) {
    $rForSplit = "$wd/R1split";
    $rRevSplit = "$wd/R2split";

    # $r4      = "$wd/$project" . "_" . "$run" . "_" . "R2.fastq";
} else {
    $rForSplit = "$wd/R1split";
    $rRevSplit = "$wd/R4split";

    # $r4      = "$wd/$project" . "_" . "$run" . "_" . "R4.fastq";
}

my $r1seqs = "$rForSplit/split_by_sample_out";
my $r4seqs = "$rRevSplit/split_by_sample_out";
my $cmd;

my $split_log = "$rForSplit/split_library_log.txt";
my @split;
my $newSamNo;

my @errors;
my $qiime = "$wd/$project" . "_" . $run . "_" . "qiime_config.txt";

###### BEGIN MAKE WORKING DIRECTORIES / VALIDATE INPUT FILES ##########
################################################

if ( !-e $pd ) {
    mkdir $pd;
    print "$pd did not exist -> making $pd\n";
}
if ( !-e $wd ) {
    mkdir $wd;
    print "$wd did not exist -> making $wd\n";
}
if ( !-e $error_log ) {
    mkdir $error_log;
    print "$error_log did not exist -> making $error_log\n";
}
if ( !-e $stdout_log ) {
    mkdir $stdout_log;
    print "$stdout_log did not exist -> making $stdout_log\n";
}

my $time = strftime( "%Y-%m-%d %H:%M:%S", localtime(time) );
my $log  = "$wd/$project" . "_" . $run . "_16S_pipeline_log.txt";

my $perlScript = catfile( $pipelineDir, "scripts", "log_version.pl" );
system( $^X, $perlScript, $log );

open my $logFH, ">>$log" or die "Cannot open $log for writing: $OS_ERROR";
print $logFH "$time\n";

my $index1Input;    # In one-step runs, these are the SAME files pointed to by
my $index2Input;    # $readsForInput and $readsRevInput

my $readsForInput;
my $readsRevInput;

if ($inDir) {       # Search for raw fastq's (or fastq.gz's)
    if ($oneStep) {

        # These are the only input files needed for 1step???
        my @r1s = glob("$inDir/*R1.fastq $inDir/*R1.fastq.gz");
        my @r2s = glob("$inDir/*R2.fastq $inDir/*R2.fastq.gz");

        # my @r3s = glob("$inDir/*R3.fastq $inDir/*R3.fastq.gz");
        # my @r4s = glob("$inDir/*R4.fastq $inDir/*R4.fastq.gz");

        if (
               scalar @r1s == 1
            && scalar @r2s == 1

            # && scalar @r3s == 1
            # && scalar @r4s == 1
          )
        {
            $index1Input = $readsForInput = $r1s[0];
            $index2Input = $readsRevInput = $r2s[0];
        } else {
            print $logFH "Couldn't find input files.\n";
            print $logFH "Files found in $inDir:\n";

            my $printme = join "\n", @r1s;
            print $logFH "$printme\n" if $printme;
            $printme = join "\n", @r2s;
            print $logFH "$printme\n" if $printme;

            # $printme = join "\n", @r3s;
            # print $logFH "$printme\n" if $printme;
            # $printme = join "\n", @r4s;
            # print $logFH "$printme\n" if $printme;
            die
"Could not find a complete and exclusive set of raw files. Since --1step given, input directory must"
              . "have exactly one R1, and R2 file. See $log for a list of the files"
              . " found.";
        }
    } else {
        my @i1s = glob("$inDir/*I1.fastq $inDir/*I1.fastq.gz");
        my @i2s = glob("$inDir/*I2.fastq $inDir/*I2.fastq.gz");
        my @r1s = glob("$inDir/*R1.fastq $inDir/*R1.fastq.gz");
        my @r2s = glob("$inDir/*R2.fastq $inDir/*R2.fastq.gz");
        my @r3s = glob("$inDir/*R3.fastq $inDir/*R3.fastq.gz");
        my @r4s = glob("$inDir/*R4.fastq $inDir/*R4.fastq.gz");

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
            print $logFH "Couldn't find input files.\n";
            print $logFH "Files found in $inDir:\n";
            my $printme = join "\n", @i1s;
            print $logFH "$printme\n" if $printme;
            $printme = join "\n", @i2s;
            print $logFH "$printme\n" if $printme;
            $printme = join "\n", @r1s;
            print $logFH "$printme\n" if $printme;
            $printme = join "\n", @r2s;
            print $logFH "$printme\n" if $printme;
            $printme = join "\n", @r3s;
            print $logFH "$printme\n" if $printme;
            $printme = join "\n", @r4s;
            print $logFH "$printme\n" if $printme;
            die
"Could not find a complete and exclusive set of raw files. Input directory must"
              . " contain exactly one set of I1, I2, R1, and R2, OR R1, R2, R3, and R4.\n"
              . "See $log for a list of the files found.";
        }
    }
} else {
    if ($oneStep) {

  # Correct??? Documentation say that r1 -> R1.fastq and r2, r3, r4 -> R2.fastq.
  # What actually happened was, R1 and R2 files autofound
  # from $inDir were used for barcodes, and silently nothing happened for
  # split_libraries (because the two files were already in the run directory)
  # pipeline was never coded to run 1-step WITHOUT -i
        $index1Input = $readsForInput = $r1file;
        $index2Input = $readsRevInput = $r4file;
    } else {
        $readsForInput = $r1file;
        $index1Input   = $r2file;
        $index2Input   = $r3file;
        $readsRevInput = $r4file;
    }
}

# In the multidbg branch, combine this section with qiime_and_validate

my @cmds;

# Regardless of one-step/two-step, get forward and reverse reads as .fastq
if ( $readsForInput !~ /.gz$/ ) { # if the index files aren't .gz, read in place
} else {                          # otherwise...
                                  # discard any existing file in working dir
                                  # (we can't test if it's complete)

    # Rename *[R|I][1|2].fastq.gz to <WD>/<PROJ>_<RUN>_[R|I][1|2].fastq
    my $dest = File::basename($readsForInput);
    $dest = substr( $dest, ( length($dest) - 10 ), 8 );  # get suffix (sans .gz)
    $dest = "$wd/$project" . "_" . "$run" . "_" . $dest;
    push( @cmds, "gzip --decompress $readsForInput > $dest" )
      ;    #copy to our working dir
    $readsForInput = $dest;
    if ($oneStep) {
        $index1Input = $dest;
    }
}
if ( $readsRevInput !~ /.gz$/ ) {    # same for reverse reads
} else {
    my $dest = File::basename($readsRevInput);
    $dest = substr( $dest, ( length($dest) - 10 ), 8 );
    $dest = "$wd/$project" . "_" . "$run" . "_" . $dest;
    push( @cmds, "gzip --decompress $readsRevInput > $dest" );
    $readsRevInput = $dest;
    if ($oneStep) {
        $index2Input = $dest;
    }
}

if ( !$oneStep ) {

    # For two step, we need to do the same for the index files
    if ( $index1Input !~ /.gz$/ ) {
    } else {
        my $dest = File::basename($index1Input);
        $dest = substr( $dest, ( length($dest) - 10 ), 8 );
        $dest = "$wd/$project" . "_" . "$run" . "_" . $dest;
        push( @cmds, "gzip --decompress $index1Input > $dest" );
        $index1Input = $dest;
    }
    if ( $index2Input !~ /.gz$/ ) {    # same for reverse reads
    } else {
        my $dest = File::basename($index2Input);
        $dest = substr( $dest, ( length($dest) - 10 ), 8 );
        $dest = "$wd/$project" . "_" . "$run" . "_" . $dest;
        push( @cmds, "gzip --decompress $index2Input > $dest" );
        $index2Input = $dest;
    }
}

# Execute the commands queued above
if ( scalar @cmds > 0 ) {
    print "---Uncompressing raw file(s) to $wd. Commands:\n";
    execute_and_log( @cmds, 0, $dryRun );
    @cmds = ();
}
if ($oneStep) {
    print $logFH
      "$project raw reads files:\n\t$readsForInput\n\t$readsRevInput\n";
} else {
    print $logFH
      "$project raw reads files:\n\t$readsForInput\n\t$readsRevInput\n"
      . "$project raw index files:\n\t$index1Input\n\t$index2Input\n";
}

if ( ( !$dbg ) || $dbg eq "qiime_and_validation" ) {

    ###### BEGIN CHECK OF QIIME CONFIGURATION ###########
    #####################################################

    $cmd = "print_qiime_config.py > $qiime";
    print "\tcmd=$cmd\n" if $dbg;
    system($cmd) == 0
      or die "system($cmd) failed with exit code: $?"
      if !$dryRun;

    ###### BEGIN VALIDATION OF MAPPING FILE ###########
    ################################################

    @errors = glob("$error_log/*.log");
    if (@errors) {
        foreach my $error (@errors) {
            $cmd = "rm $error";
            print "\tcmd=$cmd\n" if $dbg;
            system($cmd) == 0
              or die "system($cmd) failed with exit code: $?"
              if !$dryRun;
        }
    }

    print "--Validating $map\n";
    $cmd = "validate_mapping_file.py -m $map -s -o $error_log";
    print "\tcmd=$cmd\n" if $dbg;
    system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

    my $mappingError = glob("$error_log/*.log");
    if ($mappingError) {
        open MAPERROR, "<$mappingError"
          or die "Cannot open $mappingError for " . "reading: $OS_ERROR";
        $_ = <MAPERROR>;
        close MAPERROR;
        chomp;
        if ( $_ =~ /No errors or warnings found in mapping file./ ) {
            print "---Map passed validation.\n";
        } else {
            die "***Error in mapping file. See $mappingError for "
              . "details. Exiting.\n";
        }
    } else {
        die "validate_mapping_file.py did not produce an error log";
    }
}

if ( $dbg eq "qiime_and_validation" ) {
    die "Finished printing QIIME configuration and validating mapping file.";
}

open MAP, "<$map" or die "Cannot open $map for reading: $OS_ERROR";
my $nSamples    = 0;
my $extctrl     = 0;
my $pcrpos      = 0;
my $pcrneg      = 0;
my $projSamples = 0;
my $linecount   = 0;
my $null        = 0;
while (<MAP>) {
    chomp;
    $linecount++;
    if ( $linecount > 1 )    ## don't count header as sample
    {
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
    }
}
close MAP;

$nSamples = $projSamples + $extctrl + $pcrpos + $pcrneg + $null;

my $barcodes = "$wd/barcodes.fastq";

if ( ( !$dbg ) || $dbg eq "extract_barcodes" ) {

    ###### BEGIN EVALUATION OF SAMPLES VIA MAPPING FILE ###########
    ###############################################################
    print $logFH "PROJECT: $project\nVARIABLE REGION: $var\n"
      . "NO. SAMPLES: $nSamples\nNO. NULLS: $null\n"
      . "NO. EXTRACTION NEGATIVE CONTROLS: $extctrl\n"
      . "NO. PCR POSITIVE CONTROLS: $pcrpos\nNO. PCR NEGATIVE CONTROLS: $pcrneg\n"
      . "R VERSION: $R\nPECAN MODELS: $models\nQIIME CONFIGURATION DETAILS: \n"
      . "see $qiime\nMAPPING FILE: $map";
    if ($oneStep) {
        print $logFH "\nPCR PREPARATION METHOD: 1-Step\n\n";
    } else {
        print $logFH "\nPCR PREPARATION METHOD: 2-Step\n\n";
    }

    ## change to full path to full barcodes (flag)
    my $count = 0;

    print "--Checking for existence of $barcodes\n";

    my $step1;
    my $step2;
    my $step3;

    ###### BEGIN BARCODES ##########
    #######################################

    if ( !-e $barcodes ) {    # FIXME test if barcodes is the right size
        print "---Barcode files not found or were empty\n";
        my @cmds;

        my $start = time;
        print "---Extracting barcodes\n";
        $step1 = "extract_barcodes.py";

        my $mapOpt = "";
        my $bcLen  = 8;       # 2-step pcr
        if ($oneStep) {

            # Was in original code, but maybe unnecessary?
            $mapOpt = " -m \$map";

            $bcLen = 12;      # 1-step pcr
        }
        $cmd =
"qsub -b y -l mem_free=1G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log extract_barcodes.py -f $index1Input -r $index2Input -c barcode_paired_end --bc1_len $bcLen --bc2_len $bcLen $mapOpt -o $wd";
        print "\tcmd=$cmd\n" if $dbg;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;

        check_error_log( $error_log, $step1 );

        print "---Waiting for barcode extraction to complete.\n";

        #FIXME again, test size of barcodes
        while ( !( -e $barcodes ) ) { sleep 1; }

        ## try adding else statement, if needed
        my $duration = time - $start;
        print "---Barcode extraction complete\n";
        print $logFH "---Duration of barcode extraction: $duration s\n";
        print "---Duration of barcode extraction: $duration s\n";
    } else {
        open BAR, "<$barcodes"
          or die "Cannot open $barcodes for reading: " . "$OS_ERROR";
        while (<BAR>) { }
        if ( $. = $nSamples ) {
            print "---$barcodes already exists and contains $. entries,"
              . " as expected\n";
            print $logFH "-> $barcodes already exists and contains $. "
              . "entries, as expected\n";
        } else {

            # Recover from this?
            print "---$barcodes already exists, but contains $. "
              . "entries, while there are $nSamples samples. Exiting.\n";
            print $logFH "-> $barcodes already exists, but contains $."
              . " entries, while there are $nSamples samples. Exiting.\n";
            die;
        }
        close BAR;
    }

    # Delete unneeded output of extract_barcodes.py
    @cmds = ();
    push @cmds, "rm -rf $wd/reads1.fastq";
    push @cmds, "rm -rf $wd/reads2.fastq";
    execute_and_log( @cmds, 0, $dryRun );
}

if ( $dbg eq "extract_barcodes" ) {
    die "Finished extracting barcodes\n";
}

###### BEGIN SPLIT LIBRARIES ##########
#######################################

## think of way to ensure the consistent read order in the r1 and r4 files.
## print headers of r1 and r4, comm r1 r4 - to ensure the seqIDs are the same order.
if ( !$dbg || $dbg eq "demultiplex" ) {
    my @cmds;
    my $rForSeqsFq = "$rForSplit/seqs.fastq";
    my $rRevSeqsFq = "$rRevSplit/seqs.fastq";

    # Just helps with printing
    my $revName;
    if ($oneStep) {
        $revName = "R2";
    } else {
        $revName = "R4";
    }

    if ( !-e $rForSeqsFq || !-e $rRevSeqsFq ) {
        print "---Producing $rForSeqsFq and $rRevSeqsFq\n";

#         if ($oneStep) {
#             print "---Obtaining $project-specific R1 and R2 fastq files\n";
#         } else {
#             print "---Obtaining $project-specific R1 and R4 fastq files\n";
#             if ($inDir) {
#                 $cmd =
# "zcat $readsForInput > $r1 | zcat $readsRevInput > $r4 ";
#                 print "\tcmd=$cmd\n" if $dbg;
#                 system($cmd) == 0
#                   or die "system($cmd) failed with exit code: $?"
#                   if !$dryRun;
#                 print $logFH "---$project barcode and index files copied from"
#                   . " $inDir to $r1 and $r4\n";
#             } else {
#                 $cmd = "zcat $r1file > $r1 | zcat $r4file > $r4 ";
#                 print "\tcmd=$cmd\n" if $dbg;
#                 system($cmd) == 0
#                   or die "system($cmd) failed with exit code: $?"
#                   if !$dryRun;
#                 print $logFH "---$project barcode and index files copied from "
#                   . "$r1file and $r4file to $r1 and $r4\n";
#             }
#         }

        my $start = time;
        my $step2 = "split_libraries_fastq.py";
        @errors = glob("$error_log/$step2.e*");
        if (@errors) {
            foreach my $error (@errors) {
                push @cmds, "rm $error";
            }
            execute_and_log( @cmds, $dryRun );
            @cmds = ();
        }
        ## including the stitch script before so that demultiplex happens by barcode NOT order
        my $barcodeType;
        if ($oneStep) {
            $barcodeType = "24";
        } else {
            $barcodeType = "16";
        }

        push @cmds,
"qsub -cwd -b y -l mem_free=1G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log $step2 -i $readsForInput -o $rForSplit -b $barcodes -m $map --max_barcode_errors 1 --store_demultiplexed_fastq --barcode_type $barcodeType -r 999 -n 999 -q 0 -p 0.0001";
        push @cmds,
"qsub -cwd -b y -l mem_free=1G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log $step2 -i $readsRevInput -o $rRevSplit -b $barcodes -m $map --max_barcode_errors 1 --store_demultiplexed_fastq --barcode_type $barcodeType -r 999 -n 999 -q 0 -p 0.0001";
        print $logFH "Demultiplexing commands: \n";
        execute_and_log( @cmds, $logFH, $dryRun );
        @cmds = ();
        print $logFH "\n\n";

        check_error_log( $error_log, $step2 );

        print "---Waiting for R1 seqs.fastq to complete.\n";
        while ( !( -e $rForSeqsFq ) ) { sleep 1; }
        my $duration = time - $start;
        print $logFH "---Duration of R1 seqs.fastq production: $duration s\n";
        print "---Duration of R1 seqs.fastq production: $duration s\n";

        print "---Waiting for $revName seqs.fastq to complete.\n";
        while ( !( -e $rRevSeqsFq ) ) { sleep 1; }
        $duration = time - $start;
        print $logFH "---Duration of $revName seqs.fastq production: "
          . "$duration s\n";
        print "---Duration of $revName seqs.fastq production: $duration s\n";

    } else {
        print "-> $rForSeqsFq and $rRevSeqsFq already produced. Demultiplexing "
          . "...\n";
    }

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
                print "$sample $nReads\n";
            }
        }
    }
    close SPLIT;

    $newSamNo = $nSamples - scalar @split;
    if ( scalar @split ne $nSamples ) {
        print $logFH "The following "
          . scalar @split
          . " samples returned 0 "
          . "reads, due to either a lack of barcodes or a filtering of those "
          . "reads in split_libraries_fastq.py:\n";
        foreach my $x (@split) {
            print $logFH "   $x\n";
        }
        print $logFH "Number of samples after split_libraries_fastq.py: "
          . "$newSamNo\n";
    } elsif ( scalar @split eq $nSamples ) {
        print $logFH "No samples were successfully demultiplexed. Is the "
          . "mapping file correct? Exiting.\n";
        die;
    }

    ###### BEGIN SPLIT BY SAMPLE ##########
    #######################################

    print "--Checking if $project R1 & $revName seqs.fastq files were split by "
      . "sample ID\n";

    my $n_fq1 = 0;

    my @forFilenames = glob("$r1seqs/*.fastq");
    my @revFilenames = glob("$r4seqs/*.fastq");

    if (   scalar(@forFilenames) != $newSamNo
        || scalar(@revFilenames) != $newSamNo )
    {
        print "There are "
          . scalar @forFilenames
          . " split files. Waiting for $newSamNo "
          . "files\n";
        my $step3 = "split_sequence_file_on_sample_ids.py";
        @errors = glob("$error_log/$step3.e*");
        if (@errors) {
            foreach my $error (@errors) {
                push( @cmds, "rm $error" );
            }
            execute_and_log( @cmds, 0, $dryRun );
            @cmds = ();
        }
        print "---Sample specific files not found or completed... Splitting "
          . "$project seqs.fastq files by sample ID\n";
        print $logFH "There are "
          . scalar(@forFilenames)
          . " sample specific "
          . "files found (expected $newSamNo)... Splitting $project seqs.fastq "
          . "files by sample ID\n";
        execute_and_log( "rm -rf $r1seqs; rm -rf $r4seqs", 0, $dryRun );

        while ( !( -e $rForSeqsFq ) ) { sleep 1; }
        execute_and_log(
"qsub -cwd -b y -l mem_free=1G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log -V $step3 -i $rForSeqsFq --file_type fastq -o $r1seqs",
            0, $dryRun
        );

        while ( !( -e $rRevSeqsFq ) ) { sleep 1; }
        execute_and_log(
"qsub -cwd -b y -l mem_free=1G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log -V $step3 -i $rRevSeqsFq --file_type fastq -o $r4seqs",
            0, $dryRun
        );

        check_error_log( $error_log, $step3 );

        ## the $nSamples needs to be altered if a sample has 0 reads, because the sample-specific fastq won't be produced
        my $n_fq   = 0;
        my $nLines = 0;

        # Count the number of lines in R1split/seqs.fastq
        open R1, "<$rForSplit/seqs.fastq";
        while (<R1>) { }
        my $R1lines = $.;
        close R1;

        while ( $n_fq != $newSamNo || $nLines != $R1lines ) {
            $n_fq   = 0;
            $nLines = 0;
            my @filenames = glob("$r1seqs/*.fastq");
            if (@filenames) {
                foreach my $file (@filenames) {
                    $n_fq++;    # Count the number of files in the directory

                    # Also keep a running total of lines in the
                    # split_by_sample_out files
                    open SAMPLE, "<$file";
                    while (<SAMPLE>) { }
                    $nLines += $.;
                    close SAMPLE;
                }
            }
        }

        print "--All samples ($n_fq) and reads (@{[$nLines / 4]}) accounted"
          . " for in $r1seqs\n";

        $n_fq   = 0;
        $nLines = 0;

        # Count the number of lines in R4split/seqs.fastq
        open R4, "<$rRevSplit/seqs.fastq";
        while (<R4>) { }
        my $R4lines = $.;
        close R4;

        while ( $n_fq != $newSamNo || $nLines != $R4lines ) {
            $n_fq   = 0;
            $nLines = 0;
            my @filenames = glob("$r4seqs/*.fastq");
            if (@filenames) {
                foreach my $file (@filenames) {
                    $n_fq++;    # Count the number of files in the directory

                    # Also keep a running total of lines in the
                    # split_by_sample_out files
                    open SAMPLE, "<$file";
                    while (<SAMPLE>) { }
                    $nLines += $.;
                    close SAMPLE;
                }
            }
        }

        print "--All samples ($n_fq) and reads (@{[$nLines / 4]}) accounted for"
          . " in $r4seqs\n";
    } else {
        print "--$newSamNo sample-specific files present as expected.\n";
        print $logFH "$newSamNo sample-specific files present as expected.\n";
    }

    # Remove temporarily decompressed files
    print "---Removing decompressed raw files from $wd\n";

    @cmds = ();
    push( @cmds, "rm -rf $wd/${run}_R1.fastq" );
    push( @cmds, "rm -rf $wd/${run}_R2.fastq" );
    push( @cmds, "rm -rf $wd/${run}_R3.fastq" );
    push( @cmds, "rm -rf $wd/${run}_R2.fastq" );
    push( @cmds, "rm -rf $wd/${run}_I1.fastq" );
    push( @cmds, "rm -rf $wd/${run}_I2.fastq" );
    execute_and_log( @cmds, 0, $dryRun );
}

if ( $dbg eq "demultiplex" ) {
    die "Finished demultiplexing libraries";
}

print "--Checking if target primers have been removed from $project R1 & R4"
  . " sample-specific files..\n";
print $logFH
  "Checking if target primers have been removed from $project R1 & R4"
  . " sample-specific files..\n";
my @r1tcfiles = glob("$wd/*R1_tc.fastq");
my @r4tcfiles = glob("$wd/*R2_tc.fastq");

###### BEGIN TAGCLEANING ##########
###################################

my $start = time;
if ( !$dbg || $dbg eq "tagclean" ) {
    if ( !( defined $newSamNo ) ) {
        open SPLIT, "<$split_log"
          or die "Cannot open $split_log for writing: $OS_ERROR";
        while (<SPLIT>) {
            if ( $_ =~ /\t/ ) {
                chomp;
                my ( $sample, $nReads ) = split, /\t/;
                chomp $nReads;
                if ( $nReads eq "0" ) {

                    #print $logFH "Adding sample: $sample\n";
                    push @split, $sample;
                } else {
                    print "$sample $nReads\n";
                }
            }
        }
        close SPLIT;

        $newSamNo = $nSamples - scalar @split;
    }

    if ( scalar @r1tcfiles != $newSamNo || scalar @r4tcfiles != $newSamNo ) {
        if ($oneStep) {
            if ( $var eq "V3V4" ) {
                print "Removing V3V4 primers from all sequences\n";
                my $filename;
                opendir R1, $r1seqs or die "Cannot open directory $r1seqs\n";
                while ( $filename = readdir R1 ) {
                    if ( $filename =~ /.fastq/ ) {
                        my @suffixes = ( ".fastq", ".fq" );
                        my $Prefix   = basename( $filename, @suffixes );
                        my $tc       = "$wd/$Prefix" . "_R2_tc";
                        $cmd =
"perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r1seqs/$filename -out $tc -line_width 0 -verbose -tag5 GGACTACHVGGGTWTCTAAT -mm5 2 -trim_within 50";
                        print "\tcmd=$cmd\n" if $dbg;
                        system($cmd) == 0
                          or die "system($cmd) failed with exit code: $?"
                          if !$dryRun;
                    }
                }
                close R1;

                opendir R4, $r4seqs or die "Cannot open directory $r4seqs\n";
                while ( $filename = readdir R4 ) {
                    my @suffixes = ( ".fastq", ".fq" );
                    my $Prefix   = basename( $filename, @suffixes );
                    my $tc       = "$wd/$Prefix" . "_R1_tc";
                    $cmd =
"perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r4seqs/$filename -out $tc -line_width 0 -verbose -tag5 ACTCCTACGGGAGGCAGCAG -mm5 2 -trim_within 50";
                    print "\tcmd=$cmd\n" if $dbg;
                    system($cmd) == 0
                      or die "system($cmd) failed with exit code: $?"
                      if !$dryRun;
                }
                close R4;
            }

            if ( $var eq "V4" ) {
                print "--Removing V4 primers from all sequences\n";
                my $filename;
                opendir R1, $r1seqs or die "Cannot open directory $r1seqs\n";
                while ( $filename = readdir R1 ) {
                    if ( $filename =~ /.fastq/ ) {
                        my @suffixes = ( ".fastq", ".fq" );
                        my $Prefix   = basename( $filename, @suffixes );
                        my $tc       = "$wd/$Prefix" . "_R1_tc";
                        $cmd =
"qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r1seqs/$filename -out $tc -line_width 0 -verbose -tag5 GTGCCAGCMGCCGCGGTAA -mm5 2";
                        print "\tcmd=$cmd\n" if $dbg;
                        system($cmd) == 0
                          or die "system($cmd) failed with exit code: $?"
                          if !$dryRun;
                    }
                }
                close R1;

                opendir R4, $r4seqs or die "Cannot open directory $r4seqs\n";
                while ( $filename = readdir R4 ) {
                    my @suffixes = ( ".fastq", ".fq" );
                    my $Prefix   = basename( $filename, @suffixes );
                    my $tc       = "$wd/$Prefix" . "_R2_tc";
                    $cmd =
"qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r4seqs/$filename -out $tc -line_width 0 -verbose -tag5 ACTCCTACGGGAGGCAGCAG -mm5 2";
                    print "\tcmd=$cmd\n" if $dbg;
                    system($cmd) == 0
                      or die "system($cmd) failed with exit code: $?"
                      if !$dryRun;
                }
                close R4;
            }
        } else {
            if ( $var eq "V3V4" ) {
                print $logFH "...Removing V3V4 primers from all sequences.\n";
                my $filename;
                opendir R1, $r1seqs or die "Cannot open directory $r1seqs\n";
                while ( $filename = readdir R1 ) {
                    if ( $filename =~ /.fastq/ ) {
                        my @suffixes = ( ".fastq", ".fq" );
                        my $Prefix   = basename( $filename, @suffixes );
                        my $tc       = "$wd/$Prefix" . "_R1_tc";
                        $cmd =
"cp $r1seqs/$filename /home/jolim/redirected_files/R1split/$filename";
                        print "\tcmd=$cmd\n" if $dbg;

                        system($cmd) == 0
                          or die "system($cmd) failed with exit code: $?"
                          if !$dryRun;

                        $cmd =
"qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r1seqs/$filename -out $tc -line_width 0 -verbose -tag5 ACTCCTACGGGAGGCAGCAG -mm5 2";

                        print "\tcmd=$cmd\n" if $dbg;
                        system($cmd) == 0
                          or die "system($cmd) failed with exit code: $?"
                          if !$dryRun;
                    }
                }
                close R1;

                opendir R4, $r4seqs or die "Cannot open directory $r4seqs\n";
                while ( $filename = readdir R4 ) {
                    if ( $filename =~ /.fastq/ ) {
                        my @suffixes = ( ".fastq", ".fq" );
                        my $Prefix   = basename( $filename, @suffixes );
                        my $tc       = "$wd/$Prefix" . "_R2_tc";
                        $cmd =
"cp $r4seqs/$filename /home/jolim/redirected_files/R4split/$filename";
                        print "\tcmd=$cmd\n" if $dbg;

                        system($cmd) == 0
                          or die "system($cmd) failed with exit code: $?"
                          if !$dryRun;

                        $cmd =
"qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r4seqs/$filename -out $tc -line_width 0 -verbose -tag5 GGACTACHVGGGTWTCTAAT -mm5 2";

                        print "\tcmd=$cmd\n" if $dbg;
                        system($cmd) == 0
                          or die "system($cmd) failed with exit code: $?"
                          if !$dryRun;
                    }
                }
                close R4;
            }
            if ( $var eq "V4" ) {
                print "...Removing V4 primers from all sequences.\n";
                my $filename;
                opendir R1, $r1seqs or die "Cannot open directory $r1seqs\n";
                while ( $filename = readdir R1 ) {
                    if ( $filename =~ /.fastq/ ) {
                        my @suffixes = ( ".fastq", ".fq" );
                        my $Prefix   = basename( $filename, @suffixes );
                        my $tc       = "$wd/$Prefix" . "_R1_tc";
                        $cmd =
"qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r1seqs/$filename -out $tc -line_width 0 -verbose -tag5 GTGCCAGCMGCCGCGGTAA -mm5 2";
                        print "\tcmd=$cmd\n" if $dbg;
                        system($cmd) == 0
                          or die "system($cmd) failed with exit code: $?"
                          if !$dryRun;
                    }
                }
                close R1;

                opendir R4, $r4seqs or die "Cannot open directory $r4seqs\n";
                while ( $filename = readdir R4 ) {
                    my @suffixes = ( ".fastq", ".fq" );
                    my $Prefix   = basename( $filename, @suffixes );
                    my $tc       = "$wd/$Prefix" . "_R2_tc";
                    $cmd =
"qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r4seqs/$filename -out $tc -line_width 0 -verbose -tag5 ACTCCTACGGGAGGCAGCAG -mm5 2";
                    print "\tcmd=$cmd\n" if $dbg;
                    system($cmd) == 0
                      or die "system($cmd) failed with exit code: $?"
                      if !$dryRun;
                }
                close R4;
            }
            if ( $var eq "ITS" ) {
                print "...Removing ITS primers from all sequences.\n";
                my $filename;
                opendir R1, $r1seqs or die "Cannot open directory $r1seqs\n";
                while ( $filename = readdir R1 ) {
                    if ( $filename =~ /.fastq/ ) {
                        my @suffixes = ( ".fastq", ".fq" );
                        my $Prefix   = basename( $filename, @suffixes );
                        my $tc       = "$wd/$Prefix" . "_R1_tc";
                        $cmd =
"qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r1seqs/$filename -out $tc -line_width 0 -verbose -tag5 CTGCCCTTTGTACACACCGC -mm5 2";
                        print "\tcmd=$cmd\n" if $dbg;
                        system($cmd) == 0
                          or die "system($cmd) failed with exit code: $?"
                          if !$dryRun;
                    }
                }
                close R1;

                opendir R4, $r4seqs or die "Cannot open directory $r4seqs\n";
                while ( $filename = readdir R4 ) {
                    my @suffixes = ( ".fastq", ".fq" );
                    my $Prefix   = basename( $filename, @suffixes );
                    my $tc       = "$wd/$Prefix" . "_R2_tc";
                    $cmd =
"qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r4seqs/$filename -out $tc -line_width 0 -verbose -tag5 TTTCGCTGCGTTCTTCATCG -mm5 2";
                    print "\tcmd=$cmd\n" if $dbg;
                    system($cmd) == 0
                      or die "system($cmd) failed with exit code: $?"
                      if !$dryRun;
                }
                close R4;
            }
        }

        my @files  = glob("$wd/*R1_tc.fastq");
        my $nFiles = @files;
        while ( $nFiles != $newSamNo ) {
            @files  = glob("$wd/*R1_tc.fastq");
            $nFiles = @files;
        }
        print "---All tagcleaned R1 samples accounted for in $wd\n";

        @files  = glob("$wd/*R4_tc.fastq");
        $nFiles = @files;
        while ( $nFiles != $newSamNo ) {
            @files  = glob("$wd/*R2_tc.fastq");
            $nFiles = @files;
        }
        print "---All tagcleaned R4 (R2) samples accounted for in $wd\n";

        my $duration = time - $start;
        print $logFH "...Primer sequences removed from $newSamNo samples. "
          . "Beginning DADA2.\n";
        print $logFH "--Duration of tagcleaning: $duration s\n";
    } else {
        print "--$newSamNo sample-specific, tag-cleaned files present as "
          . "expected.\n";
        print $logFH "...$newSamNo sample-specific, tag-cleaned files present "
          . "as expected. Beginning DADA2.\n";
    }
}

if ( $dbg eq "tagclean" ) { die "tagclean finished"; }

###### BEGIN DADA2 ##########
#############################
if ( ( !$dbg ) || $dbg eq "dada2" ) {
    my $dada2 = "$wd/dada2_abundance_table.rds";

    my $truncLen;
    if ( $f && $r ) {
        $truncLen = "c($f,$r)";
    }

    if ( !-e $dada2 ) {
        print "--Removing old filtered fastq files from previous runs\n";
        $cmd = "rm -rf $wd/filtered";
        print "\tcmd=$cmd\n" if $dbg;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;

        print "Running DADA2 with fastq files in $wd\n";
        print $logFH "Running DADA2 for $var region";
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
                    $phix, $maxLen,   $minLen, $minQ
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
                    $phix, $maxLen,   $minLen, $minQ
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
                    $phix, $maxLen,   $minLen, $minQ
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
                    $phix, $maxLen,   $minLen, $minQ
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
                    $phix, $maxLen,   $minLen, $minQ
                );
            }
        }
    }

    my $rt     = "$wd/dada2_part1_rTmp.R";
    my $projrt = "$wd/$project" . "_" . $run . "_dada2_part1_rTmp.R";
    if ( !-e $projrt ) {
        $cmd = "mv $rt $projrt";
        print "\tcmd=$cmd\n" if $dbg;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;
        print $logFH "\tDADA2-specific commands can be found in $projrt\n";
    }

    my $rtout     = "$wd/dada2_part1_rTmp.Rout";
    my $projrtout = "$wd/$project" . "_" . $run . "_dada2_part1_rTmp.Rout";
    if ( !-e $projrtout ) {
        $cmd = "mv $rtout $projrtout";
        print "\tcmd=$cmd\n" if $dbg;
        system($cmd) == 0
          or die "system($cmd) failed with exit code: $?"
          if !$dryRun;
        print $logFH "\tDADA2-specific commands with output can be found in "
          . "$projrtout\n";
    }

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
        print $logFH scalar @removed
          . " samples removed during dada2 filtering:\n";
        for my $x (@removed) {
            print $logFH "$x\n";
        }
    }
    print $logFH "\n";

    my $dadaTbl = "$wd/dada2_part1_stats.txt";
    if ( -e $dadaTbl ) {
        print $logFH "\nFor $var region, dada2 used the following filtering "
          . "requirements:\n$truncLen\n$maxN\n$maxEE\n$truncQ\n$phix\n";
        print $logFH "dada2 completed successfully!\nAbundance table for "
          . "$project run $run located at $wd/dada2_abundance_table.rds\n";
        print $logFH "See $dadaTbl for dada2 table of reads surviving by "
          . "step\n";
        close $logFH;
    } else {
        print "---dada2 did not complete successfully, something went wrong!\n"
          . "---Check $projrtout.\n";
    }
}
###### COMPLETING $logFH FILE ##############
#########################################

open $logFH, "<$log" or die "Cannot open $log for reading: $OS_ERROR";
while (<$logFH>) {
    if ( $_ =~ /dada2 completed successfully!/ ) {
        print "---dada2 completed successfully\n";
        print "---Run-specific abundance table written to "
          . "$wd/dada2_abundance_table.rds\n";
        print "---See $log for processing details\n";
    }
}
close $logFH;

## moving final files to directory created within /local/projects/16S_DATA/projects/
## user specified directory for final storage of files.

## if reference based chim rm requested { perform vsearch filtered.fastq}

####################################################################
##                               SUBS
####################################################################

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
  rownames(track)<-gsub("_R1_tc.fastq","",rownames(track))
  track[names(v),3]<-v
  colnames(track) <- c("input", "filtered", "merged")
  write.table(track, "dada2_part1_stats.txt", quote=FALSE, append=FALSE, sep=\t, row.names=TRUE, col.names=TRUE)
  ~;
    run_R_script($Rscript);
}

sub run_R_script {

    chdir $wd;
    my $Rscript = shift;

    my $outFile = "dada2_part1_rTmp.R";
    open OUT, ">$outFile", or die "cannot write to $outFile: $!\n";
    print OUT "$Rscript";
    close OUT;

    my $cmd = "$R CMD BATCH $outFile";
    system($cmd) == 0 or die "system($cmd) failed:$?\n";

    my $outR = $outFile . "out";
    open IN, "$outR" or die "Cannot open $outR for reading: $OS_ERROR\n";
    my $exitStatus = 1;

    foreach my $line (<IN>) {
        if (   $line eq /Error/
            && $line ne /learnErrors/
            && $line ne /error rates/
            && $line ne /errors/ )
        {
            print "R script crashed at\n$line";
            print "check $outR for details\n";
            $exitStatus = 0;
            die;
        }
    }
    close IN;
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

sub check_error_log {
    my $error_log = shift;
    my $step      = shift;
    print "---Checking for $step errors\n";
    my $error;
    ## make $error array so that goes through all errors for that step.
  RECHECK:
    $error = glob("$error_log/$step.e*");
    if ( defined $error ) {
        print "---Checking for errors in $error.\n";
        open ERROR, "<$error" or die "Can't open $error.\n";
        while (<ERROR>) {
            if ( $_ =~ /error/ || $_ =~ /Error/ || $_ =~ /ERROR/ ) {
                print "Error in $step. See $error.\n";
                print $logFH "Error in $step. See $error.\n";
                die;
            } else {
                print "---No errors found in $step.\n";
                print $logFH "---No errors found in $step.\n";
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
# the second to last argument must evaluate false to prevent loggin to log file.
# the last argument is 1 if this is a dry run
sub execute_and_log {
    my $dryRun = pop @_;
    my $logFH  = pop @_;
    for (@_) {
        print "\t$_\n" if $dbg;
        if ($logFH) {
            print $logFH "\t$_\n";
        }
        system($_) == 0
          or die "system($_) failed with exit code: $?"
          if !$dryRun;
    }
}

exit 0;
