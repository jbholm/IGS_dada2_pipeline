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

=item Z<>* ./libraries/fwd/seqs.fastq

=item Z<>* ./libraries/rev/seqs.fastq

=back

=item 4. Performs tag-cleaning of each sample-specific file. Required inputs:

=over

=item Z<>* ./demultiplexed/<sample_id>_*.fastq

=item Z<>* ./demultiplexed/<sample_id>_*.fastq

=back

=item 5. Runs the forward and reverse reads through the dada2 pipeline for the 
V3V4 16S rRNA gene region. Alternatively, analysis of the V4 or ITS region may
be specified.

=over

=item Z<>* ./<sample_id>_R1_tc.fastq[.gz]

=item Z<>* ./<sample_id>_R2_tc.fastq[.gz]

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

=item B<-h>, B<--help>

Print help message and exit successfully.

=item B<--working-dir>=PATH, B<-wd> PATH

Indicate an existing directory in which to place output, and from which the
run name will be parsed. The last directory on PATH must be named after the run.
(Many of the result files will be named after the run.)

=item B<--qsub-project>=space, B<-qp> space

Indicate which qsub-project space should be used for all qsubmissions.

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

=item B<--barcodes>=file

The full path to the Qiime-formatted mapping file. Default is in the 
config; use "./data/barcodes/UDI_ALL_INDEX_MAP_corrected.txt" when the ticket or gel 
images use phrases like "set A" "set B" "set C" "set D" "XT", or 
"normal APJR4 indices".

=item B<--1Step>

Use this flag if the data are prepared by 1-Step PCR (only forward & reverse read files
available). This processes the input files correctly and activates appropriate
parameters during adapter trimming, quality trimming/filtering, and denoising.

=back

=head2 BARCODE EXTRACTION AND DEMULTIPLEXING

=over

=item B<--extract_barcodes_rev_comp_args> "ARGS"

Specify the arguments --rev_comp_bc1 and/or --rev_comp_bc2 to QIIME 
extract_barcodes.py. This will override the pipeline's configured operation 
based on the instrument ID found in the FASTQ headers. (See config.json key
"instrument_ID_to_config"). Giving an empty string will force the pipeline to
reverse-complement neither index.

=item B<--extract_barcodes_append_args> "ARGS"

These arguments will be appended to the command calling QIIME 
extract_barcodes.py. This can be used to override some arguments, such as 
--bc1_len and --bc2_len. See http://qiime.org/scripts/extract_barcodes.html.

Without --onestep, the pipeline sets the length of both barcodes equal to the 
length of the first index.

=item B<--troubleshoot_barcodes>

Try all four permutations of reverse-complementing and switching the 
concatenation order of the indexes. Whichever of the four permutations yields a
successful demux is used. Note: Performing these transformations on the indexes 
may coincidentally yield a barcode that seems to be correct, even though the 
overall demux is incorrect. 

=item B<--split_libraries_args> "ARGS"

Command line options to append to every call to QIIME split_libraries_fastq.py.
See http://qiime.org/scripts/split_libraries_fastq.html.

=back

=head2 TRIMMING, FILTERING, AND DENOISING

=over

This pipeline stores configurations for various hypervariable regions in its 
config file, located at <pipeline_path>/config.json. Some configurations are 
specific to one-step PCR library prep too. When the targeted variable region is 
not described in the config file, use the following options to specify the 
processing parameters. These options will also override any parameters specified 
by the config file.

=item B<--var-reg>={V3V4, V4, ITS, OMPA}, B<-v> {V3V4, V4, ITS, OMPA}

The targeted variable region.

=item B<--fwd_primer>="SEQUENCE"

The primer(s) to trim from the forward reads. This value is passed on to 
bbduk.sh's "literal" option, so IUPAC ambiguity codes and comma-separated 
sequences are allowed.

=item B<--rev_primer>="SEQUENCE"

The primer(s) to trim from the reverse reads. See the note about B<--fwd_primer>.

=item B<--trim-maxlength>=LENGTH

After trimming primers, filter out reads longer than LENGTH. This value is
passed on to bbduk.sh's maxlen option. When primer trimming could yield various 
output read lengths, set this parameter to the maximum expected output read 
length.

=item B<--trim_poly_g>, B<--notrim_poly_g>

Trim poly G from the 5' end (or don't). Without either of these options, the
fastq header of one of the reads is examined before trimming; the starting letter
of the instrument ID determines poly-G trimming. See config.json for the exact
letters that trigger.

=item B<--dada2-truncQ>=SCORE

Trim reads at the first instance of a quality score less than or equal to SCORE.

=item B<--dada2-truncLen-f>=LENGTH

Forward reads will be trimmed to LENGTH. Reads shorter than LENGTH will be
removed.

=item B<--dada2-truncLen-r>=LENGTH

Reverse reads will be trimmed to LENGTH. Reads shorter than LENGTH will be 
removed.

=item B<--amplicon_length>=LENGTH

The expected amplicon length. When determining the ideal trimming lengths, any
combination that prevents mate pair merging (due to excessive trimming) will be 
avoided. Amplicon length is not needed when B<--dada2-truncLen-f> and 
B<--dada2-truncLen-R> are both given.

=item B<--dada2-minLen>=LENGTH

Remove reads with length less than LENGTH. This filter is enforced after 
trimming.

=item B<--dada2-minQ>=SCORE

After trimming, reads contain a quality score less than SCORE will be discarded.

=item B<--dada2-maxEE>=VALUE

After trimming, reads with higher than VALUE "expected errors" will be 
discarded. 

=item B<--dada2-rmPhix>, B<--no-dada2-rmPhix>

Remove reads aligning to the Phi.X genome. (Give --no-dada2-rmPhix to override
a default value in the config file.)

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

BEGIN {
    use File::Spec::Functions;
    use File::Basename;
    $pipelineDir = dirname(__FILE__);
    $scriptsDir  = catdir( $pipelineDir, "scripts" );

}
use lib $scriptsDir;    # .pm files in ./scripts/ can be loaded
use lib catdir( $pipelineDir, "lib", "perl5" );
use lib catdir( $pipelineDir, "lib" );

require Version;

use File::Find;
use Pod::Usage;
use English      qw( -no_match_vars );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Cwd          qw(abs_path getcwd);
use File::Temp   qw/ tempfile /;
use File::Temp   qw/ :mktemp  /;
use POSIX;
require File::Copy;
require List::Util;
use JSON;
require IO::Tee;
use File::Path            qw( remove_tree );
use File::Copy::Recursive qw(dircopy );
use Data::Dumper;
use Hash::Merge qw( merge );
use Path::Tiny  qw(path);
use Schedule::SGE;
use Capture::Tiny qw(capture_stderr);
use PathFinder;
use File::Copy::Recursive qw(rcopy_glob);
use Data::Dumper;
use Clone 'clone';
$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################
my @dbg;
my $oneStep       = 0;
my $dada2mem      = "30G";
my $tbshtBarcodes = 0;
my $delete        = 1;
my $trySkip       = 1;

GetOptions(
    "raw-path|i=s"                      => \my $raw_dir,
    "r1=s"                              => \my $r1,
    "r2=s"                              => \my $r2,
    "i1=s"                              => \my $i1,
    "i2=s"                              => \my $i2,
    "barcodes=s"                        => \my $map,
    "var-reg|v=s"                       => \my $var,
    "help|h!"                           => \my $help,
    "d|debug=s"                         => \@dbg,
    "verbose!"                          => \my $verbose,
    "dry-run!"                          => \my $dryRun,
    "skip!"                             => \$trySkip,
    "extract_barcodes_rev_comp_args=s"  => \my $extract_barcodes_revcomp_args,
    "extract_barcodes_append_args=s"    => \my $extract_barcodes_append_args,
    "split_libraries_args=s"            => \my $split_libraries_args,
    "fwd_primer=s"                      => \my $primer_L,
    "rev_primer=s"                      => \my $primer_R,
    "trim-maxlength=i"                  => \my $trimMaxLen,
    "trim_poly_g!"                      => \my $trim_poly_g,
    "amplicon_length=i"                 => \my $amplicon_length,
    "dada2-truncLen-f|for=i"            => \my $truncLenL,
    "dada2-truncLen-r|rev=i"            => \my $truncLenR,
    "dada2-truncQ=s"                    => \my $truncQ,
    "dada2-minLen=s"                    => \my $minLen,
    "dada2-minQ=s"                      => \my $minQ,
    "dada2-maxEE=s"                     => \my $maxEE,
    "dada2-rmPhix!"                     => \my $phix,
    "dada2_error_estimation_function=s" => \my $error_estimation_function,
    "dada2-mem:s"                       => \$dada2mem,
    "1Step!"                            => \$oneStep,
    "working-dir|wd=s"                  => \my $wd,
    "troubleshoot_barcodes!"            => \$tbshtBarcodes,
    "delete!"                           => \$delete,
    )
##add option for final resting place of important data

    or pod2usage( verbose => 0, exitstatus => 1 );

if ($help) {
    pod2usage( verbose => 2, exitstatus => 0 );
    exit 1;
}

####################################################################
## PARAMETER CHECKING PT 1
####################################################################

# Check existence/compatibility of required parameters
if (!$raw_dir
    && !(
           ( !$oneStep && $r1 && $r2 && $i1 && $i2 )
        || ( $oneStep && $r1 && $r2 )
    )
    )
{
    my $parameters = $oneStep ? "-r1, -r2" : "-r1, -r2, -i1, -i2";
    die "\n\tPlease provide the location of the raw sequencing files "
        . "(single directory => -i)\n\t\tOR \n\tFull paths to each raw file => "
        . "$parameters)\n\n";
}
if ( $raw_dir && ( $r1 || $r2 || $i1 || $i2 ) ) {
    die
        "\n\tInput directory (-i) and raw files (-r*|-i*) were both specified. Please "
        . "provide one or the other.\n\n";
}
if ( $oneStep && ( $i1 || $i2 ) ) {
    die "\n\t--1Step is incompatible with -i1 and -i2.\n\n";
}
if ( !$wd ) {
    die "\n***Please choose a working directory (-wd).";
}

if ( $truncLenL && !$truncLenR ) {
    die "***\nPlease provide truncation lengths for forward and reverse "
        . "reads\n";
}

# Refine and validate all variables that refer to the filesystem
my (@paths) = (
    { path => \$wd,      name => "Working directory" },
    { path => \$raw_dir, name => "Raw file directory" },
    { path => \$r1,      name => "Raw forward reads file" },
    { path => \$r2,      name => "Raw reverse reads file" },
    { path => \$i1,      name => "Raw forward index file" },
    { path => \$i2,      name => "Raw reverse index file" },
    { path => \$map,     name => "Mapping file" }
);
foreach (@paths) {
    if ( ${ $$_{path} } ) {    # If the variable is a non-empty string...
        my $copy = ${ $$_{path} };
        $copy =~ s/\/$//;              # remove any trailing slash
        my $relPath = $copy;
        $copy = abs_path($relPath);    # get abs path
          # At the same time, abs_path returns undef if the path doesn't exist
          # so we can verify the existence of each file and directory
        if ( !defined $copy ) {
            die $$_{name}
                . " not found. Looked for "
                . $relPath . ".\n"
                . "Current working directory: "
                . getcwd() . "\n";
        }
        ${ $$_{path} } = $copy
            ;    # external variable referenced by path has now been edited.
    }
}

# After shaving off potential trailing slash,
# Split on the last path separator, then check if the last string is a word
$wd =~ s/\/$//;
if ( ( File::Spec->splitpath("$wd") )[2] !~ /[^\/]/ ) {
    die "Working directory (-wd) must be a path to a directory";
}

# now we can start logging to the log file
my $GLOBAL_PATHS = PathFinder->new( wd => $wd, pipeline_dir => $pipelineDir );

my ( $readsForInput, $readsRevInput, $index1Input, $index2Input );
if ( !@dbg || grep( /^barcodes$/, @dbg ) || grep( /^demux$/, @dbg ) ) {

  # Find the input files, and find their putative locations if the index files
  # were to be decompressed to our working directory.
  # In one-step runs, $index1Input and $index2Input are the SAME files
  # pointed to by $readsForInput and $readsRevInput
    ( $readsForInput, $readsRevInput, $index1Input, $index2Input )
        = $raw_dir ? find_raw_files(
        dir          => $raw_dir,
        pf           => $GLOBAL_PATHS,
        pcr_one_step => $oneStep
        )
        : $oneStep ? ( $r1, $r2, $r1, $r2 )
        :            ( $r1, $r2, $i1, $i2 );
}

####################################################################
## INITIALIZATION
###################################################################
$GLOBAL_PATHS->{'raw_dir'} = $raw_dir;
$GLOBAL_PATHS->{'r1'}      = $r1;
$GLOBAL_PATHS->{'r2'}      = $r2;
$GLOBAL_PATHS->{'i1'}      = $i1;
$GLOBAL_PATHS->{'i2'}      = $i2;
$GLOBAL_PATHS->set_config( config() )
    ;    # could even put sub config() in GLOBAL_PATHS???

chdir $GLOBAL_PATHS->{'wd'};

# if we're in the global run directory, share all files with group
my $global_run_storage
    = abs_path( $GLOBAL_PATHS->get_config()->{"run_storage_path"} );
if ( $GLOBAL_PATHS->{'wd'} =~ /^${global_run_storage}/ ) {
    umask(0002);
}

my $run  = basename( $GLOBAL_PATHS->{'wd'} );
my $time = strftime( "%Y-%m-%d %H:%M:%S", localtime(time) );

local $SIG{__WARN__} = sub {
    $GLOBAL_PATHS->print("WARNING: $_[0]");
    print STDERR "WARNING: $_[0]";
};    # Warnings go to log file, stdout, and stderr

# Now that we're done checking parameters, I want the more serious error
# messages to be printed to the log file.
$SIG{__DIE__} = sub {
    die @_
        if $^S; # If die was called inside an eval, just allow that to happen
                # otherwise print to the log file and then allow die to happen
    $GLOBAL_PATHS->print( @_ . "\n" );
    return 1;
};

print "Logging to: " . $GLOBAL_PATHS->part1_log() . "\n";
$GLOBAL_PATHS->print( "PIPELINE VERSION: " . Version::version() . "\n" );
$GLOBAL_PATHS->print("$time\n");

my @dropouts;

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

if ( !-e $GLOBAL_PATHS->part1_error_log() ) {
    mkdir $GLOBAL_PATHS->part1_error_log();
    print $GLOBAL_PATHS->part1_error_log()
        . " did not exist -> making "
        . $GLOBAL_PATHS->part1_error_log() . "\n";
}
else {
    # Remove old error logs
    my @oldLogs = glob( $GLOBAL_PATHS->part1_error_log() . "/*.*" );
    my @cmds;
    foreach my $log (@oldLogs) {
        my @dontDelete = (
            $GLOBAL_PATHS->part1_error_log() . "/illumina_dada2.pl.stderr", );
        if (List::Util::none { $_ eq $log }
            @dontDelete
            )
        {
            push @cmds, "rm $log";
        }
    }
    if (@cmds) {
        execute_and_log(
            cmds    => \@cmds,
            dry_run => $dryRun,
            msg     => "Removing stale error logs."
        );
    }
}

if ( !-e $GLOBAL_PATHS->part1_stdout_log() ) {
    mkdir $GLOBAL_PATHS->part1_stdout_log();
    print $GLOBAL_PATHS->part1_stdout_log()
        . " did not exist -> making "
        . $GLOBAL_PATHS->part1_stdout_log() . "\n";
}

my $metadata = project_metadata();
$metadata = project_metadata(
    metadata => { "params" => { "platform" => "ILLUMINA" } } );

my $qiime = "$GLOBAL_PATHS->{'wd'}/${run}_qiime_config.txt";

if ( $tbshtBarcodes && @dbg ) {
    $GLOBAL_PATHS->print(
        "Disabling --debug sections because --troubleshoot_barcodes was given.\n"
    );
    @dbg = ();
}

if (@dbg) {
    $GLOBAL_PATHS->print("DBG FLAGS: ");
    for (@dbg) {
        $GLOBAL_PATHS->print("$_ ");
    }

    my $e  = grep( /^barcodes$/,     @dbg );
    my $de = grep( /^demux$/,        @dbg );
    my $s  = grep( /^splitsamples$/, @dbg );
    my $t  = grep( /^tagclean$/,     @dbg );
    my $da = grep( /^dada2$/,        @dbg );
    if ( $e + $de + $s + $t + $da != scalar @dbg ) {
        die "Illegal debug option. Legal debug options are "
            . "barcodes, demux, splitsamples, tagclean, and dada2.";
    }
    if ($trySkip) {
        $trySkip = 0;
        $GLOBAL_PATHS->print(
            "\nDisabling skipping because --debug sections were selected.\n");
    }

    $GLOBAL_PATHS->print("\n");
}
$GLOBAL_PATHS->print("RUN: $run\n");
if ($oneStep) {
    $GLOBAL_PATHS->print("PCR PREPARATION METHOD: 1-Step\n");
}
else {
    $GLOBAL_PATHS->print("PCR PREPARATION METHOD: 2-Step\n");
}

# Resolve all parameters based on instrument model
my $first_fastq;
if ( !@dbg || grep( /^barcodes$/, @dbg ) || grep( /^demux$/, @dbg ) ) {
    $first_fastq = $index1Input;
}
elsif ( grep( /^splitsamples$/, @dbg ) ) {
    $first_fastq = $GLOBAL_PATHS->fwd_library();
}
elsif ( grep( /^tagclean$/, @dbg ) ) {
    $first_fastq = ( $GLOBAL_PATHS->splitsamples_output_fwd() )[0];
}
elsif ( grep( /^dada2$/, @dbg ) ) {
    $first_fastq = ( $GLOBAL_PATHS->trimmed_files() )[0];
}

my $instrument_id = get_instrument_id( fastq => $first_fastq );
$GLOBAL_PATHS->print("INSTRUMENT ID: ${instrument_id}");
my $instrument_profile = instrument_id_to_profile(
    config => $GLOBAL_PATHS->get_config(),
    id     => $instrument_id
);
if ( $GLOBAL_PATHS->load_config_profile("$instrument_profile") ) {
    $GLOBAL_PATHS->print(
        " (loading config profile \"$instrument_profile\")\n");
}
else {
    $GLOBAL_PATHS->print(" (no stored config profile)\n");
}

$GLOBAL_PATHS->override_default_param(
    value => $extract_barcodes_revcomp_args,
    name  => "extract_barcodes_rev_comp_args",
);
$GLOBAL_PATHS->override_default_param( value => $map, name => "barcodes", );
$GLOBAL_PATHS->override_default_param(
    value => $split_libraries_args,
    name  => "split_libraries_args",
);

# Resolve all parameters for the variable region
$var = defined $var ? uc $var : "";
$GLOBAL_PATHS->print("VARIABLE REGION: $var");
my $pcrstep = $oneStep ? "onestep" : "twostep";
if ( $GLOBAL_PATHS->load_config_profile("${pcrstep}_${var}") ) {
    $GLOBAL_PATHS->print(" (loading config profile \"${pcrstep}_${var}\")\n");
}
else {
    $GLOBAL_PATHS->print(" (requires parameters given on command line)\n");
}

# Resolve primer trimming parameters
$GLOBAL_PATHS->override_default_param(
    value => $primer_L,
    name  => "fwd primer"
);
$GLOBAL_PATHS->override_default_param(
    value => $primer_R,
    name  => "rev primer"
);

$GLOBAL_PATHS->print(
    "R VERSION: " . $GLOBAL_PATHS->get_config()->{'R'} . "\n" );

# Resolve DADA2 quality-trim and denoise parameters
# might need to put these next two in eval{ } blocks
$GLOBAL_PATHS->override_default_param(
    value    => $truncLenL,
    name     => "fwd trim length",
    optional => 1
);
$GLOBAL_PATHS->override_default_param(
    value    => $truncLenR,
    name     => "rev trim length",
    optional => 1
);
if (   !$GLOBAL_PATHS->get_param("fwd trim length")
    || !$GLOBAL_PATHS->get_param("rev trim length") )
{
    # if either of the trim lengths missing...
    if ( defined($amplicon_length) ) {
        $GLOBAL_PATHS->print("AMPLICON LENGTH: $amplicon_length\n");
        $GLOBAL_PATHS->print("TRIMMING PARAMETERS WILL BE OPTIMIZED");
        if (   $GLOBAL_PATHS->get_param("fwd trim length")
            || $GLOBAL_PATHS->get_param("rev trim length") )
        {
            $GLOBAL_PATHS->print(
                " (Any trim length parameters given will be ignored)\n");
        }
        else {
            $GLOBAL_PATHS->print("\n");
        }
    }
    else {
        die
            "TRIMMING PARAMETERS NOT PROVIDED. CANNOT OPTIMIZE TRIM LENGTH BECAUSE AMPLICON LENGTH NOT KNOWN\n";
    }
}
elsif ( defined($amplicon_length) ) {
    $GLOBAL_PATHS->print(
        "Amplicon length not needed because both trim lengths are known. Ignoring.\n"
    );
}

# process the rest of configurable params
my $key_values = {
    "trim quality"   => $truncQ,
    "min length"     => $minLen,
    "filter quality" => $minQ,
    "max EE"         => $maxEE,
    "remove phix"    => $phix
};
foreach ( keys %{$key_values} ) {
    $GLOBAL_PATHS->override_default_param(
        value => $key_values->{$_},
        name  => "$_",
    );
}

$GLOBAL_PATHS->localize_file_params();
my @params = (
    "bc_len",
    "extract_barcodes_rev_comp_args",
    "barcodes",
    "split_libraries_args",
    "fastqc_limits",
    "fwd primer",
    "rev primer",
    "primer_allowed_mm",
    "trim_poly_g",
    "max length",
    "fwd trim length",
    "rev trim length",
    "remove phix",
    "trim quality",
    "filter quality",
    "max EE",
    "min length",
    "dada2_error_estimation_function"
);
$GLOBAL_PATHS->print("\nFINAL PARAMS:\n");
foreach (@params) {
    $GLOBAL_PATHS->print( uc($_) . ": " );
    if ( defined $GLOBAL_PATHS->get_param($_) ) {
        $GLOBAL_PATHS->print( $GLOBAL_PATHS->get_param($_) . "\n" );
    }
    else {
        $GLOBAL_PATHS->print("UNKNOWN (non-terminal)\n");
    }
}

my $skipMe;

if (   !@dbg
    || grep( /^barcodes$/,     @dbg )
    || grep( /^demux$/,        @dbg )
    || grep( /^splitsamples$/, @dbg ) )
{
    ###### BEGIN CHECK OF QIIME CONFIGURATION ###########
    #####################################################
    $GLOBAL_PATHS->print("\n");
    my $cmd
        = qiime_cmd( 'print_qiime_config.py', $GLOBAL_PATHS->get_config() )
        . " > $qiime";
    execute_and_log(
        cmds    => [$cmd],
        logger  => $GLOBAL_PATHS,
        dry_run => $dryRun,
        msg     => "Printing QIIME configuration details\n"
    );
    $GLOBAL_PATHS->print("See: $qiime\n");

    my $map_log = catfile( $GLOBAL_PATHS->part1_error_log(),
        basename( $GLOBAL_PATHS->get_param("barcodes"), ".txt" ) . ".log" );

    if ($trySkip) {
        $skipMe = skippable(
            inputs        => [ $GLOBAL_PATHS->get_param("barcodes") ],
            outputs       => [ $GLOBAL_PATHS->part1_local_map(), $map_log ],
            checksums_in  => $metadata->{"checkpoints"}{"map"},
            checksums_out => $metadata->{"checkpoints"}{"meta"},
            step_name     => "map validation",
            previous_step_skipped => 1,
            pf                    => $GLOBAL_PATHS
        );
    }

    if ( !$skipMe ) {
        $GLOBAL_PATHS->print( "Copying "
                . $GLOBAL_PATHS->get_param("barcodes") . " to "
                . $GLOBAL_PATHS->part1_local_map()
                . "\n" );
        File::Copy::copy(
            $GLOBAL_PATHS->get_param("barcodes"),
            $GLOBAL_PATHS->part1_local_map()
        ) or die "Copy failed: $!";
        preprocess_map( $GLOBAL_PATHS->part1_local_map(), $run );

        ###### BEGIN VALIDATION OF MAPPING FILE ###########
        ################################################
        $GLOBAL_PATHS->print(
            "MAPPING FILE: " . $GLOBAL_PATHS->part1_local_map() . "\n" );

        my $cmd
            = qiime_cmd( 'validate_mapping_file.py',
            $GLOBAL_PATHS->get_config() )
            . " -m "
            . $GLOBAL_PATHS->part1_local_map()
            . " -s -o "
            . $GLOBAL_PATHS->part1_error_log();
        execute_and_log(
            cmds    => [$cmd],
            logger  => $GLOBAL_PATHS,
            dry_run => $dryRun,
            msg     => "Validating map from "
                . $GLOBAL_PATHS->get_param("barcodes") . "\n"
        );
        my $mappingError
            = glob( $GLOBAL_PATHS->part1_error_log() . "/*.log" );
        if ($mappingError) {
            open MAPERROR, "<$mappingError"
                or die "Cannot open $mappingError for "
                . "reading: $OS_ERROR";
            $_ = <MAPERROR>;    # gets first line

            chomp;
            if ( $_ =~ /No errors or warnings found in mapping file./ ) {
                $GLOBAL_PATHS->print("---Map passed validation.\n");

            }
            else {
                while (<MAPERROR>) {
                    if ( $_ =~ m/^Errors -----------------------------$/ ) {
                        my $nextLine = <MAPERROR>;
                        if ( $nextLine
                            =~ m/^Warnings ---------------------------$/ )
                        {
                            $GLOBAL_PATHS->print(
                                "---Warnings during map validation. See "
                                    . "$mappingError for details.\n" );
                        }
                        else {
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
        }
        else {
            die
                "validate_mapping_file.py terminated but did not signal success. Normally a success message is printed in its error file.";
        }
        File::Copy::move( $GLOBAL_PATHS->part1_local_map() . "_corrected.txt",
            $GLOBAL_PATHS->part1_local_map() );
        project_metadata( metadata =>
                { "map" => { "file" => $GLOBAL_PATHS->part1_local_map() } } );

        ###### BEGIN EVALUATION OF SAMPLES VIA MAPPING FILE ###########
        ###############################################################
        open MAP, "<" . $GLOBAL_PATHS->part1_local_map()
            or die "Cannot open "
            . $GLOBAL_PATHS->part1_local_map()
            . " for reading: $OS_ERROR";
        my $nSamples = 0;
        while (<MAP>) {
            chomp;
            ## don't count header as sample
            if ( $. > 1 ) {

                #  don't count any line if it doesn't start with four
                # tab-separated fields; the first three must contain
                # non-whitespace chars
                if ( $_ =~ /^(\S+\t){3}/ ) {
                    $nSamples++;
                }
                elsif ( $_ =~ /\S/ ) {

                  # QIIME's validate_mapping_file seems to already check this:
                    die
                        "In mapping file the line $. does not have four tab-separated fields.";
                }
            }
        }
        close MAP;

        $GLOBAL_PATHS->print("NO. SAMPLES: $nSamples\n\n");

        if ( !@dbg ) {
            cacheChecksums(
                files     => [ $GLOBAL_PATHS->get_param("barcodes") ],
                step_name => "map",
                pf        => $GLOBAL_PATHS
            );
            cacheChecksums(
                files => [ $mappingError, $GLOBAL_PATHS->part1_local_map() ],
                step_name => "meta",
                pf        => $GLOBAL_PATHS
            );
        }
    }
    else {

   # It is possible the first run had a map error, but the user failed to fix.
   # skippable would return true, but the log file would still inform us of
   # the error.
        my $mappingError
            = glob( catfile( $GLOBAL_PATHS->part1_error_log(), "*.log" ) );

        open MAPERROR, "<", $mappingError
            or die "Cannot open $mappingError for " . "reading: $OS_ERROR";
        $_ = <MAPERROR>;
        close MAPERROR;
        chomp;
        if ( !$_ =~ /No errors or warnings found in mapping file./ ) {
            die "***Unresolved error in mapping file. See $mappingError for "
                . "details. Exiting.\n";
        }
        $GLOBAL_PATHS->print(
            "Mapping file has been validated already. Moving on.\n");
    }

}
else {
    $GLOBAL_PATHS->print("NOT USING QIIME\n");
}

###### BEGIN BARCODES ##########
#######################################
my $newSamNo;
if ($tbshtBarcodes) {

    # make three other directories
    my %dirs_params = (
        "$GLOBAL_PATHS->{'wd'}_rc_fwd"     => "--rev_comp_bc1",
        "$GLOBAL_PATHS->{'wd'}_rc_rev"     => "--rev_comp_bc2",
        "$GLOBAL_PATHS->{'wd'}_rc_fwd_rev" => "--rev_comp_bc1 --rev_comp_bc2",
        "$GLOBAL_PATHS->{'wd'}"            => ""
    );

    # initialize directories by copying all contents from cwd
    for my $dir ( keys %dirs_params ) {

        sub init {
            my $dir = shift;
            if ( $dir ne $GLOBAL_PATHS->{'wd'} ) {
                remove_tree($dir) if ( -e $dir );
                mkdir $dir or die "Unable to create $dir\n";
                rcopy_glob( "$GLOBAL_PATHS->{'wd'}/*", $dir )
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
            . "\n" );
    my $summary = "";

    for my $dir ( keys %dirs_params ) {

        sub barcodes_demux {
            my %arg             = @_;
            my $dir             = delete $arg{wd};
            my $barcode_ori     = delete $arg{barcode_ori};
            my $switch_barcodes = delete $arg{switch_barcodes} // 0;

     # clone this and replace wd because demux can take place in any directory
            my $paths = $GLOBAL_PATHS->clone();
            $paths->set_wd($dir);
            $GLOBAL_PATHS->print("Now working in: $dir.\n");
            $GLOBAL_PATHS->print( "Please refer to " . $paths->part1_log() );

            # Extract barcodes
            if ( !$switch_barcodes ) {
                barcodes(
                    oriParams   => $barcode_ori,
                    pathfinder  => $paths,
                    index1      => $index1Input,
                    index2      => $index2Input,
                    reads1      => $readsForInput,
                    reads2      => $readsRevInput,
                    bc_len      => $GLOBAL_PATHS->get_param("bc_len"),
                    append_args => $extract_barcodes_append_args
                );
            }
            else {
                barcodes(
                    oriParams   => $barcode_ori,
                    pathfinder  => $paths,
                    index1      => $index2Input,
                    index2      => $index1Input,
                    reads1      => $readsForInput,
                    reads2      => $readsRevInput,
                    bc_len      => $GLOBAL_PATHS->get_param("bc_len"),
                    append_args => $extract_barcodes_append_args
                );
            }

            # demux, and check for success
            $success = demux(
                pathfinder  => $paths,
                die_on_fail => 0,
                append_args =>
                    $GLOBAL_PATHS->get_param("split_libraries_args")
            );
            if ($success) {
                my $msg
                    = "Demux may have succeeded when extract_barcodes.py called with "
                    . "command-line options \"${barcode_ori}\"";
                $msg .= " and switched indexes" if $switch_barcodes;
                $msg .= "\n";
                $GLOBAL_PATHS->print("\n$msg\n");
                $summary .= $msg;
            }
            else {
                my $msg = "Demux failed when extract_barcodes.py called with "
                    . "command-line options \"${barcode_ori}\"";
                $msg     .= " and switched indexes" if $switch_barcodes;
                $msg     .= "\n";
                $summary .= $msg;
            }
            $paths->close();
        }
        barcodes_demux( wd => $dir, barcode_ori => $dirs_params{$dir} );
        barcodes_demux(
            wd              => "${dir}_switched",
            barcode_ori     => $dirs_params{$dir},
            switch_barcodes => 1
        );

    }

    my $msg
        = "\nFinished troubleshooting barcodes.\n"
        . "Please inspect the split library logs in each trial directory to determine which demux was correct.\n"
        . "Correct files can be moved back to this run directory, and the pipeline can be continued using '--debug splitsamples --debug tagclean --debug dada2'.\n\n";
    $GLOBAL_PATHS->print($msg);
    $GLOBAL_PATHS->print("SUMMARY:\n$summary\n");
    exit 0;

}
else {
    if ( ( !@dbg ) || grep( /^barcodes$/, @dbg ) ) {
        barcodes(
            oriParams =>
                $GLOBAL_PATHS->get_param("extract_barcodes_rev_comp_args"),
            pathfinder  => $GLOBAL_PATHS,
            index1      => $index1Input,
            index2      => $index2Input,
            reads1      => $readsForInput,
            reads2      => $readsRevInput,
            bc_len      => $GLOBAL_PATHS->get_param("bc_len"),
            append_args => $extract_barcodes_append_args
        );
    }
    if ( !@dbg || grep( /^demux$/, @dbg ) ) {
        demux(
            pathfinder  => $GLOBAL_PATHS,
            die_on_file => 1,
            append_args => $GLOBAL_PATHS->get_param("split_libraries_args")
        );
    }

    # Then continue with code below the two subroutines
}

sub barcodes {
    my %arg         = @_;
    my $oriParams   = delete $arg{oriParams}  // '';
    my $paths       = delete $arg{pathfinder} // $GLOBAL_PATHS;
    my $index1      = delete $arg{index1};
    my $index2      = delete $arg{index2};
    my $reads1      = delete $arg{reads1};
    my $reads2      = delete $arg{reads2};
    my $bcLen       = delete $arg{bc_len}  // '';
    my $oneStep     = delete $arg{onestep} // 0;  # try to handle this outside
    my $append_args = delete $arg{append_args} // "";

    my $outputDir = $paths->{'wd'};
    my $barcodes  = catfile( $outputDir, "barcodes.fastq" );

    my $nSamples = count_samples( $paths->part1_local_map() );

    ## change to full path to full barcodes (flag)
    my $count = 0;

    $paths->print( "Forward and reverse reads will be obtained from:\n"
            . "\t$reads1\n\t$reads2\n" );
    $paths->print( "Index 1 and Index 2 will be obtained from:\n"
            . "\t$index1\n\t$index2\n" );
    my @files = ( $reads1, $reads2, $index1, $index2 );
    my @names = ( "RF", "RR", "I1", "I2" );

    my $input  = \@files;
    my $output = [$barcodes];
    if ($trySkip) {
        $skipMe = skippable(
            inputs                => $input,
            outputs               => $output,
            checksums_in          => $metadata->{"checkpoints"}{"raw"},
            checksums_out         => $metadata->{"checkpoints"}{"barcodes"},
            step_name             => "barcode extraction",
            in_names              => \@names,
            previous_step_skipped => 1,
            pf                    => $GLOBAL_PATHS
        );
    }

    my $step1;

    if ( !$skipMe ) {

        my $start = time;
        $step1 = "extract_barcodes.py";

        # Get index files as .fastq
        my @cmds = ();
        if ( $index1 !~ /\.gz$/ )    # please change this to for-loop...
        {

            # if the index files aren't .gz, just read in place.
        }
        else {                       # Otherwise...
                                     # decompress to our run directory
            my $index1_decompressed = catfile( $outputDir, "I1.fastq" );
            push( @cmds,
                "gzip --decompress --force < $index1 > $index1_decompressed"
            );
            $paths->print(
                "---Decompressing $run barcode and index file from $index1 to $index1_decompressed\n"
            );
            $index1 = $index1_decompressed;
        }
        if ( $index2 !~ /.gz$/ ) {    # same for reverse reads
        }
        else {
            my $index2_decompressed = catfile( $outputDir, "I2.fastq" );
            push( @cmds,
                "gzip --decompress --force < $index2 > $index2_decompressed"
            );
            $paths->print(
                "---Decompressing $run barcode and index file from $index2 to $index2_decompressed\n"
            );
            $index2 = $index2_decompressed;
        }

        # Execute the commands queued above
        if ( scalar @cmds > 0 ) {
            execute_and_log(
                cmds    => \@cmds,
                dry_run => $dryRun,
                logger  => $paths,
                msg     =>
                    "Decompressing raw files containing indexes to run directory...\n",
                verbose => 1
            );
            $paths->print("---Raw index file(s) decompressed.\n");
            @cmds = ();
        }

        # Sanity check
        my $rawLnsFwd = count_lines($index1);
        my $rawLnsRev = count_lines($index2);
        if ( $rawLnsFwd == $rawLnsRev ) {
            $paths->print(
                ( $rawLnsFwd / 4 . " barcoded reads to process.\n" ) );
        }
        else {
            $paths->print(
                "Warning: raw index files have differing numbers of lines.\n"
            );
        }

        if ( !$bcLen ) {
            $bcLen = first_seq_length($index1);
            $paths->print("Detected barcode length as $bcLen \n");
        }

        my $cmd
            = qiime_cmd( 'extract_barcodes.py', $GLOBAL_PATHS->get_config() )
            . " -f $index1 -r $index2 -c barcode_paired_end --bc1_len $bcLen "
            . "--bc2_len $bcLen -o $outputDir $oriParams $append_args";

        execute_and_log(
            cmds    => [$cmd],
            logger  => $paths,
            dry_run => $dryRun,
            msg     => "Waiting for barcode extraction to complete...\n",
            verbose => 1
        );

        # Sanity check
        my $barcodeLns = count_lines("$barcodes");
        if ( $rawLnsFwd != $barcodeLns ) {
            die
                "Forward index file had $rawLnsFwd lines but barcodes.fastq had $barcodeLns lines.";
        }

        my $duration = time - $start;
        $paths->print( "---Barcode extraction complete... "
                . ( $barcodeLns / 4 )
                . " barcodes extracted\n" );
        $paths->print("---Duration of barcode extraction: $duration s\n");

        # Record hash of barcodes.fastq, and inputs
        if ( !@dbg ) {
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

      # Delete unneeded output of barcodes.py. extract_barcodes.py outputs
      # each of the supplied FASTQ's with the barcodes removed, and names them
      # 'reads1.fastq' and 'reads2.fastq'
        @cmds = ();
        push @cmds, "rm -rf $outputDir/reads1.fastq";
        push @cmds, "rm -rf $outputDir/reads2.fastq";

        # If two-step, delete the temporarily decompressed index files now
        if ( !$oneStep ) {
            if ( -e catfile( $outputDir, "I1.fastq" ) ) {
                push @cmds, "rm -rf " . catfile( $outputDir, "I1.fastq" );
            }
            if ( -e catfile( $outputDir, "I2.fastq" ) ) {
                push @cmds, "rm -rf " . catfile( $outputDir, "I2.fastq" );
            }
        }
        execute_and_log(
            cmds    => \@cmds,
            dry_run => $dryRun,
            msg     => "Cleaning up after extract_barcodes.py...\n",
            logger  => $paths
        );

        if ( @dbg && !grep( /^demux$/, @dbg ) ) {
            $paths->print( "Finished extracting barcodes. Terminated "
                    . "because --debug demux was not specified.\n" );
            exit 0;
        }
    }
    else {
        $paths->print( "Barcodes already extracted as "
                . basename($barcodes)
                . ". Moving on.\n" );
    }
}

###### BEGIN SPLIT LIBRARIES ##########
#######################################

sub demux {
    my %arg         = @_;
    my $paths       = delete $arg{pathfinder}  // $GLOBAL_PATHS;
    my $append_args = delete $arg{append_args} // "";
    my $die_on_fail = delete $arg{die_on_fail} // 0;

    $paths->print("Demuxing in: $paths->{'wd'}\n");

    my $barcodes = "$paths->{'wd'}/barcodes.fastq";
    my $nSamples = count_samples( $paths->part1_local_map() );
    my $step2;
    my $step3;

    my $split_log
        = catfile( $paths->fwd_demux_dir(), "split_library_log.txt" );

    my @cmds;

   # split_libraries_fastq.py accepts fastq.gz, so no need to
   # convert_to_local_if_gz... Unless it's one-step, in which case we can save
   # time by using the already-decompressed raw files.
    if ($oneStep) {
        ( $readsForInput, $readsRevInput )
            = convert_to_local_if_gz( $paths->{'wd'}, $readsForInput,
            $readsRevInput );
    }

    if ($trySkip) {
        $skipMe = skippable(
            inputs => [
                $readsForInput, $readsRevInput,

                # FIXME need to copy local map to other directory?
                $paths->part1_local_map(), $barcodes
            ],
            outputs      => [ $paths->fwd_library(), $paths->rev_library() ],
            checksums_in => {
                "RF"       => $metadata->{"checkpoints"}{"raw"}->{"RF"},
                "RR"       => $metadata->{"checkpoints"}{"raw"}->{"RR"},
                "map"      => $metadata->{"checkpoints"}{"map"},
                "barcodes" => $metadata->{"checkpoints"}{"barcodes"}
            },
            checksums_out         => $metadata->{"checkpoints"}{"library"},
            step_name             => "demultiplexing",
            in_names              => [ "RF", "RR", "map", "barcodes" ],
            previous_step_skipped => 1,
            pf                    => $paths
        );
    }

    if ( !$skipMe ) {
        my $start  = time;
        my $step2  = "split_libraries_fastq.py";
        my $script = qiime_cmd( $step2, $GLOBAL_PATHS->get_config() );

 # qiime's split_libraries_fastq accepts some keywords too, such as "golay_12"
        my $barcodeType = first_seq_length($barcodes);

        @cmds = ();
        my $cmd
            = $GLOBAL_PATHS->get_config()->{'executor'}
            . " -N $step2 -e "
            . $paths->part1_error_log() . " -o "
            . $paths->part1_stdout_log()
            . " $script -i $readsForInput -o "
            . $paths->fwd_demux_dir()
            . " -b $barcodes -m "
            . $paths->part1_local_map()
            . " --max_barcode_errors 1 --store_demultiplexed_fastq --barcode_type $barcodeType -r 999 -n 999 -q 0 -p 0.0001";
        $cmd = join( " ", $cmd, $append_args );
        push @cmds, $cmd;
        $cmd
            = $GLOBAL_PATHS->get_config()->{'executor'}
            . " -N $step2 -e "
            . $paths->part1_error_log() . " -o "
            . $paths->part1_stdout_log()
            . " $script -i $readsRevInput -o "
            . $paths->rev_demux_dir()
            . " -b $barcodes -m "
            . $paths->part1_local_map()
            . " --max_barcode_errors 1 --store_demultiplexed_fastq --barcode_type $barcodeType -r 999 -n 999 -q 0 -p 0.0001";
        $cmd = join( " ", $cmd, $append_args );
        push @cmds, $cmd;

        execute_and_log(
            cmds    => \@cmds,
            logger  => $paths,
            dry_run => $dryRun,
            msg     => "Demultiplexing to get the run library...\n"
        );

        $paths->print(
            "---Waiting for fwd and rev seqs.fastq to complete....\n");
        $paths->print("---Monitoring $step2 error logs....\n");
        $paths->check_error_log( prefix => $step2 );
        while (!( -e $paths->fwd_library() )
            || !( -e $paths->rev_library() ) )
        {
            sleep 1;

            # check_error_log allows pipeline to terminate if the qsubbed
            # split_libraries_fastq.py prints error
            $paths->check_error_log( prefix => $step2 );
        }
        my $duration = time - $start;
        $paths->print(
            "---Duration of fwd and rev seqs.fastq production: $duration s\n"
        );

        ##Check split_library log for 0's
        my @dropouts = readSplitLog(
            file    => $split_log,
            verbose => 1,
            logger  => $paths
        );

        $newSamNo = $nSamples - scalar @dropouts;
        if ( scalar @dropouts > 0 && scalar @dropouts ne $nSamples ) {
            $paths->print( "---The following "
                    . scalar @dropouts
                    . " samples returned 0 "
                    . "reads, due to either a lack of barcodes or a filtering of those "
                    . "reads in split_libraries_fastq.py:\n" );
            foreach my $x (@dropouts) {
                $paths->print("   $x\n");
            }
            $paths->print(
                      "---Number of samples after split_libraries_fastq.py: "
                    . "$newSamNo\n" );
        }
        elsif ( scalar @dropouts == $nSamples ) {
            if ($die_on_fail) {
                $paths->print(
                    "---No samples were successfully demultiplexed. Is the "
                        . "mapping file correct? Exiting.\n" );
                die;
            }
            else {
                return 0;
            }
        }
        elsif ( scalar @dropouts == 0 ) {
            $paths->print("---Reads from all samples were demultiplexed.\n");
        }

        @cmds = ();
        push @cmds,
              catfile( $scriptsDir, "get_split_library_stats.sh" ) . " < "
            . catfile( $paths->fwd_demux_dir(), "split_library_log.txt" )
            . " > "
            . catfile( $paths->fwd_demux_dir(), "split_library_stats.txt" );
        push @cmds,
              catfile( $scriptsDir, "get_split_library_stats.sh" ) . " < "
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

        if ( $dada2mem eq "30G" ) {    # if it's the default
            open( my $sll, "<",
                ${paths}->fwd_demux_dir() . "/split_library_log.txt" )
                or die "Can't open fwd split library log: $!";
            while ( my $line = <$sll> ) {
                if ( $line =~ m/^Total number seqs written\s+(\d+)/ ) {
                    $dada2mem = List::Util::max(
                        POSIX::ceil( 4.2 * log($1) - 4.2 ) / 10, 30 );
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
              $GLOBAL_PATHS->get_config()->{'executor'}
            . " -l mem_free=300M -e "
            . $paths->part1_error_log() . " -o "
            . $paths->part1_stdout_log()
            . " $binary --limits "
            . $GLOBAL_PATHS->get_param("fastqc_limits")
            . " --outdir "
            . $paths->fwd_demux_dir() . " "
            . $paths->fwd_library();
        push @cmds,
              $GLOBAL_PATHS->get_config()->{'executor'}
            . " -l mem_free=300M -e "
            . $paths->part1_error_log() . " -o "
            . $paths->part1_stdout_log()
            . " $binary --limits "
            . $GLOBAL_PATHS->get_param("fastqc_limits")
            . " --outdir "
            . $paths->rev_demux_dir() . " "
            . $paths->rev_library();
        execute_and_log(
            cmds    => \@cmds,
            logger  => $paths,
            dry_run => $dryRun,
            msg     => "Assessing library quality...\n"
        );

        # Remove temporarily decompressed files
        if ($oneStep) {
            my @cmds = ();
            push( @cmds, "rm -rf $readsForInput" );
            push( @cmds, "rm -rf $readsRevInput" );
            execute_and_log(
                cmds    => \@cmds,
                dry_run => $dryRun,
                logger  => $paths,
                msg     =>
                    "---Removing decompressed raw files from run directory\n"
            );
        }

        if ( !@dbg ) {
            cacheChecksums(
                files     => [ $paths->fwd_library(), $paths->rev_library() ],
                step_name => "library",
                pf        => $paths
            );
        }
        project_metadata(
            metadata => { "map" => { "file" => $paths->part1_local_map() } }
        );

    }
    else {
        $paths->print( "Library already produced as "
                . basename( $paths->fwd_library() ) . " and "
                . basename( $paths->rev_library() )
                . ". Moving on.\n" );
        my @dropouts = readSplitLog(
            file    => $split_log,
            verbose => 0,
            logger  => $paths
        );

# Still need number of demuxed samples for later, even though we didn't actually demux this time
        $newSamNo = $nSamples - scalar @dropouts;
    }

    if ( @dbg && !grep( /^splitsamples$/, @dbg ) ) {
        die "Finished demultiplexing libaries. Terminated "
            . "because -d splitsamples was not specified.";
    }
    return 1;
}

if ( !@dbg || grep( /^splitsamples$/, @dbg ) ) {
    ###### BEGIN SPLIT BY SAMPLE ##########
    #######################################

    # need to look for both fastqs and fastq.gz's because after tagcleaning
    # all raw demuxed files are gzipped. Look for either
    {
        my @fwdOutput = $GLOBAL_PATHS->splitsamples_output_fwd();
        my @revOutput = $GLOBAL_PATHS->splitsamples_output_rev();

        if ($trySkip) {
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
                outputs       => [ @fwdOutput, @revOutput ],
                checksums_in  => $metadata->{"checkpoints"}{"library"},
                checksums_out => $metadata->{"checkpoints"}{"samples"},
                step_name     => "sample splitting",
                previous_step_skipped => 1,
                pf                    => $GLOBAL_PATHS
            );
        }
    }

    if ( !$skipMe ) {
        my $step3  = "split_sequence_file_on_sample_ids.py";
        my $script = qiime_cmd( $step3, $GLOBAL_PATHS->get_config() );
        my @cmds   = ();

        my $fwd_dir = mkdtemp("demultiplex_XXXXXX");
        my $rev_dir = mkdtemp("demultiplex_XXXXXX");

        push @cmds,
              "$script -i "
            . $GLOBAL_PATHS->fwd_library()
            . " --file_type fastq -o $fwd_dir";
        push @cmds,
              " $script -i "
            . $GLOBAL_PATHS->rev_library()
            . " --file_type fastq -o $rev_dir";

        execute_and_log(
            cmds    => \@cmds,
            logger  => $GLOBAL_PATHS,
            dry_run => $dryRun,
            msg     => "Splitting $run seqs.fastq " . "files by sample ID\n",
            qsub    => 1,
            mem_GB  => "5",
            name    => $step3,
            paths   => $GLOBAL_PATHS,
            config  => $GLOBAL_PATHS->get_config(),
        );

        $GLOBAL_PATHS->print("Temp fwd dir: $fwd_dir\n") if $verbose;
        $GLOBAL_PATHS->print("Temp rev dir: $rev_dir\n") if $verbose;
        my @fwdOutput   = glob( catfile( $fwd_dir, "*.fastq" ) );
        my @revOutput   = glob( catfile( $rev_dir, "*.fastq" ) );
        my $msg_printed = 0;
        while ( scalar @fwdOutput != scalar @revOutput ) {
            if ( !$msg_printed ) {
                $GLOBAL_PATHS->print(
                    "Waiting for all forward and reverse sample FASTQ's to appear in $fwd_dir and $rev_dir.\n"
                );
                $msg_printed = 1;
            }
            @fwdOutput = glob( catfile( $fwd_dir, "*.fastq" ) );
            @revOutput = glob( catfile( $rev_dir, "*.fastq" ) );
        }

        # Append _R1 or _R2 to each filename
        $GLOBAL_PATHS->print(
            "Appending _R1 or _R2 to each demuxed FASTQ filename.\n");
        remove_tree( $GLOBAL_PATHS->sample_dir() )
            if ( -e $GLOBAL_PATHS->sample_dir() );
        mkdir $GLOBAL_PATHS->sample_dir()
            or die( "Can't create directory \""
                . $GLOBAL_PATHS->sample_dir()
                . "\": $!\n" );
        foreach my $oldname (@fwdOutput) {
            my ( $name, $path, $suffix ) = fileparse( $oldname, ".fastq" );
            my $newname
                = catfile( $GLOBAL_PATHS->sample_dir(), "${name}_R1.fastq" );
            File::Copy::move( $oldname, $newname )
                or die("Can't move $oldname to $newname: $!\n");
            $oldname = $newname;
        }
        foreach my $oldname (@revOutput) {
            my ( $name, $path, $suffix ) = fileparse( $oldname, ".fastq" );
            my $newname
                = catfile( $GLOBAL_PATHS->sample_dir(), "${name}_R2.fastq" );
            File::Copy::move( $oldname, $newname )
                or die("Can't move $oldname to $newname: $!\n");
            $oldname = $newname;
        }
        remove_tree($fwd_dir);
        remove_tree($rev_dir);

        gzip_sge(
            pattern => catfile( $GLOBAL_PATHS->sample_dir(), "*.fastq" ),
            paths   => $GLOBAL_PATHS,
            config  => $GLOBAL_PATHS->get_config()
        );

        $GLOBAL_PATHS->print("---Demultiplexed FASTQ's compressed.\n");

        $GLOBAL_PATHS->print( "---Samples demultiplexed in "
                . $GLOBAL_PATHS->sample_dir()
                . "\n" );
        my @finalFwdOutput = $GLOBAL_PATHS->splitsamples_output_fwd();
        my @finalRevOutput = $GLOBAL_PATHS->splitsamples_output_rev();

        if ( !@dbg ) {
            cacheChecksums(
                files     => [ @finalFwdOutput, @finalRevOutput ],
                step_name => "samples",
                pf        => $GLOBAL_PATHS
            );
        }
    }
    else {
        $GLOBAL_PATHS->print(
            "Library already split to $newSamNo sample-specific files. Moving on.\n"
        );
    }

    # Remove temporarily decompressed files
    if ($oneStep) {
        my ( $readsForInput, $readsRevInput, $index1Input, $index2Input )
            = $GLOBAL_PATHS->{'raw_dir'} ? find_raw_files(
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
        push( @cmds, "rm -rf $readsForInput" );
        push( @cmds, "rm -rf $readsRevInput" );
        execute_and_log(
            cmds    => \@cmds,
            dry_run => $dryRun,
            msg     =>
                "---Removing decompressed raw files from $GLOBAL_PATHS->{'wd'}\n"
        );
    }

    if ( @dbg && !grep( /^tagclean$/, @dbg ) ) {
        $GLOBAL_PATHS->print(
            "Finished splitting library to sample-specific FASTQs. Terminated "
                . "because --debug tagclean was not specified.\n" );
        exit 0;
    }
}

###### BEGIN TAGCLEANING ##########
###################################

my $start = time;
if ( !@dbg || grep( /^tagclean$/, @dbg ) ) {
    my $do = 1;

    my @inputsF = $GLOBAL_PATHS->splitsamples_output_fwd();
    my @inputsR = $GLOBAL_PATHS->splitsamples_output_rev();

    if ( scalar @inputsF != scalar @inputsR ) {
        die "Unequal numbers of forward and reverse sample FASTQ's.\n";
    }
    my $newSamNo = scalar @inputsF;

    my @fwdTcFiles = glob("$GLOBAL_PATHS->{'wd'}/*R1_tc.fastq.gz");
    my @revTcFiles = glob("$GLOBAL_PATHS->{'wd'}/*R2_tc.fastq.gz");

    if ($trySkip) {
        $skipMe = skippable(
            inputs                => [ ( @inputsF, @inputsR ) ],
            outputs               => [ @fwdTcFiles, @revTcFiles ],
            checksums_in          => $metadata->{"checkpoints"}{"samples"},
            checksums_out         => $metadata->{"checkpoints"}{"tagcleaned"},
            step_name             => "primer removal",
            previous_step_skipped => 1,
            pf                    => $GLOBAL_PATHS
        );
    }

    if ( !$skipMe ) {
        my @cmds;
        my $bbduk = $GLOBAL_PATHS->get_config()->{"part1"}->{"bbduk.sh"};
        foreach my $in_F (@inputsF) {
            my ( $name, $dir, $ext ) = fileparse( $in_F, qr{\.gz} );
            my @suffixes = ("_R1.fastq");
            my $prefix   = basename( $name, @suffixes );
            my $in_R     = catfile( $GLOBAL_PATHS->sample_dir(),
                "${prefix}_R2.fastq${ext}" );
            my $out_F = catfile( $GLOBAL_PATHS->{'wd'},
                "${prefix}_R1_tc.fastq.gz" );
            my $out_R = catfile( $GLOBAL_PATHS->{'wd'},
                "${prefix}_R2_tc.fastq.gz" );
            my $stats_F = catfile( $GLOBAL_PATHS->part1_stdout_log(),
                "${prefix}_R1_tc.stats" );
            my $stats_R = catfile( $GLOBAL_PATHS->part1_stdout_log(),
                "${prefix}_R2_tc.stats" );
            my $stats_poly_g = catfile(
                $GLOBAL_PATHS->part1_stdout_log(),
                "${prefix}_poly_g_tc.stats"
            );

            sub primer_length_range {
                my $primers        = shift;
                my @primers        = split ",", $primers;
                my @primer_lengths = map length, @primers;
                return (
                    List::Util::min(@primer_lengths),
                    List::Util::max(@primer_lengths)
                );
            }

            sub get_read_length {
                my $fastq_gz = shift;
                my $read     = `gunzip < $fastq_gz | head -n 2 | tail -n 1`;
                chomp $read;
                return length($read);
            }

            my ( $k_F, $max_fwd_primer )
                = primer_length_range(
                $GLOBAL_PATHS->get_param("fwd primer") );
            my ( $k_R, $max_rev_primer )
                = primer_length_range(
                $GLOBAL_PATHS->get_param("rev primer") );

            # We use maxlen purely to filter out untrimmed reads, so just
            # set this to readlength - 1 by default, since all reads are
            # same length
            my $readlen = get_read_length($in_F);
            my $maxlen  = $trimMaxLen ? $trimMaxLen : ( $readlen - 1 );

            my $hdist = $GLOBAL_PATHS->get_param("primer_allowed_mm");
            my $mink_F
                = List::Util::min( ( $max_fwd_primer - $hdist, $k_F ) );
            my $mink_R
                = List::Util::min( ( $max_rev_primer - $hdist, $k_R ) );

            # trim fwd reads first (skipr2) then trim rev reads (skipr1)
            my $cmd
                = "$bbduk -Xmx4915m -Xms4327m in=$in_F in2=$in_R literal="
                . $GLOBAL_PATHS->get_param("fwd primer")
                . " copyundefined out=stdout.fastq stats=$stats_F overwrite=t ziplevel=9 ktrim=l k=$k_F rcomp=f hdist=$hdist mink=$mink_F hdist2=0 skipr2=t restrictleft=30 minlen=60";
            $cmd
                .= " 2>>"
                . catfile( $GLOBAL_PATHS->part1_error_log(),
                "bbduk.sh.stderr" );
            $cmd
                .= " | $bbduk -Xmx4915m -Xms4327m in=stdin.fastq interleaved=t literal="
                . $GLOBAL_PATHS->get_param("rev primer")
                . " copyundefined out=$out_F out2=$out_R stats=$stats_R overwrite=t ziplevel=9 ktrim=l k=$k_R rcomp=f hdist=$hdist mink=$mink_R hdist2=0 skipr1=t minlen=60 restrictleft=30 maxlen=$maxlen";
            $cmd
                .= " 2>>"
                . catfile( $GLOBAL_PATHS->part1_error_log(),
                "bbduk.sh.stderr" );
            if ( $GLOBAL_PATHS->get_param("trim_poly_g") ) {
                $cmd =~ s/out=${out_F} out2=${out_R}/out=stdout\.fastq/g;
                $cmd
                    .= " | $bbduk -Xmx4915m -Xms4327m in=stdin.fastq interleaved=t literal=GGGGGGGGGGGGGGGGGGGG copyundefined out=$out_F out2=$out_R stats=$stats_poly_g overwrite=t ziplevel=9 ktrim=r k=20 rcomp=f hdist=0 mink=1 hdist2=0 minlen=60";
                $cmd .= " 2>>"
                    . catfile( $GLOBAL_PATHS->part1_error_log(),
                    "bbduk.sh.stderr" );
            }
            push @cmds, $cmd;
        }

        # }

       # Do not submit bbduk.sh to grid. Not sure why. One thing that needs to
       # to be changed before submitting to grid is pre-pending the executor
       # string.
        execute_and_log(
            cmds    => \@cmds,
            logger  => $GLOBAL_PATHS,
            dry_run => $dryRun,
            msg     => "Removing $var primers from all sequences.\n",
            qsub    => 0
        );

        my @fwdFiles = glob("$GLOBAL_PATHS->{'wd'}/*R1_tc.fastq.gz");
        for ( my $i = 0; $i < scalar @fwdFiles; $i++ ) {
            my $revFile = $fwdFiles[$i] =~ s/R1_tc.fastq.gz$/R2_tc.fastq.gz/r;
            if ( -z $fwdFiles[$i] || -z $revFile ) {
                if ( -z $fwdFiles[$i] ) {
                    $GLOBAL_PATHS->print(
                        "Primer trimming filtered all reads from $fwdFiles[$i]. Deleting both sample FASTQs.\n"
                    );
                }
                else {
                    $GLOBAL_PATHS->print(
                        "Primer trimming filtered all reads from $revFile. Deleting both sample FASTQs.\n"
                    );
                }
                unlink $fwdFiles[$i];
                unlink $revFile;
            }
        }

        if ( !@dbg ) {
            cacheChecksums(
                files =>
                    [ glob( catfile( $GLOBAL_PATHS->sample_dir(), "*" ) ) ],
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
    }
    else {
        $GLOBAL_PATHS->print(
            "---Primers already removed from forward and reverse reads of $newSamNo samples. Moving on.\n"
        );
    }

    if ( @dbg && !grep( /^dada2$/, @dbg ) ) {
        $GLOBAL_PATHS->print(
            "Finished removing primers. Terminated because --debug dada2 was not "
                . "specified.\n" );
        exit 0;
    }
}

my $stats_file  = "dada2_part1_stats.txt";
my $counts_file = "dada2_abundance_table.rds";

###### BEGIN DADA2 ##########
#############################
if ( ( !@dbg ) || grep( /^dada2$/, @dbg ) ) {
    my @outputs = ( $stats_file, $counts_file );
    my @inputs  = glob("$GLOBAL_PATHS->{'wd'}/*R[1|2]_tc.fastq.gz");

    if ($trySkip) {
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

    if ( !$skipMe ) {
        my $r_out = dada2(
            truncLenL => $GLOBAL_PATHS->get_param("fwd trim length"),
            truncLenR => $GLOBAL_PATHS->get_param("rev trim length"),
            maxEE     => $GLOBAL_PATHS->get_param("max EE"),
            truncQ    => $GLOBAL_PATHS->get_param("trim quality"),
            "rm.phix" => $GLOBAL_PATHS->get_param("remove phix"),
            minLen    => $GLOBAL_PATHS->get_param("min length"),
            minQ      => $GLOBAL_PATHS->get_param("filter quality"),
            error_estimation_function =>
                $GLOBAL_PATHS->get_param("dada2_error_estimation_function"),
            amplicon_length => $amplicon_length,
            dada2mem        => $dada2mem,
            config          => $GLOBAL_PATHS->get_config(),
            paths           => $GLOBAL_PATHS,
            multithread => 32
        );

###### EVALUATING DADA2 OUTPUT ##########
#########################################
        my @removed;
        ## do a scan of this file for those samples not passing filtering (0 sequences surviving - and log this)

        open ROUT, "<${r_out}"
            or die "Cannot open ${r_out} for reading: $OS_ERROR";
        while (<ROUT>) {
            chomp;
            if ( $_ =~ /The filter removed all reads/ ) {
                push @removed, $_;
            }
        }
        close ROUT;

        if ( scalar @removed > 0 ) {
            $GLOBAL_PATHS->print(
                scalar @removed
                    . " samples removed during dada2 filtering:\n" );
            for my $x (@removed) {
                $GLOBAL_PATHS->print("$x\n");
            }
        }
        $GLOBAL_PATHS->print("\n");

        if ( List::Util::all { -e $_ } @outputs ) {

            # "dada2 completed successfully" is a key phrase that causes
            # an appropriate message printed to STDOUT
            $GLOBAL_PATHS->print(
                      "dada2 completed successfully!\nAbundance table for "
                    . "run $run located at $GLOBAL_PATHS->{'wd'}/dada2_abundance_table.rds\n"
            );
            $GLOBAL_PATHS->print(
                      "See $outputs[0] for dada2 table of reads surviving by "
                    . "step\n\n" );

            if ( !@dbg ) {
                cacheChecksums(
                    files     => \@outputs,
                    step_name => "dada2part1out",
                    pf        => $GLOBAL_PATHS
                );
            }
        }
        else {
            $GLOBAL_PATHS->print(
                "---dada2 did not complete successfully, something went wrong!\n"
                    . "---Check R output file (usually R.stderr in SGE stderr directory).\n"
            );
            die;
        }

    }
    else {
        $GLOBAL_PATHS->print(
            "DADA2 already processed the same tagcleaned files during a previous run. Doing nothing.\n"
        );
    }
}

if ( $delete && -e $stats_file && -e $counts_file ) {
    my @delete_files = (
        "$GLOBAL_PATHS->{'wd'}/barcodes.fastq",
        catfile( $GLOBAL_PATHS->fwd_demux_dir(), "seqs.fna" ),
        catfile( $GLOBAL_PATHS->rev_demux_dir(), "seqs.fna" ),
        $GLOBAL_PATHS->fwd_library(),
        $GLOBAL_PATHS->rev_library(),
        tagcleaned_files(),
        glob(
            catfile( $GLOBAL_PATHS->filtered_dir(), "*_filt\\.fastq\\.gz" )
        )
    );

    for my $file (@delete_files) {
        unlink $file;
    }
    rmdir $GLOBAL_PATHS->filtered_dir();
}

###### COMPLETING $logFH FILE ##############
#########################################

END {
    if ( defined ${global_run_storage} ) {
        if ( $GLOBAL_PATHS->{'wd'} =~ /^${global_run_storage}/ ) {
            find(
                sub {
                    chmod( 0664, $File::Find::name ) unless -d;
                    chmod( 0775, $File::Find::name ) if -d;
                },
                $GLOBAL_PATHS->{'wd'}
            );
            chmod( 0775, $GLOBAL_PATHS->{'wd'} );
        }
    }

    if ( defined $GLOBAL_PATHS ) {
        $GLOBAL_PATHS->print("\n\n\n");
        $GLOBAL_PATHS->close();
    }

}

####################################################################
##                               SUBS
####################################################################
sub read_json {
    my $file = shift;
    my $mode = shift;

    my $hash = {};

    if ( -e $file ) {

        # get existing metadata on filesystem
        my $json;
        {
            local $/;    #Enable 'slurp' mode
            open( my $FH, $mode, $file );
            seek $FH, 0, 0 or die;
            $json = <$FH>;
            close $FH;
        }

        eval { $hash = decode_json($json); };
    }

    return $hash;
}

sub config {
    my %arg             = @_;
    my $new_params      = delete $arg{new_params} // {};
    my $orig_param_file = "$pipelineDir/config.json";

    my $params = read_json( $orig_param_file, "<" );

    return $params;
}

sub project_metadata {
    my %arg          = @_;
    my $any_args     = scalar %arg;
    my $new_metadata = delete %arg{"metadata"} // {};
    my $replace      = delete $arg{"replace"}  // 0;
    my $pf           = delete %arg{"pf"}       // $GLOBAL_PATHS;

    my $metadataFile = "$pf->{'wd'}/.meta.json";
    my $old_metadata = read_json( $metadataFile, "+<" );

    # If no arguments, user wants to GET existing metadata
    if ( !$any_args ) {
        return $old_metadata;
    }
    else {
        # If arguments, user wants to UPDATE AND GET existing metadata
        my $combined_metadata;

        if ( !$replace ) {

            # recursively merge params and checkpoints with $project_metadata
            Hash::Merge::set_behavior('LEFT_PRECEDENT');
            $combined_metadata = merge( $old_metadata, $new_metadata );
        }
        else {
            # replace each present element of new_metadata
            $combined_metadata = $old_metadata;
            foreach my $step ( keys %$new_metadata ) {
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

sub preprocess_map {
    my $map_filepath = shift;
    my $run_name     = shift;

    $run_name =~ s/[^a-zA-Z0-9\.]/\./g;

    rename( $map_filepath, $map_filepath . '.bak' );
    open( IN,  '<' . $map_filepath . '.bak' ) or die $!;
    open( OUT, '>' . $map_filepath )          or die $!;
    while (<IN>) {
        s/^([a-zA-Z0-9\.]+)\t/$run_name\.$1\t/g;
        print OUT $_;
    }
    close(IN);
    close(OUT);
    unlink( $map_filepath . '.bak' );

    return;
}

sub qiime_cmd {
    my $script = shift;
    my $config = shift;

    my $python
        = catfile( abs_path( $config->{'part1'}->{"qiime1/bin"} ), "python" );
    $script
        = catfile( abs_path( $config->{'part1'}->{"qiime1/bin"} ), $script );
    return "$python -s $script";
}

#' @_[0] Filenames of inputs to checksum
#' @_[1] Filenames of outputs to checksum (not guaranteed to be existing files)
#' @_[2] A hash pointer containing the checksums for inputs
#' @_[3] A hash pointer containing the checksums for outputs
#' @_[4] The name of the current step, for log messages.
#' @_[5] The names of the inputs (should match keys of @_[2])
#' @_[6] The names of the outputs (should match keys of @_[3])
sub skippable {
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
        for grep { not defined $checksumsIn{$_} } keys %checksumsIn;
    delete $checksumsOut{$_}
        for grep { not defined $checksumsOut{$_} } keys %checksumsOut;

    if ($verbose) {
        $pf->print( "Input files: " . Dumper($inputsR) . "\n" );
        $pf->print( "Cached files: " . Dumper(%checksumsIn) . "\n" );
        $pf->print( "Output files: " . Dumper($outputsR) . "\n" );
        $pf->print( "Cached files: " . Dumper(%checksumsOut) . "\n" );
    }

    if (   scalar @$outputsR == 0
        || List::Util::none { defined $_ } @$outputsR
        || List::Util::none { -e $_ } @$outputsR )
    {
        $pf->print( ucfirst("$step has not been run yet. Running now...\n") );
        return 0;
    }    # no output files found by caller

    if ( !List::Util::all { -e $_ } @$outputsR ) {

        # Caller named output files but not all existed
        $pf->print(
            ucfirst(
                "$step output is incomplete. Removing outputs files and executing $step...\n"
            )
        );
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
        $inNames = [ map { File::Spec->abs2rel( $_, $GLOBAL_PATHS->{'wd'} ) }
                @$inputsR ];
    }
    if ( !defined $outNames ) {
        $outNames = [ map { File::Spec->abs2rel( $_, $GLOBAL_PATHS->{'wd'} ) }
                @$outputsR ];
    }

    if ( !setIdent( $inNames, [ keys %checksumsIn ] ) ) {

        # The input files are named differently from last time (or pipeline
        # provided inconsistent file nicknames)
        $pf->print(
            "Input to $step does not match previous records. Discarding $step ",
            "outputs and repeating $step.\n"
        );
        unlink @$outputsR;
        return 0;
    }
    if ( !setIdent( $outNames, [ keys %checksumsOut ] ) ) {

# The outputs files are named differently from last time (or pipeline provided inconsistent file nicknames)
        $pf->print(
            ucfirst($step)
                . " output does not match previous records. Discarding $step outputs and repeating $step...\n"
        );
        unlink @$outputsR;
        return 0;
    }

    $pf->print(
        ucfirst("$step") . " output exists. Trying to validate...\n" );

    if ( !$prevSkip ) {

        # If the previous step was not skipped, the inputs must be checksummed
        if ($verbose) {
            $pf->print( "Checksumming input files:\n" . Dumper($inputsR) );
        }
        my $i          = -1;
        my $inputValid = List::Util::all {
            if ($verbose) {
                $pf->print("Checksumming input file: $_\n");
            }
            $i++;
            my @stdout   = split( /\s/, `sha512sum $_` );
            my $checksum = $stdout[0];

            $checksum eq $checksumsIn{ $inNames->[$i] };
        }
        @$inputsR;

        # If all inputs are/nt valid give appropriate messages
        if ($inputValid) {
            $pf->print("Input to $step matches previous records.\n");
        }
        else {
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
    if ($verbose) {
        $pf->print( "Checksumming output files:\n" . Dumper($outputsR) );
    }
    my $outputValid = List::Util::all {
        $i++;
        if ($verbose) {

            $pf->print("Checksumming output file:$_\n");
        }
        my @stdout   = split( /\s/, `sha512sum $_` );
        my $checksum = $stdout[0];

        $checksum eq $checksumsOut{ $outNames->[$i] };
    }
    @$outputsR;

    # If all outputs are/nt valid give appropriate messages
    if ($outputValid) {
        $pf->print("Output to $step matches previous records.\n");
        return 1;
    }
    else {
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
sub cacheChecksums {
    my %arg   = @_;
    my $files = delete %arg{"files"};
    my $step  = delete %arg{"step_name"};
    my $pf    = delete %arg{"pf"};
    my $names = delete %arg{"file_names"}
        // [ map { File::Spec->abs2rel( $_, $pf->{'wd'} ) } @$files ];

    $pf->print("Saving progress in run folder...\n");

    my $destroyNextStep = sub {
        if ( $step2number{$step} + 1 < scalar keys %number2step ) {

           # Reset the next set of checksums so the next step is forced to run
            $metadata->{"checkpoints"}
                { $number2step{ $step2number{$step} + 1 } } = {};
        }
    };
    my $update = 0;
    my @keys   = keys %{ $metadata->{"checkpoints"}{$step} };
    if ( !setIdent( \@keys, $names ) ) {
        $update = 1;
        $destroyNextStep->()
            ;    # MUST DO THIS BEFORE ANY NEW CHECKSUMS ARE WRITTEN
         # PREVENTS "CHIMERIC" CHECKSUM FILE IN CASE OF PIPELINE INTERRUPTION HERE
    }
    my %newChecksums;

    for my $i ( 0 .. scalar @$files - 1 ) {
        my $file = $files->[$i];
        if ( !-e -f $file ) {
            die "Non-existant file passed to cacheChecksums(): $file\n";
        }
        my $name     = $names->[$i];
        my @output   = split( /\s/, `sha512sum $file` );
        my $checksum = $output[0];
        if ( !$update )
        { # !$update means that $checksums{$step}->{ $file } is defined for all files
            if ( $metadata->{"checkpoints"}{$step}->{$name} ne $checksum ) {
                $destroyNextStep->();
            }
        }

        # Store the checksum of the files just produced.
        $newChecksums{$name} = $checksum;
    }

    $metadata->{"checkpoints"}{$step} = \%newChecksums;
    if ( List::Util::all { defined $_ } $metadata->{"checkpoints"} ) {
        project_metadata(
            metadata => { "checkpoints" => $metadata->{"checkpoints"} },
            pf       => $pf,
            replace  => 1
        );
        $pf->print("\n");
    }
}

sub get_instrument_id {
    my %arg   = @_;
    my $fastq = delete %arg{"fastq"};

    my $id;
    if ( $fastq =~ /\.gz$/ ) {
        $id = `gunzip < $fastq | head -n 1`;
    }
    else {
        open my $seqs, $fastq or die "Could not open $fastq: $!";
        while (<$seqs>) {
            $id = $_;
            last;
        }
        close $seqs;
    }
    chomp $id;

    $id =~ s/^@//;      # Remove the first @
    $id =~ s/:.*$//;    # Remove everything after the first colon separator
     # QIIME split_libraries_fastq.py inserts a sample_seq field before the instrument ID
    $id =~ s/^.* //
        ; # if there are any spaces, remove everything up to and including the last one
    return ($id);
}

sub instrument_id_to_profile {
    my %arg    = @_;
    my $config = delete %arg{"config"};
    my $id     = delete %arg{"id"};

    my $instrument_id_to_profile = $config->{"instrument_ID_to_profile"};
    foreach my $id_pattern ( keys %{$instrument_id_to_profile} ) {
        if ( $id =~ $id_pattern ) {
            return ( $instrument_id_to_profile->{$id_pattern} );
        }
    }

    die "Unrecognized instrument ID: $id. Recognized patterns: "
        . ( join ", ", keys %{$instrument_id_to_profile} )
        . ". Please update config.json.\n";
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
    my $file = shift;
    my $ans  = count_lines($file) - 1;
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

    if ($oneStep) {
        if ( scalar @r1s == 1 && scalar @r2s == 1 ) {
            $index1Input = $readsForInput = $r1s[0];
            $index2Input = $readsRevInput = $r2s[0];
        }
        else {
            my $files_found = "";
            my @file_arrays = ( \@r1s, \@r2s );
            foreach my $file_array_ref (@file_arrays) {
                my @file_array = @$file_array_ref;
                if ( scalar @file_array ) {
                    $files_found .= ( join "\n", @file_array );
                }
            }

            my $message
                = "Could not find a complete and exclusive set of raw files."
                . " Since --1step given, input directory must have exactly one"
                . " R1 and R2 file. Files found:\n"
                . $files_found;

            $pf->print($message);

            die $message;

        }
    }
    else {
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
            && scalar @r4s == 0 )
        {
            $index1Input   = $i1s[0];
            $index2Input   = $i2s[0];
            $readsForInput = $r1s[0];
            $readsRevInput = $r2s[0];
        }
        elsif (scalar @i1s == 0
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
        }
        else {
            my $files_found = "";
            my @file_arrays = ( \@i1s, \@i2s, \@r1s, \@r2s, \@r3s, \@r4s );
            foreach my $file_array_ref (@file_arrays) {
                my @file_array = @$file_array_ref;
                if ( scalar @file_array ) {
                    $files_found .= ( join "\n", @file_array );
                }
            }
            my $message
                = "Could not find a complete and exclusive set of raw files. Input directory must"
                . " contain exactly one set of I1, I2, R1, and R2, OR R1, R2, R3, and R4.\n"
                . "File found in $dir:\n$files_found";
            $pf->print($message);
            die $message;
        }
    }
    return ( $readsForInput, $readsRevInput, $index1Input, $index2Input );
}

# Given *R1.fastq(.gz), if it's .gz, then removes the extension and gives a
# filepath in the local directory
# @param 0 The run directory
# @params 1... The full path to the original file
sub convert_to_local_if_gz {
    my $wd    = shift;
    my @files = @_;
    my @ans;
    foreach my $file (@files) {
        if ( $file =~ /\.gz$/ ) {

            # Rename *[R|I][1|2].fastq.gz to <WD>/<RUN>_[R|I][1|2].fastq
            my $dest = basename($file);

            my $suffix = substr( $dest, ( length($dest) - 11 ), 8 )
                ;    # get suffix (sans .gz)

            $wd =~ s/\/$//;    # remove any trailing slash
            my @dirs = File::Spec->splitdir($wd);
            $file = "$wd/$dirs[scalar(@dirs) - 1]" . "_$suffix";
        }
        push( @ans, $file );
    }
    return @ans;
}

sub first_seq_length {
    my $file = shift;
    open my $reads, $file or die "Could not open $file: $!";
    my $first_index;
    while (<$reads>) {
        chomp;
        $first_index = $_ if $. == 2;
        last              if $. == 2;
    }
    close $reads;

    return ( length $first_index );
}

sub readSplitLog {
    my %arg       = @_;
    my $pf        = delete %arg{"logger"};
    my $split_log = delete %arg{"file"};
    my $verbose   = delete %arg{"verbose"} // undef;

    ##Check split_library log for 0's
    open SPLIT, "<$split_log"
        or die "Cannot open $split_log for writing: " . "$OS_ERROR";
    my @dropouts;
    while (<SPLIT>) {
        if ( $_ =~ /\t/ ) {
            chomp;
            my ( $sample, $nReads ) = split, /\t/;
            chomp $nReads;
            if ( $nReads eq "0" ) {
                push @dropouts, $sample;
            }
            else {
                if ($verbose) {
                    $pf->print("$_\n");
                }
            }
        }
    }
    close SPLIT;
    return @dropouts;
}

sub gzip_sge {
    my %arg    = @_;
    my $patt   = delete %arg{"pattern"};
    my $paths  = delete %arg{"paths"};
    my $config = delete %arg{"config"};

    # Make file list and script for SGE array job
    # NB demuxed_seqs.lst is the name required by gzip_array_template.sh
    my $file_list = catfile( $paths->commands_dir(), "demuxed_seqs.lst" );
    open( my $FH, '>', $file_list ) or die $!;
    my $cmd = "find $patt -maxdepth 1 -type f -printf \"\%p\\n\"";
    $FH->print(`$cmd`);
    $FH->close();

    my $array_script
        = catfile( $paths->commands_dir(), "gzip_array_template.sh" );
    File::Copy::copy(
        catfile( $paths->pipeline_dir(), "data", "gzip_array_template.sh" ),
        $array_script )
        or die "Copy failed: $!";
    my $file = path($array_script);

    my $data   = $file->slurp_utf8;
    my $nFiles = count_lines($file_list);
    $data =~ s/^#\$ -t 1-.*$/#\$ -t 1-$nFiles/gm;
    $data =~ s/demuxed_seqs.lst/$file_list/g;
    $file->spew_utf8($data);

    $config = clone($config);
    $config->{"qsub_b"} = "n";

    execute_and_log(
        cmds    => [$array_script],
        logger  => $GLOBAL_PATHS,
        dry_run => $dryRun,
        msg     => "Compressing demultiplexed FASTQ's...\n",
        qsub    => 1,
        mem_GB  => "1",
        name    => "gzip",
        paths   => $GLOBAL_PATHS,
        config  => $config,
    );
}

sub dada2 {
    my %arg                       = @_;
    my $truncLenL                 = delete %arg{"truncLenL"};
    my $truncLenR                 = delete %arg{"truncLenR"};
    my $maxEE                     = delete %arg{"maxEE"};
    my $truncQ                    = delete %arg{"truncQ"};
    my $rm_phix                   = delete %arg{"rm.phix"};
    my $minLen                    = delete %arg{"minLen"};
    my $minQ                      = delete %arg{"minQ"};
    my $error_estimation_function = delete %arg{"error_estimation_function"};
    my $amplicon_length           = delete %arg{"amplicon_length"};
    my $dada2mem                  = delete %arg{"dada2mem"};
    my $config                    = delete %arg{"config"};
    my $paths                     = delete %arg{"paths"};
    my $multithread               = delete %arg{"multithread"} // undef;

    $dada2mem =~ s/[a-zA-Z]//g;

    my $Rscript
        = catfile( $pipelineDir, "scripts", "filter_and_denoise_illumina.R" );
    my $args = "--maxN=0 --maxEE=$maxEE --truncQ=$truncQ --memory=$dada2mem";
    $args .= " --rm.phix" if ($rm_phix);
    $args .= " --multithread $multithread" if ($multithread);
    $args .= " --error_estimation_function=$error_estimation_function"
        if ($error_estimation_function);
    if ( $truncLenL && $truncLenR ) {
        $args .= " --truncLenL=$truncLenL --truncLenR=$truncLenR";
    }
    else {
        $args
            .= " --optimize_trim --amplicon_length=$amplicon_length --verbose";
    }

    chdir $GLOBAL_PATHS->{'wd'};

    my $exitStatus   = 1;
    my $previous_mem = 0;
    my $try          = 1;
    my $R_out = catfile( $GLOBAL_PATHS->part1_error_log(), "R.stderr" );

    while ( $exitStatus == 1 && ( $dada2mem != $previous_mem || $try < 3 ) ) {
        my $cmd
            = "rm -rf "
            . $GLOBAL_PATHS->filtered_dir() . " "
            . catfile( $GLOBAL_PATHS->{'wd'}, "dada2_part1_stats.txt " )
            . catfile( $GLOBAL_PATHS->{'wd'}, "dada2_abundance_table.rds " )
            . catfile( $GLOBAL_PATHS->{'wd'}, ".RData" );
        execute_and_log(
            cmds    => [$cmd],
            dry_run => $dryRun,
            msg     =>
                "Removing old filtered fastq files, stats, and Rout files from previous runs\n"
        );

        $cmd
            = "-e ${R_out} -o "
            . $GLOBAL_PATHS->part1_stdout_log()
            . " -N Rscript \""
            . $GLOBAL_PATHS->get_config()->{R}
            . "script $Rscript $args\"";
        execute_and_log(
            cmds    => [$cmd],
            logger  => $GLOBAL_PATHS,
            dry_run => $dryRun,
            msg     =>
                "Running DADA2 with fastq files in $GLOBAL_PATHS->{'wd'} for $var region...\n",
            qsub   => 1,
            mem_GB => "${dada2mem}",
            config => $config,
            paths  => $paths,
            name   => "R_dada2"
        );

        if ( -e "dada2_part1_stats.txt" ) {    # sign of success
            $GLOBAL_PATHS->print("---R script completed without errors.\n");
            $GLOBAL_PATHS->print(
                      "---DADA2-specific commands can be found in "
                    . "$Rscript\n" );
            $exitStatus = 0;                   # Move on from DADA2
        }
        else {
            open IN, "<${R_out}"
                or die
                "Could not open ${R_out} to read DADA2 R log: $OS_ERROR\n";
            my $line = <IN>;
            while ($line) {
                if (    # signs of bad termination that might be fixed
                        # with more memory
                    (      $line =~ /Error in/
                        || $line =~ /Execution halted/
                        || $line =~ /encountered errors/
                        || $line =~ /Error:/
                        || $line =~ /Traceback:/
                    )
                    && $line !~ /learnErrors/   # False positives; don't match
                    && $line !~ /error rates/
                    && $line !~ /errors/
                    )
                {
                    # Don't change $exitStatus, so DADA2 is restarted

                    $GLOBAL_PATHS->print("---R script crashed at: $line\n");

                   # Preserve the last R log file that errored. Get rid of the
                   # old R output file, then run R again.
                    File::Copy::move( "${R_out}", "${R_out}.old" );
                    $GLOBAL_PATHS->print(
                        "---See ${R_out}.old for details.\n");
                    $GLOBAL_PATHS->print("---Attempting to restart R...\n");

                    $previous_mem = $dada2mem;
                    $dada2mem     = List::Util::min( $dada2mem * 2, 160 );
                    $args =~ s/--memory=[0-9]*/--memory=$dada2mem/g;

                    last;
                }
                elsif ( $line =~ /(E|e)rror/ || $line =~ /ERROR/ ) {
                    my $error_msg
                        = "Error in DADA2. See ${R_out} for details.";
                    $GLOBAL_PATHS->print($error_msg);
                    die $error_msg;
                }
                $line = <IN>;
            }
            close IN;

            $try += 1;
        }
    }
    if ( $try > 2 ) {
        my $msg
            = "Won't attempt to restart R, as the same command has been attempted twice.";
        $GLOBAL_PATHS->print($msg);
        die $msg;
    }
    if ( $exitStatus == 1 ) {
        my $msg
            = "R output doesn't indicate error, but analysis files were not produced.";
        $GLOBAL_PATHS->print($msg);
        die $msg;
    }

    return ${R_out};
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

################################################################################
# Execute the given commands;
#  par 1: array ref containing commands
# par 2: Lexical filehandle that overrides printing to STDOUT. Pass 0 for no override.
# par 3: $dryRun parameter
# par 4: Human-friendly message to be printed to STDOUT.
# par 5: true if commands are to be qsubbed (causes return value to be array of SGE job IDs)
# Remaining args: cmds to be run through system()
sub execute_and_log {
    my %arg      = @_;
    my $cmds     = delete %arg{"cmds"};
    my @commands = @$cmds;
    my $fh       = delete %arg{"logger"}  // $GLOBAL_PATHS;
    my $dryRun   = delete %arg{"dry_run"} // 0;
    my $cuteMsg  = delete %arg{"msg"}     // 0;
    my $qsub     = delete %arg{"qsub"}    // 0;
    my $mem_GB   = delete %arg{"mem_GB"}  // "1";
    my $config   = delete %arg{"config"}  // undef;
    my $paths    = delete %arg{"paths"}   // undef;
    my $name     = delete %arg{"name"}    // "";

    if ( !@commands ) {
        warn "execute_and_log() called with no commands\n";
    }

    $fh->print("$cuteMsg\n");

    if ($qsub) {
        my $l                   = "mem_free=${mem_GB}G";
        my @allowed_exit_status = (0);
        my $exit_status         = 1;
        my $try                 = 0;
        my $success             = 0;

        while ( !grep( /^${exit_status}$/, @allowed_exit_status )
            && $try < $config->{"qsub_allowed_tries"} )
        {
            $try += 1;
            print STDERR "Try: $try\n";
            print STDERR "Max: $config->{'qsub_allowed_tries'}\n";

            my @pids;
            my $sge = Schedule::SGE->new(
                -executable => {
                    qsub =>
                        '/usr/local/packages/sge-root/bin/lx24-amd64/qsub',
                    qstat =>
                        '/usr/local/packages/sge-root/bin/lx24-amd64/qstat',
                    qacct =>
                        '/usr/local/packages/sge-root/bin/lx24-amd64/qacct',
                    qdel => '/usr/local/packages/sge-root/bin/lx24-amd64/qdel'
                },
                -l           => "$l",
                -use_cwd     => 1,
                -project     => $config->{"qsub_P"},
                -name        => $name,
                -output_file => $paths->part1_stdout_log(),
                -error_file  => $paths->part1_error_log(),
                -verbose     => $verbose,
                -w           => $config->{"qsub_w"},
                -b           => $config->{"qsub_b"}
            );

            foreach my $cmd (@commands) {
                if ($verbose) {    # print each command
                    $fh->print("\$ $cmd\n");
                }
                else {
                    STDOUT->print("\$ $cmd\n");
                }

                $sge->command($cmd);
                push @pids, $sge->execute();
            }

            my $done;
            while ( !$done ) {
                $done = 0;
                $fh->print("Checking PIDs for job status...\n") if $verbose;
                $sge->status();
                foreach (@pids) {
                    $fh->print("PID: $_\n") if $verbose;
                    my @job_stats = @{ $sge->brief_job_stats($_) };
                    $fh->print( "Job stats:\n" . Dumper(@job_stats) )
                        if $verbose;
                    if (@job_stats) {
                        last;
                    }
                    else {
                        $done = 1;
                        my $qacct = $sge->qacct($_);
                        $exit_status = $qacct->{"exit_status"};

                        if (   $qacct->{"failed"} != "0"
                            || $exit_status != "0" )
                        {
                            print STDERR
                                "Job $_ failed with exit_status $exit_status. Was try #${try} of the same step.\n";
                            if ( $exit_status == "139" ) {
                                $mem_GB *= 2;
                                $l
                                    = "$l hostname='!$qacct->{'hostname'}' mem_free=${mem_GB}G";
                            }
                        }
                        else {
                            $success = 1;
                        }
                    }
                }
                sleep 1;
            }
        }

        if ( !$success ) {
            die "Reached max retries. Quitting.\n";
        }
    }
    else {

        foreach my $cmd (@commands) {
            if ($verbose) {    # print each command
                $fh->print("\$ $cmd\n");
            }
            else {
                STDOUT->print("\$ $cmd\n");
            }
            $fh->print( system($cmd) );
        }

    }

    $fh->print( "Finished " . lcfirst($cuteMsg) . "\n" );
}

sub tagcleaned_files {
    return glob( catdir( $GLOBAL_PATHS->{'wd'}, "*R[1|2]_tc.fastq.gz" ) );
}

exit 0;
