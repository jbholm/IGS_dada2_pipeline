#!/usr/bin/env perl
#. /usr/local/packages/usepackage/share/usepackage/use.bsh
use File::Basename;
use File::Spec;
my $pipelineDir = dirname(__FILE__);

=head1 NAME

  illumina_dada2.pl

=head1 DESCRIPTION

  The script can be launched from any location on the IGS server, it automatically 
  produces a directory in /local/groupshare/ravel named after the project and run ID provided. 

  Beginning with the path to four raw Illumina sequencing files (R1, R2, R3, R4),
  a mapping file, a project ID, a run ID, and specifying the targeted variable region,
  thi script:
    1. Produces individual .fastq files for each sample listed in the mapping file
    2. Performs tag-cleaning of each file
    3. Runs the R1 (forward) and R4 (later called R2, reverse) files through the dada2 
       pipeline for either the 16S rRNA gene V3V4 or V4 regions.

=head1 SYNOPSIS

  THE FOLLOWING STEPS ARE REQUIRED BEFORE THIS SCRIPT WILL RUN:
  
  1. Qlogin to RHEL7: 
        qlogin -P jravel-lab -l mem_free=500M -q interactive.q
  2. Then enter: 
        export LD_LIBRARY_PATH=/usr/local/packages/gcc/lib64
        source /usr/local/packages/usepackage/share/usepackage/use.bsh
        use python-2.7
        use qiime
  3. then run script as below: 
  
  FOR 2-STEP
  illumina_dada2.pl -i <input directory> -p <project name> -r <run ID> -m <mapping file> -v <variable region> -sd <storage directory>

  OR

  illumina_dada2.pl -r1 <full path to raw R1 file> -r2 <full path to raw R2 file>
                    -r3 <full path to raw R2 file> -r4 <full path to raw R2 file>
                    -p <project name> -r <run ID> -m <mapping file> -v <variable region>
                    -sd <storage directory>
  
  or qsub:
  
  qsub -cwd -b y -l mem_free=1G -P jravel-lab -q threaded.q -pe thread 4 -V -e <path_to_logs> -o <path_to_logs> /home/jholm/bin/illumina_dada2_1_step_optional.pl -i <path to raw files> -p <project name> -r <run id> -v <variable region> -m <full path to mapping file> -sd <storage directory>

  FOR 1-STEP
  illumina_dada2.pl -i <input directory> -p <project name> -r <run ID> -m <mapping file> -v <variable region> -sd <storage directory> --1Step

  OR

  illumina_dada2.pl -r1 <full path to raw R1 file> -r2 <full path to raw R2 file>
                    -r3 <full path to raw R2 file> -r4 <full path to raw R2 file>
                    -p <project name> -r <run ID> -m <mapping file> -v <variable region>
                    -sd <storage directory> --1Step
  or qsub:
  
  qsub -cwd -b y -l mem_free=1G -P jravel-lab -q threaded.q -pe thread 4 -V -e <path_to_logs>
  -o <path_to_logs> /home/jholm/bin/illumina_dada2.pl -i <path to raw files> -p <project name> 
  -r <run id> -v <variable region> -m <full path to mapping file> -sd <storage directory> --1Step

=head1 OPTIONS

=over

=item B<--input-Directory, -i>
  Single full path to directory containing raw R1, R2, R3, and R4 files

=item B<--raw r1, -r1>
  Full path to raw R1 read file (forward read file, or r1) (.gz).

=item B<--raw r2, -r2>
  Full path to raw R2 read file (barcode file, or i1) (.gz).

=item B<--raw r3, -r3>
  Full path to raw R3 read file (barcode file, or i2) (.gz).

=item B<--raw r4, -r4>
  Full path to raw R4 read file (reverse read file, r2 or r4)(.gz).

=item B<--Project-ID, -p>
  The project ID.

=item B<--Sequencing-Run-ID, -r>
  A specific run ID. 

=item B<--Mapping-File, -m>
  The full path to the Qiime-formatted mapping file.

=item B<--input-Directory, -v>
  The targeted variable region current options:
  V3V4 or V4 or ITS

=item B<--1Step, --1Step>
  Use this flag if the data are prepared by 1-Step PCR 
  (only r1 & r2 raw files available)

=item B<--storage-dir, -sd>
  Indicate if you would like the pipeline to run within
  scratch ( /local/scratch/)
    OR
  groupshare (/local/groupshare/ravel)

=item B<-h|--help>
  Print help message and exit successfully.

=item B<--qsub-project, -qp>
  Indicate which qsub-project space should be used for 
  all qsubmissions

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

$OUTPUT_AUTOFLUSH = 1;

####################################################################
##                             OPTIONS
####################################################################

GetOptions(
  "raw-path|i=s"      => \my $inDir,
  "r1-path|r1=s"      => \my $r1file,
  "r2-path|r2=s"      => \my $r2file,
  "r3-path|r3=s"      => \my $r3file,
  "r4-path|r4=s"      => \my $r4file,
  "project-name|p=s"  => \my $project,
  "run-ID|r=s"        => \my $run,
  "map|m=s"           => \my $map,
  "var-reg|v=s"       => \my $var,
  "help|h!"           => \my $help,
  "debug"             => \my $debug,
  "dry-run"           => \my $dryRun,
  "skip-err-thld"     => \my $skipErrThldStr,
  "dada2-truncLen-f|for=i"    => \my $f,
  "dada2-truncLen-r|rev=i"    => \my $r,
  "dada2-maxN"        => \my $maxN,
  "dada2-maxEE"       => \my $maxEE,
  "dada2-truncQ"      => \my $truncQ,
  "dada2-rmPhix"      => \my $phix,
  "dada2-maxLen"      => \my $maxLen, 
  "dada2-minLen"      => \my $minLen,
  "dada2-minQ"        => \my $minQ,
  "1Step"             => \my $oneStep,
  "storage-dir|sd=s"  => \my $sd,
  "qsub-project|qp=s" => \my $qp,
  )
##add option for final resting place of important data

  or pod2usage(verbose => 0,exitstatus => 1);

if ($help)
{
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if (!$inDir && !$r1file && !$r2file && !$r3file && !$r4file)
{
  print "\n\tPlease provide the location of the raw sequencing files (single directory => -i)\n\t\tOR \n\tFull paths to each R1, R2, R3, and R4 file => -r1, -r2, -r3, -r4)\n\n";
  pod2usage(verbose => 2,exitstatus => 0);
  exit 1;
}

if (!$map)
{
  print "\n\tPlease provide a full path to the project mapping file (-m)\n\n";
  exit 1;
}
if (!$var)
{
  print "\n\tPlease indicate the targeted variable region (-v V3V4 or -v V4 or -v ITS)\n\n";
  exit 1;
}
my $models;
if ($var eq 'V3V4')
{
  $models = "/local/projects-t2/jholm/PECAN/v1.0/V3V4/merged_models/";
}
if ($var eq 'V4' || $var eq 'ITS')
{
  $models = "PECAN not used.\n";
}
if (!$sd)
{
  print "\n***Please choose a storage directory (-sd), either 'scratch' or 'groupshare'.";
  print "\nOR provide a full path to an existing directory.";
  pod2usage(verbose => 2,exitstatus => 0);
  exit;
} 

##"dada2-truncLen-f|f=s"    => \my $truncLenF,
##"dada2-truncLen-r|r=s"    => \my $truncLenR,
my $truncLen;

if($f && $r)
{
	$truncLen = "c($f,$r)";
}

if ($f && !$r)
{
	print "***\nPlease provide truncation lengths for forward and reverse reads\n";
  	exit;
}


my $qproj;

if ($qp)
{
  $qproj = $qp;
}
else
{
  print "\nqsub-project ID (--qp) not provided. Using jravel-lab as default\n";
  $qproj = "jravel-lab";
}

#local $ENV{use} = "/usr/local/packages/usepackage/share/usepackage/use.bsh";
#local $ENV{LD_LIBRARY_PATH} = "/usr/local/packages/gcc/lib64";
my $R = "/usr/local/packages/r-3.4.0/bin/R";


####################################################################
##                               MAIN
####################################################################
my $pd;
my $wd;

if ($sd eq "scratch")
{
  $pd = "/local/scratch/$project";
  $wd = "/local/scratch/$project/$run";
}
elsif ($sd eq "groupshare")
{
  $pd = "/local/groupshare/ravel/$project";
  $wd = "/local/groupshare/ravel/$project/$run";
}
else
{
  $pd = "$sd/$project";
  $wd = "$sd/$project/$run";
}

my $error_log = "$wd/qsub_error_logs";
my $stdout_log = "$wd/qsub_stdout_logs";
## Instead of working directory, give option for provided directory (or current working directory)

###### BEGIN MAKE WORKING DIRECTORIES ##########
################################################

if (!-e $pd)
{
  mkdir $pd;
  print "$pd did not exist -> making $pd\n";
}
if (!-e $wd)
{
  mkdir $wd;
  print "$wd did not exist -> making $wd\n";
}
if (!-e $error_log)
{
  mkdir $error_log;
  print "$error_log did not exist -> making $error_log\n";
}
if (!-e $stdout_log)
{
  mkdir $stdout_log;
  print "$stdout_log did not exist -> making $stdout_log\n";
} 

my $time = strftime("%Y-%m-%d %H:%M:%S", localtime(time));
my $log = "$wd/$project"."_".$run."_16S_pipeline_log.txt";
truncate $log, 0;

my $perlScript = File::Spec->catfile($pipelineDir, "scripts", "log_version.pl");
system($^X, $perlScript, $log);

open my $logFh, ">>$log" or die "Cannot open $log for writing: $OS_ERROR";
print $logFh "$time\n";

###### BEGIN CHECK OF QIIME CONFIGURATION ###########
#####################################################
my $qiime = "$wd/$project"."_".$run."_"."qiime_config.txt";
my $cmd = "print_qiime_config.py > $qiime";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

###### BEGIN VALIDATION OF MAPPING FILE ###########
################################################
my @errors;
#my $cmd;
@errors = glob("$error_log/*.log");
if (@errors)
{
  foreach my $error (@errors)
  {
    $cmd = "rm $error";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun; 
  }
}

print "--Validating $map\n";
$cmd = "qsub -b y -l mem_free=1G -P $qproj -V -e $error_log -o $stdout_log validate_mapping_file.py -m $map -s -o $error_log";
print "\tcmd=$cmd\n" if $dryRun || $debug;
system($cmd) == 0 or die "system($cmd) failed:$?\n" if !$dryRun;

RETURN0:
my $mappingError = glob("$error_log/*.log");
if ($mappingError)
{
  open MAPERROR, "<$mappingError" or die "Cannot open $mappingError for reading: $OS_ERROR";
  while ($_ = <MAPERROR>)
  {
    #print $line;
    chomp;
    if ($_ =~ /No errors or warnings found in mapping file./)
    {
      print $logFh "---Map passed validation.\n";
    }
    else
    {
      print $logFh "***Error in mapping file. See $mappingError for details. Exiting.\n";
      exit;
    }
  }
}
else
{
  goto RETURN0;
}
close MAPERROR;

open MAP, "<$map" or die "Cannot open $map for reading: $OS_ERROR";
my $nSamples = 0;
my $extctrl = 0;
my $pcrpos = 0;
my $pcrneg = 0;
my $projSamples = 0;
my $linecount = 0;
my $null = 0;
while (<MAP>)
{
  chomp;
  $linecount++;
  if ($linecount > 1) ## don't count header as sample
  {
    if ($_ =~ "EXTNTC" || $_ =~ "NTC.EXT" || $_ =~ "NTCEXT" || $_ =~ "EXT.NTC" || $_ =~ "NTC" || $_ =~ "EXTNEG")
    {
      $extctrl++;
    }
    elsif ($_ =~ "PCRPOS" || $_ =~ "PCR.pos" || $_ =~ "pCR.pos" || $_ =~ "POSCTRL" || $_ =~ "POS.CTRL" || $_ =~ "POSCON")
    {
      $pcrpos++;
    }
    elsif ($_ =~ "PCRNTC" || $_ =~ "PCR.NEG" || $_ =~ "PCR.NTC"|| $_ =~ "PCRNEG" || $_ =~ "PCRNEGCTRL")
    {
      $pcrneg++;
    }
    elsif ($_ =~ /NULL/ )
    {
      $null++;
    }
    else
    {
      $projSamples++;
    }
  }
}
close MAP;


$nSamples = $projSamples + $extctrl + $pcrpos + $pcrneg + $null;



###### BEGIN EVALUATION OF SAMPLES VIA MAPPING FILE ###########
###############################################################
print $logFh "PROJECT: $project\nVARIABLE REGION: $var\nNO. SAMPLES: $nSamples\nNO. NULLS: $null\nNO. EXTRACTION NEGATIVE CONTROLS: $extctrl\nNO. PCR POSITIVE CONTROLS: $pcrpos\nNO. PCR NEGATIVE CONTROLS: $pcrneg\nR VERSION: $R\nPECAN MODELS: $models\nQIIME CONFIGURATION DETAILS: \nsee $qiime\nMAPPING FILE: $map";
if ($oneStep)
{
  print $logFh "\nPCR PREPARATION METHOD: 1-Step\n\n";
}
else
{
  print $logFh "\nPCR PREPARATION METHOD: 2-Step\n\n";
}

my $r1 = "$wd/$project"."_"."$run"."_"."R1.fastq"; ## change to full path (flag)
my $r2 = "$wd/$project"."_"."$run"."_"."R2.fastq"; ## change to full path (flag)
my $r3 = "$wd/$project"."_"."$run"."_"."R3.fastq"; ## change to full path - usr provides full path (flag)

my $barcodes = "$wd/barcodes.fastq"; ## change to full path to full barcodes (flag)
my $count = 0;

print "--Checking for existence of $barcodes\n";


my $step1;
my $step2;
my $step3;

###### BEGIN 1STEP BARCODES ##########
#######################################

if ($oneStep)
{
  my $r1file = glob("$inDir/*R1.fastq*");
  my $r2file = glob("$inDir/*R2.fastq*");
  if (!-e $barcodes )
  {
    print "---Barcode files not found or were empty\n";
    if ($inDir && $r1file != /.gz/)
    {
      $r1 = $r1file;
      $r2 = $r2file;
    }
    else
    {
      print "---Copying barcode and index files to $wd\n";
      $cmd = "zcat $r1file > $r1 | zcat $r2file > $r2 ";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
      print $logFh "$project barcode and index files copied from $r2file and $r3file to $r1 and $r2\n";
    }
    my $start = time;
    print "---Extracting barcodes and index files\n";
    ## add in (possible replace) awk version of concatenation 
    $step1 = "extract_barcodes.py";
    @errors = glob("$error_log/$step1.e*");
    if (@errors)
    {
      foreach my $error (@errors)
      {
        $cmd = "rm $error";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun; 
      }
    }
    $cmd = "qsub -b y -l mem_free=1G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log extract_barcodes.py -f $r1 -r $r2 -c barcode_paired_end --bc1_len 12 --bc2_len 12 -m $map -o $wd";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    
    check_error_log($error_log, $step1);

    print "---Waiting for barcode extraction to complete.\n";

    RETURN1:
    if(!(-e $barcodes))
    {
      goto RETURN1;
    }
    else
    {
      ## try adding else statement, if needed
      my $duration = time - $start;
      print "---Barcode extraction complete\n";
      print $logFh "Duration of barcode extraction: $duration s\n";
      print "---Duration of barcode extraction: $duration s\n";
    }
  }
  else
  {
    open BAR, "<$barcodes" or die "Cannot open $barcodes for reading: $OS_ERROR";
    while (<BAR>)
    {
      $count++;
    }
    close BAR;
    if ($count = $nSamples)
    {
      print "---$barcodes already exists and contains $count entries, as expected\n";
      print $logFh "-> $barcodes already exists and contains $count entries, as expected\n";
    }
    else
    {
      print "---$barcodes already exists, but contains $count entries, while there are $nSamples samples\n";
      print $logFh "-> $barcodes already exists, but contains $count entries, while there are $nSamples samples\n";
      exit 0;
    }
  }
}

###### BEGIN 2STEP BARCODES ##########
#######################################
else
{
  ## user specified path to barcodes.fastq --- if present, use, if not, then check for existence in cwd, if not, produce
  if (!-e $barcodes )
  {
    print "---Barcode files not found or were empty\n";
    if (!-e $r2 || !-e $r3 )
    {
      my $r2 = "$wd/$project"."_"."$run"."_"."R2.fastq"; ## change to full path (flag)
      my $r3 = "$wd/$project"."_"."$run"."_"."R3.fastq";

      if ($inDir)
      {
        print "---Copying barcode and index files to $wd\n";
        $cmd = "zcat $inDir/*R2.fastq.gz > $r2 | zcat $inDir/*R3.fastq.gz > $r3 ";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
        print $logFh "$project barcode and index files copied from $inDir to $r2 and $r3\n";
      }
      else
      {
        print "---Copying barcode and index files to $wd\n";
        $cmd = "zcat $r2file > $r2 | zcat $r3file > $r3 ";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
      }
    }

    my $start = time;
    print "---Extracting barcodes and index files\n";
    ## add in (possible replace) awk version of concatenation 
    $step1 = "extract_barcodes.py";
    @errors = glob("$error_log/$step1.e*");
    if (@errors)
    {
      foreach my $error (@errors)
      {
        $cmd = "rm $error";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun; 
      }
    }
    $cmd = "qsub -b y -l mem_free=1G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log extract_barcodes.py --input_type barcode_paired_end -f $r2 -r $r3 --bc1_len 8 --bc2_len 8 -o $wd";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    
    check_error_log($error_log, $step1);

    print "---Waiting for barcode extraction to complete.\n";

    RETURN1:
    if(!(-e $barcodes))
    {
    	goto RETURN1;
    }
    else
    {
      ## try adding else statement, if needed
      my $duration = time - $start;
      print "---Barcode extraction complete\n";
      print $logFh "---Duration of barcode extraction: $duration s\n";
      print "---Duration of barcode extraction: $duration s\n";
    }
  }
  else
  {
    open BAR, "<$barcodes" or die "Cannot open $barcodes for reading: $OS_ERROR";
    while (<BAR>)
    {
      $count++;
    }
    close BAR;
    if ($count = $nSamples)
    {
      print "---$barcodes already exists and contains $count entries, as expected\n";
      print $logFh "-> $barcodes already exists and contains $count entries, as expected\n";
    }
    else
    {
      print "---$barcodes already exists, but contains $count entries, while there are $nSamples samples. Exiting.\n";
      print $logFh "-> $barcodes already exists, but contains $count entries, while there are $nSamples samples. Exiting.\n";
      exit 0;
    }
  }

  my $reads1 = "$wd/reads1.fastq";
  my $reads2 = "$wd/reads2.fastq";
  $cmd = "rm -rf $reads1";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  $cmd = "rm -rf $reads2";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
}
###### BEGIN SPLIT LIBRARIES ##########
#######################################

my $r1split;
my $r4split;
my $r4;

if ($oneStep)
{
  $r1split = "$wd/R1split";
  $r4split = "$wd/R2split";
  $r4 = "$wd/$project"."_"."$run"."_"."R2.fastq";
}
else
{
  $r1split = "$wd/R1split";
  $r4split = "$wd/R4split";
  $r4 = "$wd/$project"."_"."$run"."_"."R4.fastq";
}

my $r1fq = "$r1split/seqs.fastq";
my $r4fq = "$r4split/seqs.fastq";

## think of way to ensure the consistent read order in the r1 and r4 files.
## print headers of r1 and r4, comm r1 r4 - to ensure the seqIDs are the same order.


print "--Checking for existence of $r1fq and $r4fq\n";
if (!-e $r1fq || !-e $r4fq)
{
  if ($oneStep)
  {
    print "---Obtaining $project-specific R1 and R2 fastq files\n";
  }
  else
  {
    print "---Obtaining $project-specific R1 and R4 fastq files\n";
    if ($inDir)
    {
      $cmd = "zcat $inDir/*R1.fastq.gz > $r1 | zcat $inDir/*R4.fastq.gz > $r4 ";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
      print $logFh "---$project barcode and index files copied from $inDir to $r1 and $r4\n";
    }
    else
    {
      $cmd = "zcat $r1file > $r1 | zcat $r4file > $r4 ";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
      print $logFh "---$project barcode and index files copied from $r1file and $r4file to $r1 and $r4\n";
    }
  }
    
  print "---Producing $r1fq and $r4fq\n";
  my $start = time;
  $step2 = "split_libraries_fastq.py";
  @errors = glob("$error_log/$step2.e*");
  if (@errors)
  {
    foreach my $error (@errors)
    {
      $cmd = "rm $error";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun; 
    }
  }
  ## including the stitch script before so that demultiplex happens by barcode NOT order
  if ($oneStep)
  {
    $cmd = "qsub -cwd -b y -l mem_free=1G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log split_libraries_fastq.py -i $r1 -o $r1split -b $barcodes -m $map --max_barcode_errors 1 --store_demultiplexed_fastq --barcode_type 24 -r 999 -n 999 -q 0 -p 0.0001";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    print $logFh "Demultiplexing command F: \n\t$cmd\n\n";
    $cmd = "qsub -cwd -b y -l mem_free=1G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log split_libraries_fastq.py -i $r4 -o $r4split -b $barcodes -m $map --max_barcode_errors 1 --store_demultiplexed_fastq --barcode_type 24 -r 999 -n 999 -q 0 -p 0.0001";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    print $logFh "Demultiplexing command R: \n\t$cmd\n\n";
  }
  else
  {
    $cmd = "qsub -cwd -b y -l mem_free=1G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log split_libraries_fastq.py -i $r1 -o $r1split -b $barcodes -m $map --max_barcode_errors 1 --store_demultiplexed_fastq --barcode_type 16 -r 999 -n 999 -q 0 -p 0.0001";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    print $logFh "Demultiplexing command F: \n\t$cmd\n\n";
    $cmd = "qsub -cwd -b y -l mem_free=1G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log split_libraries_fastq.py -i $r4 -o $r4split -b $barcodes -m $map --max_barcode_errors 1 --store_demultiplexed_fastq --barcode_type 16 -r 999 -n 999 -q 0 -p 0.0001";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    print $logFh "Demultiplexing command R: \n\t$cmd\n\n";
  }

  check_error_log($error_log, $step2);

  $r1fq = "$r1split/seqs.fastq";
  $r4fq = "$r4split/seqs.fastq";

  print "---Waiting for R1 seqs.fastq to complete.\n";
  RETURN2:
  if(!(-e $r1fq))
  {
  	goto RETURN2;
  }
  else
  {

  }
  my $duration = time - $start;
  print $logFh "---Duration of R1 seqs.fastq production: $duration s\n";
  print "---Duration of R1 seqs.fastq production: $duration s\n";

  if ($oneStep)
  {
    print "---Waiting for R2 seqs.fastq to complete.\n";
  }
  else
  {
    print "---Waiting for R4 seqs.fastq to complete.\n";
  }
  RETURN3:
  if(!(-e $r4fq))
  {
  	goto RETURN3;
  }
  $duration = time - $start;
  if ($oneStep)
  {
    print $logFh "---Duration of R2 seqs.fastq production: $duration s\n";
    print "---Duration of R2 seqs.fastq production: $duration s\n";
  }
  else
  {
    print $logFh "Duration of R4 seqs.fastq production: $duration s\n";
    print "---Duration of R4 seqs.fastq production: $duration s\n";
  }
}
else
{
  print "-> $r1fq and $r4fq fastq files already produced. Demultiplexing ...\n";
}

##Check split_library log for 0's
my $split_log = "$r1split/split_library_log.txt";
my @split;
my $newSamNo;
open SPLIT, "<$split_log" or die "Cannot open $split_log for writing: $OS_ERROR";
while (<SPLIT>)
{
	if ($_ =~ /\t/)
	{
		chomp;
		my ($sample, $nReads) = split, /\t/;
    	chomp $nReads;
		if ( $nReads eq "0")
		{
			#print $logFh "Adding sample: $sample\n";
	  		push @split, $sample;
		}
		else
		{
			print "$sample $nReads\n";
		}
	}
}
close SPLIT;
#pop @split; ## to remove the final line of "Total number seqs written" ---commented out on 01/21/19  because it was deleting the last element in the array of samples 
$newSamNo = $nSamples - scalar @split;
if (scalar @split ne $nSamples)
{
  print $logFh "The following ". scalar @split." samples returned 0 reads, due to either a lack of barcodes or a filtering of those reads in split_libraries_fastq.py:\n";
  foreach my $x (@split)
  {
    print $logFh "   $x\n";
  }
  print $logFh "Number of samples after split_libraries_fastq.py: $newSamNo\n";
}
elsif (scalar @split eq $nSamples)
{
  print $logFh "No samples were successfully demultiplexed. Is the mapping file correct? Exiting.\n";
  exit;
}

###### BEGIN SPLIT BY SAMPLE ##########
#######################################

if ($oneStep)
  {
    print "--Checking if $project R1 & R2 seqs.fastq files were split by sample ID\n";
  }
else
  {
    print "--Checking if $project R1 & R4 seqs.fastq files were split by sample ID\n";
  }

my $r1seqs = "$r1split/split_by_sample_out";
my $r4seqs = "$r4split/split_by_sample_out";

COUNT0:
my $n_fq1 = 0;

my @filenames = glob("$r1seqs/*.fastq");
my @r4filenames = glob("$r4seqs/*.fastq");

if (scalar(@filenames) != $newSamNo || scalar(@r4filenames) != $newSamNo)
{
  print "There are ".scalar @filenames." split files. Waiting for $newSamNo files\n";
  my $step3 = "split_sequence_file_on_sample_ids.py";
  @errors = glob("$error_log/$step3.e*");
  if (@errors)
  {
    foreach my $error (@errors)
    {
      $cmd = "rm $error";
      print "\tcmd=$cmd\n" if $dryRun || $debug;
      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun; 
    }
  }
  print "---Sample specific files not found or completed... Splitting $project seqs.fastq files by sample ID\n";
  print $logFh "There are ". scalar(@filenames)." sample specific files found (expected $newSamNo)... Splitting $project seqs.fastq files by sample ID\n";
  $cmd = "rm -rf $r1seqs; rm -rf $r4seqs";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  RETURN4:
  if(!(-e $r1fq))
  {
  	goto RETURN4;
  }
  $cmd = "qsub -cwd -b y -l mem_free=1G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log -V split_sequence_file_on_sample_ids.py -i $r1fq --file_type fastq -o $r1seqs";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  print $logFh "$cmd\n";

  RETURN5:
  if(!(-e $r4fq))
  {
  	goto RETURN5;
  }
  $cmd = "qsub -cwd -b y -l mem_free=1G -P $qproj -q threaded.q -pe thread 4 -V -e $error_log -o $stdout_log -V split_sequence_file_on_sample_ids.py -i $r4fq --file_type fastq -o $r4seqs";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  print $logFh "$cmd\n";

  check_error_log($error_log, $step3);

  ## the $nSamples needs to be altered if a sample has 0 reads, because the sample-specific fastq won't be produced

 
  COUNT0:
  my $n_fq = 0;
  my @filenames = glob("$r1seqs/*.fastq");
  {
    if (@filenames)
    {
      foreach my $file (@filenames)
      {
        $n_fq++;
      }
      if ($n_fq != $newSamNo)
      {
        goto COUNT0;
      }
      else
      {
        print "--All samples accounted for in $r1seqs\n"
      }
    }
    else 
    {
       goto COUNT0;
    }
  }
  COUNT1:
  $n_fq = 0;
  my @r4filenames = glob("$r4seqs/*.fastq");
  {
    if (@r4filenames)
    {
      foreach my $file (@r4filenames)
      {
        $n_fq++;
      }
      if ($n_fq != $newSamNo)
      {
        goto COUNT0;
      }
      else
      {
        print "--All samples accounted for in $r4seqs\n"
      }
    }
    else 
    {
       goto COUNT0;
    }
  }
}
else
{
  print "--$newSamNo sample-specific files present as expected.\n";
  print $logFh "$newSamNo sample-specific files present as expected.\n";
}

print "--Checking if target primers have been removed from $project R1 & R4 sample-specific files..\n";
print $logFh "Checking if target primers have been removed from $project R1 & R4 sample-specific files..\n";
my @r1tcfiles = glob("$wd/*R1_tc.fastq");
my @r4tcfiles = glob("$wd/*R2_tc.fastq");

###### BEGIN TAGCLEANING ##########
###################################

my $start = time;
if (scalar(@r1tcfiles) != $newSamNo || scalar(@r4tcfiles) != $newSamNo)
{
  if ($oneStep)
  {
    if ($var eq "V3V4")
    {
      print "Removing V3V4 primers from all sequences\n";
      my $filename;
      opendir R1, $r1seqs or die "Cannot open directory $r1seqs\n";
      while ( $filename = readdir R1) 
      {
        if ($filename =~ /.fastq/)
        {
          my @suffixes = (".fastq",".fq");
          my $Prefix = basename($filename, @suffixes);
          my $tc = "$wd/$Prefix"."_R2_tc";
          $cmd = "qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log/tagcleanerR2_$filename -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r1seqs/$filename -out $tc -line_width 0 -verbose -tag5 GGACTACHVGGGTWTCTAAT -mm5 2 -trim_within 50";
          print "\tcmd=$cmd\n" if $dryRun || $debug;
          system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
          ##print $logFh "$cmd\n";
        }
      }
      close R1;

      opendir R4, $r4seqs or die "Cannot open directory $r4seqs\n";
      while ( $filename = readdir R4) 
      {
        my @suffixes = (".fastq",".fq");
        my $Prefix = basename($filename, @suffixes);
        my $tc = "$wd/$Prefix"."_R1_tc";
        $cmd = "qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log/tagcleanerR1_$filename -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r4seqs/$filename -out $tc -line_width 0 -verbose -tag5 ACTCCTACGGGAGGCAGCAG -mm5 2 -trim_within 50";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
        ##print $logFh "$cmd\n";
      }
      close R4;
    }

    if ($var eq "V4")
    {
      print "--Removing V4 primers from all sequences\n";
      my $filename;
      opendir R1, $r1seqs or die "Cannot open directory $r1seqs\n";
      while ( $filename = readdir R1) 
      {
        if ($filename =~ /.fastq/)
        {
          my @suffixes = (".fastq",".fq");
          my $Prefix = basename($filename, @suffixes);
          my $tc = "$wd/$Prefix"."_R1_tc";
          $cmd = "qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r1seqs/$filename -out $tc -line_width 0 -verbose -tag5 GTGCCAGCMGCCGCGGTAA -mm5 2";
          print "\tcmd=$cmd\n" if $dryRun || $debug;
          system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
        }
      }
      close R1;

      opendir R4, $r4seqs or die "Cannot open directory $r4seqs\n";
      while ( $filename = readdir R4) 
      {
        my @suffixes = (".fastq",".fq");
        my $Prefix = basename($filename, @suffixes);
        my $tc = "$wd/$Prefix"."_R2_tc";
        $cmd = "qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r4seqs/$filename -out $tc -line_width 0 -verbose -tag5 ACTCCTACGGGAGGCAGCAG -mm5 2";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
        #print $logFh "$cmd\n";
      }
      close R4;
    }
  }
	else
  {
  	if ($var eq "V3V4")
  	{
  	  print $logFh "...Removing V3V4 primers from all sequences.\n";
  	  my $filename;
  	  opendir R1, $r1seqs or die "Cannot open directory $r1seqs\n";
  	  while ( $filename = readdir R1) 
  	  {
  	    if ($filename =~ /.fastq/)
  	    {
  	      my @suffixes = (".fastq",".fq");
  	      my $Prefix = basename($filename, @suffixes);
  	      my $tc = "$wd/$Prefix"."_R1_tc";
  	      $cmd = "qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r1seqs/$filename -out $tc -line_width 0 -verbose -tag5 ACTCCTACGGGAGGCAGCAG -mm5 2";
  	      print "\tcmd=$cmd\n" if $dryRun || $debug;
  	      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  	      ##print $logFh "$cmd\n";
  	    }
  	  }
  	  close R1;

  	  opendir R4, $r4seqs or die "Cannot open directory $r4seqs\n";
  	  while ( $filename = readdir R4) 
  	  {
  	    my @suffixes = (".fastq",".fq");
  	    my $Prefix = basename($filename, @suffixes);
  	    my $tc = "$wd/$Prefix"."_R2_tc";
  	    $cmd = "qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r4seqs/$filename -out $tc -line_width 0 -verbose -tag5 GGACTACHVGGGTWTCTAAT -mm5 2";
  	    print "\tcmd=$cmd\n" if $dryRun || $debug;
  	    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  	    ##print $logFh "$cmd\n";
  	  }
  	  close R4;
  	}
  	if ($var eq "V4")
  	{
  	  print "...Removing V4 primers from all sequences.\n";
  	  my $filename;
  	  opendir R1, $r1seqs or die "Cannot open directory $r1seqs\n";
  	  while ( $filename = readdir R1) 
  	  {
  	    if ($filename =~ /.fastq/)
  	    {
  	      my @suffixes = (".fastq",".fq");
  	      my $Prefix = basename($filename, @suffixes);
  	      my $tc = "$wd/$Prefix"."_R1_tc";
  	      $cmd = "qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r1seqs/$filename -out $tc -line_width 0 -verbose -tag5 GTGCCAGCMGCCGCGGTAA -mm5 2";
  	      print "\tcmd=$cmd\n" if $dryRun || $debug;
  	      system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  	      #print $logFh "$cmd\n";
  	    }
  	  }
  	  close R1;

  	  opendir R4, $r4seqs or die "Cannot open directory $r4seqs\n";
  	  while ( $filename = readdir R4) 
  	  {
  	    my @suffixes = (".fastq",".fq");
  	    my $Prefix = basename($filename, @suffixes);
  	    my $tc = "$wd/$Prefix"."_R2_tc";
  	    $cmd = "qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r4seqs/$filename -out $tc -line_width 0 -verbose -tag5 ACTCCTACGGGAGGCAGCAG -mm5 2";
  	    print "\tcmd=$cmd\n" if $dryRun || $debug;
  	    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  	    #print $logFh "$cmd\n";
  	  }
  	  close R4;
  	}
    if ($var eq "ITS")
    {
      print "...Removing ITS primers from all sequences.\n";
      my $filename;
      opendir R1, $r1seqs or die "Cannot open directory $r1seqs\n";
      while ( $filename = readdir R1) 
      {
        if ($filename =~ /.fastq/)
        {
          my @suffixes = (".fastq",".fq");
          my $Prefix = basename($filename, @suffixes);
          my $tc = "$wd/$Prefix"."_R1_tc";
          $cmd = "qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r1seqs/$filename -out $tc -line_width 0 -verbose -tag5 CTGCCCTTTGTACACACCGC -mm5 2";
          print "\tcmd=$cmd\n" if $dryRun || $debug;
          system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
          #print $logFh "$cmd\n";

          # $cmd = "qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $tc -out $tc -line_width 0 -verbose -tag5 CTGCCCTTTGTACACACCGC -mm5 2";
          # print "\tcmd=$cmd\n" if $dryRun || $debug;
          # system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
          # print $logFh "$cmd\n";

          # $cmd = "qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $tc -out $tc -line_width 0 -verbose -tag5 CTGCCCTTTGTACACACCGC -mm5 2";
          # print "\tcmd=$cmd\n" if $dryRun || $debug;
          # system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
          # print $logFh "$cmd\n";
        }
      }
      close R1;

      opendir R4, $r4seqs or die "Cannot open directory $r4seqs\n";
      while ( $filename = readdir R4) 
      {
        my @suffixes = (".fastq",".fq");
        my $Prefix = basename($filename, @suffixes);
        my $tc = "$wd/$Prefix"."_R2_tc";
        $cmd = "qsub -cwd -b y -l mem_free=200M -P $qproj -V -e $error_log -o $stdout_log perl /usr/local/packages/tagcleaner-0.16/bin/tagcleaner.pl -fastq $r4seqs/$filename -out $tc -line_width 0 -verbose -tag5 TTTCGCTGCGTTCTTCATCG -mm5 2";
        print "\tcmd=$cmd\n" if $dryRun || $debug;
        system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
        #print $logFh "$cmd\n";
      }
      close R4;
    }
  }

  COUNT2:
  my $n_tc = 0;
  @r1tcfiles = glob("$wd/*R1_tc.fastq");
  {
    if (@r1tcfiles)
    {
      foreach my $file (@r1tcfiles)
      {
        $n_tc++;
      }
      if ($n_tc != $newSamNo)
      {
        goto COUNT2;
      }
      else
      {
        print "---All tagcleaned R1 samples accounted for in $wd\n"
      }
    }
    else 
    {
       goto COUNT2;
    }
  }
  COUNT3:
  my $n_tcr4 = 0;
  @r4tcfiles = glob("$wd/*R2_tc.fastq");
  {
    if (@r4tcfiles)
    {
      foreach my $file (@r4tcfiles)
      {
        $n_tcr4++;
      }
      if ($n_tcr4 != $newSamNo)
      {
        goto COUNT3;
      }
      else
      {
        print "---All tagcleaned R4 (R2) samples accounted for in $wd\n"
      }
    }
    else 
    {
       goto COUNT3;
    }
  }
  my $duration = time - $start;
  print $logFh "...Primer sequences removed from $newSamNo samples. Beginning DADA2.\n";
  print $logFh "--Duration of tagcleaning: $duration s\n";
}
else
{
  print "--$newSamNo sample-specific, tag-cleaned files present as expected.\n";
  print $logFh "...$newSamNo sample-specific, tag-cleaned files present as expected. Beginning DADA2.\n";
}

###### BEGIN DADA2 ##########
#############################

my $dada2 = "$wd/dada2_abundance_table.rds";
if (!-e $dada2)
{
  print "--Removing old filtered fastq files from previous runs\n";
  $cmd = "rm -rf $wd/filtered";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;

  print "Running DADA2 with fastq files in $wd\n";
  print $logFh "Running DADA2 for $var region";
  chdir $pd;  
  if ($oneStep)
  {
    if ($var eq "V3V4")
    {
      if (!$truncLen)
      {
        $truncLen = "c(255,255)";
      }
      if (!$maxN)
      {
        $maxN = 0;
      }
      if (!$maxEE)
      {
        $maxEE = "c(2,2)";
      }
       if (!$truncQ)
      {
        $truncQ = 2;
      }
      if (!$phix)
      {
        $phix = "TRUE";
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
      dada2($run, $truncLen, $maxN, $maxEE, $truncQ, $phix, $maxLen, $minLen, $minQ);
    }

    if ($var eq "V4")
    { 
      if (!$truncLen)
      {
        $truncLen = 200;
      }
      if (!$maxN)
      {
        $maxN = 0;
      }
      if (!$maxEE)
      {
        $maxEE = "c(2,2)";
      }
       if (!$truncQ)
      {
        $truncQ = 2;
      }
      if (!$phix)
      {
        $phix = "TRUE";
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
      dada2($run, $truncLen, $maxN, $maxEE, $truncQ, $phix, $maxLen, $minLen, $minQ);
    }
  }
  else
  {
    if ($var eq "V3V4")
    {
      if (!$truncLen)
      {
        $truncLen = "c(255,225)";
      }
      if (!$maxN)
      {
        $maxN = 0;
      }
      if (!$maxEE)
      {
        $maxEE = "c(2,2)";
      }
       if (!$truncQ)
      {
        $truncQ = 2;
      }
      if (!$phix)
      {
        $phix = "TRUE";
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
      dada2($run, $truncLen, $maxN, $maxEE, $truncQ, $phix, $maxLen, $minLen, $minQ);
    }

    if ($var eq "V4")
    { 
      if (!$truncLen)
      {
        $truncLen = 200;
      }
      if (!$maxN)
      {
        $maxN = 0;
      }
      if (!$maxEE)
      {
        $maxEE = "c(2,2)";
      }
       if (!$truncQ)
      {
        $truncQ = 2;
      }
      if (!$phix)
      {
        $phix = "TRUE";
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
      dada2($run, $truncLen, $maxN, $maxEE, $truncQ, $phix, $maxLen, $minLen, $minQ);
    }

    if ($var eq "ITS")
    { 
      if (!$truncLen)
      {
        $truncLen = 0;
      }
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
        $phix = "TRUE";
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
      dada2($run, $truncLen, $maxN, $maxEE, $truncQ, $phix, $maxLen, $minLen, $minQ);
    }
  }
}

my $rt = "$wd/dada2_part1_rTmp.R";
my $projrt = "$wd/$project"."_".$run."_dada2_part1_rTmp.R";
if (!-e $projrt)
{
  $cmd = "mv $rt $projrt";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  print $logFh "\tDADA2-specific commands can be found in $projrt\n";
}

my $rtout = "$wd/dada2_part1_rTmp.Rout";
my $projrtout = "$wd/$project"."_".$run."_dada2_part1_rTmp.Rout";
if (!-e $projrtout)
{
  $cmd = "mv $rtout $projrtout";
  print "\tcmd=$cmd\n" if $dryRun || $debug;
  system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  print $logFh "\tDADA2-specific commands with output can be found in $projrtout\n";
}

###### EVALUATING DADA2 OUTPUT ##########
#########################################
my @removed;
## do a scan of this file for those samples not passing filtering (0 sequences surviving - and log this)
open RTOUT, "<$projrtout" or die "Cannot open $rtout for reading: $OS_ERROR";
while (<RTOUT>)
{
  chomp;
  if ($_ =~ /filterAndTrim/)
  {
    #print "Splitting $_\n";
    my @filtTrim = split(/\s/, $_);
    for my $x (@filtTrim)
    {
      #print "Searching $x\n";
      if ($x =~ /truncLen/)
      {
        chomp;
        $truncLen = $x;
      }
      if ($x =~ /maxN/)
      {
        chomp;
        $maxN = $x;
      }
      if ($x =~ /maxEE/)
      {
        chomp;
        $maxEE = $x;
      }
      if ($x =~ /truncQ/)
      {
        chomp;
        $truncQ = $x;
      }
      if ($x =~ /phix/)
      {
        chomp;
        $phix = $x;
      }
    }
  }
  if ($_ =~ /The filter removed all reads/)
  {
  	push @removed, $_;
  }
}
close RTOUT;

if (scalar @removed > 0)
{
	print $logFh scalar @removed . " samples removed during dada2 filtering:\n";
	for my $x (@removed)
	{
		print $logFh "$x\n";
	}
}
print $logFh "\n";

my $dadaTbl = "$wd/dada2_part1_stats.txt";
if (-e $dadaTbl)
{
  print $logFh "\nFor $var region, dada2 used the following filtering requirements:\n$truncLen\n$maxN\n$maxEE\n$truncQ\n$phix\n";
  print $logFh "dada2 completed successfully!\nAbundance table for $project run $run located at $wd/dada2_abundance_table.rds\n";
  print $logFh "See $dadaTbl for dada2 table of reads surviving by step\n";
  close $logFh;
}
else
{
  print "---dada2 did not complete successfully, something went wrong!\n---Check $projrtout.\n";
}

###### COMPLETING $logFh FILE ##############
#########################################

open $logFh, "<$log" or die "Cannot open $log for reading: $OS_ERROR";
while (<$logFh>)
{
  if ($_=~ /dada2 completed successfully!/)
  {
    print "---dada2 completed successfully\n";
    print "---Run-specific abundance table written to $wd/dada2_abundance_table.rds\n";
    print "---See $log for processing details\n";
    print "---Removing original R1, R2, R3, and R4 files from $wd\n";
    $cmd = "rm -rf $r1";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    $cmd = "rm -rf $r2";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    $cmd = "rm -rf $r3";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
    $cmd = "rm -rf $r4";
    print "\tcmd=$cmd\n" if $dryRun || $debug;
    system($cmd) == 0 or die "system($cmd) failed with exit code: $?" if !$dryRun;
  }
}
close $logFh;

## moving final files to directory created within /local/projects/16S_DATA/projects/
## user specified directory for final storage of files. 

## if reference based chim rm requested { perform vsearch filtered.fastq}

####################################################################
##                               SUBS
####################################################################


sub dada2
{
  my $run = shift;
  my $truncLen = shift;
  my $maxN = shift;
  my $maxEE = shift;
  my $truncQ = shift;
  my $phix = shift;
  my $maxLen = shift;
  my $minLen = shift;
  my $minQ = shift;

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
  run_R_script( $Rscript );
}

sub run_R_script {

  chdir $wd;
  my $Rscript = shift;

  my $outFile = "dada2_part1_rTmp.R";
  open OUT, ">$outFile",  or die "cannot write to $outFile: $!\n";
  print OUT "$Rscript";
  close OUT;

  my $cmd = "$R CMD BATCH $wd/$outFile";
  system($cmd) == 0 or die "system($cmd) failed:$?\n";

  my $outR = $outFile . "out";
  open IN, "$outR" or die "Cannot open $outR for reading: $OS_ERROR\n";
  my $exitStatus = 1;

  foreach my $line (<IN>)
  {
    if ($line eq /Error/ && $line ne /learnErrors/ && $line ne /error rates/ && $line ne /errors/)
    {
      print "R script crashed at\n$line";
      print "check $outR for details\n";
      $exitStatus = 0;
      exit;
    }
  }
  close IN;
}

sub source 
{
    my $name = shift;

    open my $fh, "<$name" or die "could not open $name: $!";

    while (<$fh>) {
        chomp;
        my ($k, $v) = split /=/, $_, 2;
        $v =~ s/^(['"])(.*)\1/$2/; #' fix highlighter
        $v =~ s/\$([a-zA-Z]\w*)/$ENV{$1}/g;
        $v =~ s/`(.*?)`/`$1`/ge; #dangerous
        $ENV{$k} = $v;
    }
}

sub check_error_log
{
  my $error_log = shift;
  my $step = shift;
  print "---Checking for $step errors\n";
  my $error;
  ## make $error array so that goes through all errors for that step.
  RECHECK:
  $error = glob("$error_log/$step.e*");
  if (defined $error)
  {
    print "---Checking for errors in $error.\n";
    open ERROR, "<$error" or die "Can't open $error.\n";
    while (<ERROR>)
    {
      if ($_ =~ /error/ || $_ =~ /Error/ || $_ =~ /ERROR/)
      {
        print "Error in $step. See $error.\n";
        print $logFh "Error in $step. See $error.\n";
        exit;
      }
      else
      {
        print "---No errors found in $step.\n";
        print $logFh "---No errors found in $step.\n";
      }
    }
  }
  else
  {
    ##print "---Re-checking for $step error file\n";
    goto RECHECK;
  }

  close ERROR;
}

exit 0;
