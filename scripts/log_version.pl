#!/usr/bin/env perl
use File::Basename;
my $scriptsDir = dirname(__FILE__);
my $pipelineDir = dirname($scriptsDir);
my $debug = $ARGV[1];

# Retrieve and log pipeline version
open my $versionF, "<", "${pipelineDir}/VERSION" or die "Cannot open".
  "${pipelineDir}/VERSION to read pipeline version: ${OS_ERROR}\n";
chomp(my $version = <$versionF>);
if (!$debug) { #Remove pre-release number if debug mode
	my @versionSplit = split(/-/, $version);
	$version = $versionSplit[0];
}
close $versionF or warn "Closing ${pipelineDir}/VERSION failed\n";

open (my $fh, ">>", "${ARGV[0]}") or die "Cannot open $ARGV[0] for append";
print $fh "PIPELINE VERSION: ${version}\n";
close $fh;
