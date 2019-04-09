#!/usr/bin/env perl
use File::Basename;
my $scriptsDir = dirname(__FILE__);
my $pipelineDir = dirname($scriptsDir);

# Retrieve and log pipeline version
open my $versionF, "<", "${pipelineDir}/VERSION" or die "Cannot open".
  "${pipelineDir}/VERSION to read pipeline version: ${OS_ERROR}\n";
chomp(my $version = <$versionF>);
close $versionF or warn "Closing ${pipelineDir}/VERSION failed\n";

print "Logging to: ${ARGV[0]}\n";
open (my $fh, ">>", "${ARGV[0]}") or die "Cannot open $ARGV[0] for append";
print $fh "PIPELINE VERSION: ${version}\n";
close $fh;
