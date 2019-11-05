#!/usr/bin/env perl

package Version;
use Exporter;
@ISA    = ('Exporter');
@EXPORT = ('version');
use File::Basename;
my $scriptsDir = dirname(__FILE__);
my $pipelineDir = dirname($scriptsDir);

sub version {

  # Retrieve and log pipeline version
  open my $versionF, "<", "${pipelineDir}/VERSION" or die "Cannot open".
    "${pipelineDir}/VERSION to read pipeline version: ${OS_ERROR}\n";
  chomp(my $version = <$versionF>);
  close $versionF;

  return($version);
}

1;
