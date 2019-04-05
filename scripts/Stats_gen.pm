#!/usr/bin/perl

package Stats_gen;
use Exporter;
@ISA = ('Exporter');
@EXPORT = ('combine_dada2_stats' );
use strict;
use warnings;
use File::Spec::Functions;
use File::Spec::Functions qw(splitdir);

##############################################################################################################################################################
######################  SCRIPT TO GENERATE A CONCATENATED DADA2 STATS FILE WHICH INCLUDES BOTH PART1 & PART2 STATS INFO  ###################################################################################################################################################################################################

########################Provide the path to the parent directory along with the script
#####Usage: /home/jgeorge/scripts/16s/stats_gen.pl /local/scratch/jgeorge/HMP

sub combine_dada2_stats {
my $path = $_[0] or die "Please specify which directory to search";
$path =~ s/\/$//;    # Remove possible trailing slash

########get an array of path to the folders for each run within the project
my @list_of_dirs = grep { -d } glob "$path/*";

  ####to get the project name from the path specified
my @words = splitdir( $path );
my $project_name = $words[ ( scalar(@words) - 1 ) ]; 

print "$project_name is the name of the project\n";

  ####path to the concatenated dada2 part1 stats file for every runs within the project#########################################
my $out_file = catfile(@words, "${project_name}_dada2_part1_concatenated_stats.txt");
print "$out_file is the concatenated part 1 stats file\n";
open( OUTFILE, '>', $out_file ) or die "Could not open file $!";
print OUTFILE "sample_id\tinput\tfiltered\tmerged\n";

########Go inside each run folders and concatenate the dada2+part1_stats file of each run to out_file
foreach (@list_of_dirs) {
    opendir( SUBDIR, $_ ) || die "Could not open file $!";
    my @folders = readdir(SUBDIR);
    closedir(SUBDIR);
    foreach my $file (@folders) {
        if ( $file =~ m/stats.txt/ ) {
            my $file_name = "$_/" . "dada2_part1_stats.txt";

            open( INFILE, '<', $file_name ) or die "Could not open file $!";
            while (my $read = <INFILE> ) {
                chomp $read;
                if ( $read !~ m/input/ ) {
                    print OUTFILE "$read\n";
                }
            }
        }
    }
}
close OUTFILE;
close INFILE;

###########Read the contents of dada2_part2 stats file for the project to an array

my $filename2 = catfile(($path), "dada2_part2_stats.txt");
open( my $inputFH, '<', $filename2 ) or die "Could not open file $!";
my @stats;
while (<$inputFH>) {
    chomp $_;
    push( @stats, $_ );
}

############Open the concatenated part1 stats file and generate a stats_cmp file which includes part1 and part2 info for all samples in the project both in the same file and include the ones that got omitted in part2 stats file with the number 0#########################################

open( INFILE3, '<', $out_file )
  or die "Could not open file $!";    ####open concatenated part1  stats file

my $outfile2 = "$path" . "/stats_file_cmp.txt";
open( OUTFILE2, ">", $outfile2 ) or die "Could not open file $!";

print OUTFILE2 "sample_id\tinput\tfiltered\tmerged\tnonchimeric\n";
while ( my $line = <INFILE3> ) {
    if ( $line !~ m/^\s*input/ ) {
        chomp $line;
        $line =~ s/\s/\t/g;
        my @info     = split( '\t', $line );
        my $sampleid = $info[0];
        print OUTFILE2 "$line\t";
        my $count = 0;
        foreach my $check (@stats) {
            my @split = split( '\t', $check );

            if ( $split[0] =~ m/^$sampleid$/ ) {
                print OUTFILE2 "$split[1]";
                $count++;
            }
        }
        if ( $count == 0
          ) ########### to include the ones that got omitted in dada2 part2 stats.txt file by adding 0 in the non-chimeric column
        {
            print OUTFILE2 "0\n";
        } else {
            print OUTFILE2 "\n";
        }
    }
}

close INFILE3;
close OUTFILE2;
}

1;
