package PathFinder;

use 5.010; 
use strict;
use warnings;
use File::Basename;
use File::Spec::Functions;
use English qw( -no_match_vars );
use Cwd qw(getcwd);


sub new {
    my ($class, %args) = @_;
    my $self = bless { %args }, $class;
    $self->open();
    return $self;
}

sub set_wd {
    my ($self, $wd) = @_;
    #$self->close();
    $self->{'wd'} = $wd;
    $self->open();
}

sub open {
    my ($self) = @_;
    open my $logFH, ">>".$self->part1_log()
        or die "Cannot open " . $self->part1_log() . " for writing: $OS_ERROR";
    my $old = select $logFH;
    $OUTPUT_AUTOFLUSH = 1;
    select $old;
    $self->{'log_tee'} = new IO::Tee(\*STDOUT, $logFH);
}

sub close {
    my ($self) = @_;
    my @handles = $self->{'log_tee'}->handles();
    foreach(@handles) {
        if ($_ ne \*STDOUT) {
            $_->close();
        }
    }
}

sub print {
    my $self = shift;
    my $msg      = shift;
    $self->{'log_tee'}->print($msg);
}

sub sample_method {
    my ($self) = @_;
    print $self->{sample_data};
}

sub clone {
    my $self = shift;
    my $copy = bless { %$self }, ref $self;
    $self->open();
    #$register{$copy} = localtime; # Or whatever else you need to do with a new object.
    
    return $copy;
}

sub part1_log {
    my $self = shift;
    my $run = basename($self->{'wd'});
    return catfile($self->{'wd'}, "${run}_16S_pipeline_log.txt");
}

sub part1_error_log {
    my $self = shift;
    return catdir($self->{'wd'}, "qsub_error_logs");
}

sub part1_stdout_log {
    my $self = shift;
    return catdir($self->{'wd'}, "qsub_stdout_logs");
}

sub part1_local_map {
    my $self = shift;
    return catdir($self->{'wd'}, basename($self->{'map'}));
}

sub fwd_demux_dir {
    my $self = shift;
    return catdir($self->{'wd'}, "fwdSplit");
}

sub rev_demux_dir {
    my $self = shift;
    return catdir($self->{'wd'}, "revSplit");
}

sub fwd_library {
    my $self = shift;
    print(catdir($self->fwd_demux_dir(), "seqs.fastq"));
    return catdir($self->fwd_demux_dir(), "seqs.fastq");
}

sub rev_library {
    my $self = shift;
    return catdir($self->rev_demux_dir(), "seqs.fastq");
}

sub fwd_sample_dir {
    my $self = shift;
    return catdir($self->fwd_demux_dir(), "split_by_sample_out");
}

sub rev_sample_dir {
    my $self = shift;
    return catdir($self->rev_demux_dir(), "split_by_sample_out");
}

sub splitsamples_output_fwd {
    my $self = shift;
    my %arg      = @_;
    my $exts     = delete %arg{"exts"};
    my @extensions = @$exts;

    my @globs = map {
        catfile($self->fwd_sample_dir(), "*.$_")
        } @extensions;
    my @result = glob(join " ", @globs);
    return @result;
}

sub splitsamples_output_rev {
    my $self = shift;
    my %arg      = @_;
    my $exts     = delete %arg{"exts"};
    my @extensions = @$exts;

    my @globs = map {
        catfile($self->rev_sample_dir(), "*.$_")
        } @extensions;
    my @result = glob(join " ", @globs);
    return @result;
}

sub filtered_dir {
    my $self = shift;
    return catdir($self->{'wd'}, "filtered");
}


# Check for qsub errors.
# @_[0] Prefix of error files to search.
sub check_error_log
{
    my $self = shift;
    my %arg      = @_;
    my $step     = delete %arg{"prefix"} // "";
    my $error;
    
    ## make $error array so that goes through all errors for that step.
    RECHECK:
    $error = glob(catfile($self->part1_error_log(), "$step.e*"));
    if (defined $error)
    {
        CORE::open my $error_fh, "<$error" or die "Can't open $error.\n";
        while (<$error_fh>)
        {
            if ($_ =~ /error/ || $_ =~ /Error/ || $_ =~ /ERROR/)
            {
                print "Error in $step. See $error.\n";
                self->log("Error in $step. See $error.\n");
                seek($error_fh, 0, 0);
                my $errorMessage = do {local $/; <$error_fh>};
                die $errorMessage;
            }
        }
        $error_fh->close();
    } else
    {
        ##print "---Re-checking for $step error file\n";
        goto RECHECK;
    }
}

1;