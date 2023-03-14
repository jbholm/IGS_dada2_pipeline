package PathFinder;

use 5.010; 
use strict;
use warnings;
use File::Basename;
use File::Spec::Functions;
use English qw( -no_match_vars );
use Cwd qw(getcwd);
use Hash::Merge qw( merge );

sub new {
    my ($class, %args) = @_;
    my $self = bless { %args }, $class;
    $self->open();
    mkdir $self->commands_dir();
    return $self;
}

sub set_wd {
    my ($self, $wd) = @_;
    #$self->close();
    $self->{'wd'} = $wd;
    $self->open();
}

sub set_config {
    my ($self, $config) = @_;
    $self->{'config'} = $config;
}

sub get_config {
    my $self = shift;
    return $self->{'config'};
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

sub pipeline_dir {
    my $self = shift;
    return $self->{'pipeline_dir'};
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

sub load_config_profile {
    my $self = shift;
    my $profile  = shift;

    if(exists $self->{"config"}->{'part1 param profiles'}->{"$profile"}) {
        # 'part1 params' contains universal params; 
        # "part1 param profiles" contains instrument-, PCR-, and region-specific
        # param profiles.
        my $profile = $self->{"config"}->{'part1 param profiles'}->{"$profile"};
        Hash::Merge::set_behavior('RIGHT_PRECEDENT');
        $self->{"config"}->{'part1 params'} = merge($self->{"config"}->{'part1 params'}, $profile);
        return 1;
    } else {
        return 0;
    }
}

sub override_default_param {
    my $self = shift;
    my %arg      = @_;
    my $value = delete %arg{value};
    my $name = delete %arg{name};
    my $optional = delete %arg{optional} // 0; 


    my $ans;
    if(exists $self->{'config'}->{'part1 params'}->{$name} && ! defined($value)) {
        $ans = $self->{'config'}->{'part1 params'}->{$name};
    } else {
        if (defined($value)) {
            if(exists $self->{'config'}->{'part1 params'}->{$name}) {
                $self->print("(Overriding default configured value for \"$name\" with value \"$value\" from CLI)\n");
            }
            $ans = $value;
        } else {
            if (! $optional) {
                $self->print("Config has no default value for \"$name\" and value was not provided on the command line.\n");
                die "Config has no default value for \"$name\" and value was not provided on the command line.";
            }
        }
    }

    $self->{'config'}->{'part1 params'}->{"$name"} = $ans;
    return($ans);
}

sub get_param {
    my $self = shift;
    my $param = shift;

    return $self->{"config"}->{"part1 params"}->{"$param"};
}

sub set_param {
    my $self = shift;
    my %arg  = @_;
    my $arg = \%arg;
    
    my $name = (keys %arg)[0];
    $self->{"config"}->{"part1 params"}->{"$name"} = $arg->{"$name"};

    return;
}

sub localize_file_params {
    my $self = shift;
    foreach(keys %{$self->{"config"}->{"part1 params"}}) {
        my $val = $self->get_param($_);
        my $patt = "^\\." . catfile("", "");
        if($val =~ /$patt/) {
            $self->set_param($_, catfile($self->pipeline_dir(), $val));
        }
    }
}

sub commands_dir {
    my $self = shift;
    return catdir($self->{'wd'}, "commands");
}

sub part1_local_map {
    my $self = shift;
    return catdir($self->{'wd'}, basename($self->get_param('demux_map')));
}

sub fwd_demux_dir {
    my $self = shift;
    return catdir($self->{'wd'}, "libraries", "fwd");
}

sub rev_demux_dir {
    my $self = shift;
    return catdir($self->{'wd'}, "libraries", "rev");
}

sub fwd_library {
    my $self = shift;
    return catdir($self->fwd_demux_dir(), "seqs.fastq");
}

sub rev_library {
    my $self = shift;
    return catdir($self->rev_demux_dir(), "seqs.fastq");
}

sub sample_dir {
    my $self = shift;
    return catdir($self->{'wd'}, "demultiplexed");
}

sub splitsamples_output_fwd {
    my $self = shift;
    my %arg      = @_;
    my $exts     = [ "fastq", "fastq.gz" ];
    my @extensions = @$exts;

    my @globs = map {
        catfile($self->sample_dir(), "*_R1.$_")
        } @extensions;
    my @result = glob(join " ", @globs);
    return @result;
}
sub splitsamples_output_rev {
    my $self = shift;
    my %arg      = @_;
    my $exts     = [ "fastq", "fastq.gz" ];
    my @extensions = @$exts;

    my @globs = map {
        catfile($self->sample_dir(), "*_R2.$_")
        } @extensions;
    my @result = glob(join " ", @globs);
    return @result;
}

sub trimmed_files {
    my $self = shift;
    
    my @exts = ( "fastq", "fastq.gz" );
    my @globs_f = map {
        catfile($self->{'wd'},  "R1_tc.$_")
    } @exts;
    my @globs_r = map {
        catfile($self->{'wd'},  "R2_tc.$_")
    } @exts;

    my @result = glob(join " ", (@globs_f, @globs_r));
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
                $self->print("Error in $step. See $error.\n");
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