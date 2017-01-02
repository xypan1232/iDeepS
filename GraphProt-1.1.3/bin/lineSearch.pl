#!/usr/bin/perl
use feature ':5.10';
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/ min max /;
use POSIX qw/ceil floor/;
use File::Temp qw(tempdir);
use File::Basename;
use File::Copy;
use Cwd;
use Scalar::Util qw(looks_like_number);

=head1 NAME

lineSearch.pl

=head1 SYNOPSIS

lineSearch.pl -fa fasta_for_graphs

Options:

    -fa         optimize parameters for this fasta
    -param      parameter definition for linesearch
    -mf         makefile to use for crossvalidation
    -bindir     expect RBPaffinity-specific binaries to be in this directory
    -of         write optimal parameters to this file
    -sgdopt     utilize sgd-internal optimization
                this has several effects:
                * current parameters are written to intermediate parameter files
                * intermediate parameter files contain default and extreme values
                  for parameters optimized directly by sgd: R, D, EPOCHS, LAMBDA
    -tmp				prefix for temporary directory (default: /var/tmp/)
    -debug      enable debug output
    -help       brief help message
    -man        full documentation

=head1 DESCRIPTION

=cut

###############################################################################
# create temporary directory
# adds an error handler that deletes the directory in case of error
# SIGUSR{1/2} are sent by the sge prior to the uncatchable SIGKILL if the
# option -notify was set
###############################################################################
my $tmp_template = 'lineSearch-XXXXXX';
my $tmp_prefix   = '/var/tmp/';
$SIG{'INT'}  = 'end_handler';
$SIG{'TERM'} = 'end_handler';
$SIG{'ABRT'} = 'end_handler';
$SIG{'USR1'} = 'end_handler';
$SIG{'USR2'} = 'end_handler';

sub end_handler {
  print STDERR "signal ", $_[0], " caught, cleaning up temporary files\n";

  # change into home directory. deletion of the temporary directory will
  # fail if it is the current working directory
  chdir();
  File::Temp::cleanup();
  die();
}

###############################################################################
# parse command line options
###############################################################################
my $fa;
my $param;
my $mf;
my $bindir;
my $of;
my $sgdopt;
my $help;
my $man;
my $tmpopt;
my $debug;
my $result = GetOptions(
  "fa=s"     => \$fa,
  "param=s"  => \$param,
  "mf=s"     => \$mf,
  "bindir=s" => \$bindir,
  "of=s"     => \$of,
  "sgdopt"   => \$sgdopt,
  "tmp=s"		 => \$tmpopt,
  "help"     => \$help,
  "man"      => \$man,
  "debug"    => \$debug
);
pod2usage( -exitstatus => 1, -verbose => 1 ) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;
($result)       or pod2usage(2);
( defined $fa ) or pod2usage("error: -fa parameter mandatory");
( -f $fa )      or die "error: could not find file '$fa'";
( -f $mf )      or die "error: could not find file '$mf'";
( -f $param )   or die "error: could not find file '$param'";

###############################################################################
# main
###############################################################################

$tmp_prefix = $tmpopt if (defined $tmpopt);

# global variables
my $CURRDIR = cwd();
($debug) and say STDERR "cwd: $CURRDIR";
my $top_correlation = 0;
my @top_rounds      = 0;
my $n_rounds        = 1;
my $basename        = $fa;
$basename =~ s/.fa//;
my %result_cache;

# these parameters can be optimized internally by eden
my %sgd_internal_optimization =
  ( 'R' => 1, 'D' => 1, 'EPOCHS' => 1, 'LAMBDA' => 1 );

# gspan creation depends on these parameters
my %gspan_params =
  ( 'ABSTRACTION' => 1, 'STACK' => 1, 'CUE' => 1, 'VIEWPOINT' => 1 );

# binaries
my $libsvm         = '~/src/libsvm-3.0/svm-train';
my $libsvm_options = ' -v 5 -s 3 -t 0';

# we optimize these parameters
my @parameters;

# valid values for parameters
my %parameters;
open PARAM, $param;
while (<PARAM>) {
  next if /^\s*#/;
  my ( $id, $default, @values ) = split ' ';
  push @parameters, $id;
  $parameters{$id}{default} = $default;
  $parameters{$id}{current} = $default;
  $parameters{$id}{values}  = \@values;

  # if doing the sge-internal optimization, we just use the rightmost value
  # specified in the parameter file
  if ( defined $sgdopt and defined $sgd_internal_optimization{$id} ) {
    $parameters{$id}{default} = $values[-1];
    $parameters{$id}{current} = $values[-1];
  }
}

# print important variables
if ($debug) {
  say STDERR 'parameters to optimize: ', join( ', ', @parameters );
  say STDERR 'keys in hash: ', join( ', ', keys %{ $parameters{'epsilon'} } );
  while ( my ( $param, $param_h ) = each %parameters ) {
    while ( my ( $key, $values ) = each %{$param_h} ) {
      if ( $key eq 'values' ) {
        say STDERR join( '/', $param, $key ), ': ', join( ', ', @{$values} );
      } else {
        say STDERR join( '/', $param, $key ), ': ', $values;
      }
    }
  }
}

################################################################################
# given the current parameters,
# determine the filename to be used for caching purposes.
# the filename should encode all neccessary parameters
sub get_cached_gspan_filename {
  my ( $basename, @parameters ) = @_;
  my $parameter_string;

  for my $par (@parameters) {

    # use parameter if gspan depends on it
    if ( defined $gspan_params{$par} ) {
      $parameter_string .= $par . $parameters{$par}{current};
    }
  }

  return $basename . '.' . $parameter_string;
}

# create temporary directory to cache gspan files
my $gspan_cache_dir_obj = File::Temp->newdir(
  "GspanCache_" . $tmp_template,
  DIR     => $tmp_prefix,
  CLEANUP => 1
);
my $gspan_cache_dir = $gspan_cache_dir_obj->dirname;
say STDERR "CACHE DIRECTORY: $gspan_cache_dir";

# main loop: do until finished
my $optimization_finished = 0;
do {
  $debug and say STDERR '';
  $debug and say STDERR 'starting round ' . $n_rounds;

  # use this to track if we have to do an evaluation even though there
  # is always only one choice per parameter
  my $started_at_least_one_evaluation = 0;

  # optimize each parameter
  # if doing sgd-internal optimization, do this first
  if ( defined $sgdopt ) {

    # do sgd-internal optimization and save optimal values under current key
    my $tmpdir_obj =
      File::Temp->newdir( $tmp_template, DIR => $tmp_prefix, CLEANUP => 1 );
    my $tmpdir     = $tmpdir_obj->dirname;
    my $param_file = $of;

    say STDERR "\n*** optimizing sgd-internal parameters";

    # rewrite to fit to line search input fasta filename
    $param_file =~ s/.param/.ls_sgdopt.param/;
    $debug and say STDERR "writing temporaray parameters to '$param_file'";
    my $cv_file = $basename . '.cv';

    # copy relevant files to tmp
    copy( $fa, $tmpdir );

    # get cached gspan filename
    my $cached_gspan_fname =
      get_cached_gspan_filename( $basename, @parameters );
    my $cached_gspan_fname_full = "$gspan_cache_dir/$cached_gspan_fname";

    # if exists, copy to temp dir
    if ( -f $cached_gspan_fname_full ) {
      $debug
        and say STDERR "reusing cached gspan $cached_gspan_fname"
        . " -> $tmpdir/${basename}_sgdopt.gspan.gz";
      copy( $cached_gspan_fname_full, "$tmpdir/$basename.gspan.gz" );
    }
    copy( 'PARAMETERS',                $tmpdir );
    copy( 'EXPERIMENT_SPECIFIC_RULES', $tmpdir );

    # test parameter combination / get value from previous run
    # create parameter file
    chdir($tmpdir);
    open PARS, '>', $param_file or die "error: can't open $param_file for writing";
    print STDERR 'parameters: ';

    # distinguish parameter source by parameter
    for my $par (@parameters) {
      my $parameter_value;
      if ( defined $sgd_internal_optimization{$par} ) {

        # if the parameter is optimized by the sgd, always use value strored in
        # default field
        $parameter_value = $parameters{$par}{default};
      } else {

        # all other parameters should use the values that are currently best
        $parameter_value = $parameters{$par}{current};
      }
      print STDERR $par, ' ', $parameter_value, ";";
      say PARS $par, ' ', $parameter_value;
    }
    print STDERR "\n";
    close PARS;

    # in a first step prepare parameter file
    # this will do the sgd-internal optimization
    my $final_param_file = $of;
    $final_param_file =~ s/.param/.ls.param/;
    my $prepare_param =
      "make -f $bindir/$mf -e PWD=$bindir $final_param_file";
    $debug and say STDERR "calling '$prepare_param'";
    system("$prepare_param") == 0 or die "$prepare_param failed: $?";

    # look up new best parameters from file
    my $pfile = $of;
    $pfile =~ s/.param/.ls.param/;
    open PFILE, '<', $pfile or die "error opening file '$pfile'";
    while (<PFILE>) {
      my ( $par, $val ) = split(/\s/);

      # set new current parameters based on sgd-optimization
      $parameters{$par}{current} = $val;
    }

    # if cached gspan filename not exists, copy to cache
    if ( not -f $cached_gspan_fname_full ) {
      $debug and say STDERR "caching gspan $cached_gspan_fname";
      copy( "$tmpdir/$basename.gspan.gz", $cached_gspan_fname_full );
    }

    # exit temp directory
    chdir($CURRDIR);
  }

  # optimize other parameters
  for my $par (@parameters) {

    # some parameters are varied by sgd
    if ( defined $sgdopt and defined $sgd_internal_optimization{$par} ) {
      $debug
        and say 'skipping parameter because of sgd-internal optimization: '
        . $par;
      # exit temp directory
      chdir($CURRDIR);
      next;
    }

    # some parameters don't vary
    # anyway, we have to ensure that we get at least one result
    # so in case of a non-varying parameter, do the analysis anyway
    # if none has been done so far
    if ( @{ $parameters{$par}{values} } <= 1
      and $started_at_least_one_evaluation )
    {
      $debug and say 'skipping parameter that does not vary: ' . $par;
      # exit temp directory
      chdir($CURRDIR);
      next;
    }
    $started_at_least_one_evaluation = 1;

    say STDERR "\n*** optimizing parameter $par, "
      . "round: $n_rounds, "
      . "current best: $top_correlation";
    for my $try_this ( @{ $parameters{$par}{values} } ) {

      # set new parameter
      $parameters{$par}{current} = $try_this;

      my $error;
      my $correlation;

      # get current parameter key for hash lookup
      my $param_key = get_current_param_key();
      my $cached_value_found;

      # look for cached result
      if ($cached_value_found) {
        say STDERR 'cached value found: ', $param_key;

        # retrieve cached result and be done
        $correlation = $result_cache{$param_key};
      } else {
        my $tmpdir_obj =
          File::Temp->newdir( $tmp_template, DIR => $tmp_prefix, CLEANUP => 1 );
        my $tmpdir = $tmpdir_obj->dirname;

        my $param_file = $of;

        # rewrite to fit to line search input fasta filename
        $param_file =~ s/.param/.ls.param/;
        $debug and say STDERR "writing temporaray parameters to '$param_file'";
        my $cv_file = $basename . '.cv';

        # copy relevant files into tmp
        copy( $fa, $tmpdir );

        # get cached gspan filename
        my $cached_gspan_fname =
          get_cached_gspan_filename( $basename, @parameters );
        my $cached_gspan_fname_full = "$gspan_cache_dir/$cached_gspan_fname";

        # if exists, copy to temp dir
        if ( -f $cached_gspan_fname_full ) {
          $debug and say STDERR "reusing cached gspan $cached_gspan_fname";
          copy( $cached_gspan_fname_full, "$tmpdir/$basename.gspan.gz" );
        }
        copy( 'PARAMETERS',                $tmpdir );
        copy( 'EXPERIMENT_SPECIFIC_RULES', $tmpdir );

        # test parameter combination / get value from previous run
        # create parameter file
        chdir($tmpdir);
        open PARS, '>', $param_file;
        print STDERR 'parameters: ';
        for my $par (@parameters) {
          my $parameter_value = $parameters{$par}{current};
          print STDERR $par, ' ', $parameter_value, ";";
          say PARS $par, ' ', $parameter_value;
        }
        print STDERR "\n";
        close PARS;

        # call Makefile for cv
        my $exec =
          "make cv -f $bindir/$mf -e CV_FILES=$cv_file -e PWD=$bindir";
        $debug and say STDERR $exec;
        system("$exec") == 0 or die "$exec failed: $?";

        # parse result
        open RES, '<', $cv_file;
        $correlation = <RES>;
        close RES;
        chomp $correlation;

        # save result for later reference
        $result_cache{ get_current_param_key() } = $correlation;
        say STDERR "correlation: $correlation";

        # if cached gspan filename not exists, copy to cache
        if ( not -f $cached_gspan_fname_full ) {
          $debug and say STDERR "caching gspan $cached_gspan_fname";
          copy( "$tmpdir/$basename.gspan.gz", $cached_gspan_fname_full );
        }

        # exit temp directory
        chdir($CURRDIR);
      }

      # save state info
      if ( $correlation > $top_correlation ) {

        # if the new result is better, save these parameters
        $top_correlation = $correlation;
        for my $par (@parameters) {
          $parameters{$par}{currentbest} = $parameters{$par}{current};
        }
      }
    }

    # set current to the best parameter combination
    for my $par (@parameters) {
      $parameters{$par}{current} = $parameters{$par}{currentbest};
    }
  }

  # do a maximum of 3 rounds
  if ( $n_rounds++ > 3 ) {
    say STDERR "\n";
    say STDERR "maximum of 3 rounds reached, stopping";
    $optimization_finished = 1;
  }

  # stop if the last round improved correlation by less than 0.01
  push @top_rounds, $top_correlation;
  if ( $top_rounds[-1] - $top_rounds[-2] < 0.01 ) {
    say STDERR "\n";
    say STDERR "improvement to last round < 0.01, stopping";
    $optimization_finished = 1;
  }
} while ( not $optimization_finished );

say STDERR "top values from rounds: ", join( '; ', @top_rounds );

# write final parameters
open OUT, '>', $of;
for my $par (@parameters) {
  print STDERR $par, ' ', $parameters{$par}{current}, ";";
  say OUT $par, ' ', $parameters{$par}{current};
}
print STDERR "\n";
close OUT;

chdir();

sub get_current_param_key {
  my $current_param_key;
  for my $par (@parameters) {
    my $value = $parameters{$par}{current};
    $current_param_key .= $par . ';' . $value . ';';
  }
  return $current_param_key;
}
