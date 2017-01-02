#!/usr/local/perl/bin/perl
use feature ':5.10';
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/ min max /;
use POSIX qw/ceil floor/;
use File::Temp qw(tempdir);
use File::Basename;

=head1 NAME

=head1 SYNOPSIS

Options:

    -shrep      shrep gspan file
    -acc        accessibility gspan file
    -debug      enable debug output
    -help       brief help message
    -man        full documentation

=head1 DESCRIPTION

=cut

###############################################################################
# parse command line options
###############################################################################
my $acc_gspan;
my $shrep_gspan;
my $help;
my $man;
my $debug;
my $result = GetOptions( "help" => \$help,
  "man"     => \$man,
  "debug"   => \$debug,
  "shrep=s" => \$shrep_gspan,
  "acc=s"   => \$acc_gspan );
pod2usage( -exitstatus => 1, -verbose => 1 ) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;
( defined $shrep_gspan and defined $shrep_gspan ) or pod2usage("error: shrep and acc parameters are mandatory");
( -f $shrep_gspan ) or die "error: file '$shrep_gspan' not found";
( -f $acc_gspan )   or die "error: file '$acc_gspan' not found";
($result)           or pod2usage(2);

###############################################################################
# main
###############################################################################
open SHREP, '<', $shrep_gspan;
open ACC,   '<', $acc_gspan;

my $read_from_shrep      = 1;
my $seen_graph_header    = 0;
my $skipped_shrep_header = '';
my $skipped_acc_header   = '';

while (1) {
  if ( $read_from_shrep and defined $skipped_shrep_header ) {
    print $skipped_shrep_header;
    $seen_graph_header    = 1;
    $skipped_shrep_header = '';
  }

  if ( not $read_from_shrep and defined $skipped_acc_header ) {
    print $skipped_acc_header;
    $seen_graph_header  = 1;
    $skipped_acc_header = '';
  }

  my $line = $read_from_shrep ? <SHREP> : <ACC>;

  # case: EOF
  if ( not defined $line ) {

    # when shrep is finished, we have to finish acc as well
    # otherwise we end
    if ($read_from_shrep) {
      $seen_graph_header = 0;
      $read_from_shrep   = 0;
      $debug and say STDERR "finished reading shreps";
      next;
    } else {
      $debug and say STDERR "finished reading accs";
      last;
    }
  }

  # case: start of graph
  if ( $line =~ /^t # / ) {
    if ($seen_graph_header) {

      # reset, read line. switch input file
      $seen_graph_header = 0;
      if ($read_from_shrep) {
        $skipped_shrep_header = $line;
      } else {
        $line =~ s/^t/a/;
        $skipped_acc_header = $line;
      }
      $read_from_shrep = 1 - $read_from_shrep;
      next;
    } else {

      # mark acc graphs
      if ( not $read_from_shrep ) {
        $line =~ s/^t/a/;
      }
      print $line;
    }
  }

  # the usual: just print stuff
  print $line;
}

close SHREP;
close ACC;
