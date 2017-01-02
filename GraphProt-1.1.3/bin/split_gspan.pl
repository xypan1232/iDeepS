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
use File::Copy;

=head1 NAME

split_gspan.pl -param_file $*.param -gspan_file $*.gspan.gz -feature_dir $@.FEATURE_DIR

=head1 SYNOPSIS

Options:

    -gspan_file  split this gspan
    -feature_dir put split gspans into this directory
    -group_size  split into groups of group_size graphs (default: 10)
    -debug       enable debug output
    -help        brief help message
    -man         full documentation

=head1 DESCRIPTION

=cut

###############################################################################
# parse command line options
###############################################################################
my $help;
my $man;
my $debug;
my $param_file;
my $gspan_file;
my $feature_dir;
my $group_size;
my $result = GetOptions(
  "help"          => \$help,
  "man"           => \$man,
  "debug"         => \$debug,
  "feature_dir=s" => \$feature_dir,
  "gspan_file=s"  => \$gspan_file,
  "group_size=i"  => \$group_size
);
pod2usage( -exitstatus => 1, -verbose => 1 ) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;
( defined $group_size ) or $group_size = 10;
( $group_size > 0 )
  or
  pod2usage("error: choose positive integer as group_size (is '$group_size')");
($result) or pod2usage(2);

###############################################################################
# main
###############################################################################
open my $gspan, "zcat $gspan_file | ";

# count graphs
my $graph_id = 0;

# number used for split gspan file
my $gspan_id = 0;

# filehandle for split gspans
my $splitfile;

# filename for split gspans
my $outfname;

while ( my $gspan_line = <$gspan> ) {

  # at new graph start
  if ( $gspan_line =~ /^t/ ) {

    # write groupsize graphs per file
    if ( $graph_id % $group_size == 0 ) {

      # close old file
      if ( defined $splitfile ) {
        close $splitfile;
      }

      # create new file for writing
      $gspan_id++;
      $outfname = $feature_dir . "/" . $gspan_id . ".gspan.gz";
      open $splitfile, "| gzip > $outfname";

    }

    # count new graph
    $graph_id++;

  }

  # print to split gspan file
  print $splitfile $gspan_line;
}
