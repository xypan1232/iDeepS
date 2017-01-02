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

combineFeatures.pl

=head1 SYNOPSIS

combineFeatures.pl [features] [probabilities] [groups]

combines feature vectors by summing up. information about weights and group
membership has to be supplied

Options:

    -debug      enable debug output
    -help       brief help message
    -man        full documentation

=head1 DESCRIPTION

=cut

###############################################################################
# write features
###############################################################################
sub writeFeatures {
  my ($fh_ref) = @_;
  for my $key ( sort { $a <=> $b } keys %{$fh_ref} ) {
    print ' ', $key, ':', $fh_ref->{$key};
  }
  print "\n"
}

###############################################################################
# combine features
###############################################################################
sub combineFeatures {
  my ( $vecs_h, $probs_h ) = @_;

  # get total prob
  my $totalprob;
  for my $p ( @{$probs_h} ) {
    $totalprob += $p;
  }

  my %fh;
  for my $feats ( @{$vecs_h} ) {
    my $p = shift @{$probs_h};
    my @pairs = split( ' ', $feats );
    for my $pair (@pairs) {
      my ( $feature, $weight ) = split( ':', $pair );
      $fh{$feature} += ( $p / $totalprob ) * $weight
    }
  }

  return \%fh;
}

###############################################################################
# parse command line options
###############################################################################
my $help;
my $man;
my $debug;
my $result = GetOptions( "help" => \$help,
  "man"   => \$man,
  "debug" => \$debug );
pod2usage( -exitstatus => 1, -verbose => 1 ) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;
($result) or pod2usage(2);
( @ARGV == 3 ) or pod2usage(2);
my $feat_f = shift @ARGV;
( -f $feat_f ) or pod2usage("error reading $feat_f");
my $prob_f = shift @ARGV;
( -f $prob_f ) or pod2usage("error reading $prob_f");
my $group_f = shift @ARGV;
( -f $group_f ) or pod2usage("error reading $group_f");

###############################################################################
# main
###############################################################################
open FEAT,  $feat_f;
open PROB,  $prob_f;
open GROUP, $group_f;

my $current_gid = 1;
my @vecs;
my @probs;
while ( defined( my $feat = <FEAT> ) and defined( my $prob = <PROB> ) and defined( my $group = <GROUP> ) ) {
  if ( $current_gid == $group ) {

    # add to queue
    push @vecs,  $feat;
    push @probs, $prob;
  } else {

    # main part
    my $fh_h = combineFeatures( \@vecs, \@probs );
    writeFeatures($fh_h);

    # clean up
    undef(@vecs);
    undef(@probs);

    # redo
    $current_gid = $group;
    redo;
  }
}

# one last time
my $fh_h = combineFeatures( \@vecs, \@probs );
writeFeatures($fh_h);
