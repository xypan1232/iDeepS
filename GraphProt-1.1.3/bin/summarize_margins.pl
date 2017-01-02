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

summarize_margins.pl -W WINDOW_SIZE

reads nucleotide-wise margins from stdin or file,
computes min, max, mean, median, sum for centered window of size WINDOW_SIZE

input format (tab separated): sequence id, sequence position, margin
output is input plus: min, max, mean, median, sum

=head1 SYNOPSIS

Options:

    -W          window size
    -debug      enable debug output
    -help       brief help message
    -man        full documentation

=head1 DESCRIPTION

=cut

###############################################################################
# parse command line options
###############################################################################
my $winsize;
my $help;
my $man;
my $debug;
my $result = GetOptions( "help" => \$help,
  "man"   => \$man,
  "debug" => \$debug,
  "W=i"   => \$winsize );
pod2usage( -exitstatus => 1, -verbose => 1 ) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;
($result) or pod2usage(2);

( defined $winsize ) or pod2usage("error: parameter W mandatory");
( $winsize > 0 ) or pod2usage("error: window size must be larger 0");

###############################################################################
# functions
###############################################################################

###############################################################################
# meanArray
# compute mean of values in array, ignore NA values
###############################################################################
sub meanArray {
  my ($aref) = @_;
  my $mean   = 0;
  my $nvals  = 0;
  foreach my $val (@$aref) {
    if ( $val ne "NA" ) {
      $mean += $val;
      $nvals++;
    }
  }
  if ( $nvals > 0 ) {
    $mean = $mean / $nvals;
  } else {
    $mean = "NA";
  }
  return $mean;
}

###############################################################################
# medianArray
# compute median of values in array
###############################################################################
sub medianArray {
  my ($aref) = @_;
  my @values = @$aref;

  my $median;
  my $mid = int @values / 2;
  my @sorted_values = sort { $a <=> $b } @values;
  if ( @values % 2 ) {
    $median = $sorted_values[$mid];
  } else {
    $median = ( $sorted_values[ $mid - 1 ] + $sorted_values[$mid] ) / 2;
  }

  return $median;
}

###############################################################################
# sumArray
# compute sum of values in array, ignore NA values
###############################################################################
sub sumArray {
  my ($aref) = @_;
  my $sum = 0;
  foreach my $val (@$aref) {
    if ( $val ne "NA" ) {
      $sum += $val;
    }
  }
  return $sum;
}

###############################################################################
# minArray
# compute min of values in array, ignore NA values
###############################################################################
sub minArray {
  my ($aref) = @_;
  my @vals = ();
  foreach my $val (@$aref) {
    if ( $val ne "NA" ) {
      push @vals, $val;
    }
  }
  return min @vals;
}

###############################################################################
# maxArray
# compute min of values in array, ignore NA values
###############################################################################
sub maxArray {
  my ($aref) = @_;
  my @vals = ();
  foreach my $val (@$aref) {
    if ( $val ne "NA" ) {
      push @vals, $val;
    }
  }
  return max @vals;
}

sub summarize {
  my ($aref) = @_;

  my $min    = minArray($aref);
  my $max    = maxArray($aref);
  my $mean   = meanArray($aref);
  my $median = medianArray($aref);
  my $sum    = sumArray($aref);

  return join( "\t", $min, $max, $mean, $median, $sum );
}

sub cmp_wins {

  # array of lines to summarize
  my ($linestack_aref) = @_;
  my @lines = @{$linestack_aref};

  # array of lines containing summarized values
  my @summaries;

  # queue of margins to summarize (the window)
  my @winqueue;

  # handle left border and full windows
  foreach my $line (@lines) {

    # push margin into queue
    my ( undef, undef, $margin ) = split( "\t", $line );
    push @winqueue, $margin;

    # if window size reached, remove oldest item from queue
    shift @winqueue if ( @winqueue > $winsize );

    # compute summaries if at least ceil(W/2) items in queue
    if ( @winqueue >= ceil( $winsize / 2 ) ) {

      # summarize window values
      push @summaries, summarize( \@winqueue );
    }
  }

  # handle remaining windows on right border
  while ( shift @winqueue ) {
    push @summaries, summarize( \@winqueue ) if ( @winqueue >= floor( $winsize / 2 ) + 1 )
  }

  # print output
  my $nlines     = @lines;
  my $nsummaries = @summaries;
  $debug and say STDERR "found $nsummaries summaries for $nlines lines";
  foreach my $line (@lines) {
    my $summary = shift @summaries;
    next if (not defined $summary);
    say join( "\t", $line, $summary );
  }

}

###############################################################################
# main
###############################################################################

# parse input; data for each sequence is forwarded to the windowing function
my $current_seqid = 0;    # fist sequence id is expected to be 0
my @linestack;
while ( my $line = <> ) {
  chomp $line;
  my ($seqid) = split( "\t", $line );
  if ( $seqid == $current_seqid ) {

    # if sequence id not changed, just record line
    push @linestack, $line;
  } else {

    # push lines to windowing function
    $debug and print STDERR "handling id '$current_seqid': ";
    cmp_wins( \@linestack );

    # start new linestack using current line
    @linestack = ($line);

    # update current_seqid
    $current_seqid = $seqid;
  }
}

# push lines to windowing function
$debug and print STDERR "handling id '$current_seqid': ";
cmp_wins( \@linestack )
