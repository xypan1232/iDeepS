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

vertex2ntmargins.pl -dict train.dict < NTMARGINS

=head1 SYNOPSIS

vertex2ntmargins.pl -dict train.dict -vertexoffsets train.vertex_offsets < NTMARGINS

Options:

	-dict		vertex margin dictionary
    -debug      enable debug output
    -help       brief help message
    -man        full documentation

=head1 DESCRIPTION

=cut

###############################################################################
# parse command line options
###############################################################################
my $help;
my $man;
my $debug;
my $dict;
my $vertex_offsets;
my $result = GetOptions( "help" => \$help,
  "man"    => \$man,
  "debug"  => \$debug,
  "dict=s" => \$dict );
pod2usage( -exitstatus => 1, -verbose => 1 ) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;
($result) or pod2usage(2);

###############################################################################
# main
###############################################################################

# read dictionary and fill %id2pos
# this maps sequence and vertice ids to nucleotide positions
my %id2pos;
open DICT, $dict or die $!;
while (<DICT>) {
  chomp;
  my ( $seq_id, $vertex_id, $nt, $pos ) = split(/\s/);
  $id2pos{$seq_id}{$vertex_id} = $pos;
}
close DICT;

# second step: map vertice id to nucleotide positions using dictionary
my %pos2margin;
my %pos2margin_count;
while (<>) {
  my ( $seq_id, $vertex_id, $margin ) = split(/\s/);
#  defined $id2pos{$seq_id} or die "error: nothing known about id '$seq_id' in id2pos";
	defined $id2pos{$seq_id}{$vertex_id} or next;
  my $pos = $id2pos{$seq_id}{$vertex_id};
  $pos2margin{$seq_id}{$pos} += $margin;
  $pos2margin_count{$seq_id}{$pos} += 1;
}

# average scores from different shreps
foreach my $seq_id ( keys %pos2margin ) {
    foreach my $pos ( keys %{$pos2margin{$seq_id}} ) {
        if ($pos2margin_count{$seq_id}{$pos} != 0) {
            $pos2margin{$seq_id}{$pos} /= $pos2margin_count{$seq_id}{$pos}
        }
    }
}

# print results
my @seqids = sort { $a <=> $b } ( keys %pos2margin );
for my $seq_id (@seqids) {
  my @positions = sort { $a <=> $b } ( keys %{ $pos2margin{$seq_id} } );
  for my $pos (@positions) {
    say join( "\t", $seq_id, $pos, $pos2margin{$seq_id}{$pos} );
  }
}
