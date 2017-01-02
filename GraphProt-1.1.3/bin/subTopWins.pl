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

subTopWins.pl

=head1 SYNOPSIS

subTopWins.pl --input $< --locations $*.top_wins --win_size $(MARGINS_WINDOW)

select windows of length max --win_size defined by --locations from --input

Options:
	-input      plain sequences
	-locations  which positions to cut out
	-win_size   max length of cut out substring
    -debug      enable debug output
    -help       brief help message
    -man        full documentation

=head1 DESCRIPTION

=cut

###############################################################################
# parse command line options
###############################################################################
my $input;
my $location;
my $win_size;
my $help;
my $man;
my $debug;
my $result = GetOptions (	"input=s"	=> \$input,
							"locations=s" => \$location,
							"win_size=i" => \$win_size,
							"help"	=> \$help,
							"man"	=> \$man,
							"debug" => \$debug);
pod2usage(-exitstatus => 1, -verbose => 1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
($result) or pod2usage(2);

###############################################################################
# main
###############################################################################

# read input
open INPUT, '<', $input;
my @input = <INPUT>;
#say STDERR "first input seq: $input[0]";

# read locations, print substrings in order as defined by the locations file
open LOC, '<', $location;
while (my $loc = <LOC>) {
	my ($seq_id, $centeredpos, $margin) = split("\t", $loc);
#	say STDERR $seq_id;
	# get sequence
	my $seq = $input[$seq_id];
	chomp $seq;
	# extract subsequence
	# sequence position read from dictionary is 0-based
#	say STDERR $seq;
#	say STDERR $pos;
#	say STDERR length($seq) - ($pos-2);
	my $pos = max(0, $centeredpos - floor($win_size/2));
	my $subseq = substr($seq, $pos, min($win_size, length($seq) - $pos));
	say $subseq;
}
