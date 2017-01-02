#!/usr/local/perl/bin/perl
use feature ':5.10';
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/ min max /;
use POSIX qw/ceil floor/;

=head1 NAME

catTable.pl

=head1 SYNOPSIS

catTables.pl [LIST OF FILES]

append tab separated tables, add row containing source file name

Options:

    -head       assume that every input file has a one-line header;
                expand header of first input file; skip headers of remaining
                files
    -nskip      skip the first n lines of each input file
    -debug      enable debug output
    -help       brief help message
    -man        full documentation

=head1 DESCRIPTION

=cut


###############################################################################
# parse command line options
###############################################################################
my $head;
my $nskip;
my $help;
my $man;
my $result = GetOptions (	"head"	=> \$head,
							"nskip=i" => \$nskip,
							"help"	=> \$help,
							"man"	=> \$man);
pod2usage(-exitstatus => 1, -verbose => 1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
($result) or pod2usage(2);
($nskip) or ($nskip = 0);

###############################################################################
# main
###############################################################################
my $print_head = 1;
while (my $fname = shift @ARGV) {
	open IN, $fname;
	
	# skip nskip lines
	for (my $i=0; $i<$nskip; $i++) { <IN>; }
	
	if ($head) {
		my $header = <IN>;
		chomp $header;
		if ($print_head) {
			say join("\t", "fname", $header);
		}
		$print_head = 0; # only do this once
	}
	
	while (my $line = <IN>) {
		chomp $line;
		say join("\t", $fname, $line);
	}
	
	close IN;
}

###############################################################################
# stub
# in:
# out:
###############################################################################
sub stub {
    my ($val) = @_
}
