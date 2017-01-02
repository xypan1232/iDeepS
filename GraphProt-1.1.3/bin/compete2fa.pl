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

createFasta.pl

=head1 SYNOPSIS

usage: createFasta.pl --target-dir=DIR [List of files]

takes RNAcompete text files with affinity and writes accompanying fasta
file with same basename including numbered id and affinity within second
field after id

Options:

	-target-dir write fasta output to this directory
    -debug      enable debug output
    -help       brief help message
    -man        full documentation

=head1 DESCRIPTION

=cut

###############################################################################
# get filename parts of input files
###############################################################################
#my ( $prefix, $path, $suffix ) = fileparse( $infname, "\.[^.]*" );

###############################################################################
# create temporary directory
###############################################################################
#my $tmp_template = 'template-XXXXXX';
#my $tmp_prefix = '/var/tmp/';
#my $tmpdir = tempdir($tmp_template, DIR => $tmp_prefix, CLEANUP => 1);

###############################################################################
# parse command line options
###############################################################################
my $help;
my $man;
my $dir;
my $result = GetOptions( "target-dir=s" => \$dir,
  "help" => \$help,
  "man"  => \$man );
pod2usage( -exitstatus => 1, -verbose => 1 ) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;
($result) or pod2usage(2);

###############################################################################
# main
###############################################################################

while ( my $fname = shift @ARGV ) {
  my ( $prefix, $path, $suffix ) = fileparse( $fname, "\.[^.]*" );
  my $ofname = "$prefix.fa";
  my $i      = 1;

  open IN, $fname or die "error: $!";
  open OUT, ">", "$dir/$ofname" or die "error: $!";

  while ( my $line = <IN> ) {
    chomp $line;
    my ( $affinity, $sequence ) = split( /\s/, $line );
    say OUT ">${prefix}_$i $affinity";
    say OUT "$sequence";
    $i++;
  }

  close IN;
  close OUT;
}

###############################################################################
# stub
# in:
# out:
###############################################################################
sub stub {
  my ($val) = @_
}
