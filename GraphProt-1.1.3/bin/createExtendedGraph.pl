#!/usr/bin/perl
use File::Basename;
use lib dirname($0);    # search skript directory for modules
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

createExtendedGraph.pl

=head1 SYNOPSIS

createExtendedGraph.pl -fasta=FASTA

takes sequences of fasta file and generates graph output for
fabrizio. context accessibilities are computed using a modified version
of RNAplfold that prints out a table of the context accessibilities for u=1

Options:

    -W          RNAplfold window size (default: 150)
    -L          RNAplfold bp-span (default: 100)
    -nocontext  only use accessibilities, no context probabilities
    -cutoff     do not use probabilities below or at this value (default: 0)
    -nostruct   do not compute structure part
    -vp         enable special viewpoint handling of nucleotides (default: off)
    -fasta      fasta file including input sequences
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
my $fasta;
my $path;
my $nostruct;
my $vp;
my $cutoff;
my $nocontext;
my $W;
my $L;
my $result = GetOptions( "help" => \$help,
  "man"       => \$man,
  "fasta=s"   => \$fasta,
  "nostruct"  => \$nostruct,
  "vp"        => \$vp,
  "cutoff=f"  => \$cutoff,
  "nocontext" => \$nocontext,
  "W=i"       => \$W,
  "L=i"       => \$L );
pod2usage( -exitstatus => 1, -verbose => 1 ) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;
($result) or pod2usage(2);
( defined $cutoff ) or $cutoff = 0;
( defined $W )      or $W      = 150;
( defined $L )      or $L      = 100;

##################################################################################
# This method parses a fasta file and is useful if the header ID is not-unique!!
# It returns the header lines in an array and the sequence lines in the same
# order. It is then up to the user to extract the parts from the header that is
# necessary for the script.
# Furthermore, the method deals with multiple lines, and returns a single sequence
# without line breaks.
# Input:
#       file        The name of the fasta file
# Output:
#   (1) An array reference for each header line
#   (2) An array reference for each sequence in same order as the headers
##################################################################################
sub read_fasta_with_nonunique_headers {
  my ($file) = @_;
  my $FUNCTION = "read_fasta_file in Sequences.pm";

  my $header    = "";
  my $seqstring = "";
  my @headers   = ();
  my @sequences = ();
  open( IN_HANDLE, "<$file" ) || die "ERROR in $FUNCTION:\n" .
    "Couldn't open the following file in package Tool," .
    " sub read_fasta_file: $file\n";

  while ( my $line = <IN_HANDLE> ) {
    chomp($line);

    # header (can contain one space after > symbol)
    if ( $line =~ /^\>(.*)/ ) {
      if ($header) {
        $seqstring =~ s/\s*//g;    ## do not allow spaces in sequence
        push( @headers,   $header );
        push( @sequences, $seqstring );
        $seqstring = "";
      }
      $header = $1;
    } else {
      $seqstring .= $line if ($header);
    }
  }

  if ($header) {
    $seqstring =~ s/\s*//g;        ## do not allow spaces in sequence
    push( @headers,   $header );
    push( @sequences, $seqstring );
    $seqstring = "";
  }
  return ( \@headers, \@sequences );
}

###############################################################################
# computeAccessibilities
# in: RNA sequence
# out: reference to array of accessibilities [P,E,H,I,M]
###############################################################################
sub computeAccessibilities {
  my ($seq) = @_;
  chomp $seq;
  my $len    = length($seq);
  my $W_call = min( $len, $W );
  my $L_call = min( $len, $L );

  # compute accessibilities
  open ACCS, "echo $seq | /home/maticzkd/src/ViennaRNA-1.8.4-context_stdout/Progs/RNAplfold " .
    "-u 1 -W ${W_call} -L ${L_call} |";
  my @accs = <ACCS>;
  close ACCS;

  # parse accs
  my @P;
  my @E;
  my @H;
  my @I;
  my @M;
  my $parsedseq = $accs[0];
  chomp $parsedseq;

  if ( $parsedseq ne $seq ) {
    die "error parsing accessibilities: expecting sequence '$seq', got '$parsedseq'"
  }
  for ( my $i = 1 ; $i <= $len ; $i++ ) {
    my ( $pos, $P, $E, $H, $I, $M ) = split( "\t", $accs[$i] );
    if ( $pos != $i ) {
      die "error parsing accessibilities: expecting position '$i', got '" . $accs[$i] . "'";
    }
    push @P, $P;
    push @E, $E;
    push @H, $H;
    push @I, $I;
    push @M, $M;
  }

  return [ \@P, \@E, \@H, \@I, \@M ];
}

###############################################################################
# main
###############################################################################
($fasta) or die "error: specify fasta file";
( -f $fasta ) or die "error: no such file '$fasta'";

# my ($fasta_ref, $order_aref, $header_ref) = read_fasta_file($fasta);
my ( $headers_aref, $seqs_aref ) = read_fasta_with_nonunique_headers($fasta);
my $n = 0;
foreach my $header ( @{$headers_aref} ) {
  my $seq = shift @{$seqs_aref};
  my ( $id, $affinity ) = split( /\s/, $header );

  # print out how many stuff we already computed
  $n++;
  if ( $n % 100 == 0 ) {
    print STDERR '.';
    if ( $n % 1000 == 0 ) {
      print STDERR $n / 1000, 'K';
    }
  }

  say join( ' ', 't', '#', $id, $affinity );
  my $graph_id = 0;

  # create backbone edges and vertices
  my @seq = split( //, $seq );
  my $pos = 1;    # position in sequence, 1-based
  foreach my $nt (@seq) {
    my $vertice_label;
    if ( defined $vp ) {

      # with viewpoint option enabled, set backbone vertices according to
      # uppercase/lowercase nucleotide
      $vertice_label = ( $nt =~ /[a-z]/ ) ? 'V' : 'v';
    } else {

      # otherwise just use the regular 'v'
      $vertice_label = 'v';
    }
    say join( ' ', $vertice_label, $graph_id, $nt, $pos++ );
    if ( $graph_id > 0 ) {
      say join( ' ', 'e', $graph_id - 1, $graph_id, 'b' );    # b like backbone
    }
    $graph_id++;
  }

  # skip accessibility part if requested
  next if ($nostruct);

  # compute accessibilities on the fly
  my @accs = @{ computeAccessibilities($seq) };

  # for each seq position sort context acc labels
  # and generate graph
  for ( my $pos = 0 ; $pos < length($seq) ; $pos++ ) {

    # at first create the structure context subgraph
    my %accs;
    if ($nocontext) {
      %accs = (
        'U' => $accs[0][$pos] );
    } else {
      %accs = (
        'U' => 1,              # ensure that label for unpaired is always on top
        'E' => $accs[1][$pos],
        'H' => $accs[2][$pos],
        'I' => $accs[3][$pos],
        'M' => $accs[4][$pos] );
    }

    foreach my $key ( keys %accs ) {
      if ( $accs{$key} <= $cutoff ) {
        delete $accs{$key};
      }
    }

    my $num_vertices  = keys %accs;
    my $unpaired_prob = $accs[0][$pos];
    if ( $num_vertices > 1 or $unpaired_prob > $cutoff ) {

      # if unpaired values above cutoff exist do the whole shebang

      my @keys = sort { $accs{$b} <=> $accs{$a} } keys %accs;
      my @node_ids = ( $pos, $graph_id .. ( $graph_id + $num_vertices - 1 ) );

      # save id of unpaired vertice for later
      my $unpaired_id = $graph_id;
      $graph_id += $num_vertices;

      while (@keys) {
        say join( ' ', 'v', $node_ids[1], $keys[0], '0' );

        # don't link the first vertice
        if ( @keys != $num_vertices ) {
          say join( ' ', 'e', $node_ids[0], $node_ids[1], 'a' ); # a like accessibility
        }
        shift @keys;
        shift @node_ids;
      }

      # and connect the subgraphs depending on the relation
      # of the paired and unpaired probabilities
      my $paired_prob = 1 - $unpaired_prob;
      if ( $paired_prob > $cutoff ) {

        # do it the normal way

        # create the paired probability vertice
        my $paired_id = $graph_id;
        say join( ' ', 'v', $paired_id, 'P', '0' );
        $graph_id++;

        if ( $unpaired_prob > $paired_prob ) {
          say join( ' ', 'e', $pos, $unpaired_id, 'a' );  # a like accessibility
          say join( ' ', 'e', $unpaired_id, $paired_id, 'a' ); # a like accessibility
        }
        else {
          say join( ' ', 'e', $pos, $paired_id, 'a' );    # a like accessibility
          say join( ' ', 'e', $paired_id, $unpaired_id, 'a' ); # a like accessibility
        }
      } else {

        # leave out paired probability
        say join( ' ', 'e', $pos, $unpaired_id, 'a' );    # a like accessibility
      }
    } else {

      # else just print out the paired probability

      # create the paired probability vertice
      my $paired_id = $graph_id;
      say join( ' ', 'v', $paired_id, 'P', '0' );
      $graph_id++;

      say join( ' ', 'e', $pos, $paired_id, 'a' );        # a like accessibility
    }
  }
}
