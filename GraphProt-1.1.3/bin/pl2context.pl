my @U, @E, @H, @I, @M;
while (<>) {
  if ( /^>/ and @U > 0 ) {

    # print stuff
    print '>seq_', $i++, "\n";
    my @P = map( 1 - $_, @U );
    print join( "\t", @P ), "\n";
    print join( "\t", @H ), "\n";
    print join( "\t", @I ), "\n";
    print join( "\t", @M ), "\n";
    print join( "\t", @E ), "\n";

    # delete all arrays
    @U = ();
    @E = ();
    @H = ();
    @I = ();
    @M = ();
  } elsif (/^[ACGU]+/) {

    # do nothing
  } elsif (/^\d+/) {

    # save values
    chomp;
    my ( undef, $U, $H, $I, $M, $E ) = split( "\t", $_ );
    push @U, $U;
    push @H, $H;
    push @I, $I;
    push @M, $M;
    push @E, $E;
  }
}

# print again
print '>seq_', $i++, "\n";
my @P = map( 1 - $_, @U );
print join( "\t", @P ), "\n";
print join( "\t", @H ), "\n";
print join( "\t", @I ), "\n";
print join( "\t", @M ), "\n";
print join( "\t", @E ), "\n";
