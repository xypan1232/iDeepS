#!/usr/bin/perl
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/ min max /;
use POSIX qw/ceil floor/;
use File::Temp qw(tempdir);
use File::Basename;
use Cwd qw/abs_path getcwd/;
use File::Copy;

=head1 NAME

GraphProt.pl

=head1 SYNOPSIS

GraphProt.pl -mode {regression,classification} -action {ls,cv,train,predict,predict_profile,predict_has,motif}

Options:

    -mode        'regression' or 'classification'
                     default: classification
    -action      what should GraphProt do?
                     ls: optimize parameters
                     cv: run a crossvalidation
                     train: train a model
                     predict: predict binding for a whole site
                     predict_profile: predict binding profiles
                     predict_has: predict high-affinity sites
                     motif: create sequence and structure motifs given a model
    -onlyseq     use GraphProt sequence models
    -prefix      this prefix is used for all results
                     default: GraphProt
    -model       GraphProt model
    -fasta       fasta file containing binding sites
    -help        brief help message
    -man         full documentation

Specify parameters as determined by -action ls:

    -params      parameter config file created by -action ls, alternatively
                 use the following options to manually enter the desired
                 parameters

Graph and Feature options:

    -abstraction RNAshapes abstraction level [RNA structure graphs]
                     default: 3
    -R           GraphProt radius
                     default: 1
    -D           GraphProt distance
                     default: 4
    -bitsize     GraphProt bitsize used for feature encoding
                     default: 14

Classification options:

    -negfasta    fasta file containing negative class sequences
    -lambda      SGD parameter lambda  [classification]
                     default: 10e-4
    -epochs      SGD parameter epochs  [classification]
                     default: 10

Regression options:

    -affinities  list of affinities
                     one value per line, same order as binding sites (fasta)
    -c           SVR parameter c       [regression]
                     default: 1
    -epsilon     SVR parameter epsilon [regression]
                     default: 0.1

Prediction options:

    -percentile  keep only regions with average score above a percentile
                 as high-affinity sites
                     default: 99

Motif options:

    -motif_len   set length of motifs
                     (default: 12)
    -motif_top_n use use n top-scoring subsequences for motif creation
                     (default: 1000)

=head1 DESCRIPTION

=cut

#Advanced options:
#
#    -keep-tmp    retain temporary files


###############################################################################
# create temporary directory
# adds an error handler that deletes the directory in case of error
# SIGUSR{1/2} are sent by the sge prior to the uncatchable SIGKILL if the
# option -notify was set
###############################################################################
$SIG{'INT'}  = 'end_handler';
$SIG{'TERM'} = 'end_handler';
$SIG{'ABRT'} = 'end_handler';
$SIG{'USR1'} = 'end_handler';
$SIG{'USR2'} = 'end_handler';

sub end_handler {
  print STDERR "signal '", $_[0], "' caught, cleaning up temporary files\n";

  # change into home directory. deletion of the temporary directory will
  # fail if it is the current working directory
  chdir();
  File::Temp::cleanup();
  die();
}

###############################################################################
# defaults
###############################################################################
my $def_abstraction = 3;
my $def_R           = 1;
my $def_D           = 4;
my $def_bitsize     = 14;
my $def_epochs      = 10;
my $def_lambda      = 10e-4;
my $def_epsilon     = 0.1;
my $def_c           = 1;
my $def_prefix      = 'GraphProt';
my $def_percentile  = 99;
my $def_motif_len   = 12;
my $def_motif_top_n = 1000;

###############################################################################
# parse command line options
###############################################################################

my $mode;
my $action;
my $onlyseq;
my $prefix;
my $fasta;
my $model;
my $params_fn;
my $negfasta;
my $affys;
my $R;
my $D;
my $c;
my $epsilon;
my $epochs;
my $lambda;
my $abstraction;
my $bitsize;
my $percentile;
my $motif_len;
my $motif_top_n;
my $help;
my $man;
my $debug;
my $keep_tmp;
my $result = GetOptions(
  "mode=s"        => \$mode,
  "action=s"      => \$action,
  "prefix=s"      => \$prefix,
  "onlyseq"       => \$onlyseq,
  "fasta=s"       => \$fasta,
  "negfasta=s"    => \$negfasta,
  "affinities=s"  => \$affys,
  "model=s"       => \$model,
  "params=s"      => \$params_fn,
  "R=i"           => \$R,
  "D=i"           => \$D,
  "c=f"           => \$c,
  "epsilon=f"     => \$epsilon,
  "epochs=i"      => \$epochs,
  "lambda=f"      => \$lambda,
  "abstraction=i" => \$abstraction,
  "bitsize=i"     => \$bitsize,
  "percentile=f"  => \$percentile,
  "motif_len=i"   => \$motif_len,
  "motif_top_n=i" => \$motif_top_n,
  "help"          => \$help,
  "man"           => \$man,
  "debug"         => \$debug,
  "keep-tmp"      => \$keep_tmp
);
pod2usage( -exitstatus => 1, -verbose => 1 ) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

###############################################################################
# check program paths
###############################################################################

my $scriptdir = abs_path( dirname($0) );

# check fastapl
`perl $scriptdir/bin/fastapl --version`;
if ( $? != 256 ) {
	print STDERR "\n\nfastapl call returned: '$?'\n";
  print STDERR "error running fastapl. perhaps your perl version is very old?\n";
  exit;
}

# check RNAshapes
`RNAshapes -h`;
if ( $? != 0 ) {
  print STDERR "please check if RNAshapes is installed and in your PATH.\n";
  exit;
}

# check EDeN
`$scriptdir/EDeN/EDeN -h`;
if ( $? != 0 ) {
  print STDERR
    "please check if the EDeN binary is executable on your system.\n";
  exit;
}

# check make
`make -h`;
if ( $? != 0 ) {
  print STDERR "please check if GNU make is installed and in your PATH.\n";
  exit;
}

# check perf
`perf`;
if ( $? != 256 ) {
  print STDERR "please check if perf is installed and in your PATH.\n";
  exit;
}

if ( defined $mode and $mode eq 'regression' ) {

  # check svm-train
  `svm-train --help`;
  if ( $? != 256 ) {
    print STDERR
"please check if libsvm is installed and in your PATH (can't call smv-train).\n";
    exit;
  }
  `svm-predict`;
  if ( $? != 256 ) {
    print STDERR
"please check if libsvm is installed and in your PATH (can't call svm-predict).\n";
    exit;
  }
}

if ( defined $action and $action eq 'motif' ) {
  `R --help`;
  if ( $? != 0 ) {
    print STDERR
      "please check if R is installed and in your path (can't call R).\n";
    exit;
  }
  `echo 'library(plyr)' | R --slave &> /dev/null`;
  if ( $? != 0 ) {
    print STDERR
      "please make sure that the R package 'plyr' is installed.\nfrom within R, call: install.packages('plyr')";
    exit;
  }
  `echo 'library(stats)' | R --slave &> /dev/null`;
  if ( $? != 0 ) {
    print STDERR
      "please make sure that the R package 'stats' is installed.\nfrom within R, call: install.packages('stats')";
    exit;
  }
  my $weblogo_version = `weblogo --version`;
  chomp $weblogo_version;
  if ( $? != 0 ) {
    print STDERR
      "please make sure that WebLogo is installed and in your PATH (can't call weblogo).\n";
    exit;
  }
}

###############################################################################
# check parameters
###############################################################################

($result) or pod2usage(2);

( defined $mode ) or $mode = 'classification';
( defined $action ) or pod2usage("please specify -action\n");

if (defined $params_fn) {
	(-f $params_fn) or pod2usage("error: file '$params_fn' not found");
}

my $quit_with_help = 0;

sub check_param_fasta {
  if ( not defined $fasta ) {
    print STDERR
      "missing parameter: please specify a set of binding sites (fasta)\n";
    $quit_with_help = 1;
    return 0;
  }
  if ( not -f $fasta ) {
    pod2usage("error: can't find file '$fasta'");
    $quit_with_help = 1;
    return 0;
  }
}

sub check_param_negfasta {
  if ( not defined $negfasta ) {
    print STDERR "missing parameter: please specify a set of unbound sites\n";
    $quit_with_help = 1;
    return 0;
  }
  if ( not -f $negfasta ) {
    pod2usage("error: can't find file '$negfasta'");
    $quit_with_help = 1;
    return 0;
  }
}

sub check_param_affys {
  if ( not defined $affys ) {
    print STDERR
      "missing parameter: please specify the list of affinities (affinities)\n";
    $quit_with_help = 1;
    return 0;
  }
  if ( not -f $affys ) {
    pod2usage("error: can't find file '$affys'");
    $quit_with_help = 1;
    return 0;
  }
}

sub check_param_model {
  if ( not defined $model ) {
    print STDERR
      "missing parameter: please specify the GraphProt model to use\n";
    $quit_with_help = 1;
    return 0;
  }
  if ( not -f $model ) {
    pod2usage("error: can't find file '$model'");
    $quit_with_help = 1;
    return 0;
  }
}

sub check_param_R {

  #    if (not defined $R) {
  #        print STDERR "missing parameter: please specify the radius (R)\n";
  #        $quit_with_help = 1;
  #        return 0;
  #    }
  defined $R or $R = $def_R;
  if ( $R < 0 ) {
    pod2usage("error: please specify a positive radius (R)");
    $quit_with_help = 1;
    return 0;
  }

  print "R: $R\n";
}

sub check_param_D {

  #    if (not defined $D) {
  #        print STDERR "missing parameter: please specify the distance (D)\n";
  #        $quit_with_help = 1;
  #        return 0;
  #    }
  defined $D or $D = $def_D;
  if ( $D < 0 ) {
    pod2usage("error: please specify a positive distance (D)");
    $quit_with_help = 1;
    return 0;
  }

  print "D: $D\n";
}

sub check_param_bitsize {

#    if (not defined $bitsize) {
#        print STDERR "missing parameter: please specify the bitsize (bitsize)\n";
#        $quit_with_help = 1;
#        return 0;
#    }
  defined $bitsize or $bitsize = $def_bitsize;
  if ( $bitsize < 8 ) {
    pod2usage(
      "error: please specify a positive bitsize larger than 8 (bitsize)");
    $quit_with_help = 1;
    return 0;
  }

  print "bitsize: $bitsize\n";
}

sub check_param_abstraction {

#    if (not defined $abstraction) {
#        print STDERR "missing parameter: please specify the RNAshapes abstraction level (abstraction)\n";
#        $quit_with_help = 1;
#        return 0;
#    }
  defined $abstraction or $abstraction = $def_abstraction;
  if ( $abstraction < 1 or $abstraction > 5 ) {
    pod2usage(
"error: please specify an RNAshapes abstraction level between 1 and 5 (abstraction)"
    );
    $quit_with_help = 1;
    return 0;
  }

  print "abstraction: $abstraction\n";
}

sub check_param_c {

#    if (not defined $c) {
#        print STDERR "missing parameter: please specify SVR parameter c (c)\n";
#        $quit_with_help = 1;
#        return 0;
#    }
  defined $c or $c = $def_c;
  if ( $c <= 0 ) {
    pod2usage("error: please specify the SVR parameter c (c)");
    $quit_with_help = 1;
    return 0;
  }

  print "c: $c\n";
}

sub check_param_epsilon {

#    if (not defined $epsilon) {
#        print STDERR "missing parameter: please specify SVR parameter epsilon (epsilon)\n";
#        $quit_with_help = 1;
#        return 0;
#    }
  defined $epsilon or $epsilon = $def_epsilon;
  if ( $c <= 0 ) {
    pod2usage("error: please specify a positive epsilon (epsilon)");
    $quit_with_help = 1;
    return 0;
  }

  print "epsilon: $epsilon\n";
}

sub check_param_epochs {

#    if (not defined $epochs) {
#        print STDERR "missing parameter: please specify SGD parameter epochs (epochs)\n";
#        $quit_with_help = 1;
#        return 0;
#    }
  defined $epochs or $epochs = $def_epochs;
  if ( $epochs <= 0 ) {
    pod2usage("error: please specify a value larger 0 (epochs)");
    $quit_with_help = 1;
    return 0;
  }

  print "epochs: $epochs\n";
}

sub check_param_lambda {
  defined $lambda or $lambda = $def_lambda;

#    if (not defined $lambda) {
#        print STDERR "missing parameter: please specify SGD parameter lambda (lambda)\n";
#        $quit_with_help = 1;
#        return 0;
#    }
  if ( $lambda <= 0 ) {
    pod2usage("error: please specify a positive lambda (lambda)");
    $quit_with_help = 1;
    return 0;
  }

  print "lambda: $lambda\n";
}

sub check_param_percentile {
  defined $percentile or $percentile = $def_percentile;

  if ( $percentile <= 0 or $percentile >= 100 ) {
    pod2usage("error: please specify a percentile in (0,100)");
    $quit_with_help = 1;
    return 0;
  }
}

sub check_param_motif_len {
    defined $motif_len or $motif_len = $def_motif_len;

    if ($motif_len < 4) {
        pod2usage("error: please specify a motif length >= 4");
        $quit_with_help = 1;
        return 0;
    }

    print "creating motif of length $motif_len\n";
}

sub check_param_motif_top_n {
    defined $motif_top_n or $motif_top_n = $def_motif_top_n;

    if ($motif_top_n < 1) {
        pod2usage("error: please specify a number > 0 for parameter motif_top_n");
        $quit_with_help = 1;
        return 0;
    }

    print "using top $motif_top_n scoring subsequences for motif creation\n";
}

sub check_params_motif {
    check_param_motif_len;
    check_param_motif_top_n;
}

sub parse_params_regression {
	if (defined $params_fn) {

		# read parameters file
		my %params;
		open PARAMS, '<', $params_fn;
		while (my $line = <PARAMS>) {
			my ($key, $value) = split ' ', $line;
			$params{$key} = $value;
		}
		close PARAMS;

		# parse parameters
		$R = $params{'R:'};
		$D = $params{'D:'};
		$bitsize = $params{'bitsize:'};
		$epsilon = $params{'epsilon:'};
		$c = $params{'c:'};
		if (not defined $onlyseq) {
			$abstraction = $params{'abstraction:'};
		}
	}
}

sub check_params_regression {
  print "Using parameters:\n";

  check_param_R;
  check_param_D;
  check_param_bitsize;
  check_param_c;
  check_param_epsilon;
  # only check when creating a structure model
  defined $onlyseq or check_param_abstraction;

  print "\n";
}

sub parse_params_classification {
	if (defined $params_fn) {

		# read parameters file
		my %params;
		open PARAMS, '<', $params_fn;
		while (my $line = <PARAMS>) {
			my ($key, $value) = split ' ', $line;
			$params{$key} = $value;
		}
		close PARAMS;

		# parse parameters
		$R = $params{'R:'};
		$D = $params{'D:'};
		$bitsize = $params{'bitsize:'};
		$epochs = $params{'epochs:'};
		$lambda = $params{'lambda:'};
		if (not defined $onlyseq) {
			$abstraction = $params{'abstraction:'};
		}
	}
}

sub check_params_classification {

	print "Using parameters:\n";

  check_param_R;
  check_param_D;
  check_param_bitsize;
  check_param_epochs;
  check_param_lambda;
  # only check when creating a structure model
  defined $onlyseq or check_param_abstraction;

  print "\n";
}

defined $action or pod2usage("please specify the GraphProt action\n");
defined $prefix or $prefix = $def_prefix;

if ( $mode eq 'regression' ) {
  if ( $action eq 'ls' ) {
    check_param_fasta;
    check_param_affys;
  } elsif ( $action eq 'cv' ) {
    check_param_fasta;
    check_param_affys;
    parse_params_regression;
    check_params_regression;
  } elsif ( $action eq 'train' ) {
    check_param_fasta;
    check_param_affys;
    parse_params_regression;
    check_params_regression;
  } elsif ( $action eq 'predict' ) {
    check_param_fasta;
    check_param_model;
    parse_params_regression;
    check_params_regression;
  } elsif ( $action eq 'predict_profile' ) {
    print STDERR "sorry, invalid action in regression setting\n";
    exit 2;
  } elsif ( $action eq 'predict_has' ) {
    print STDERR "sorry, invalid action in regression setting\n";
    exit 2;
  } elsif ( $action eq 'motif' ) {
    print STDERR "sorry, invalid action in regression setting\n";
    exit 2;
  } else {
    pod2usage("error: unknown action '$action'\n");
  }
} elsif ( $mode eq 'classification' ) {
  if ( $action eq 'ls' ) {
    check_param_fasta;
    check_param_negfasta;
  } elsif ( $action eq 'cv' ) {
    check_param_fasta;
    check_param_negfasta;
    parse_params_classification;
    check_params_classification;
  } elsif ( $action eq 'train' ) {
    check_param_fasta;
    check_param_negfasta;
    parse_params_classification;
    check_params_classification;
  } elsif ( $action eq 'predict' ) {
    check_param_fasta;
    check_param_model;
    parse_params_classification;
    check_params_classification;
  } elsif ( $action eq 'predict_profile' ) {
    check_param_fasta;
    check_param_model;
    parse_params_classification;
    check_params_classification;
  } elsif ( $action eq 'predict_has' ) {
    check_param_fasta;
    check_param_model;
    parse_params_classification;
    check_params_classification;
    check_param_percentile;
  } elsif ( $action eq 'motif' ) {
    check_param_fasta;
    check_param_model;
    parse_params_classification;
    check_params_classification;
    check_params_motif;
  } else {
    pod2usage("error: unknown action '$action'\n");
  }
} else {
  pod2usage("error: unknown mode '$mode'\n");
}

$quit_with_help and pod2usage();

# TODO: check input files

###############################################################################
# main
###############################################################################

sub parse_param_file {
  my ($param_fname) = @_;
  my $parsed_params = "";

  open LSPARS, "<", $param_fname;
  while ( my $paramline = <LSPARS> ) {
    my ( $param_name, $param_val ) = split( ' ', $paramline );

    # prettyprint parameters
    $param_name =~ s/^b/bitsize/;
    $param_name =~ s/^e/epsilon/;
    $param_name =~ s/^ABSTRACTION/abstraction/;
    $param_name =~ s/^EPOCHS/epochs/;
    $param_name =~ s/^LAMBDA/lambda/;

    # skip defaults
    next if ( $param_name =~ /DIRECTED/ );
    next if ( $param_name =~ /CUE/ );
    next if ( $param_name =~ /STACK/ );
    next if ( $param_name =~ /VIEWPOINT/ );

    $parsed_params .= "$param_name: $param_val\n";
  }
  close LSPARS;

  return $parsed_params;
}

# set up temporary directory
my $tmp_template = 'GraphProt_tmp-XXXXXX';
my $tmp_prefix   = abs_path( getcwd() );
my $tmpdir       = tempdir( $tmp_template, DIR => $tmp_prefix, CLEANUP => 1 );

# write parameters
if ( $action ne "ls" ) {
  open PARAMETERS, ">", "$tmpdir.param";
  defined $R           and print PARAMETERS "R $R\n";
  defined $D           and print PARAMETERS "D $D\n";
  defined $c           and print PARAMETERS "c $c\n";
  defined $epsilon     and print PARAMETERS "e $epsilon\n";
  defined $epochs      and print PARAMETERS "EPOCHS $epochs\n";
  defined $lambda      and print PARAMETERS "LAMBDA $lambda\n";
  defined $abstraction and print PARAMETERS "ABSTRACTION $abstraction\n";
  defined $bitsize     and print PARAMETERS "b $bitsize\n";
  print PARAMETERS "CUE nil\n";
  print PARAMETERS "STACK nil\n";
  print PARAMETERS "VIEWPOINT --vp\n";
  print PARAMETERS "DIRECTED DIRECTED\n";
  close PARAMETERS;
}

# collect make call
my $makecall = "make -C $scriptdir -e PWD=$scriptdir -e TMPDIR=${tmp_prefix}";

# use sequence graphs
if ( defined $onlyseq ) {
  $makecall .= " -e GRAPH_TYPE=SEQUENCE";
} else {

  # use structure graphs
  $makecall .= " -e GRAPH_TYPE=CONTEXTSHREP";
}

if ( $mode eq 'regression' ) {
  if ( $action eq 'ls' ) {

    # copy input files
    copy( $fasta, "$tmpdir.ls.fa" );
    copy( $affys, "$tmpdir.ls.affy" );

    # add parameters
    $makecall .= " -e SVM=SVR -e DO_LINESEARCH=YES -e DO_SGDOPT=NO";

    # add targets
    $makecall .= " $tmpdir.param";
    system("$makecall");

    # parse and report parameters
    my $lspars = parse_param_file("$tmpdir.param");
    print "optimized parameters:\n";
    print $lspars;
    open NICEPARAMS, '>', "$prefix.params";
    print NICEPARAMS $lspars;
    close NICEPARAMS;
    print "parameters written to file $prefix.params\n";
  } elsif ( $action eq 'cv' ) {

    # copy input files
    copy( $fasta, "$tmpdir.train.fa" );
    copy( $affys, "$tmpdir.train.affy" );

    # add parameters
    $makecall .= " -e SVM=SVR";

    # add targets
    $makecall .= " $tmpdir.train.cv_svr";
    system("$makecall");
    move( "$tmpdir.train.cv_svr", "$prefix.cv_results" );
    print "crossvalidation results written to file $prefix.cv_results\n";
  } elsif ( $action eq 'train' ) {

    # copy input files
    copy( $fasta, "$tmpdir.train.fa" );
    copy( $affys, "$tmpdir.train.affy" );

    # add parameters
    $makecall .= " -e SVM=SVR";

    # add targets
    $makecall .= " $tmpdir.train.model";
    system("$makecall");
    move( "$tmpdir.train.model", "$prefix.model" );
    print "GraphProt model written to file $prefix.model\n";
  } elsif ( $action eq 'predict' ) {

    # copy input files
    copy( $fasta, "$tmpdir.test.fa" );
    if ( defined $affys ) {
      copy( $affys, "$tmpdir.test.affy" );
    } else {
      open FASTA, "<", $fasta;
      open DUMMY, ">", "$tmpdir.test.affy";
      while ( my $line = <FASTA> ) {
        if ( $line =~ /^>/ ) {
          print DUMMY "1\n";
        }
      }
      close DUMMY;
      close FASTA;
    }
    copy( $model, "$tmpdir.train.model" );

    # add parameters
    $makecall .= " -e SVM=SVR";

    # add targets
    $makecall .= " $tmpdir.test.predictions_affy";
    print "$makecall\n";
    system("$makecall");

    # prepare results
    system("$scriptdir/bin/fastapl_printid.pl < ${tmpdir}.test.fa > ${tmpdir}.test.fa.id");
    system('awk "{if (\$1 > 0) {\$1=1}; if(\$1 < 0) {\$1=-1}; print \$1}"' . " < ${tmpdir}.test.predictions_svr > $tmpdir.test.predicted_class");
    system('paste -d "\t" ' . "${tmpdir}.test.fa.id $tmpdir.test.predicted_class ${tmpdir}.test.predictions_svr > $prefix.predictions");
    print "GraphProt predictions written to file $prefix.predictions\n";

  } elsif ( $action eq 'predict_profile' ) {
    print STDERR "sorry, invalid action in regression setting\n";
    exit 2;
  } elsif ( $action eq 'motif' ) {
    print STDERR "sorry, invalid action in regression setting\n";
    exit 2;
  } else {
    pod2usage("error: unknown action '$action'\n");
  }
} elsif ( $mode eq 'classification' ) {
  if ( $action eq 'ls' ) {

    # copy files
    copy( $fasta,    "$tmpdir.ls.positives.fa" );
    copy( $negfasta, "$tmpdir.ls.negatives.fa" );

    # add parameters
    $makecall .= " -e SVM=SGD  -e DO_LINESEARCH=YES";

    # add targets
    $makecall .= " $tmpdir.param";
    system("$makecall");

    # parse and report parameters
    my $lspars = parse_param_file("$tmpdir.param");
    print "optimized parameters:\n";
    print $lspars;
    open NICEPARAMS, '>', "$prefix.params";
    print NICEPARAMS $lspars;
    close NICEPARAMS;
    print "parameters written to file $prefix.params\n";
  } elsif ( $action eq 'cv' ) {

    # copy files
    copy( $fasta,    "$tmpdir.train.positives.fa" );
    copy( $negfasta, "$tmpdir.train.negatives.fa" );

    # add parameters
    $makecall .= " -e SVM=SGD ";

    # add targets
    $makecall .= " $tmpdir.train.cv.perf";
    system("$makecall");

    # prepare output
    system("$scriptdir/bin/fastapl_printid.pl < ${tmpdir}.train.fa > ${tmpdir}.train.fa.id");
    system("cat $tmpdir.train.cv.predictions_class | awk '\$1!=0' | sed 's/^-1/0/g' | " .
    	'paste -d "\t" ' . "${tmpdir}.train.fa.id - > $prefix.cv_predictions");
    move( "$tmpdir.train.cv.perf", "$prefix.cv_results" );
    print
      "GraphProt crossvalidation results written to files $prefix.cv_predictions and $prefix.cv_results\n";
  } elsif ( $action eq 'train' ) {

    # copy files
    copy( $fasta,    "$tmpdir.train.positives.fa" );
    copy( $negfasta, "$tmpdir.train.negatives.fa" );

    # add parameters
    $makecall .= " -e SVM=SGD ";

    # add targets
    $makecall .= " $tmpdir.train.model";
    system("$makecall");

    # copy model
    move( "$tmpdir.train.model", "$prefix.model" );
    print "GraphProt model written to file $prefix.model\n";
  } elsif ( $action eq 'predict' ) {

    # copy files
    copy( $model, "$tmpdir.train.model" );
    copy( $fasta, "$tmpdir.test.positives.fa" );

    # enable negfasta (create empty if not supplied)
    if ( defined $negfasta ) {
      copy( $negfasta, "$tmpdir.test.negatives.fa" );
    } else {
      system("touch $tmpdir.test.negatives.fa");
    }

    # add parameters
    $makecall .= " -e SVM=SGD ";

    # add targets
    $makecall .= " $tmpdir.test.predictions_class";
    system("$makecall");

    # prepare results
    system("$scriptdir/bin/fastapl_printid.pl < ${tmpdir}.test.fa > ${tmpdir}.test.fa.id");
    system('paste -d "\t" ' . "${tmpdir}.test.fa.id " . "${tmpdir}.test.predictions_sgd" . ' | tr " " "\t" > ' . "$prefix.predictions");
    print "GraphProt predictions written to file $prefix.predictions\n";
  } elsif ( $action eq 'predict_profile' ) {

    # copy files
    copy( $model, "$tmpdir.train.model" );
    copy( $fasta, "$tmpdir.test.positives.fa" );
    system("touch $tmpdir.test.negatives.fa");

    # add parameters
    $makecall .= " -e SVM=SGD ";

    # add targets
    $makecall .= " $tmpdir.test.nt_margins $tmpdir.test.nt_margins.summarized";
    system("$makecall");

    # copy results
    copy( "$tmpdir.test.nt_margins", "$prefix.profile" );
    print
      "GraphProt nucleotide-wise margins written to file $prefix.profile\n";
  } elsif ( $action eq 'predict_has' ) {

    # copy files
    copy( $model, "$tmpdir.train.model" );
    copy( $fasta, "$tmpdir.test.positives.fa" );
    system("touch $tmpdir.test.negatives.fa");

    # add parameters
    $makecall .= " -e SVM=SGD -e HAS_PERCENTILE=$percentile";

    # add targets
    $makecall .= " $tmpdir.test.nt_margins.perc_thresh";
    system("$makecall");

    # copy results
    copy( "$tmpdir.test.nt_margins.perc_thresh", "$prefix.has" );
    print
      "GraphProt high affinity sites written to file $prefix.has\n";
} elsif ( $action eq 'motif' ) {

    # copy files
    copy( $model, "$tmpdir.train.model" );
    copy( $fasta, "$tmpdir.test.positives.fa" );
    system("touch $tmpdir.test.negatives.fa");

    # add parameters
    $makecall .= " -e SVM=SGD ";
    $makecall .= " -e MARGINS_WINDOW=${motif_len} ";
    $makecall .= " -e TOP_WINDOWS=${motif_top_n} ";

    # add targets
    $makecall .= " $tmpdir.test.sequence_top_wins.truncated";
    $makecall .= " $tmpdir.test.sequence_top_wins.truncated.logo.png";
    system("$makecall");
    if ( not defined $onlyseq ) {
      $makecall .= " $tmpdir.test.struct_annot_top_wins.truncated";
      $makecall .= " $tmpdir.test.struct_annot_top_wins.truncated.logo.png";
      system("$makecall");
    }

    # copy results
    copy( "$tmpdir.test.sequence_top_wins.truncated",
      "$prefix.sequence_motif" );
    copy( "$tmpdir.test.sequence_top_wins.truncated.logo.png",
      "$prefix.sequence_motif.png" );
    print "GraphProt sequence motif written to file $prefix.sequence_motif.png\n";
    if ( not defined $onlyseq ) {
      copy( "$tmpdir.test.struct_annot_top_wins.truncated",
        "$prefix.structure_motif" );
      copy( "$tmpdir.test.struct_annot_top_wins.truncated.logo.png",
        "$prefix.structure_motif.png" );
      print
        "GraphProt structure motif written to file $prefix.structure_motif.png\n";
    }
  } else {
    pod2usage("error: unknown action '$action'\n");
  }
} else {
  pod2usage("error: unknown mode '$mode'\n");
}

if (defined $keep_tmp) {
	foreach my $tmpfile (glob "$tmpdir.*") {
		my $newfile = $tmpfile;
		$newfile =~ s/$tmpdir\./$prefix\./;
		$debug and print STDERR "renaming: ", $tmpfile, " -> ", $newfile, "\n";
		rename $tmpfile, $newfile;
	}
} else {
	# clean up
	unlink glob "$tmpdir.*";
	chdir();
  File::Temp::cleanup();
}
