package StructureLibrary::RNAtoolsConfig;
#require Exporter;
#@ISA = qw(Exporter);
#@EXPORT = qw(
#			);
#@EXPORT_OK = qw(
#);

use strict;
use vars qw(%CONFIG);

# you may need to change some of these variables to fit to your system
our %CONFIG = (
    # All the configurations go here
    VERSION => '0.01',
    TMPDIR => '/var/tmp/',
    LIB_DIR => '/usr/local/user/RNAtools/StructureLibrary/',
    RNAFOLD => '/usr/local/vrna/1.8.5/bin/RNAfold',
    RNAFOLDPARAMS => ' -p -d2 -noLP ',
    LOCALFOLD =>  '~/bin/RNAtools/centred_hack.pl',
    LOCALFOLDPARAMS => ' -W 300 -L 150 -noLP -skipbordernts 15',
    ALL_LOCALFOLDS => '~/bin/RNAtools/local-fold.pl',
    RNAPLFOLD => '/usr/local/vrna/1.8.5/bin/RNAplfold',
    RNAPLFOLDPARAMS => ' -W 200 -L 150 -noLP -c 0',
    RNASUBOPT => '/usr/local/vrna/1.8.5/bin/RNAsubopt',
    RNASUBOPTPARAMS => '-s -p -d2 -noLP',
    RFOLD => '/home/maticzkd/cluster/src/rfold-0.1-1/src/run_rfold',
    RNASHAPES => '/usr/local/rnashapes/2.1.6/bin/RNAshapes',
    DOTVIEW => '/home/sita/bin/RNAtools/StructureLibrary/dotview.pl',
    RT	=>	0.61632, # RT = 1.98717 cal/mol/K * 310.15 = 0.61632 kcal/mol
    AUTHORS  => ["Sita Lange"],
    STR_FORMATS   => {
                   DOTBRACKET => '().-',
                   CONSTRAINT => '().|x<>',
                  },
);
1;
