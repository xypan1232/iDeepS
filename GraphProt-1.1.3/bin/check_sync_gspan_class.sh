#! /bin/bash
[ $# -ne 2 ] &&  echo "Usage:  $(basename $0) GSPAN.gz class" && echo "Instead of: $(basename $0) $*"  && exit 1

# this skript will check if the class file and the (gzipped) gspan are in sync
# iE if the number of lines in the class file matches the number of graphs
# (regexp '^t') in the gspan file

GSPAN_FILE=$1
CLASS_FILE=$2

N_CLASS=`cat  ${CLASS_FILE} | wc -l`
N_GRAPHS=`zcat ${GSPAN_FILE} | grep -c '^t'`

#echo $N_GRAPHS;
#echo $N_CLASS;

[ $N_GRAPHS -ne $N_CLASS ] && echo "error: gspan and class not in sync [${N_GRAPHS} graphs vs. ${N_CLASS} class entries]" && exit 1

echo "OK: gspan and class in sync" && exit 0
