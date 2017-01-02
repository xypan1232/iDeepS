## generate usage message
################################################################################

usage:
	@echo
	@echo Plese run GraphProt using the command \"GraphProt.pl\".
	@echo \"GraphProt.pl --help\" will display a short summary of the required parameters.

## user defined  parameters are located here
################################################################################
include PARAMETERS


## general behaviour
################################################################################
SHELL:=$(BASH)
.DELETE_ON_ERROR:
ifeq ($(SECONDARY),YES)
# don't delete intermediate files
.SECONDARY:
endif


## paths
################################################################################
# expect binaries to reside in pwd/bin, otherwise this variable must be overwritten
PWD:=$(shell pwd)
# overwritten by GraphProt.pl to put it into same directory
TMPDIR=$(PWD)
BINDIR:=$(PWD)/bin
DATADIR:=$(PWD)/data


## project-internal tools
################################################################################
LINESEARCH:=$(PERL) $(BINDIR)/lineSearch.pl
COMBINEFEATURES:=$(PERL) $(BINDIR)/combineFeatures.pl
CREATE_EXTENDED_ACC_GRAPH:=$(PERL) $(BINDIR)/createExtendedGraph.pl
MERGE_GSPAN:=$(PERL) $(BINDIR)/merge_gspan.pl
FILTER_FEATURES:=$(PERL) $(BINDIR)/filter_features.pl
SUMMARIZE_MARGINS:=$(PERL) $(BINDIR)/summarize_margins.pl
MARGINS2BG:=$(PERL) $(BINDIR)/margins2bg.pl
VERTEX2NTMARGINS:=$(PERL) $(BINDIR)/vertex2ntmargins.pl
PLOTLC:=$(BASH) $(BINDIR)/plotlc.sh
CHECK_SYNC_GSPAN_CLASS:=$(BASH) $(BINDIR)/check_sync_gspan_class.sh
SELECT_TOP_WIN_SHREPS:=$(BINDIR)/selectTopWinShreps.R
SUBSET_NT_BY_THRESHOLD:=$(PERL) $(BINDIR)/subsetNTsbythreshold.pl
SUBSET_TOP_WINS:=$(PERL) $(BINDIR)/subTopWins.pl
MEDIAN_AWK:=$(BINDIR)/median.awk
GSPAN_SPLIT_GRAPHS:=$(BINDIR)/gspan_split_shreps.awk
SPLIT_GSPAN:=$(PERL) $(BINDIR)/split_gspan.pl
GENERICSUBMITTOCLUSTER:=$(BINDIR)/generic_submit_to_cluster.sh
MOREGENERICSGESUBMITSCRIPT:=$(BINDIR)/more_generic_submit_to_cluster.sh
GSPAN2VDICT:=$(BINDIR)/gspan2vertex_dict.awk
FASTAPL:=$(PERL) $(BINDIR)/fastapl
FASTA2GSPAN:=$(PERL) $(PWD)/fasta2shrep_gspan.pl
CAT_TABLES:=$(PERL) /home/maticzkd/co/MiscScripts/catTables.pl
SVMSGDNSPDK:=$(PWD)/EDeN/EDeN
THRESH_MARGINS:=$(PWD)/bin/subset_margins_by_quant.R

## set appropriate id (used to determine which parameter sets to use)
################################################################################
ifeq ($(SVM),SVR)
METHOD_ID=svr
endif
ifeq ($(SVM),TOPSVR)
METHOD_ID=svr
endif
ifeq ($(SVM),SGD)
METHOD_ID=sgd
endif


## set targets for RNAcompete evaluation
################################################################################
ifeq ($(EVAL_TYPE),RNACOMPETE)
# filenames for full data sets
FULL_BASENAMES:=$(patsubst %,%_data_full_A,$(PROTEINS)) \
			$(patsubst %,%_data_full_B,$(PROTEINS))

# filenames of data sets containing only weakly structured sequences
BRUIJN_BASENAMES:=$(patsubst %,%_data_bruijn_A,$(PROTEINS)) \
			$(patsubst %,%_data_bruijn_B,$(PROTEINS))

# extract prefixes for further assembling of target filenames
ifeq ($(TRAINING_SETS),FULL)
BASENAMES:=$(FULL_BASENAMES)
else
ifeq ($(TRAINING_SETS),WEAK)
BASENAMES:=$(BRUIJN_BASENAMES)
else
BASENAMES:=$(FULL_BASENAMES) $(BRUIJN_BASENAMES)
endif
endif

# general class statistics
CSTAT_FILES:=$(patsubst %,%.cstats,$(FULL_BASENAMES))

# generate staticstics on positive/negative composition
classstats : summary.cstats $(CSTAT_FILES)

endif


## set targets for CLIP-seq evaluation
################################################################################
ifeq ($(EVAL_TYPE),CLIP)
BASENAMES=$(PROTEINS)
endif


## set targets common to all evaluations
################################################################################
# parameter files (from linesearch or default values)
PARAM_FILES:=$(patsubst %,%.param,$(BASENAMES))
# results of crossvalidations
CV_FILES:=$(patsubst %,%.train.cv,$(BASENAMES))
# plots crossvalidation result curves
CV_PLOTS:=\
	$(patsubst %,%.train.cv.prplot.svg,$(BASENAMES)) \
	$(patsubst %,%.train.cv.rocplot.svg,$(BASENAMES)) \
	$(patsubst %,%.train.cv.accplot.svg,$(BASENAMES))
# models
MODEL_FILES:=$(patsubst %,%.train.model,$(BASENAMES))
# final results spearman correlation
CORRELATION_FILES:=$(patsubst %,%.test.correlation,$(BASENAMES))
# final results from perf
PERF_FILES:=$(patsubst %,%.test.perf,$(BASENAMES))
# plot prediction result curves
PERF_PLOTS:=\
	$(patsubst %,%.test.prplot.svg,$(BASENAMES)) \
	$(patsubst %,%.test.rocplot.svg,$(BASENAMES)) \
	$(patsubst %,%.test.accplot.svg,$(BASENAMES))
# nucleotide-wise margins
TESTPART_FILES:=$(patsubst %,%.test.nt_margins.summarized,$(BASENAMES))
# nucleotide-wise margins as bigWig
TESTPART_BIGWIG:=\
	$(patsubst %,%.test.nt_margins.summarized.plus.bw,$(BASENAMES)) \
	$(patsubst %,%.test.nt_margins.summarized.minus.bw,$(BASENAMES))
# files for learningcurve
LC_FILES:=$(patsubst %,%.lc.png,$(BASENAMES))

### motif stuff
# sequence motifs
SEQMOTIFS:=$(BASENAMES:%=%.test.sequence_top_wins.truncated.logo.png)
SEQMOTIFS_BIT:=$(BASENAMES:%=%.test.sequence_top_wins.truncated.logo_bit.png)
# structural context motifs
STRUCTMOTIFS:=$(BASENAMES:%=%.test.struct_annot_top_wins.truncated.logo.png)
# simplified accessibility motifs
ACCMOTIFS:=$(BASENAMES:%=%.test.struct_annot_top_wins.truncated.pup.logo.png)


## general feature and affinity creation (overridden where apropriate)
################################################################################
ifeq ($(SGEARRAY),NO)
%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan.gz %.affy | %.param
	$(CHECK_SYNC_GSPAN_CLASS) $*.gspan.gz $*.affy && \
	OMP_NUM_THREADS=1 $(SVMSGDNSPDK) -a FEATURE -i $< -r $(RADIUS) -d $(DISTANCE) -b $(BITSIZE) -g $(DIRECTED)
	cat $<.feature | grep -v \"^\$\" | paste -d' ' $*.affy - > $@
	-rm -f $<.feature
endif

ifeq ($(SGEARRAY),YES)
%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan.gz %.affy | %.param
	-rm -rf $@.FEATURE_DIR
	mkdir -p $@.FEATURE_DIR
	mkdir -p $@.FEATURE_DIR/SGEOUT
	$(CHECK_SYNC_GSPAN_CLASS) $*.gspan.gz $*.affy;
	$(SPLIT_GSPAN) -gspan_file $*.gspan.gz -feature_dir $@.FEATURE_DIR -group_size=100;
	ls -l $@.FEATURE_DIR/*.gspan.gz | wc -l > $@.NSEQS;
	echo "$(SVMSGDNSPDK) -a FEATURE -r $(RADIUS) -d $(DISTANCE) -b $(BITSIZE) -g $(DIRECTED) -i " > $@.FEATURE_DIR/edencall
	ssh `whoami`@$(SGE_SUBMIT_HOST) \
	'$(SGE_EXPORT); cd $(PWD); $(SGE_BIN_PATH)/qsub -N `echo feature_$@ | tr "/" "_"` -t 1-`cat $@.NSEQS` -o $@.FEATURE_DIR/SGEOUT -e $@.FEATURE_DIR/SGEOUT $(MOREGENERICSGESUBMITSCRIPT) $@.FEATURE_DIR $@.FEATURE_DIR/edencall'
	( for i in `seq 1 \`cat $@.NSEQS\``; \
	do cat $@.FEATURE_DIR/$$i.gspan.gz.feature; \
	done ) | \
	grep -v \"^\$\" | paste -d' ' $*.affy - > $@
	-rm -rf $@.NSEQS $@.FEATURE_DIR
endif

%.feature.gz : %.feature
	gzip -f $<

# extract affinities from fasta
# expected to reside in last field of fasta header
%.affy : %.fa
	$(FASTAPL) -e '@ry = split(/\s/,$$head); print $$ry[-1], "\n"' < $< > $@
#	$(FASTAPL) -e 'print $$head[-1], "\n"' < $< > $@

## receipes specific to graph type
################################################################################
ifeq ($(GRAPH_TYPE),SEQUENCE)
# line search parameters
LSPAR:=$(DATADIR)/ls.$(METHOD_ID).onlyseq.parameters

%.gspan.gz : VIEWPOINT=$(subst nil,,$(shell grep '^VIEWPOINT ' $*.param | cut -f 2 -d' '))
%.gspan.gz : %.fa | %.param
	$(FASTA2GSPAN) $(VIEWPOINT) --seq-graph-t -nostr -stdout -fasta $< | \
	gzip > $@; exit $${PIPESTATUS[0]}

%.sequence : %.gspan.gz
	zcat $< | awk '/^u/{print $$NF}' | tr 'tT' 'uU' > $@

%.sequence.nt_subset : %.sequence %.nt_margins
	$(SUBSET_NT_BY_THRESHOLD) -input $< -locations <(cat $*.nt_margins | \
	awk -v thresh=`cat $*.nt_margins | \
	cut -f 3 | sort -nr | \
	awk -f $(MEDIAN_AWK)` '$$3 > thresh') | \
	tr 'ACGU' '_' > $@

%.top_wins : %.nt_margins.summarized $(SELECT_TOP_WIN_SHREPS)
	$(RBIN) --slave --no-save --args $< < $(SELECT_TOP_WIN_SHREPS) | \
	sort -k3,3nr | \
	head -n $(TOP_WINDOWS) | \
	sort -k1,1n > $@
endif

################################################################################
ifeq ($(GRAPH_TYPE),SHREP)
# line search parameters
LSPAR:=$(DATADIR)/ls.$(METHOD_ID).shrep.parameters

%.gspan.gz : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.gspan.gz : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.gspan.gz : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.gspan.gz : VIEWPOINT=$(subst nil,,$(shell grep '^VIEWPOINT ' $*.param | cut -f 2 -d' '))
%.gspan.gz : %.fa | %.param
	$(FASTA2GSPAN) -fasta $< -stdout \
	--seq-graph-t --seq-graph-alph \
	$(STACK) $(CUE) $(VIEWPOINT) \
	-t $(ABSTRACTION) \
	-M $(SHREPS_MAX) \
	-wins '$(SHAPES_WINS)' \
	-shift '$(SHAPES_SHIFT)' | \
	gzip > $@; exit $${PIPESTATUS[0]}
endif

################################################################################
ifeq ($(GRAPH_TYPE),CONTEXTSHREP)
# line search parameters
LSPAR:=$(DATADIR)/ls.$(METHOD_ID).shrep.parameters

ifeq ($(SGEARRAY),YES)
%.gspan.gz : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.gspan.gz : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.gspan.gz : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.gspan.gz : VIEWPOINT=$(subst nil,,$(shell grep '^VIEWPOINT ' $*.param | cut -f 2 -d' '))
%.gspan.gz : %.fa | %.param
	-rm -rf $@.GSPAN_DIR
	$(FASTA2GSPAN) -fasta $< \
	-sge -group 100 \
	-sge-submit-user=maticzkd \
	-sge-submit-host=$(SGE_SUBMIT_HOST) \
	-sge-script $(GENERICSUBMITTOCLUSTER) \
	-o $@.GSPAN_DIR \
	--seq-graph-t --seq-graph-alph -abstr \
	$(STACK) $(CUE) $(VIEWPOINT) \
	-t $(ABSTRACTION) \
	-M $(SHREPS_MAX) \
	-wins '$(SHAPES_WINS)' \
	-shift '$(SHAPES_SHIFT)'
	NSEQS=`ls -l $@.GSPAN_DIR/*.gspan.bz2 | wc -l`; \
	(for i in `seq 1 $$NSEQS`; \
	do bzcat $@.GSPAN_DIR/$$i.group.gspan.bz2;\
	done ) | gzip > $@
	-rm -rf $@.GSPAN_DIR
endif

ifeq ($(SGEARRAY),NO)
%.gspan.gz : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.gspan.gz : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.gspan.gz : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.gspan.gz : VIEWPOINT=$(subst nil,,$(shell grep '^VIEWPOINT ' $*.param | cut -f 2 -d' '))
%.gspan.gz : %.fa | %.param
	$(FASTA2GSPAN) -fasta $< -stdout \
	--seq-graph-t --seq-graph-alph -abstr \
	$(STACK) $(CUE) $(VIEWPOINT) \
	-t $(ABSTRACTION) \
	-M $(SHREPS_MAX) \
	-wins '$(SHAPES_WINS)' \
	-shift '$(SHAPES_SHIFT)' | \
	gzip > $@; exit $${PIPESTATUS[0]}
endif

# for motif creation with contextshreps we evaluate each shrep as a single graph
%.motif.gspan.gz : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.motif.gspan.gz : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.motif.gspan.gz : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.motif.gspan.gz : VIEWPOINT=$(subst nil,,$(shell grep '^VIEWPOINT ' $*.param | cut -f 2 -d' '))
%.motif.gspan.gz : %.test.fa | %.param
	$(FASTA2GSPAN) -fasta $< -stdout \
	$(STACK) $(CUE) $(VIEWPOINT) \
	--seq-graph-t --seq-graph-alph -abstr \
	--abstr-out $*.test.struct_annot \
	-t $(ABSTRACTION) \
	-M $(SHREPS_MAX) \
	-wins '$(SHAPES_WINS)' \
	-shift '$(SHAPES_SHIFT)' | \
	awk -f $(GSPAN_SPLIT_GRAPHS) | \
	awk '/^t/{ntpos=1; print}/^[vV]/&&NF==4{print $$1,$$2,$$3,ntpos++}/^[vV]/&&NF!=4{print}/^e/{print}' | \
	gzip > $@; exit $${PIPESTATUS[0]}



# different filenames for motif creation
# compute margins of graph vertices
# vertex_margins format: seqid verticeid margin
%.motif.vertex_margins : EPOCHS=$(shell grep '^EPOCHS ' $*.param | cut -f 2 -d' ')
%.motif.vertex_margins : LAMBDA=$(shell grep '^LAMBDA ' $*.param | cut -f 2 -d' ')
%.motif.vertex_margins : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.motif.vertex_margins : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.motif.vertex_margins : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.motif.vertex_margins : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.motif.vertex_margins : %.motif.gspan.gz %.test.class %.train.model | %.param
	$(SVMSGDNSPDK) \
	-a TEST_PART \
	-m $*.train.model \
	-i $*.motif.gspan.gz \
	-t $*.test.class \
	-g $(DIRECTED) \
	-r $(RADIUS) \
	-d $(DISTANCE) \
	-b $(BITSIZE) \
	-e $(EPOCHS) \
	-l $(LAMBDA)
	mv $<.prediction_part $@

%.test.sequence : %.motif.gspan.gz
	zcat $< | awk '/^t/{print $$NF}' > $@

# compute nucleotide-wise margins from vertice margins
%.motif.nt_margins : %.motif.vertex_margins %.motif.vertex_dict
	cat $< | \
	$(VERTEX2NTMARGINS) -dict $*.motif.vertex_dict > $@

# removed this because it won't work with sequence only
# I think I added this for the viewpoint stuff
# check using testclip
#	awk '$$2!=0' > $@

%.sequence.nt_subset : %.sequence %.motif.nt_margins
	$(SUBSET_NT_BY_THRESHOLD) \
	-input $< -locations <(cat $*.motif.nt_margins | \
	awk -v thresh=`cat $*.nt_margins | \
	cut -f 3 | \
	sort -nr | \
	awk -f $(MEDIAN_AWK)` '$$3 > thresh') | \
	tr 'ACGU' '_' > $@

%.test.top_wins : %.motif.nt_margins.summarized $(SELECT_TOP_WIN_SHREPS)
	$(RBIN) --slave --no-save --args $< < $(SELECT_TOP_WIN_SHREPS) | \
	sort -k3,3nr | \
	head -n $(TOP_WINDOWS) | \
	sort -k1,1n > $@

%.struct_annot.nt_subset : %.struct_annot %.motif.nt_margins
	$(SUBSET_NT_BY_THRESHOLD) \
	-input $< \
	-locations <(cat $*.nt_margins | \
	awk -v thresh=`cat $*.nt_margins | cut -f 3 | sort -nr | awk -f $(MEDIAN_AWK)` '$$3 > 10') | \
	tr 'EHSIMB' '_' > $@

%.struct_annot_top_wins : %.struct_annot %.top_wins
	$(SUBSET_TOP_WINS) \
	--input $< \
	--locations $*.top_wins \
	--win_size $(MARGINS_WINDOW) > $@

%.pup : %
	cat $< | tr 'HBIEM' 'U' | tr 'S' 'P' > $@

%.struct_annot_top_wins.truncated.logo.png : %.struct_annot_top_wins.truncated
	cat $< | awk '{print ">"i++"\n"$$0}' | \
	$(WEBLOGO) -F png_print -o $@ \
	--alphabet 'BEHIMS' \
	--errorbars NO \
	--fineprint '' \
	--units probability \
	--color 'red' 'S' 'Stacking' \
	--color blue E External \
	--color green M Multiloop \
	--color black H Hairpin \
	--color 'dark orange' I InternalLoop \
	--color purple B Bulge \
	--show-yaxis NO \
	--show-xaxis NO

%.struct_annot_top_wins.truncated.pup.logo.png : %.struct_annot_top_wins.truncated.pup
	cat $< | awk '{print ">"i++"\n"$$0}' | \
	$(WEBLOGO) -F png_print -o $@ \
	--color-scheme classic \
	--alphabet 'UP' \
	--color red P 'Paired' \
	--color black U 'Unpaired' \
	--errorbars NO --fineprint '' \
	--units probability \
	--show-yaxis NO \
	--show-xaxis NO

# these only work with contextshreps

# calculate all motifs
motif: seqmotif structmotif accmotif

# calculate structural context motifs
structmotif: $(STRUCTMOTIFS)

# calculate simplified accessibility motifs
accmotif: $(ACCMOTIFS)
endif


## receipes specific to SVM type
################################################################################
# support vector regression
ifeq ($(SVM),SVR)
# results from cross validation
%.cv_svr : C=$(shell grep '^c ' $*.param | cut -f 2 -d' ')
%.cv_svr : EPSILON=$(shell grep '^e ' $*.param | cut -f 2 -d' ')
%.cv_svr : %.feature | %.param
	$(SVRTRAIN) -s 3 -t 0 -m $(SVR_CACHE) -c $(C) -p $(EPSILON) -h 0 -v $(CV_FOLD) $< > $@

# final result of cross validation: squared correlation coefficient
%.cv : %.cv_svr
	cat $< | \
	grep 'Cross Validation Squared correlation coefficient' | \
	$(PERL) -ne 'print /(\d+.\d+)/' > $@

# SVR model
%.model : C=$(shell grep '^c' $*.param | cut -f 2 -d' ')
%.model : EPSILON=$(shell grep '^e' $*.param | cut -f 2 -d' ')
%.model : %.feature | %.param
	$(SVRTRAIN) -s 3 -t 0 -m $(SVR_CACHE) -c $(C) -p $(EPSILON) $< $@

# SVR predictions
%.test.predictions_svr : %.train.model %.test.feature
	$(SVRPREDICT) $*.test.feature $< $@

# affinities and predictions default format
%.predictions_affy : %.predictions_svr %.affy
	paste $*.affy $< > $@

# class membership and predictions default format
%.predictions_class : %.predictions_svr %.class
	paste $*.class $< > $@

endif

## stochastic gradient descent
################################################################################
ifeq ($(SVM),SGD)
# extract single performance measure, used for linesearch decisions
%.cv : %.cv.perf
	cat $< | grep 'APR' | awk '{print $$NF}' > $@

# train model; this one directly works on gspans
%.model : EPOCHS=$(shell grep '^EPOCHS ' $*.param | cut -f 2 -d' ')
%.model : LAMBDA=$(shell grep '^LAMBDA ' $*.param | cut -f 2 -d' ')
%.model : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.model : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.model : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.model : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.model : %.feature.gz %.gspan.gz %.class | %.param
	$(CHECK_SYNC_GSPAN_CLASS) $*.gspan.gz $*.class && \
	$(SVMSGDNSPDK) -a TRAIN \
	-i $< -f SPARSE_VECTOR -t $*.class -m $@ \
	-r $(RADIUS) -d $(DISTANCE) -g $(DIRECTED) -b $(BITSIZE) \
	-e $(EPOCHS) -l $(LAMBDA)

# evaluate model
%.test.predictions_sgd : EPOCHS=$(shell grep '^EPOCHS ' $*.param | cut -f 2 -d' ')
%.test.predictions_sgd : LAMBDA=$(shell grep '^LAMBDA ' $*.param | cut -f 2 -d' ')
%.test.predictions_sgd : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.test.predictions_sgd : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.test.predictions_sgd : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.test.predictions_sgd : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.test.predictions_sgd : %.train.model %.test.feature.gz %.test.gspan.gz | %.param
	$(SVMSGDNSPDK) -a TEST \
	-m $< -i $*.test.feature.gz -f SPARSE_VECTOR -t $*.test.class \
	-r $(RADIUS) -d $(DISTANCE) -b $(BITSIZE) -g $(DIRECTED) \
	-e $(EPOCHS) -l $(LAMBDA)
	mv $*.test.feature.gz.prediction $*.test.predictions_sgd

# affinities and predictions default format
%.predictions_affy : %.predictions_sgd %.affy
	cat $< | awk '{print $$2}' | paste $*.affy - > $@

# class membership and predictions default format: class{-1,1}, prediction
%.predictions_class : %.predictions_sgd %.class
	cat $< | awk '{print $$2}' | paste $*.class - > $@

# results from crossvalidation cast into default format: class{-1,1}, prediction
%.cv.predictions_class : EPOCHS=$(shell grep '^EPOCHS ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : LAMBDA=$(shell grep '^LAMBDA ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : %.feature.gz %.gspan.gz %.class | %.param
	$(CHECK_SYNC_GSPAN_CLASS) $*.gspan.gz $*.class && \
	$(SVMSGDNSPDK) -a CROSS_VALIDATION -c $(CV_FOLD) \
	-m $*.model -i $< -f SPARSE_VECTOR -t $*.class \
	-r $(RADIUS) -d $(DISTANCE) -b $(BITSIZE) -g $(DIRECTED) \
	-e $(EPOCHS) -l $(LAMBDA) &> $@.log
	cat $<.cv_predictions | awk '{print $$2==1?1:-1, $$4}' > $@
	-rm  -f $<.cv_predictions $*.model_*

# compute margins of graph vertices
# vertex_margins format: seqid verticeid margin
%.test.vertex_margins : EPOCHS=$(shell grep '^EPOCHS ' $*.param | cut -f 2 -d' ')
%.test.vertex_margins : LAMBDA=$(shell grep '^LAMBDA ' $*.param | cut -f 2 -d' ')
%.test.vertex_margins : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.test.vertex_margins : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.test.vertex_margins : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.test.vertex_margins : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.test.vertex_margins : %.test.gspan.gz %.test.class %.train.model | %.param
	$(CHECK_SYNC_GSPAN_CLASS) $*.test.gspan.gz $*.test.class && \
	$(SVMSGDNSPDK) -a TEST_PART \
	-m $*.train.model -i $< -t $*.test.class \
	-r $(RADIUS) -d $(DISTANCE) -b $(BITSIZE) -g $(DIRECTED) \
	-e $(EPOCHS) -l $(LAMBDA) &> $@.log
	mv $<.prediction_part $@

# compute margins of graph vertices
# vertex_margins format: seqid verticeid margin
%.train.vertex_margins : EPOCHS=$(shell grep '^EPOCHS ' $*.param | cut -f 2 -d' ')
%.train.vertex_margins : LAMBDA=$(shell grep '^LAMBDA ' $*.param | cut -f 2 -d' ')
%.train.vertex_margins : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.train.vertex_margins : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.train.vertex_margins : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.train.vertex_margins : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.train.vertex_margins : %.train.gspan.gz %.train.class %.train.model | %.param
	$(CHECK_SYNC_GSPAN_CLASS) $*.test.gspan.gz $*.test.class && \
	$(SVMSGDNSPDK) -a TEST_PART \
	-m $*.train.model -i $< -t $*.train.class \
	-r $(RADIUS) -d $(DISTANCE) -b $(BITSIZE) -g $(DIRECTED) \
	-e $(EPOCHS) -l $(LAMBDA) &> $@.log
	mv $<.prediction_part $@

# dictionary of all graph vertices
# dict file format: seqid verticeid nt pos
%.vertex_dict : %.gspan.gz
	zcat $< | \
	awk -f $(GSPAN2VDICT) > $@

# compute nucleotide-wise margins from vertice margins
%.nt_margins : %.vertex_margins %.vertex_dict
	cat $< | \
	$(VERTEX2NTMARGINS) -dict $*.vertex_dict > $@

# format (tab separated): sequence id, sequence position, margin,
#                         min, max, mean, median, sum
%.nt_margins.summarized : %.nt_margins
	@echo ""
	@echo "summarizing nucleotide-wise margins:"
	$(SUMMARIZE_MARGINS) -W $(MARGINS_WINDOW) < $< > $@

%.nt_margins.perc_thresh : %.nt_margins.summarized
	@echo ""
	@echo "selecting high affinity sites ($(PERCENTILE)th percentile)"
	cat $< | \
	cut -f 1,2,6 > $@.tmp
	$(RBIN) --slave --args $(HAS_PERCENTILE) $@.tmp $@ < $(THRESH_MARGINS)
	@rm $@.tmp

%.nt_margins.summarized.plus.bedGraph : %.nt_margins.summarized %.bed
	@echo ""
	@echo "converting margins to bedGraph"
	cat $< | \
	$(MARGINS2BG) \
	-strand plus \
	-bed $*.bed \
	--aggregate $(MARGINS_MEASURE) | \
	$(BEDTOOLS) sort > $@; \
	exit $${PIPESTATUS[0]}

%.nt_margins.summarized.minus.bedGraph : %.nt_margins.summarized %.bed
	@echo ""
	@echo "converting margins to bedGraph"
	cat $< | \
	$(MARGINS2BG) \
	-strand minus \
	-bed $*.bed \
	--aggregate $(MARGINS_MEASURE) | \
	$(BEDTOOLS) sort > $@; \
	exit $${PIPESTATUS[0]}

# compute learningcurve
# svmsgdnspdk creates LEARNINGCURVE_SPLITS many files of the format
# output.lc_predictions_{test,train}_fold{1..LEARNINGCURVE_SPLITS.ID}
# format: id, class, prediction, margin
# we evaluate each one and summarize the probabilities in the following format:
# SPLIT train_performance test_performance
# this is done using svmsgdnspdk default parameters
%.lc.perf : EPOCHS=$(shell grep '^EPOCHS ' $*.train.param | cut -f 2 -d' ')
%.lc.perf : LAMBDA=$(shell grep '^LAMBDA ' $*.train.param | cut -f 2 -d' ')
%.lc.perf : RADIUS=$(shell grep '^R ' $*.train.param | cut -f 2 -d' ')
%.lc.perf : DISTANCE=$(shell grep '^D ' $*.train.param | cut -f 2 -d' ')
%.lc.perf : BITSIZE=$(shell grep '^b ' $*.train.param | cut -f 2 -d' ')
%.lc.perf : DIRECTED=$(shell grep '^DIRECTED ' $*.train.param | cut -f 2 -d' ')
%.lc.perf : %.train.feature.gz %.train.class %.train.gspan.gz
	-rm -f $*.train.dat_lc;
	LC=10; \
	NUM_REP=10; \
	lcn=$$((LC+1)); \
	for r in $$(seq 1 $$NUM_REP); \
	do $(SVMSGDNSPDK) -a LEARNING_CURVE -g $(DIRECTED) -r $(RADIUS) -d $(DISTANCE) -b $(BITSIZE) -e $(EPOCHS) -l $(LAMBDA) -i $< -f SPARSE_VECTOR -t $*.train.class -p $$lcn -R $$r > /dev/null; \
	for i in $$(seq 1 $$LC); \
	do \
	dim=$$(cat  $*.train.gspan.gz.lc_predictions_train_fold_$$i | wc -l); \
	echo -n "calculating learningcurve iteration $$i"; \
	echo -n "$$dim " >> $*.train.dat_lc; \
	cat $*.train.gspan.gz.lc_predictions_train_fold_$$i | \
	awk '{print $$2,$$4}' | $(PERF) -APR -ROC -ACC -t 0 -PRF 2> /dev/null | \
	awk '{printf("%s %s ",$$1,$$2)}END{printf("\n")}' >> $*.train.dat_lc; \
	cat $*.train.gspan.gz.lc_predictions_test_fold_$$i | awk '{print $$2,$$4}' | \
	$(PERF) -APR -ROC -ACC -t 0 -PRF 2> /dev/null | \
	awk '{printf("%s %s ",$$1,$$2)}END{printf("\n")}' >> $*.train.dat_lc; \
	done; \
	done; \
	cat $*.train.dat_lc | awk 'NR%2==1{printf("%s ",$$0)}NR%2==0{print $$0}' | column -t > $@
	-rm -f $*.train.gspan.gz.lc_predictions_t*_fold_* $*.train.dat_lc

%.lc.png : %.lc.perf
	$(PLOTLC) $< $@
	cat $@.fit_log | \
	grep 'ats' | \
	grep '=' | \
	tail -n 1 | \
	awk '{print $$3}' > $@.train_limit

endif


## evaluations specific to RNAcompete analysis
################################################################################
ifeq ($(EVAL_TYPE),RNACOMPETE)

# class memberships {-1,0,1}
%.class : BASENAME=$(firstword $(subst _, ,$<))
%.class : HT=$(shell grep $(BASENAME) $(DATA_DIR)/RNAcompete_positive_thresholds.txt | cut -f 2 -d' ')
%.class : LT=$(shell grep $(BASENAME) $(DATA_DIR)/RNAcompete_negative_thresholds.txt | cut -f 2 -d' ')
%.class : %.affy
	cat $< | \
	awk '{ if ($$1 > $(HT)) {print 1} else { if ($$1 < $(LT)) {print -1} else {print 0} } }' > $@

# some statistics about class distribution
%.cstats : BASENAME=$(firstword $(subst _, ,$<))
%.cstats : TYPE=$(word 3,$(subst _, ,$<))
%.cstats : SET=$(word 4,$(subst ., ,$(subst _, ,$<)))
%.cstats : HT=$(shell grep $(BASENAME) $(DATA_DIR)/RNAcompete_positive_thresholds.txt | cut -f 2 -d' ')
%.cstats : LT=$(shell grep $(BASENAME) $(DATA_DIR)/RNAcompete_negative_thresholds.txt | cut -f 2 -d' ')
%.cstats : HN=$(shell cat $< | grep '^>' | awk '$$NF > $(HT)' | wc -l)
%.cstats : LN=$(shell cat $< | grep '^>' | awk '$$NF < $(LT)' | wc -l)
%.cstats : %.fa
	$(PERL) -e 'print join("\t", "$(BASENAME)", "$(SET)", "$(LT)", "$(LN)", "$(HT)", "$(HN)"),"\n"' > $@

# final class summary
summary.cstats : $(CSTAT_FILES)
	( $(PERL) -e 'print join("\t", "protein", "set", "negative threshold", "negative instances", "positive threshold", "positive instances"),"\n"'; \
	cat $^ | sort -k1,2 ) > $@
endif


## evaluations specific to CLIP analysis
################################################################################
ifeq ($(EVAL_TYPE),CLIP)
# combine input sequences
%.fa : %.positives.fa %.negatives.fa %.unknowns.fa
	( $(FASTAPL) -p -1 -e '$$head .= " 1";' < $*.positives.fa; \
	  $(FASTAPL) -p -1 -e '$$head .= " -1";' < $*.negatives.fa; \
	  $(FASTAPL) -p -1 -e '$$head .= " 0";' < $*.unknowns.fa ) > $@

%.bed : %.positives.bed %.negatives.bed %.unknowns.bed
	cat $^ > $@

# for clip data, affinities are actually the class
%.class : %.affy
	cp $< $@
endif

## general motif generation
################################################################################

%.sequence_top_wins : %.sequence %.top_wins
	$(SUBSET_TOP_WINS) --input $< --locations $*.top_wins --win_size $(MARGINS_WINDOW) > $@

%.truncated : %
	cat $< | awk 'length($$0)==$(MARGINS_WINDOW)' > $@

%.sequence_top_wins.truncated.logo.png : %.sequence_top_wins.truncated
	cat $< | awk '{print ">"i++"\n"$$0}' | \
	$(WEBLOGO) -F png_print -o $@ \
	--color-scheme classic \
	--sequence-type rna \
	--errorbars NO \
	--fineprint '' \
	--units probability \
	--show-yaxis NO \
	--show-xaxis NO

%.sequence_top_wins.truncated.logo_bit.png : %.sequence_top_wins.truncated
	cat $< | awk '{print ">"i++"\n"$$0}' | \
	$(WEBLOGO) -F png_print -o $@ \
	--color-scheme classic \
	--sequence-type rna \
	--fineprint ''

## misc helper receipes
################################################################################

# if needed, unzip gspans (for svr)
%.gspan : %.gspan.gz
	zcat $< > $@

#ifeq ($(DO_LINESEARCH),NO)
## just use defaults instead of doing line search
#%.param : $(LSPAR)
#	cut -f 1,2 -d' ' < $< > $@
#else
ifeq ($(DO_SGDOPT),YES)
# do parameter optimization by line search but also use sgd-internal optimization
%.param : %.ls.fa $(LSPAR)
	$(LINESEARCH) \
	-fa $< \
	-param $(LSPAR) \
	-mf Makefile \
	-of $@ \
	-bindir $(PWD) \
	-tmp $(TMPDIR) \
	-sgdopt &> $@.log

# call sgdsvmnspdk optimization and write file containing optimized parameters
%.ls.param : EPOCHS=$(shell grep '^EPOCHS ' $*.ls_sgdopt.param | cut -f 2 -d' ')
%.ls.param : LAMBDA=$(shell grep '^LAMBDA ' $*.ls_sgdopt.param | cut -f 2 -d' ')
%.ls.param : RADIUS=$(shell grep '^R ' $*.ls_sgdopt.param | cut -f 2 -d' ')
%.ls.param : DISTANCE=$(shell grep '^D ' $*.ls_sgdopt.param | cut -f 2 -d' ')
%.ls.param : BITSIZE=$(shell grep '^b ' $*.ls_sgdopt.param | cut -f 2 -d' ')
%.ls.param : DIRECTED=$(shell grep '^DIRECTED ' $*.ls_sgdopt.param | cut -f 2 -d' ')
%.ls.param : %.ls_sgdopt.param %.ls_sgdopt.gspan.gz %.ls_sgdopt.class
	$(SVMSGDNSPDK) -a PARAMETERS_OPTIMIZATION \
	-i $*.ls_sgdopt.gspan.gz \
	-t $*.ls_sgdopt.class \
	-m $@ \
	-p $(SGDOPT_STEPS) \
	-r $(RADIUS) \
	-d $(DISTANCE) \
	-b $(BITSIZE) \
	-g $(DIRECTED) \
	-e $(EPOCHS) \
	-l $(LAMBDA) \
	&> $@.log
	( cat $< | grep -v -e '^D ' -e '^R ' -e '^EPOCHS ' -e '^LAMBDA '; \
	cat $*.ls_sgdopt.gspan.gz.opt_param | awk '{print "R",$$2,"\nD",$$4,"\nEPOCHS",$$6,"\nLAMBDA",$$8}' \
	) > $@
	rm $*.ls_sgdopt.gspan.gz.opt_param
else
# do parameter optimization by line search
%.param : %.ls.fa $(LSPAR)
	$(LINESEARCH) \
	-fa $< \
	-param $(LSPAR) \
	-mf Makefile \
	-of $@ \
	-tmp $(TMPDIR) \
	-bindir $(PWD) &> $@.log
endif
#endif

ifeq ($(LS_DO_INNER_CV),YES)
# subset fastas prior to line search
%.ls.fa : %.train.fa
	cat $< | \
	$(FASTAPL) -e 'print ">", $$head, "\t", $$seq, "\n"' | \
	$(SHUF) -n $(LINESEARCH_INPUT_SIZE) | \
	$(PERL) -ane \
	'$$seq = pop @F; $$head = join(" ", @F); print $$head, "\n", $$seq, "\n";' > \
	$@
endif

%.ls_sgdopt.fa : %.ls.fa
	cp $< $@

%.ls_sgdopt.affy : %.ls.affy
	cp $< $@

# link parameter files
%.train.param : %.param
	cp $< $@

%.test.param : %.param
	cp $< $@

# create empty input
%.unknowns.fa :
	@echo ""
	@echo "using empty set of unknowns!"
	touch $@

%.unknowns.bed :
	@echo ""
	@echo "using empty set of unknowns!"
	touch $@

# compute performance measures
# remove unknowns, set negative class to 0 for perf
%.perf : %.predictions_class
	cat $< | \
	awk '$$1!=0' | \
	sed 's/^-1/0/g' | \
	$(PERF) -confusion 2> $@.log > $@

# plot precision-recall
%.prplot : %.predictions_class
	cat $< | \
	sed 's/^-1/0/g' | \
	$(PERF) -plot pr | \
	awk 'BEGIN{p=1}/ACC/{p=0}{if (p) {print}}' > $@

%.prplot.svg : %.prplot
	cat $< | \
	$(GNUPLOT) -e "set ylabel 'precision'; set xlabel 'recall'; set terminal svg; set style line 1 linecolor rgb 'black'; plot [0:1] [0:1] '-' using 1:2 with lines;" > $@

# plot ROC curve
%.rocplot : %.predictions_class
	cat $< | \
	sed 's/^-1/0/g' | \
	$(PERF) -plot roc | \
	awk 'BEGIN{p=1}/ACC/{p=0}{if (p) {print}}' > $@

%.rocplot.svg : %.rocplot
	cat $< | \
	$(GNUPLOT) -e "set ylabel 'true positive rate'; set xlabel 'false positive rate'; set terminal svg; set style line 1 linecolor rgb 'black'; plot [0:1] [0:1] '-' using 1:2 with lines;" > $@

# plot accuracy
%.accplot : %.predictions_class
	cat $< | \
	sed 's/^-1/0/g' | \
	$(PERF) -plot roc | \
	awk 'BEGIN{p=1}/ACC/{p=0}{if (p) {print}}' > $@

%.accplot.svg : %.accplot
	cat $< | \
	$(GNUPLOT) -e "set ylabel 'accuracy'; set xlabel 'threshold'; set terminal svg; set style line 1 linecolor rgb 'black'; plot [0:1] [0:1] '-' using 1:2 with lines;" > $@


# compute correlation: correlation \t pvalue
%.correlation : %.predictions_affy
	cat $< | \
	$(RBIN) --slave \
	-e 'library(stats); data=read.table("$<", col.names=c("prediction","measurement")); t <- cor.test(data$$measurement, data$$prediction, method="spearman", alternative="greater"); write.table(cbind(t$$estimate, t$$p.value), file="$@", col.names=F, row.names=F, quote=F, sep="\t")'

results_aucpr.csv : $(PERF_FILES)
	grep -H -e APR -e ROC $^ | \
	tr ':' "\t" | \
	$(RBIN) --slave -e 'library(reshape); d<-read.table("stdin", col.names=c("id","variable","value")); write.table( cast(d), file="", row.names=F, quote=F, sep="\t")' > $@

results_correlation.csv : $(CORRELATION_FILES)
	$(CAT_TABLES) $(CORRELATION_FILES) > $@

# convert bedGraph to bigWig
%.bw : %.bedGraph $(GENOME_BOUNDS)
	$(BEDGRAPH2BIGWIG) $*.bedGraph $(GENOME_BOUNDS) $@

# # do need genome bounds
# $(GENOME_BOUNDS) :
# 	@echo ""
# 	@echo "error: require genome boundaries $@" && exit 255

## phony target section
################################################################################
.PHONY: all ls cv classstats test clean distclean learningcurve \
	motif seqmotif structmotif accmotif

# do predictions and tests for all PROTEINS, summarize results
all: $(PERF_FILES) $(CORRELATION_FILES) results_aucpr.csv results_correlation.csv

# do parameter line search for all PROTEINS
ls : $(PARAM_FILES)

# do crossvalidation
cv : $(CV_FILES)

# plot pr, roc and accuracies of crossvalidation
cvplot : $(CV_PLOTS)

# train target
train : $(MODEL_FILES)

# test target
test : $(PERF_FILES) $(CORRELATION_FILES)

# plot pr, roc and accuracies of predictions
testplot : $(PERF_PLOTS)

# compute nucleotide-wise margins
testpart : $(TESTPART_FILES)

# compute nucleotide-wise margins for genome-browser
testpart_coords : $(TESTPART_BIGWIG)

# see if additional data will help improve classification
learningcurve: $(LC_FILES)

# keep fastasd, models, predictions and results
clean:
	-rm -rf *.gspan *.gspan.gz *.threshold* *.feature *.affy \
	*.feature_filtered *.filter *.class *top_wins *.sequence *.truncated \
	*.prediction_part *_annot *.pup *.vertex_margins *.vertex_dict *.log \
	*.rocplot* *.accplot* *.prplot*

# calculate all motifs
motif: seqmotif

# calculate sequence motifs
seqmotif: $(SEQMOTIFS) $(SEQMOTIFS_BIT)

debug:
	echo $(LSPAR)

# copy all files for the GraphProt distribution
dist:
	#-rm -rf $(DIST_DIR)
	-mkdir -p $(DIST_DIR)
	# copy base progs
	cp -v GraphProt.pl fasta2shrep_gspan.pl Makefile $(DIST_DIR)
	# copy bin
	rsync -avP bin/ $(DIST_DIR)/bin --exclude=unused
	# copy data
	rsync -avP EDeN data StructureLibrary recompile_EDeN.sh $(DIST_DIR)/ --exclude='*.o'
	# copy distribution parameters
	cp PARAMETERS_dist $(DIST_DIR)/PARAMETERS
	# markdown
	cp README.md $(DIST_DIR)/README.md

# # delete all files except fastas
distclean: clean
	-rm -rf *.param *.perf *.predictions_class *.predictions_affy \
	*.predictions *.predictions_sgd *.ls.fa *.csv *model \
	*.sgeout *.class *.correlation *.cv *.cv.predictions \
	*.cv_svr *.model_* *.prplot *.prplot.svg $(LC_FILES) *.nt_margins* \
	*.logo.png *.logo_bit.png

ifeq ($(EVAL_TYPE),CLIP)
# test various stuff
runtests: testclip.param testclip.train.cv \
	testclip.test.perf testclip.test.correlation \
	testclip.test.prplot.svg
endif
ifeq ($(EVAL_TYPE),RNACOMPETE)
# test various stuff
runtests: test_data_full_A.param test_data_full_A.train.cv \
	test_data_full_A.test.perf test_data_full_A.test.correlation \
	test_data_full_A.test.prplot.svg
endif

## miscellaneous rules
################################################################################

# load genome sizes from ucsc
%.tab :
	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
	"select chrom, size from $*.chromInfo" | grep -v size > $@

# get sequence from bed using twoBitToFa
%.positives.fa : %.positives.bed
	 $(TWOBITTOFA) -bed=$< $(GENOME) $@

# get sequence from bed using twoBitToFa
%.negatives.fa : %.negatives.bed
	 $(TWOBITTOFA) -bed=$< $(GENOME) $@

# get sequence from bed using twoBitToFa
%.unknowns.fa : %.unknowns.bed
	 $(TWOBITTOFA) -bed=$< $(GENOME) $@

# get test data
testclip.train.positives.fa : $(DATADIR)/testclip.train.positives.fa
	cp -f $< $@

testclip.train.negatives.fa : $(DATADIR)/testclip.train.negatives.fa
	cp -f $< $@

testclip.test.positives.fa : $(DATADIR)/testclip.test.positives.fa
	cp -f $< $@

testclip.test.negatives.fa : $(DATADIR)/testclip.test.negatives.fa
	cp -f $< $@

test_data_full_A.test.fa : $(DATADIR)/test_data_full_A.test.fa
	cp -f $< $@

test_data_full_A.train.fa : $(DATADIR)/test_data_full_A.train.fa
	cp -f $< $@
