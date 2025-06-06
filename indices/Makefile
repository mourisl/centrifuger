#
# Makefile
# fbreitwieser, 2016-01-29 13:00
# Modified by mourisl for Centrifuger.
#

SHELL := /bin/bash

THREADS?=1
KEEP_FILES?=0

get_ref_file_names = $(addprefix $(REFERENCE_SEQUENCES_DIR)/, $(addsuffix $(1), \
	$(addprefix all-,$(COMPLETE_GENOMES)) \
	$(addprefix all-,$(addsuffix -chromosome_level,$(CHROMOSOME_LEVEL_GENOMES))) \
	$(addprefix all-,$(addsuffix -any_level,$(ANY_LEVEL_GENOMES))) \
	$(addprefix mammalian-reference-,$(MAMMALIAN_TAXIDS)) \
	$(addprefix all-compressed-,$(COMPLETE_GENOMES_COMPRESSED)) \
	$(if $(INCLUDE_CONTAMINANTS),contaminants)))

DL_DIR=downloaded-seq
TMP_DIR?=tmp_$(IDX_NAME)
TAXID_SUFFIX:=.map
REFERENCE_SEQUENCES_DIR:=reference-sequences

.PHONY: index index-name index-size .path-ok .dustmasker-ok

define USAGE

Makefile to create common indices to use with Centrifuge.

  make [OPTIONS] TARGET

OPTIONS:
    THREADS=n          Number of threads for downloading, compression and
                       index building

STANDARD TARGETS:

    p						        Download all prokaryotic (archaea+bacteria) genomes from RefSeq,
                        and compresses them at the species level

    p+h+v             	 

		p+v

		v

Alternatively, a IDX_NAME and one or more genomes may be specified as
options to build a custom database.

EXTENDED OPTIONS:
	COMPLETE_GENOMES=s
	COMPLETE_GENOMES_COMPRESSED=s
	MAMMALIAN_TAXIDS=i
	INCLUDE_CONAMINANTS=1
	DONT_DUSTMASK=1
	IDX_NAME=s

EXAMPLES:
	# Make an index with all complete bacterial and archaeal genomes, and compress
	# the bacterial genomes to the species level
	make p_compressed

	# same as:
	make COMPLETE_GENOMES=archaea COMPLETE_GENOMES_COMPRESSED=bacteria IDX_NAME=p_compressed

	# Make an index with just the human genome
	make IDX_NAME=h MAMMALIAN_TAXIDS=9606

	# All archaeal genomes and contaminant sequences from UniVec and EmVec
	make IDX_NAME=a COMPLETE_GENOMES=archaea  INCLUDE_CONTAMINANTS=1

endef
export USAGE

###################################################################################################
ifndef IDX_NAME

all:
	@echo "$$USAGE"

IDX_NAME?=$(shell basename $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST)))))

INDICES=p+h+v p+m+v p+v v p refseq_microbial refseq_full nt

p+h+v: export ANY_LEVEL_GENOMES:=viral
p+h+v: export COMPLETE_GENOMES:=archaea bacteria
p+h+v: export MAMMALIAN_TAXIDS:=9606
p+h+v: export INCLUDE_CONTAMINANTS:=1
p+h+v: export IDX_NAME:=p+h+v

p+m+v: export ANY_LEVEL_GENOMES:=viral
p+m+v: export COMPLETE_GENOMES:=archaea bacteria
p+m+v: export MAMMALIAN_TAXIDS:=10090
p+m+v: export INCLUDE_CONTAMINANTS:=1
p+m+v: export IDX_NAME:=p+m+v

p+v: export ANY_LEVEL_GENOMES:=viral
p+v: export COMPLETE_GENOMES:=archaea bacteria
p+v: export INCLUDE_CONTAMINANTS:=1
p+v: export IDX_NAME:=p+v

v: export ANY_LEVEL_GENOMES:=viral
v: export IDX_NAME:=v

p: export COMPLETE_GENOMES:=archaea bacteria
p: export IDX_NAME:=p

p_compressed: export COMPLETE_GENOMES_COMPRESSED:=archaea bacteria
p_compressed: export IDX_NAME:=p_compressed

p_compressed+h+v: export ANY_LEVEL_GENOMES:=viral
p_compressed+h+v: export COMPLETE_GENOMES_COMPRESSED:=archaea bacteria
p_compressed+h+v: export MAMMALIAN_TAXIDS:=9606
p_compressed+h+v: export INCLUDE_CONTAMINANTS:=1
p_compressed+h+v: export IDX_NAME:=p_compressed+h+v

refseq_microbial: export COMPLETE_GENOMES:=archaea bacteria fungi protozoa
refseq_microbial: export CHROMOSOME_LEVEL_GENOMES:=$(COMPLETE_GENOMES)
refseq_microbial: export ANY_LEVEL_GENOMES:=viral
##refseq_microbial: export SMALL_GENOMES:=mitochondrion plasmid plastid # TODO
refseq_microbial: export MAMMALIAN_TAXIDS:=9606 10090
refseq_microbial: export INCLUDE_CONTAMINANTS:=1
refseq_microbial: export IDX_NAME:=refseq_microbial
refseq_microbial: export CFR_BUILD_OPTS+=--ftabchars 14

refseq_full: export COMPLETE_GENOMES:=archaea bacteria fungi invertebrate plant protozoa vertebrate_mammalian vertebrate_other viral
refseq_full: export CHROMOSOME_LEVEL_GENOMES:=$(COMPLETE_GENOMES)
refseq_full: export ANY_LEVEL_GENOMES:=viral
refseq_full: export SMALL_GENOMES:=mitochondrion plasmid plastid
refseq_full: export MAMMALIAN_TAXIDS:=9606 10090
refseq_full: export INCLUDE_CONTAMINANTS:=1
refseq_full: export IDX_NAME:=refseq_full


nt: export IDX_NAME:=nt

$(INDICES):
	@echo Making: $@: $(IDX_NAME)
	$(MAKE) -f $(THIS_FILE) IDX_NAME=$(IDX_NAME)

####################################################################################################
else ## IDX_NAME is defined

DONT_DUSTMASK=1
TAXONOMY_DOWNLOAD_OPTS?=
REFERENCE_SEQUENCES=$(call get_ref_file_names,.fna)
TAXID_MAPS=$(call get_ref_file_names,$(TAXID_SUFFIX))
CFR_BUILD_OPTS?=

ifeq (nt,$(IDX_NAME))
ifeq ($(strip $(DONT_DUSTMASK)),)
REFERENCE_SEQUENCES+=nt-dusted.fna
else
REFERENCE_SEQUENCES+=nt-sorted.fna
endif
TAXID_MAPS+=nt.map
CFR_BUILD_OPTS+=--ftabchars=14
endif


ifeq ($(strip $(REFERENCE_SEQUENCES)),)
$(error REFERENCE_SEQUENCES is not set - specify at lease one of COMPLETE_GENOMES, \
COMPLETE_GENOMES_COMPRESSED, or MAMMALIAN_TAXIDS with the IDX_NAME ($(IDX_NAME)))
endif

SIZE_TABLES=$(addprefix $(REFERENCE_SEQUENCES_DIR)/all-compressed-,$(addsuffix .size,$(COMPLETE_GENOMES_COMPRESSED)))
ifneq ($(strip $(COMPLETE_GENOMES_COMPRESSED)),)
#CFR_BUILD_OPTS+=--size-table <(cat $(SIZE_TABLES))
endif

CF_DOWNLOAD_OPTS?=
CF_COMPRESS_OPTS?=
ifeq ($(strip $(DONT_DUSTMASK)),)
CF_DOWNLOAD_OPTS+=-m
else
CF_COMPRESS_OPTS+=--noDustmasker
endif

all: $(IDX_NAME).1.cfr

# vim:ft=make
endif ## ifndef IDX_NAME

$(REFERENCE_SEQUENCES_DIR):
	mkdir -p $(REFERENCE_SEQUENCES_DIR)

#$(TAXID_MAPS): | $(REFERENCE_SEQUENCES_DIR)
#	rm $(patsubst %$(TAXID_SUFFIX),%.fna, $@)
#	$(MAKE) -f $(THIS_FILE) $(patsubst %$(TAXID_SUFFIX),%.fna, $@)

nt.gz:
	curl -o nt.gz ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz

nt.fna: nt.gz
	gunzip -c nt.gz > nt.fna

accession2taxid/nucl_gb.accession2taxid.gz:
	mkdir -p accession2taxid
	curl ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz > accession2taxid/nucl_gb.accession2taxid.gz

accession2taxid/nucl_wgs.accession2taxid.gz:
	mkdir -p accession2taxid
	curl ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz > accession2taxid/nucl_wgs.accession2taxid.gz

nt.map: nt-sorted.fna

nt-sorted.fna: nt.fna accession2taxid/nucl_gb.accession2taxid.gz accession2taxid/nucl_wgs.accession2taxid.gz
	centrifuge-sort-nt.pl -m nt.map -a nt-acs-wo-mapping.txt \
		nt.fna accession2taxid/nucl_gb.accession2taxid.gz accession2taxid/nucl_wgs.accession2taxid.gz \
		> nt-sorted.fna

nt-dusted.fna: nt-sorted.fna | .dustmasker-ok
	 dustmasker -infmt fasta -in nt-sorted.fna -level 20 -outfmt fasta | sed '/^>/! s/[^AGCT]/N/g' > nt-dusted.fna

$(REFERENCE_SEQUENCES_DIR)/mammalian-reference-%.fna: | $(REFERENCE_SEQUENCES_DIR)
	@[[ -d $(TMP_DIR) ]] && rm -rf $(TMP_DIR); mkdir -p $(TMP_DIR)
	centrifuger-download -o $(TMP_DIR) -d "vertebrate_mammalian" -a "Chromosome" -t $* -c 'reference genome' -P $(THREADS) refseq > \
		$(TMP_DIR)/$(patsubst %.fna,%$(TAXID_SUFFIX), $(notdir $@))
	find $(TMP_DIR)/vertebrate_mammalian -name "*.fna.gz" | xargs zcat > $@.tmp && mv $@.tmp $@
	mv $(TMP_DIR)/$(patsubst %.fna,%$(TAXID_SUFFIX),$(notdir $@)) $(patsubst %.fna,%$(TAXID_SUFFIX),$@)
ifeq (1,$(KEEP_FILES))
	[[ -d $(DL_DIR)/vertebrate_mammalian ]] || mkdir -p $(DL_DIR)/vertebrate_mammalian
	mv $(TMP_DIR)/vertebrate_mammalian/* $(DL_DIR)/vertebrate_mammalian
else
	rm -rf $(TMP_DIR)
endif

$(REFERENCE_SEQUENCES_DIR)/all-compressed-%.fna: | $(REFERENCE_SEQUENCES_DIR) taxonomy/nodes.dmp taxonomy/names.dmp .dustmasker-ok
	[[ -d $(TMP_DIR) ]] && rm -rf $(TMP_DIR); mkdir -p $(TMP_DIR)
	centrifuger-download -o $(TMP_DIR) -d "$*" -P $(THREADS) refseq > $(TMP_DIR)/all-$*.map
	time centrifuge-compress.pl $(TMP_DIR)/$* taxonomy $(CF_COMPRESS_OPTS) -map $(TMP_DIR)/all-$*.map \
		-o $@.tmp -t $(THREADS) -maxG 50000000 2>&1 | tee centrifuge-compress-$(IDX_NAME).log && \
	mv $@.tmp.fa $@ && mv $@.tmp.size $(patsubst %.fna,%.size,$@) && \
	mv $@.tmp.map $(patsubst %.fna,%$(TAXID_SUFFIX),$@)
ifeq (1,$(KEEP_FILES))
	[[ -d $(DL_DIR)/$* ]] || mkdir -p $(DL_DIR)/$*
	mv $(TMP_DIR)/$*/* $(DL_DIR)/$*
else
	rm -rf $(TMP_DIR)
endif

$(REFERENCE_SEQUENCES_DIR)/all-%-chromosome_level.fna: | $(REFERENCE_SEQUENCES_DIR) .dustmasker-ok
	[[ -d $(TMP_DIR) ]] && rm -rf $(TMP_DIR); mkdir -p $(TMP_DIR)
	@echo Downloading and dust-masking $*
	centrifuger-download -o $(TMP_DIR) $(CF_DOWNLOAD_OPTS) -a "Chromosome" -d "$*" -P $(THREADS) refseq > \
		$(TMP_DIR)/$(patsubst %.fna,%$(TAXID_SUFFIX),$(notdir $@))
	find $(TMP_DIR)/$* -name "*.fna.gz" | xargs zcat > $@.tmp && mv $@.tmp $@
	mv $(TMP_DIR)/$(patsubst %.fna,%$(TAXID_SUFFIX),$(notdir $@)) $(patsubst %.fna,%$(TAXID_SUFFIX),$@)
ifeq (1,$(KEEP_FILES))
	[[ -d $(DL_DIR)/$* ]] || mkdir -p $(DL_DIR)/$*
	mv $(TMP_DIR)/$*/* $(DL_DIR)/$*
else
	rm -rf $(TMP_DIR)
endif

$(REFERENCE_SEQUENCES_DIR)/all-%-any_level.fna: | $(REFERENCE_SEQUENCES_DIR) .dustmasker-ok
	[[ -d $(TMP_DIR) ]] && rm -rf $(TMP_DIR); mkdir -p $(TMP_DIR)
	@echo Downloading and dust-masking $*
	centrifuger-download -o $(TMP_DIR) $(CF_DOWNLOAD_OPTS) -a "Any" -d "$*" -P $(THREADS) refseq > \
		$(TMP_DIR)/$(patsubst %.fna,%$(TAXID_SUFFIX),$(notdir $@))
	find $(TMP_DIR)/$* -name "*.fna.gz" | xargs zcat > $@.tmp && mv $@.tmp $@
	mv $(TMP_DIR)/$(patsubst %.fna,%$(TAXID_SUFFIX),$(notdir $@)) $(patsubst %.fna,%$(TAXID_SUFFIX),$@)
ifeq (1,$(KEEP_FILES))
	[[ -d $(DL_DIR)/$* ]] || mkdir -p $(DL_DIR)/$*
	mv $(TMP_DIR)/$*/* $(DL_DIR)/$*
else
	rm -rf $(TMP_DIR)
endif

$(REFERENCE_SEQUENCES_DIR)/all-%.fna: | $(REFERENCE_SEQUENCES_DIR) .dustmasker-ok
	[[ -d $(TMP_DIR) ]] && rm -rf $(TMP_DIR); mkdir -p $(TMP_DIR)
	@echo Downloading and dust-masking $*
	centrifuger-download -o $(TMP_DIR) $(CF_DOWNLOAD_OPTS) -d "$*" -P $(THREADS) refseq > \
		$(TMP_DIR)/$(patsubst %.fna,%$(TAXID_SUFFIX),$(notdir $@))
	find $(TMP_DIR)/$* -name "*.fna.gz" | xargs zcat > $@.tmp && mv $@.tmp $@
	mv $(TMP_DIR)/$(patsubst %.fna,%$(TAXID_SUFFIX),$(notdir $@)) $(patsubst %.fna,%$(TAXID_SUFFIX),$@)
ifeq (1,$(KEEP_FILES))
	[[ -d $(DL_DIR)/$* ]] || mkdir -p $(DL_DIR)/$*
	mv $(TMP_DIR)/$*/* $(DL_DIR)/$*
else
	rm -rf $(TMP_DIR)
endif

$(REFERENCE_SEQUENCES_DIR)/contaminants.fna: | $(REFERENCE_SEQUENCES_DIR)
	[[ -d $(TMP_DIR) ]] && rm -rf $(TMP_DIR); mkdir -p $(TMP_DIR)
	centrifuger-download -o $(TMP_DIR) contaminants > $(TMP_DIR)/$(patsubst %.fna,%$(TAXID_SUFFIX),$(notdir $@))
	#find $(TMP_DIR)/contaminants -name "*.fna.gz" | xargs zcat > $@.tmp && mv $@.tmp $@
	find $(TMP_DIR)/contaminants -name "*.fna" | xargs cat > $@.tmp && mv $@.tmp $@
	mv $(TMP_DIR)/$(patsubst %.fna,%$(TAXID_SUFFIX),$(notdir $@)) $(patsubst %.fna,%$(TAXID_SUFFIX),$@)
ifeq (1,$(KEEP_FILES))
	[[ -d $(DL_DIR)/contaminants ]] || mkdir -p $(DL_DIR)/contaminants
	mv $(TMP_DIR)/contaminants/* $(DL_DIR)/$*
else
	rm -rf $(TMP_DIR)
endif

DUSTMASKER_EXISTS := $(shell command -v dustmasker)
.dustmasker-ok:
ifndef DUSTMASKER_EXISTS
ifeq ($(strip $(DONT_DUSTMASK)),)
	$(error dustmasker program does not exist. Install NCBI blast+, or set option DONT_DUSTMASK=1)
endif
endif


taxonomy/names.dmp: taxonomy/nodes.dmp
taxonomy/nodes.dmp: | .path-ok
	[[ -d $(TMP_DIR) ]] && rm -rf $(TMP_DIR); mkdir -p $(TMP_DIR)
	centrifuger-download $(TAXONOMY_DOWNLOAD_OPTS) -o $(TMP_DIR)/taxonomy taxonomy
	mkdir -p taxonomy
	mv $(TMP_DIR)/taxonomy/* taxonomy && rmdir $(TMP_DIR)/taxonomy && rmdir $(TMP_DIR)

$(IDX_NAME).1.cfr: $(REFERENCE_SEQUENCES) $(SIZE_TABLES) $(TAXID_MAPS) taxonomy/nodes.dmp taxonomy/names.dmp | .path-ok
	@echo Index building prerequisites: $^
	[[ -d $(TMP_DIR) ]] && rm -rf $(TMP_DIR); mkdir -p $(TMP_DIR)
	cat $(TAXID_MAPS) > $(REFERENCE_SEQUENCES_DIR)/cat_taxid.map 
	time centrifuger-build -t $(THREADS) $(CFR_BUILD_OPTS) \
		--conversion-table $(REFERENCE_SEQUENCES_DIR)/cat_taxid.map \
		--taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
		-r $(call join_w_r,$(REFERENCE_SEQUENCES)) -o $(TMP_DIR)/$(IDX_NAME) 2>&1 | tee centrifuger-build-$(IDX_NAME).log
	mv $(TMP_DIR)/$(IDX_NAME).*.cfr . && rmdir $(TMP_DIR)


clean:
	# Removing input sequences (all required information is in the index)
	rm -rf taxonomy
	rm -rf $(DL_DIR)
	rm -rf $(TMP_DIR)
	rm -rf tmp_*
	rm -rf reference-sequences
	rm -f *.map
	rm -f *.log

# Join a list with commas
COMMA:=,
EMPTY:=
SPACE:= $(EMPTY) $(EMPTY)
CFR_REFARG:=$(SPACE)-r$(SPACE)
join_w_comma = $(subst $(SPACE),$(COMMA),$(strip $1))
join_w_r = $(subst $(SPACE),$(CFR_REFARG),$(strip $1))


THIS_FILE := $(lastword $(MAKEFILE_LIST))
PATH_OK  := $(shell command -v centrifuger-build 2> /dev/null && command -v centrifuger-download 2> /dev/null )
CF_BASE_DIR := $(shell dirname $(shell dirname $(THIS_FILE)))

error_msg := centrifuger-download and centrifuger-build are not available - please make sure they are in the path.
define n


endef

TEST_PROGRAMS=centrifuger-build centrifuger-download

ifneq ("$(wildcard $(CF_BASE_DIR)/centrifuger-build)","")
error_msg := $(error_msg)$n$nThe following command may solve this problem:$n  export PATH=$$PATH:"$(CF_BASE_DIR)"$n
endif

.path-ok:
ifndef PATH_OK
    $(error $n$(error_msg))
else
	@echo Found centrifuger-download and centrifuger-build.
endif

index-name:
	echo $(IDX_NAME)

index-size:
	du -csh $(IDX_NAME).[1234].cfr

