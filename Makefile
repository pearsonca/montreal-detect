
default:
	@echo helloworld

R := /usr/bin/env Rscript

SIMPATH := ../scala-commsim
DIGPATH := ../montreal-digest
WKDIR := detection

include ../muri-overall/references.mk
include $(SIMPATH)/dims.mk
include $(DIGPATH)/dims.mk

seq2 = $(strip $(shell for i in $$(seq $(1) $(2)); do printf '%03d ' $$i; done))
seq = $(call seq2,1,$(1))

BPATH := $(INDIR)/$(WKDIR)
DATAPATH := $(INDIR)/digest
OUTSRC := $(INDIR)/simulate/covert

$(BPATH): | $(INDIR)
	mkdir -p $@

define factorial1dir
$(BPATH)/$(1): | $(BPATH)
	mkdir -p $$@

endef

define factorial2dir
$(BPATH)/$(1): | $(BPATH)/$(dir $(1))
	mkdir -p $$@

endef

DIMVARS := USAGE PWRS TIMS SZS INTERVALS WINDOWS SCORING
SAMPIDS := $(call seq,$(SAMPN))

#$(foreach d,$(DIMVARS),$(info $($(d))))

popped = $(wordlist 2,$(words $(1)),$(1))

define indim
$(call factorial2dir,$(1))

$(if $(2),$(foreach d,$($(firstword $(2))),$(call indim,$(1)/$(d),$(call popped,$(2)))))
endef

define topdim
$(call factorial1dir,$(1))

$(foreach d,$($(firstword $(2))),$(call indim,$(1)/$(d),$(call popped,$(2))))
endef

$(foreach d,$($(firstword $(DIMVARS))),$(eval $(call topdim,$(d),$(call popped,$(DIMVARS)))))

define sampdirs
$(BPATH)/$(1)/$(2) : | $(BPATH)/$(1)
	mkdir -p $$@
endef

ALLDETECTBASEPBS :=

SNAPBASEPBS :=

define detectingbase
$(BPATH)/$(1)/$(2)/%/trim.rds: trim.R $(DATAPATH)/raw/pairs.rds $(DATAPATH)/raw/location-lifetimes.rds $(OUTSRC)/$(1)/%/cc.csv $(OUTSRC)/$(1)/%/cu.csv
	mkdir -p $$(dir $$@)
	$(R) $$^ $(subst /,$(SPACE),$(2)) > $$@

# % = sample
$(BPATH)/$(1)/$(2)/%/base.rds: pre-spinglass-detect.R $(DATAPATH)/raw/pairs.rds $(DATAPATH)/background/$(2)/base $(BPATH)/$(1)/$(2)/%/trim.rds | $(BPATH)/$(1)/$(2)
	$(R) $$^ $(subst /,$(SPACE),$(2)) > $$@

$(BPATH)/$(1)/$(2)/snapFTPR.rds: base_review.R $(wildcard $(BPATH)/$(1)/$(2)/*/base.rds)
	$(R) $$< $(BPATH)/$(1)/$(2) $(DATAPATH)/background/$(2)/base > $$@

snaps-$(subst /,-,$(1))-$(subst /,-,$(2)).pbs: snap-review.sh
	./$$^ $(subst /,-,$(1))-$(subst /,-,$(2)) $(1)/$(2) > $$@

SNAPBASEPBS += snaps-$(subst /,-,$(1))-$(subst /,-,$(2)).pbs

base-$(subst /,-,$(1))-$(subst /,-,$(2)).pbs: base-detect.sh
	./$$^ $(subst /,-,$(1))-$(subst /,-,$(2)) $(1)/$(2) $(SAMPN) > $$@

ALLDETECTBASEPBS += base-$(subst /,-,$(1))-$(subst /,-,$(2)).pbs

.PRECIOUS: $(BPATH)/$(1)/$(2)/%/trim.rds $(BPATH)/$(1)/$(2)/%/base.rds

endef

ALLDETECTPBS :=

ALLFINALPCPBS :=
# # loop over covert dims, then analysis dims, then sample N
define detecting # 1 dir for covert, 2 is dir for detection
# % = sample
$(BPATH)/$(1)/$(2)/%/acc: pre-spinglass-score.R $(DATAPATH)/background/$(dir $(2))base $(BPATH)/$(1)/$(dir $(2))%/base.rds $(BPATH)/$(1)/$(dir $(2))%/trim.rds
	mkdir -p $$@
	$(R) $$^ $(subst /,$(SPACE),$(2)) $$@

# need to get sample number in here somehow, but shouldn't be an issue
detect-$(subst /,-,$(1))-$(subst /,-,$(2)).pbs: acc-detect.sh
	./$$^ $(subst /,-,$(1))-$(subst /,-,$(2)) $(1)/$(2) $(SAMPN) > $$@

ALLDETECTPBS += detect-$(subst /,-,$(1))-$(subst /,-,$(2)).pbs

.PRECIOUS: $(BPATH)/$(1)/$(2)/%/acc.rds

$(BPATH)/$(1)/$(2)/%/pc.rds: spinglass-detect.R $(DATAPATH)/background/$(2)/acc $(DATAPATH)/background/$(2)/pc $(BPATH)/$(1)/$(2)/%/acc
	mkdir -p $$(dir $$@)
	$(R) $$^ $(lastword $(subst /,$(SPACE),$(2))) > $$@

pc-$(subst /,-,$(1))-$(subst /,-,$(2)).pbs: pc-detect.sh
	./$$^ $(subst /,-,$(1))-$(subst /,-,$(2)) $(1)/$(2) $(SAMPN) > $$@

ALLFINALPCPBS += pc-$(subst /,-,$(1))-$(subst /,-,$(2)).pbs

endef

$(foreach d,$(COVERTDIMS),\
 $(foreach b,$(BG-BASE-FACTORIAL),\
$(eval $(call detectingbase,$(d),$(b)))\
))

$(foreach d,$(COVERTDIMS),\
 $(foreach b,$(BG-FACTORIAL),\
$(eval $(call detecting,$(d),$(b)))\
))

# ./input/detection/high/hi/late/20/15/15/%/base.rds: pre-spinglass-detect.R ./input/digest/raw/pairs.rds ./input/digest/raw/location-lifetimes.rds ./input/digest/background/15/15/base ./input/simulate/covert/high/hi/late/20/%/cc.csv ./input/simulate/covert/high/hi/late/20/%/cu.csv | ./input/detection/high/hi/late/20/15/15/%

alldetectbasepbs: $(ALLDETECTBASEPBS)
allsnapsreviewpbs: $(SNAPBASEPBS)
alldetectpbs: $(ALLDETECTPBS)
allfinalpbs: $(ALLFINALPCPBS)

.PHONY: submitsomebase

submitsomebase: alldetectbasepbs
	for f in $(wordlist $(s),$(e),$(ALLDETECTBASEPBS)); do sbatch $$f; done;

submitsomeacc: alldetectpbs
	for f in $(wordlist $(s),$(e),$(ALLDETECTPBS)); do sbatch $$f; done;

submitsomepc: allfinalpbs
	for f in $(wordlist $(s),$(e),$(ALLFINALPCPBS)); do sbatch $$f; done;


# $(foreach d,$(COVERTDIMS),\
#  $(foreach b,$(BG-FACTORIAL),\
# $(eval $(call detectingsecond,$(d),$(b)))\
# ))


#.SECONDEXPANSION:

# $(foreach d,high/hi/late/20,\
#  $(foreach b,15/15/censor,\
# $(info $(call detectingsecond,$(d),$(b)))\
# ))

#
# $(foreach d,$(COVERTDIMS),\
#  $(foreach b,$(BG-FACTORIAL),\
# $(info $(call detectingsecond,$(d),$(b)))\
# ))
