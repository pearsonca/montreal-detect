
default:
	@echo helloworld

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

# # loop over covert dims, then analysis dims, then sample N
define detecting # 1 dir for covert, 2 is dir for detection
$(foreach i,$(SAMPS),$(call factorial2dir,$(1)/$(2)/$(i)))

# % = sample
$(BPATH)/$(1)/$(2)/%/base.rds: pre-spinglass-detect.R\
$(DATAPATH)/raw/pairs.rds $(DATAPATH)/raw/location-lifetimes.rds\
$(DATAPATH)/background/$(2)/base $(OUTSRC)/$(1)/%/cc.csv $(OUTSRC)/$(1)/%/cu.csv\
| $(BPATH)/$(1)/$(2)/%/
	@echo do something

# % = sample
$(BPATH)/$(1)/$(2)/%/acc.rds: pre-spinglass-score.R $(DATAPATH)/background/$(2)/base $(BPATH)/$(1)/$(2)/%/base.rds
	@echo do something

# need to get sample number in here somehow, but shouldn't be an issue
detect-$(subst /,-,$(1))-$(subst /,-,$(2)).pbs:
	@echo do something

ALLDETECTPBS += detect-$(1)-$(2).pbs

endef

define detectingsecond # % = sample/increment; 1 = covert dims, 2 = analysis dims
$(BPATH)/$(1)/$(2)/%/pc.rds: spinglass-detect.R\
$(DATAPATH)/background/$(2)/pc/$$$$(lastword $$$$(subst /,$(SPACE),$$$$*)).rds\
$(BPATH)/$(1)/$(2)/$$$$(firstword $$$$(subst /,$(SPACE),$$$$*))/acc.rds | $(BPATH)/$(1)/$(2)/$$$$(firstword $$$$(subst /,$(SPACE),$$$$*))
	mkdir -p $$(dir $$@)
	@echo do something

endef

$(foreach d,$(COVERTDIMS),\
 $(foreach b,$(BG-FACTORIAL),\
$(eval $(call detecting,$(d),$(b)))\
))

alldetectpbs: $(ALLDETECTPBS)

.SECONDEXPANSION:

$(foreach d,$(COVERTDIMS),\
 $(foreach b,$(BG-FACTORIAL),\
$(eval $(call detectingsecond,$(d),$(b)))\
))

$(foreach d,$(COVERTDIMS),\
 $(foreach b,$(BG-FACTORIAL),\
$(info $(call detectingsecond,$(d),$(b)))\
))