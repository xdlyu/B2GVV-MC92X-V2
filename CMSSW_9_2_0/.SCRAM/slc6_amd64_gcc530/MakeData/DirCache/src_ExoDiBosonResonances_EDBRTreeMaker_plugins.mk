ifeq ($(strip $(ExoDiBosonResonancesEDBRTreeMakerAuto)),)
ExoDiBosonResonancesEDBRTreeMakerAuto := self/src/ExoDiBosonResonances/EDBRTreeMaker/plugins
PLUGINS:=yes
ExoDiBosonResonancesEDBRTreeMakerAuto_files := $(patsubst src/ExoDiBosonResonances/EDBRTreeMaker/plugins/%,%,$(wildcard $(foreach dir,src/ExoDiBosonResonances/EDBRTreeMaker/plugins ,$(foreach ext,$(SRC_FILES_SUFFIXES),$(dir)/*.$(ext)))))
ExoDiBosonResonancesEDBRTreeMakerAuto_BuildFile    := $(WORKINGDIR)/cache/bf/src/ExoDiBosonResonances/EDBRTreeMaker/plugins/BuildFile
ExoDiBosonResonancesEDBRTreeMakerAuto_LOC_USE := self  boost FWCore/Framework FWCore/PluginManager FWCore/ParameterSet FWCore/Utilities CommonTools/UtilAlgos DataFormats/Common DataFormats/Candidate DataFormats/PatCandidates DataFormats/Math DataFormats/HepMCCandidate CondFormats/JetMETObjects JetMETCorrections/Objects HLTrigger/HLTcore PhysicsTools/PatUtils PhysicsTools/PatExamples PhysicsTools/PatAlgos SimDataFormats/JetMatching AnalysisDataFormats/TopObjects root
ExoDiBosonResonancesEDBRTreeMakerAuto_PRE_INIT_FUNC += $$(eval $$(call edmPlugin,ExoDiBosonResonancesEDBRTreeMakerAuto,ExoDiBosonResonancesEDBRTreeMakerAuto,$(SCRAMSTORENAME_LIB),src/ExoDiBosonResonances/EDBRTreeMaker/plugins))
ExoDiBosonResonancesEDBRTreeMakerAuto_PACKAGE := self/src/ExoDiBosonResonances/EDBRTreeMaker/plugins
ALL_PRODS += ExoDiBosonResonancesEDBRTreeMakerAuto
ExoDiBosonResonances/EDBRTreeMaker_forbigobj+=ExoDiBosonResonancesEDBRTreeMakerAuto
ExoDiBosonResonancesEDBRTreeMakerAuto_INIT_FUNC        += $$(eval $$(call Library,ExoDiBosonResonancesEDBRTreeMakerAuto,src/ExoDiBosonResonances/EDBRTreeMaker/plugins,src_ExoDiBosonResonances_EDBRTreeMaker_plugins,$(SCRAMSTORENAME_BIN),,$(SCRAMSTORENAME_LIB),$(SCRAMSTORENAME_LOGS)))
ExoDiBosonResonancesEDBRTreeMakerAuto_CLASS := LIBRARY
else
$(eval $(call MultipleWarningMsg,ExoDiBosonResonancesEDBRTreeMakerAuto,src/ExoDiBosonResonances/EDBRTreeMaker/plugins))
endif
ALL_COMMONRULES += src_ExoDiBosonResonances_EDBRTreeMaker_plugins
src_ExoDiBosonResonances_EDBRTreeMaker_plugins_parent := ExoDiBosonResonances/EDBRTreeMaker
src_ExoDiBosonResonances_EDBRTreeMaker_plugins_INIT_FUNC += $$(eval $$(call CommonProductRules,src_ExoDiBosonResonances_EDBRTreeMaker_plugins,src/ExoDiBosonResonances/EDBRTreeMaker/plugins,PLUGINS))
