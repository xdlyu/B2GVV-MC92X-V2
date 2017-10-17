ifeq ($(strip $(ExoDiBosonResonancesEDBRWLeptonicProducerAuto)),)
ExoDiBosonResonancesEDBRWLeptonicProducerAuto := self/src/ExoDiBosonResonances/EDBRWLeptonicProducer/plugins
PLUGINS:=yes
ExoDiBosonResonancesEDBRWLeptonicProducerAuto_files := $(patsubst src/ExoDiBosonResonances/EDBRWLeptonicProducer/plugins/%,%,$(wildcard $(foreach dir,src/ExoDiBosonResonances/EDBRWLeptonicProducer/plugins ,$(foreach ext,$(SRC_FILES_SUFFIXES),$(dir)/*.$(ext)))))
ExoDiBosonResonancesEDBRWLeptonicProducerAuto_BuildFile    := $(WORKINGDIR)/cache/bf/src/ExoDiBosonResonances/EDBRWLeptonicProducer/plugins/BuildFile
ExoDiBosonResonancesEDBRWLeptonicProducerAuto_LOC_USE := self  boost FWCore/Framework FWCore/PluginManager FWCore/ParameterSet DataFormats/Candidate CommonTools/CandUtils CommonTools/Utils DataFormats/PatCandidates CondFormats/JetMETObjects JetMETCorrections/Objects SimDataFormats/JetMatching DataFormats/Common DataFormats/Math DataFormats/HepMCCandidate PhysicsTools/PatUtils PhysicsTools/PatExamples AnalysisDataFormats/TopObjects root
ExoDiBosonResonancesEDBRWLeptonicProducerAuto_PRE_INIT_FUNC += $$(eval $$(call edmPlugin,ExoDiBosonResonancesEDBRWLeptonicProducerAuto,ExoDiBosonResonancesEDBRWLeptonicProducerAuto,$(SCRAMSTORENAME_LIB),src/ExoDiBosonResonances/EDBRWLeptonicProducer/plugins))
ExoDiBosonResonancesEDBRWLeptonicProducerAuto_PACKAGE := self/src/ExoDiBosonResonances/EDBRWLeptonicProducer/plugins
ALL_PRODS += ExoDiBosonResonancesEDBRWLeptonicProducerAuto
ExoDiBosonResonances/EDBRWLeptonicProducer_forbigobj+=ExoDiBosonResonancesEDBRWLeptonicProducerAuto
ExoDiBosonResonancesEDBRWLeptonicProducerAuto_INIT_FUNC        += $$(eval $$(call Library,ExoDiBosonResonancesEDBRWLeptonicProducerAuto,src/ExoDiBosonResonances/EDBRWLeptonicProducer/plugins,src_ExoDiBosonResonances_EDBRWLeptonicProducer_plugins,$(SCRAMSTORENAME_BIN),,$(SCRAMSTORENAME_LIB),$(SCRAMSTORENAME_LOGS)))
ExoDiBosonResonancesEDBRWLeptonicProducerAuto_CLASS := LIBRARY
else
$(eval $(call MultipleWarningMsg,ExoDiBosonResonancesEDBRWLeptonicProducerAuto,src/ExoDiBosonResonances/EDBRWLeptonicProducer/plugins))
endif
ALL_COMMONRULES += src_ExoDiBosonResonances_EDBRWLeptonicProducer_plugins
src_ExoDiBosonResonances_EDBRWLeptonicProducer_plugins_parent := ExoDiBosonResonances/EDBRWLeptonicProducer
src_ExoDiBosonResonances_EDBRWLeptonicProducer_plugins_INIT_FUNC += $$(eval $$(call CommonProductRules,src_ExoDiBosonResonances_EDBRWLeptonicProducer_plugins,src/ExoDiBosonResonances/EDBRWLeptonicProducer/plugins,PLUGINS))
