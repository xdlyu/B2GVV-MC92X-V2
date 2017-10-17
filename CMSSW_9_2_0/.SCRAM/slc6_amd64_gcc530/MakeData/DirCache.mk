ifeq ($(strip $(ExoDiBosonResonancesEDBRCommonAuto)),)
ExoDiBosonResonancesEDBRCommonAuto := self/src/ExoDiBosonResonances/EDBRCommon/plugins
PLUGINS:=yes
ExoDiBosonResonancesEDBRCommonAuto_files := $(patsubst src/ExoDiBosonResonances/EDBRCommon/plugins/%,%,$(wildcard $(foreach dir,src/ExoDiBosonResonances/EDBRCommon/plugins ,$(foreach ext,$(SRC_FILES_SUFFIXES),$(dir)/*.$(ext)))))
ExoDiBosonResonancesEDBRCommonAuto_BuildFile    := $(WORKINGDIR)/cache/bf/src/ExoDiBosonResonances/EDBRCommon/plugins/BuildFile
ExoDiBosonResonancesEDBRCommonAuto_LOC_USE := self  FWCore/Framework FWCore/PluginManager FWCore/ParameterSet FWCore/Utilities CommonTools/UtilAlgos DataFormats/Common DataFormats/MuonReco DataFormats/EgammaCandidates DataFormats/PatCandidates RecoEgamma/EgammaTools RecoEcal/EgammaCoreTools DataFormats/Candidate DataFormats/Math root
ExoDiBosonResonancesEDBRCommonAuto_PRE_INIT_FUNC += $$(eval $$(call edmPlugin,ExoDiBosonResonancesEDBRCommonAuto,ExoDiBosonResonancesEDBRCommonAuto,$(SCRAMSTORENAME_LIB),src/ExoDiBosonResonances/EDBRCommon/plugins))
ExoDiBosonResonancesEDBRCommonAuto_PACKAGE := self/src/ExoDiBosonResonances/EDBRCommon/plugins
ALL_PRODS += ExoDiBosonResonancesEDBRCommonAuto
ExoDiBosonResonances/EDBRCommon_forbigobj+=ExoDiBosonResonancesEDBRCommonAuto
ExoDiBosonResonancesEDBRCommonAuto_INIT_FUNC        += $$(eval $$(call Library,ExoDiBosonResonancesEDBRCommonAuto,src/ExoDiBosonResonances/EDBRCommon/plugins,src_ExoDiBosonResonances_EDBRCommon_plugins,$(SCRAMSTORENAME_BIN),,$(SCRAMSTORENAME_LIB),$(SCRAMSTORENAME_LOGS)))
ExoDiBosonResonancesEDBRCommonAuto_CLASS := LIBRARY
else
$(eval $(call MultipleWarningMsg,ExoDiBosonResonancesEDBRCommonAuto,src/ExoDiBosonResonances/EDBRCommon/plugins))
endif
ALL_COMMONRULES += src_ExoDiBosonResonances_EDBRCommon_plugins
src_ExoDiBosonResonances_EDBRCommon_plugins_parent := ExoDiBosonResonances/EDBRCommon
src_ExoDiBosonResonances_EDBRCommon_plugins_INIT_FUNC += $$(eval $$(call CommonProductRules,src_ExoDiBosonResonances_EDBRCommon_plugins,src/ExoDiBosonResonances/EDBRCommon/plugins,PLUGINS))
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
