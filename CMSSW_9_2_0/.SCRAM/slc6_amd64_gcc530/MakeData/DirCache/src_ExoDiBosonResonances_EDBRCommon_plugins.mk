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
