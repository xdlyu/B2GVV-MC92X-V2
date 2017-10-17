ALL_SUBSYSTEMS+=ExoDiBosonResonances
subdirs_src_ExoDiBosonResonances = src_ExoDiBosonResonances_EDBRCommon src_ExoDiBosonResonances_EDBRJets src_ExoDiBosonResonances_EDBRTreeMaker src_ExoDiBosonResonances_EDBRWLeptonicProducer
ALL_PACKAGES += ExoDiBosonResonances/EDBRCommon
subdirs_src_ExoDiBosonResonances_EDBRCommon := src_ExoDiBosonResonances_EDBRCommon_plugins src_ExoDiBosonResonances_EDBRCommon_python
ifeq ($(strip $(PyExoDiBosonResonancesEDBRCommon)),)
PyExoDiBosonResonancesEDBRCommon := self/src/ExoDiBosonResonances/EDBRCommon/python
src_ExoDiBosonResonances_EDBRCommon_python_parent := 
ALL_PYTHON_DIRS += $(patsubst src/%,%,src/ExoDiBosonResonances/EDBRCommon/python)
PyExoDiBosonResonancesEDBRCommon_files := $(patsubst src/ExoDiBosonResonances/EDBRCommon/python/%,%,$(wildcard $(foreach dir,src/ExoDiBosonResonances/EDBRCommon/python ,$(foreach ext,$(SRC_FILES_SUFFIXES),$(dir)/*.$(ext)))))
PyExoDiBosonResonancesEDBRCommon_LOC_USE := self  
PyExoDiBosonResonancesEDBRCommon_PACKAGE := self/src/ExoDiBosonResonances/EDBRCommon/python
ALL_PRODS += PyExoDiBosonResonancesEDBRCommon
PyExoDiBosonResonancesEDBRCommon_INIT_FUNC        += $$(eval $$(call PythonProduct,PyExoDiBosonResonancesEDBRCommon,src/ExoDiBosonResonances/EDBRCommon/python,src_ExoDiBosonResonances_EDBRCommon_python,1,1,$(SCRAMSTORENAME_PYTHON),$(SCRAMSTORENAME_LIB),,))
else
$(eval $(call MultipleWarningMsg,PyExoDiBosonResonancesEDBRCommon,src/ExoDiBosonResonances/EDBRCommon/python))
endif
ALL_COMMONRULES += src_ExoDiBosonResonances_EDBRCommon_python
src_ExoDiBosonResonances_EDBRCommon_python_INIT_FUNC += $$(eval $$(call CommonProductRules,src_ExoDiBosonResonances_EDBRCommon_python,src/ExoDiBosonResonances/EDBRCommon/python,PYTHON))
ALL_PACKAGES += ExoDiBosonResonances/EDBRJets
subdirs_src_ExoDiBosonResonances_EDBRJets := src_ExoDiBosonResonances_EDBRJets_python
ifeq ($(strip $(PyExoDiBosonResonancesEDBRJets)),)
PyExoDiBosonResonancesEDBRJets := self/src/ExoDiBosonResonances/EDBRJets/python
src_ExoDiBosonResonances_EDBRJets_python_parent := 
ALL_PYTHON_DIRS += $(patsubst src/%,%,src/ExoDiBosonResonances/EDBRJets/python)
PyExoDiBosonResonancesEDBRJets_files := $(patsubst src/ExoDiBosonResonances/EDBRJets/python/%,%,$(wildcard $(foreach dir,src/ExoDiBosonResonances/EDBRJets/python ,$(foreach ext,$(SRC_FILES_SUFFIXES),$(dir)/*.$(ext)))))
PyExoDiBosonResonancesEDBRJets_LOC_USE := self  
PyExoDiBosonResonancesEDBRJets_PACKAGE := self/src/ExoDiBosonResonances/EDBRJets/python
ALL_PRODS += PyExoDiBosonResonancesEDBRJets
PyExoDiBosonResonancesEDBRJets_INIT_FUNC        += $$(eval $$(call PythonProduct,PyExoDiBosonResonancesEDBRJets,src/ExoDiBosonResonances/EDBRJets/python,src_ExoDiBosonResonances_EDBRJets_python,1,1,$(SCRAMSTORENAME_PYTHON),$(SCRAMSTORENAME_LIB),,))
else
$(eval $(call MultipleWarningMsg,PyExoDiBosonResonancesEDBRJets,src/ExoDiBosonResonances/EDBRJets/python))
endif
ALL_COMMONRULES += src_ExoDiBosonResonances_EDBRJets_python
src_ExoDiBosonResonances_EDBRJets_python_INIT_FUNC += $$(eval $$(call CommonProductRules,src_ExoDiBosonResonances_EDBRJets_python,src/ExoDiBosonResonances/EDBRJets/python,PYTHON))
ALL_PACKAGES += ExoDiBosonResonances/EDBRTreeMaker
subdirs_src_ExoDiBosonResonances_EDBRTreeMaker := src_ExoDiBosonResonances_EDBRTreeMaker_plugins src_ExoDiBosonResonances_EDBRTreeMaker_test
ALL_COMMONRULES += src_ExoDiBosonResonances_EDBRTreeMaker_test
src_ExoDiBosonResonances_EDBRTreeMaker_test_parent := ExoDiBosonResonances/EDBRTreeMaker
src_ExoDiBosonResonances_EDBRTreeMaker_test_INIT_FUNC += $$(eval $$(call CommonProductRules,src_ExoDiBosonResonances_EDBRTreeMaker_test,src/ExoDiBosonResonances/EDBRTreeMaker/test,TEST))
ALL_PACKAGES += ExoDiBosonResonances/EDBRWLeptonicProducer
subdirs_src_ExoDiBosonResonances_EDBRWLeptonicProducer := src_ExoDiBosonResonances_EDBRWLeptonicProducer_plugins
