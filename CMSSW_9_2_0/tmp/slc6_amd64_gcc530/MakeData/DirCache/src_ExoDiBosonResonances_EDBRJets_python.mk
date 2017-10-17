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
