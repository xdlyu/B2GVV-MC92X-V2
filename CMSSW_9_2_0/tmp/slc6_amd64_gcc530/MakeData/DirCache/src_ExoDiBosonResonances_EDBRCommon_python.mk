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
