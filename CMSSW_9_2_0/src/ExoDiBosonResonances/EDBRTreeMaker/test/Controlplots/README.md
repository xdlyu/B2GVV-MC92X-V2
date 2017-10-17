(也许要配合higgscombinetools使用:https://twiki.cern.ch/twiki/pub/CMS/PKURunIISamples15Fall/Notelimit-7X.txt)

https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit

export SCRAM_ARCH=slc6_amd64_gcc481
cmsrel CMSSW_7_1_15 ### must be a 7_1_X release  >= 7_1_15;  (7.0.X and 7.2.X are NOT supported either) 
cd CMSSW_7_1_15/src 
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v5.0.1
scramv1 b clean
scramv1 b 
--------------------------------------------------
cd ../../src
cmsenv
. /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
. /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.34/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh
(需要特殊版本的ROOT)

--------------------------------------------------
You need compile some c++ files first: (for some unknown reason, you have to compile them one by one)

cd ./PDFs
root
root [0] gSystem->AddIncludePath("-I$ROOFITSYS/include")
root [1] .L PdfDiagonalizer.cc++
root [2] .q

root
root [0] gSystem->AddIncludePath("-I$ROOFITSYS/include")
root [1] .L RooRelBWRunningWidth.cxx++
root [2] .q

root
root [0] gSystem->AddIncludePath("-I$ROOFITSYS/include")
root [1] .L HWWLVJRooPdfs.cxx++
root [2] .q

root
root [0] gSystem->AddIncludePath("-I$ROOFITSYS/include")
root [1] .L Util.cxx++
root [2] .q


python g1_exo_doFit_class_pad_ratio_poisson_final_pdf_unc_btag_approval_ele_Jordan_puppi_full2016_lowsb.py --control -c el
python g1_exo_doFit_class_pad_ratio_poisson_final_pdf_unc_btag_approval_muHLT_Jordan_puppi_full2016_lowsb.py --control -c mu


Full B~H Lumi:36.42/fb
