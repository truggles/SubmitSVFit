pushd $CMSSW_BASE/src


git clone git@github.com:veelken/SVFit_standalone.git TauAnalysis/SVfitStandalone
cd TauAnalysis/SVfitStandalone
git checkout dd7cf43e3f930040959f7d700cef976307d7cec3
pushd $CMSSW_BASE/src





#MVA MET
#git cms-addpkg RecoMET/METPUSubtraction
cd RecoMET/METPUSubtraction/
git cms-addpkg RecoMET/METPUSubtraction
git cms-addpkg DataFormats/METReco
git remote add -f mvamet https://github.com/rfriese/cmssw.git
git checkout MVAMET2_beta_0.6 -b mvamet

pushd $CMSSW_BASE/src

git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections 


git clone https://github.com/CMS-HTT/LeptonEff-interface.git HTT-utilities/LepEffInterface/
pushd HTT-utilities/LepEffInterface/
git clone https://github.com/CMS-HTT/LeptonEfficiencies.git data 
popd


pushd $CMSSW_BASE/src

#METSignificance
git cms-addpkg RecoMET/METProducers




