{
  gSystem->cd("src");
  gROOT->ProcessLine(".L FuncSmearingZmumu2011ChargePtCorr2.C+");
  gROOT->ProcessLine(".L FuncSmearingZmumu2012PtCorr2.C+");
  gSystem->cd("../");
  gROOT->ProcessLine(".L analyzer.C++");
  
  unsigned maxEvents = 1000000;
  TString runPeriod("8TeV");
  analyzer("/raid/raid8/jhugon/higgsSamples/stage1/8TeV/GluGlu_HToMM_M-125.root","ggHmumu.root",runPeriod,false,true,maxEvents,0.0846512128);
  analyzer("/raid/raid8/jhugon/higgsSamples/stage1/8TeV/VBF_HToMM_M-125.root","vbfHmumu.root",runPeriod,false,true,maxEvents,0.006840064);
  analyzer("/raid/raid8/jhugon/higgsSamples/stage1/8TeV/WH_HToMuMu_M-125.root","wHmumu.root",runPeriod,false,true,maxEvents,0.0027143994);
  analyzer("/raid/raid8/jhugon/higgsSamples/stage1/8TeV/ZH_HToMuMu_M-125.root","zHmumu.root",runPeriod,false,true,maxEvents,0.0017098189);

  analyzer("/raid/raid8/jhugon/higgsSamples/stage1/8TeV/DYJetsToLL_minimal.root","DY.root",runPeriod,false,false,maxEvents,2.295515052);
  analyzer("/raid/raid8/jhugon/higgsSamples/stage1/8TeV/TTJets_minimal.root","TT.root",runPeriod,false,false,maxEvents,0.6413329165);
  analyzer("/raid/raid8/jhugon/higgsSamples/stage1/8TeV/WW_minimal.root","WW.root",runPeriod,false,false,maxEvents,0.2071593204);
  analyzer("/raid/raid8/jhugon/higgsSamples/stage1/8TeV/WZ_minimal.root","WZ.root",runPeriod,false,false,maxEvents,0.0654616994);
  analyzer("/raid/raid8/jhugon/higgsSamples/stage1/8TeV/ZZ_minimal.root","ZZ.root",runPeriod,false,false,maxEvents,0.0337205085);

//  analyzer("/raid/raid8/jhugon/higgsSamples/stage1/8TeV/SingleMuRun2012A-13Jan2013-v1_minimal.root","SingleMuRun2012A.root",runPeriod,true,false,maxEvents);
//  analyzer("/raid/raid8/jhugon/higgsSamples/stage1/8TeV/SingleMuRun2012B-13Jan2013-v1_minimal.root","SingleMuRun2012B.root",runPeriod,true,false,maxEvents);
//  analyzer("/raid/raid8/jhugon/higgsSamples/stage1/8TeV/SingleMuRun2012C-13Jan2013-v1_minimal.root","SingleMuRun2012C.root",runPeriod,true,false,maxEvents);
//  analyzer("/raid/raid8/jhugon/higgsSamples/stage1/8TeV/SingleMuRun2012D-13Jan2013-v1_minimal.root","SingleMuRun2012D.root",runPeriod,true,false,maxEvents);
  
}