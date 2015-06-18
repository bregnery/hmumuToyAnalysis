{
  gSystem->cd("src");
  gROOT->ProcessLine(".L FuncSmearingZmumu2011ChargePtCorr2.C+");
  gROOT->ProcessLine(".L FuncSmearingZmumu2012PtCorr2.C+");
  gSystem->cd("../");
  gROOT->ProcessLine(".L analyzer.C++");
  
  unsigned maxEvents = 1000000;
  TString runPeriod("13TeV");

  // Scale factors are the last input variable
  // scale factor = (cross section * SM branching ratio * Luminosity)/ Number of events generated
  //// I had to use 14 TeV cross sections. Please change to 13 TeV and recalculate
  //// if they become available.
  //analyzer("/home/bregnery/13TeV/CMSSW_7_0_9/src/UserArea/UFDiMuonsAnalyzer/stage1/ggHmumu13TeV.root","ggHmumu.root",runPeriod,false,true,maxEvents,0.0846512128);
  //analyzer("/home/bregnery/13TeV/CMSSW_7_0_9/src/UserArea/UFDiMuonsAnalyzer/stage1/VBFHmumu13TeV.root","vbfHmumu.root",runPeriod,false,true,maxEvents,0.006840064);
  analyzer("/home/bregnery/13TeV/CMSSW_7_0_9/src/UserArea/UFDiMuonsAnalyzer/stage1/VHmumu13TeV.root","wHmumu.root",runPeriod,false,true,maxEvents,0.0027143994);

  analyzer("/home/bregnery/13TeV/CMSSW_7_0_9/src/UserArea/UFDiMuonsAnalyzer/stage1/DYToLL13TeV.root","DY.root",runPeriod,false,false,maxEvents,2.295515052);
  //analyzer("/cms/data/store/user/jhugon/hmumu/stage1/8TeV/TTJets_minimal.root","TT.root",runPeriod,false,false,maxEvents,0.6413329165);
  //analyzer("/cms/data/store/user/jhugon/hmumu/stage1/8TeV/WW_minimal.root","WW.root",runPeriod,false,false,maxEvents,0.2071593204);
  //analyzer("/cms/data/store/user/jhugon/hmumu/stage1/8TeV/WZ_minimal.root","WZ.root",runPeriod,false,false,maxEvents,0.0654616994);
  //analyzer("/cms/data/store/user/jhugon/hmumu/stage1/8TeV/ZZ_minimal.root","ZZ.root",runPeriod,false,false,maxEvents,0.0337205085);

//  analyzer("/cms/data/store/user/jhugon/hmumu/stage1/8TeV/SingleMuRun2012A-13Jan2013-v1_minimal.root","SingleMuRun2012A.root",runPeriod,true,false,maxEvents);
//  analyzer("/cms/data/store/user/jhugon/hmumu/stage1/8TeV/SingleMuRun2012B-13Jan2013-v1_minimal.root","SingleMuRun2012B.root",runPeriod,true,false,maxEvents);
//  analyzer("/cms/data/store/user/jhugon/hmumu/stage1/8TeV/SingleMuRun2012C-13Jan2013-v1_minimal.root","SingleMuRun2012C.root",runPeriod,true,false,maxEvents);
//  analyzer("/cms/data/store/user/jhugon/hmumu/stage1/8TeV/SingleMuRun2012D-13Jan2013-v1_minimal.root","SingleMuRun2012D.root",runPeriod,true,false,maxEvents);
  
}
