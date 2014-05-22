
#include <TSystem.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TF1.h>
#include <TRandom3.h>

#include <vector>

#include "src/DataFormats.h"
#include "src/helpers.h"
#include "src/LumiReweightingStandAlone.h"
#include "src/SmearingTool.h"
#include "src/SmearingTool2011.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer1D.h"


// Calculate Pt error
Double_t sigma(Double_t Pt){
  Double_t sigma = Pt*(TMath::Sqrt(33.44/(Pt*Pt) + 0.0016));
	return sigma;
}

Double_t calcPtJet2(Double_t ptJet1, Double_t phiJet1, Double_t phiJet2, Double_t dimuonPt, Double_t dimuonPhi){
	Double_t pxJet2 = ptJet1*(TMath::Cos(phiJet1)) + dimuonPt*(TMath::Cos(dimuonPhi));
	Double_t pyJet2 = ptJet1*(TMath::Sin(phiJet1)) + dimuonPt*(TMath::Sin(dimuonPhi));
	Double_t ptJet2 = TMath::Sqrt(pxJet2*pxJet2 + pyJet2*pyJet2);
	return ptJet2;
}

/*
  CalcChi2 c(0.,0.,0.,...)
  c(4.567);
*/
class CalcChi2{
	private:
		const Double_t ptJet1_;
		const Double_t ptJet2_;
		const Double_t phiJet1_;
		const Double_t phiJet2_;
		const Double_t sigmaJet1_;
		const Double_t sigmaJet2_;
		const Double_t dimuonPt_;
		const Double_t dimuonPhi_;
	public:
		CalcChi2();
		CalcChi2(Double_t ptJet1, Double_t ptJet2, Double_t phiJet1, Double_t phiJet2, Double_t sigmaJet1, Double_t sigmaJet2, Double_t dimuonPt, Double_t dimuonPhi);
  Double_t operator()(double ptJet1Param);
};

CalcChi2::CalcChi2():ptJet1_(0), ptJet2_(0), phiJet1_(0), phiJet2_(0), sigmaJet1_(0), sigmaJet2_(0), dimuonPt_(0), dimuonPhi_(0){

}

CalcChi2::CalcChi2(Double_t ptJet1, Double_t ptJet2, Double_t phiJet1, Double_t phiJet2, Double_t sigmaJet1, Double_t sigmaJet2, Double_t dimuonPt, Double_t dimuonPhi):
         ptJet1_(ptJet1), ptJet2_(ptJet2), phiJet1_(phiJet1), phiJet2_(phiJet2), sigmaJet1_(sigmaJet1), sigmaJet2_(sigmaJet2), dimuonPt_(dimuonPt), dimuonPhi_(dimuonPhi){
}

Double_t CalcChi2::operator()(double ptJet1Param){
	Double_t ptJet1pred = ptJet1Param;
	Double_t ptJet2pred = calcPtJet2(ptJet1pred,phiJet1_,phiJet2_,dimuonPt_,dimuonPhi_);
	Double_t Chi2 = ((ptJet1_-ptJet1pred)*(ptJet1_-ptJet1pred))/sigmaJet1_ + ((ptJet2_-ptJet2pred)*(ptJet2_-ptJet2pred))/sigmaJet2_;
	return Chi2;
}

void analyzer (TString inputFileName,TString outputFileName, TString runPeriod, bool isData, bool isSignal, unsigned maxEvents, float SF)
{
  using namespace std;

  ///////////////////
  // Configuration

  float minMmm = 110;
  float maxMmm = 160;

  //gErrorIgnoreLevel = kError;
  const unsigned ISMEAR = 2;

  ///////////////////////////
  // Output Histograms

  setStyle();

  TH1F* dimuonMassHist = new TH1F("dimuonMass","",50,110,160);
  setHistTitles(dimuonMassHist,"M(#mu#mu) [GeV/c^{2}]","Events");
  dimuonMassHist->Sumw2();
  TH1F* nJetsHist = new TH1F("nJets","",10,0,10);
  setHistTitles(nJetsHist,"N_{Jets}","Events");
  nJetsHist->SetStats(1);
  nJetsHist->Sumw2();
  TH1F* recoPtHist = new TH1F("recoPt","",100,0,100); // adding a histogram for reco pt
  setHistTitles(recoPtHist,"pt","Events");
  recoPtHist->Sumw2();
  TH1F* recoEtaHist = new TH1F("recoEta","",100,0,2.5); // adding a histogram for reco eta
  setHistTitles(recoEtaHist,"Eta","Events");
  recoEtaHist->Sumw2();
  TH1F* nJetsEtHist = new TH1F("nJetsEtHist","",100,0,100);
  setHistTitles(nJetsEtHist,"et","Events");
  nJetsEtHist->Sumw2();
  TH1F* jetMass = new TH1F("jetMass","",50,0,200);// 2 jet mass after the cuts
  setHistTitles(jetMass,"M(2Jet) [GeV/c^{2}]","Events");
  jetMass->SetStats(1);
  jetMass->Sumw2();
  TH1F* jetEta = new TH1F("jetEta","",10,0,10); //histogram for jet eta
  setHistTitles(jetEta,"Eta","Events");
  jetEta->Sumw2();
  TH1F* H2JetPhi = new TH1F("MuMu2JetPhi","",10,-1,-0.9); // histogram for the angle between the Higgs and the jets created
  setHistTitles(H2JetPhi,"cos(Phi)","Events");
  H2JetPhi->SetStats(1);
  H2JetPhi->Sumw2();
  TH1F* DiMuonPt = new TH1F("DiMuonPt","",100,0,200);
  setHistTitles(DiMuonPt,"Pt(DiMuon) [GeV/c]","Events");
  DiMuonPt->SetStats(1);
  DiMuonPt->Sumw2();
  TH1F* PtMiss = new TH1F("PtMiss","",100,0,200);
  setHistTitles(PtMiss,"Pt(miss) [GeV/c]","Events");
  PtMiss->SetStats(1);
  PtMiss->Sumw2();
  TH1F* DijetRep = new TH1F("DijetRep","",14,-7,7);
  setHistTitles(DijetRep,"Dijet Rapidity","Events");
  DijetRep->SetStats(1);
  DijetRep->Sumw2();
  TH1F* DiMuRap = new TH1F("DiMuRap","",10,-5,5);
  setHistTitles(DiMuRap,"Dimuon Rapidity","Events");
  DiMuRap->SetStats(1);
  DiMuRap->Sumw2();
  TH1F* P1MET0histo = new TH1F("P1MET0histo","",50,0,300);
  setHistTitles(P1MET0histo,"Jet 1 Pt from MET=0 [GeV/c]","Events");
  P1MET0histo->SetStats(1);
  P1MET0histo->Sumw2();
  TH1F* Jet1P = new TH1F("Jet1P","",50,0,300);
  setHistTitles(Jet1P,"Jet 1 Pt measured [GeV/c]","Events");
  Jet1P->SetStats(1);
  Jet1P->Sumw2();
  TH1F* JetMass1MET = new TH1F("JetMass1MET","",50,0,200);
  setHistTitles(JetMass1MET,"2 Jet invariant mass with Jet 1 from MET=0 [GeV/c^2]","Events");
  JetMass1MET->SetStats(1);
  JetMass1MET->Sumw2();
  TH1F* P2MET0histo = new TH1F("P2MET0histo","",50,0,300);
  setHistTitles(P2MET0histo,"Jet 2 Pt from MET=0 [GeV/c]","Events");
  P2MET0histo->SetStats(1);
  P2MET0histo->Sumw2();
  TH1F* Jet2P = new TH1F("Jet2P","",50,0,300);
  setHistTitles(Jet2P,"Jet 2 Pt measured [GeV/c]","Events");
  Jet2P->SetStats(1);
  Jet2P->Sumw2();
  TH1F* JetMass2MET = new TH1F("JetMass2MET","",50,0,200);
  setHistTitles(JetMass2MET,"2 Jet invariant mass with Jet 2 from MET=0 [GeV/c^{2}]","Events");
  JetMass2MET->SetStats(1);
  JetMass2MET->Sumw2();
  TH1F* JetMass85MET = new TH1F("JetMass85MET","",50,0,200);
  setHistTitles(JetMass85MET,"2 Jet invariant mass with mass value closest to 85 [GeV/c^{2}]","Events");
  JetMass85MET->SetStats(1);
  JetMass85MET->Sumw2();
  TH1F* jjmmPtDiff = new TH1F("jjmmPtDiff","",50,-200,200);
  setHistTitles(jjmmPtDiff,"The Differece of Dimuon Pt and DijetPt [GeV/c]","Events");
  jjmmPtDiff->SetStats(1);
  jjmmPtDiff->Sumw2();
  TH2F* DiffDimuPt = new TH2F("DiffDimuPt","",50,-200,200,50,0,200);
  DiffDimuPt->SetStats(1);
  DiffDimuPt->Sumw2();
  TH1F* DiffMF = new TH1F("DiffMF","",50,-200,200);
  setHistTitles(DiffMF,"The Difference of Measured and fitted Pt [GeV/c]","Events");
  DiffMF->SetStats(1);
  DiffMF->Sumw2();
  TH1F* DiffMF2 = new TH1F("DiffMF2","",50,-200,200);
  setHistTitles(DiffMF2,"The Difference for Jet 2 in Pt [GeV/c]","Events");
  DiffMF2->SetStats(1);
  DiffMF2->Sumw2();
  TH1F* j1fitPt = new TH1F("j1fitPt","",50,0,200);
  setHistTitles(j1fitPt,"Jet 1 fitted Pt [GeV/c]","Events");
  j1fitPt->SetStats(1);
  j1fitPt->Sumw2();
  TH1F* j2fitPt = new TH1F("j2fitPt","",50,0,200);
  setHistTitles(j2fitPt,"Jet 2 fitted Pt [GeV/c]","Events");
  j2fitPt->SetStats(1);
  j2fitPt->Sumw2();
  TH1F* jjfitM = new TH1F("jjfitM","",50,0,200);
  setHistTitles(jjfitM,"Fitted Dijet Mass [GeV/c^{2}]","Events");
  jjfitM->SetStats(1);
  jjfitM->Sumw2();
  

  ///////////////////////////
  Double_t MASS_MUON = 0.105658367;    //GeV/c2

  //////////////////////////
  // Tree Branches
  cout << "Analyzing filename: "<< inputFileName.Data() << endl;
  if (isData)
    cout << "isData\n";
  if (isSignal)
    cout << "isSignal\n";

  TChain * tree = new TChain("tree");
  tree->Add(inputFileName);


  // These are the names of the muons (See src/DataFormats.h for definitions!)
  _MuonInfo reco1, reco2;

  tree->SetBranchAddress("reco1", &reco1);
  tree->SetBranchAddress("reco2", &reco2);

  // These are the dimuon mass, pt, rapidity, and phi
  float recoCandMass, recoCandPt, recoCandY, recoCandPhi;
  float recoCandMassRes, recoCandMassResCov;

  tree->SetBranchAddress("recoCandMass",       &recoCandMass);
  tree->SetBranchAddress("recoCandPt",         &recoCandPt);
  tree->SetBranchAddress("recoCandY",          &recoCandY);
  tree->SetBranchAddress("recoCandPhi",        &recoCandPhi);
  tree->SetBranchAddress("recoCandMassRes",    &recoCandMassRes);
  tree->SetBranchAddress("recoCandMassResCov", &recoCandMassResCov);

  // MC truth info
  float trueMass=-99999.0;
  if(!isData && tree->GetBranchStatus("trueMass"))
    tree->SetBranchAddress("trueMass", &trueMass);

  /// Higgs Boson MC truth info (after FSR)
  _genPartInfo genHpostFSR;
  if(!isData && tree->GetBranchStatus("genHpostFSR"))
    tree->SetBranchAddress("genHpostFSR", &genHpostFSR);

  _TrackInfo reco1GenPostFSR;
  if(!isData && tree->GetBranchStatus("genM1HpostFSR"))
    tree->SetBranchAddress("genM1HpostFSR", &reco1GenPostFSR);

  _TrackInfo reco2GenPostFSR;
  if(!isData && tree->GetBranchStatus("genM2HpostFSR"))
    tree->SetBranchAddress("genM2HpostFSR", &reco2GenPostFSR);

  /// the jet collection
  // these 'rawJets' already have Loose Jet ID applied, and JES corrections
  // and are cross-cleaned of tight muons
  // later, jets will have JER corrections, PUID, and basic cuts applied
  _PFJetInfo rawJets;
  tree->SetBranchAddress("pfJets",&rawJets);

  float puJetFullDisc[10];
  float puJetSimpleDisc[10];
  float puJetCutDisc[10];

  tree->SetBranchAddress("puJetFullDisc",&puJetFullDisc);
  tree->SetBranchAddress("puJetSimpleDisc",&puJetSimpleDisc);
  tree->SetBranchAddress("puJetCutDisc",&puJetCutDisc);

  float puJetFullId[10];
  float puJetSimpleId[10];
  float puJetCutId[10];

  tree->SetBranchAddress("puJetFullId",&puJetFullId);
  tree->SetBranchAddress("puJetSimpleId",&puJetSimpleId);
  tree->SetBranchAddress("puJetCutId",&puJetCutId);

  int nPU=0;
  if (!isData)
    {
      tree->SetBranchAddress("nPU",&nPU);
    }
  _VertexInfo vertexInfo;
  tree->SetBranchAddress("vertexInfo",&vertexInfo);
  _EventInfo eventInfo;
  tree->SetBranchAddress("eventInfo",&eventInfo);

  // Be careful, the met has not been well validated
  _MetInfo met;
  tree->SetBranchAddress("met",&met);

  //////////////////////////
  //for PU reweighting

  reweight::LumiReWeighting lumiWeights("pileupDists/PileUpHistMC2012Summer50nsPoissonOOTPU.root","pileupDists/PileUpHist2012ABCD.root","pileup","pileup");
  if (runPeriod == "7TeV")
  {
    cout << "Using 2011AB PU reweighting\n";
    lumiWeights = reweight::LumiReWeighting("pileupDists/PileUpHistMCFall11.root","pileupDists/PileUpHist2011AB.root","pileup","pileup");
  }
  else
  {
    cout << "Using 2012ABCD PU reweighting\n";
  }

  ///////////////////////////////
  // Which Muon Selection to Use

  bool (*muonIdFuncPtr)(_MuonInfo&);
  if (runPeriod == "7TeV")
    {
      cout << "Using 2011 Tight Muon Selection\n";
      muonIdFuncPtr = &isKinTight_2011_noIso;
    }
  else
    {
      cout << "Using 2012 Tight Muon Selection\n";
      muonIdFuncPtr = &isKinTight_2012_noIso;
    }

  /////////////////////////
  // Smearing
  SmearingTool *smearPT = new SmearingTool();
  SmearingTool2011 *smearPT2011 = new SmearingTool2011();

  /////////////////////////////
  /////////////////////////////

  unsigned nEvents = tree->GetEntries();
  unsigned reportEach=100000;
  if (nEvents/100000>reportEach)
    reportEach = nEvents/100000;

  ///////////////////////////////
  ///////////////////////////////
  ///////////////////////////////
  // Event Loop

  for(unsigned i=0; i<nEvents;i++)
  {
   
    if(i >= maxEvents)
      break;

    tree->GetEvent(i);
    if (i % reportEach == 0) cout << "Event: " << i << endl;

    // Reject events with invalid muons
    if (reco1.pt < 0. || reco2.pt < 0.)
        continue;

    /////////////////////////////////////////////////
    // Muon Resolution Smearing to match MuScleFit

    if(isSignal) // smear only signal because it has muons from higgs 
    {
      if(reco1GenPostFSR.pt<0.)
        cout << "Muon 1 Post FSR not valid!\n";
      if(reco2GenPostFSR.pt<0.)
        cout << "Muon 2 Post FSR not valid!\n";
      float ptReco1 = -1.;
      float ptReco2 = -1.;
      if(runPeriod == "7TeV")
      {
        ptReco1 = smearPT2011 -> PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge, reco1.pt, ISMEAR);
        ptReco2 = smearPT2011 -> PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge, reco2.pt, ISMEAR);
      }
      else
      {
        ptReco1 = smearPT -> PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge, reco1.pt, ISMEAR);
        ptReco2 = smearPT -> PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge, reco2.pt, ISMEAR);
      }
      TLorentzVector reco1Vec;
      TLorentzVector reco2Vec;
      reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,MASS_MUON);
      reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,MASS_MUON);
      TLorentzVector diMuonVec = reco1Vec + reco2Vec;

      reco1.pt = ptReco1;
      reco2.pt = ptReco2;
      recoCandMass = diMuonVec.M();
      recoCandPt = diMuonVec.Pt();
      recoCandY = diMuonVec.Rapidity();
      recoCandPhi = diMuonVec.Phi();
    
      reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,MASS_MUON);
      reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,MASS_MUON);
      diMuonVec = reco1Vec + reco2Vec;
    }

    //////////////////////////////////////////
    // Muon-related cuts

    if (recoCandMass > maxMmm || recoCandMass < minMmm)
        continue;

    bool muon1PassId = (*muonIdFuncPtr)(reco1);
    bool muon2PassId = (*muonIdFuncPtr)(reco2);
    if ( !(muon1PassId && muon2PassId))
        continue;

    bool muon1PassIso = (getPFRelIso(reco1) <= 0.12);
    bool muon2PassIso = (getPFRelIso(reco2) <= 0.12);
    if ( !(muon1PassIso && muon2PassIso))
        continue;
// Order muons by pt
    if (reco1.pt < reco2.pt)
    {
      _MuonInfo tmpMuon = reco1;
      reco1 = reco2;
      reco1 = tmpMuon;
    }

    // PU reweight
    float weight = SF;
    if (!isData)
    {
      weight *= lumiWeights.weight(nPU);
    }
    DiMuonPt->Fill(recoCandPt,weight);

    // Jet Part
    // Do basic selection on jets and JER corrections
    std::vector<TLorentzVector> jets;
    const float jetPtCut = 25.;
    const float jetAbsEtaCut = 2.7;
    const int jetPUIDCut = 4; // >=    tight = 7, medium = 6, loose = 4. Only loose is useful!!
    for(unsigned iJet=0; (iJet < unsigned(rawJets.nJets) && iJet < 10);iJet++)
    {
      // apply jet energy resolution corrections
      if (rawJets.genPt[iJet]>0.0 && rawJets.pt[iJet]>15.)
        rawJets.pt[iJet] = jerCorr(rawJets.pt[iJet],rawJets.genPt[iJet],rawJets.eta[iJet]); 
      bool goodPt = rawJets.pt[iJet]>jetPtCut;
      bool goodEta = fabs(rawJets.eta[iJet])<jetAbsEtaCut;
      bool goodPUID = puJetFullId[iJet] >= jetPUIDCut;
      if (goodPt && goodEta && goodPUID)
      {
        TLorentzVector tmpJetVec;
        tmpJetVec.SetPtEtaPhiM(rawJets.pt[iJet],rawJets.eta[iJet],rawJets.phi[iJet],rawJets.mass[iJet]);
        jets.push_back(tmpJetVec);
      }
    }

//    /////////////////////////////////////////////
//
//    cout << "Event: "<< eventInfo.run << ":" << eventInfo.event << endl;
//
//    // print muon-related info
//    cout << "recoCandMass: " << recoCandMass << endl;
//    cout << "muon Pt1: " << reco1.pt << endl;
//    cout << "muon Pt2: " << reco2.pt << endl;
//    cout << "muon eta1: " << reco1.eta << endl;
//    cout << "muon eta2: " << reco2.eta << endl;
//    cout << "muon iso1: " << getPFRelIso(reco1) << endl;
//    cout << "muon iso2: " << getPFRelIso(reco2) << endl;
//    cout << "PU weight: " << weight << endl;
//    cout << endl;
//
//    // print jet-related info
//    cout << "nJets: " << jets.size() << endl;
//    for (unsigned iJet=0; iJet < jets.size(); iJet++)
//    {
//      cout << "Jet "<<(iJet+1)<<": pt="<< jets[iJet].Pt() << " eta="<<jets[iJet].Eta()<<endl;
//    }
//    cout << endl;
//
//    /////////////////////////////////////////////
      // Fill Histograms

      dimuonMassHist->Fill(recoCandMass,weight);
      nJetsHist->Fill(jets.size(),weight);
      //DiMuonPt->Fill(recoCandPt,weight);
      DiMuRap->Fill(recoCandY,weight);
	  recoPtHist->Fill(reco1.pt,weight); // fill new hist
	  recoPtHist->Fill(reco2.pt,weight); // fill new hist
	  recoEtaHist->Fill(abs(reco1.eta),weight);
	  recoEtaHist->Fill(abs(reco2.eta),weight);
	  // Plot jet mass and jet energy after performing more selections
	  int njetsel = 0;
	  int index1 = 0;
	  int index2 = 0;
	  float jetM = 0;
	  for (unsigned iJet=0; iJet < jets.size(); iJet++){
		jetEta->Fill(jets[iJet].Eta(), weight);
	  }
	  for (unsigned iJet=0; iJet < jets.size(); iJet++){
		nJetsEtHist->Fill(jets[iJet].Et(),weight);
	  }
	  for (unsigned iJet=0; iJet < jets.size(); iJet++){
		if (jets[iJet].Pt() > 20){
			nJetsHist->Fill(jets.size(),weight);
			njetsel++;
			if (njetsel == 1){
				index1 = iJet;
			}
			else if (njetsel == 2){
				index2 = iJet;
			}
		}
	  }
	  if (njetsel >=2){ //&& njetsel <=3){
	    if(jets[index1].Pt()>=30 && jets[index2].Pt()>=30){
	      
            //Solve for 2 Jet Invariant Mass
			Double_t E = (jets[index1].E()+jets[index2].E());
			Double_t Px = (jets[index1].Px()+jets[index2].Px());
			Double_t Py = (jets[index1].Py()+jets[index2].Py());
			Double_t Pz = (jets[index1].Pz()+jets[index2].Pz());
			jetM = TMath::Sqrt((E*E)-(Px*Px)-(Py*Py)-(Pz*Pz));
			
			//Solve for Dijet Pt
			Double_t DijetPt = TMath::Sqrt((Px*Px)+(Py*Py));
			Double_t DiPtDiff = recoCandPt - DijetPt;
			jjmmPtDiff->Fill(DiPtDiff,weight);
			DiffDimuPt->Fill(DiPtDiff,recoCandPt,weight);

			//Solve for the cos of the angle between H and V
			Double_t DiMuPx = recoCandPt*TMath::Cos(recoCandPhi);
			Double_t DiMuPy = recoCandPt*TMath::Sin(recoCandPhi);
			Double_t Pt = TMath::Sqrt((Px*Px) + (Py*Py));
			Double_t cosPhi = ((Px*DiMuPx)+(Py*DiMuPy))/(TMath::Abs(Pt*recoCandPt));
			
			//Solve for Dijet Rapidity
			Double_t DiJetRepidity = TMath::Log((E + Pz)/(E-Pz));
			DijetRep->Fill(DiJetRepidity,weight);
			
			//Using MET=0 in order to find a better value for leading Jet P
			TLorentzVector muon1P;
			muon1P.SetPtEtaPhiM(reco1.pt,reco1.eta,reco1.phi,MASS_MUON);
			TLorentzVector muon2P;
			muon2P.SetPtEtaPhiM(reco2.pt,reco2.eta,reco2.phi,MASS_MUON);
			TLorentzVector mumuj2P = muon1P + muon2P + jets[index2];
			float j1PtMET0 = mumuj2P.Pt();
			TLorentzVector jet1PMET0 = jets[index1]*(j1PtMET0/jets[index1].Pt());
			P1MET0histo->Fill(jet1PMET0.Pt(),weight);
			Jet1P->Fill(jets[index1].Pt(),weight);
			
			//Using MET=0 jet 1 to calculate 2 jet invariant mass
			Double_t E1MET0 = (jet1PMET0.E()+jets[index2].E());
			Double_t Px1MET0 = (jet1PMET0.Px()+jets[index2].Px());
			Double_t Py1MET0 = (jet1PMET0.Py()+jets[index2].Py());
			Double_t Pz1MET0 = (jet1PMET0.Pz()+jets[index2].Pz());
			Double_t jetM1MET0 = TMath::Sqrt((E1MET0*E1MET0)-(Px1MET0*Px1MET0)-(Py1MET0*Py1MET0)-(Pz1MET0*Pz1MET0));
			JetMass1MET->Fill(jetM1MET0,weight);
			jetMass->Fill(jetM,weight);
			
			//Using MET=0 in order to find a better value for Sub-leading Jet P
			TLorentzVector mumuj1P = muon1P + muon2P + jets[index1];
			float j2PtMET0 = mumuj1P.Pt();
			TLorentzVector jet2PMET0 = jets[index2]*(j2PtMET0/jets[index2].Pt());
			P2MET0histo->Fill(jet2PMET0.Pt(),weight);
			Jet2P->Fill(jets[index2].Pt(),weight);
			
			//Using MET=0 jet 1 to calculate 2 jet invariant mass
			Double_t E2MET0 = (jet2PMET0.E()+jets[index1].E());
			Double_t Px2MET0 = (jet2PMET0.Px()+jets[index1].Px());
			Double_t Py2MET0 = (jet2PMET0.Py()+jets[index1].Py());
			Double_t Pz2MET0 = (jet2PMET0.Pz()+jets[index1].Pz());
			Double_t jetM2MET0 = TMath::Sqrt((E2MET0*E2MET0)-(Px2MET0*Px2MET0)-(Py2MET0*Py2MET0)-(Pz2MET0*Pz2MET0));
			JetMass2MET->Fill(jetM2MET0,weight);
			
			//Select MET=0 Jet Pt value that gives the closest value to 85
			if(jetM2MET0 >= jetM1MET0 && jetM1MET0 >= 85){
				JetMass85MET->Fill(jetM1MET0,weight);
			}
			else if(jetM1MET0 > jetM2MET0 && jetM2MET0 >= 85){
				JetMass85MET->Fill(jetM2MET0,weight);
			}
			else if(jetM2MET0 < jetM1MET0 && jetM1MET0 <= 85){
				JetMass85MET->Fill(jetM1MET0,weight);
			}
			else if(jetM1MET0 < jetM2MET0 && jetM2MET0 <= 85){
				JetMass85MET->Fill(jetM2MET0,weight);
			}
			
			//Select Jet Pt values using chi squared minimization 
			Double_t sigmaJet1 = sigma(jets[index1].Pt());
			Double_t sigmaJet2 = sigma(jets[index2].Pt());
			CalcChi2 c(jets[index1].Pt(),jets[index2].Pt(),jets[index1].Phi(),jets[index2].Phi(),sigmaJet1,sigmaJet2,recoCandPt,recoCandPhi);
			Double_t loBound = jets[index1].Pt() - 5*sigmaJet1;
			Double_t upBound = jets[index1].Pt() + 5*sigmaJet1;
			/*
			Double_t interval = (2*jets[index1].Pt()-20)/20;
			for(int j=0;j<20;j++){
			  float ptTempj = j*interval + 20;
			  cout << "    Jet Pt:" << ptTempj << "    " << c(ptTempj)<< endl;
			}
			*/
			ROOT::Math::Functor1D cFunctor(&c,"CalcChi2","operator()");
			ROOT::Math::GSLMinimizer1D minChi2(ROOT::Math::Minim1D::kGOLDENSECTION);
			if(loBound >= 10){
			  minChi2.SetFunction(cFunctor,jets[index1].Pt(),loBound,upBound);
			}
			else{
			  minChi2.SetFunction(cFunctor,jets[index1].Pt(),10,upBound);
			}
			minChi2.Minimize(100,0.01,0.01);
			/*
			cout << "Measured: " << jets[index1].Pt() << endl;
			cout << "Found minimum: x = " << minChi2.XMinimum() << " f(x) = " << minChi2.FValMinimum() << " Status: " << minChi2.Status() << endl;
			*/
			Double_t jet1FitPt = minChi2.XMinimum();
			Double_t jet2FitPt =calcPtJet2(minChi2.XMinimum(),jets[index1].Phi(),jets[index2].Phi(),recoCandPt,recoCandPhi);
			Double_t DiffMeasFitJ1 = jets[index1].Pt() - jet1FitPt;
			Double_t DiffMeasFitJ2 = jets[index2].Pt() - jet2FitPt;
			DiffMF->Fill(DiffMeasFitJ1,weight);
			DiffMF2->Fill(DiffMeasFitJ2,weight);
			j1fitPt->Fill(jet1FitPt,weight);
			j2fitPt->Fill(jet2FitPt,weight);

			//Calculate invariant mass using jet Fit values
			TLorentzVector j1fit = jets[index1]*(jet1FitPt/jets[index1].Pt());
			TLorentzVector j2fit = jets[index2]*(jet2FitPt/jets[index2].Pt());
			Double_t jjfitE = (j1fit.E() + j2fit.E());
			Double_t jjfitPx = (j1fit.Px() + j2fit.Px());
			Double_t jjfitPy = (j1fit.Py() + j2fit.Py());
			Double_t jjfitPz = (j1fit.Pz() + j2fit.Pz());
			Double_t jjfitMass = TMath::Sqrt(jjfitE*jjfitE - jjfitPx*jjfitPx - jjfitPy*jjfitPy - jjfitPz*jjfitPz);
			jjfitM->Fill(jjfitMass,weight);
			
			
			//cout << jets[index1].Pt() << " " << sigmaJet1 << endl;
			//CalcChi2 c(1.,1.,0.,0.,1.,1.,1.,1.);
			//Double_t test = c(1.);
			//cout << "Chi^{2} value: " << test << endl;
			//Double_t test2 = calcPtJet2(jets[index1].Pt(),jets[index1].Phi(),jets[index2].Phi(),recoCandPt,recoCandPhi);
			//cout << "Jet 2 Pt calc: " << test2 << endl;
			//cout << "Jet 2 Pt MET0: " << jet2PMET0.Pt() << endl;
			 
			if(jetM >= 60 && jetM <= 110){
			  //jetMass->Fill(jetM,weight);
			  //if(cosPhi <= -0.95){
			    //jetMass->Fill(jetM,weight);
			    H2JetPhi->Fill(cosPhi,weight);
			    //DiMuonPt->Fill(recoCandPt,weight);
			    //if(recoCandPt >= 110){
			      //jetMass->Fill(jetM,weight);
			      Double_t jetsPx = 0;
			      for(unsigned iJet=0; iJet<jets.size();iJet++){
				jetsPx = jetsPx + jets[iJet].Px();
			      }
			      Double_t jetsPy = 0;
			      for(unsigned iJet=0; iJet<jets.size();iJet++){
				jetsPy = jetsPy + jets[iJet].Py();
			      }
			      Double_t missPt = TMath::Sqrt((DiMuPx+jetsPx)*(DiMuPx+jetsPx)+(DiMuPy+jetsPy)*(DiMuPy+jetsPy));
			      //if(missPt <= 40){
				PtMiss->Fill(missPt,weight);
				//if(DiJetRepidity <= 3 && DiJetRepidity >=-3){
				//jetMass->Fill(jetM,weight);
				  //}
				  //}
				  //}
				  //}
			 }
		       }
	} 
  }
  float JetMassBin=0;
  float PtMissBin=0;
  float MuMuBin=0;

  TFile* outFile = new TFile(outputFileName,"RECREATE");
  outFile->cd();
  recoPtHist->Write();
  dimuonMassHist->Write();
  DiMuonPt->Write();
  float nbinsT = DiMuonPt->GetNbinsX();
  for(int i=0; i<= nbinsT + 1; i++){
    MuMuBin = MuMuBin + DiMuonPt->GetBinContent(i);
  }
  cout << "Events after 2 muon selection: " << MuMuBin << endl;
  nJetsHist->Write();
  recoEtaHist->Write();
  nJetsEtHist->Write();
  jetMass->Write();
  float nbins = jetMass->GetNbinsX();
  for(int i=0; i<= nbins + 1; i++){
    //cout << "Bin Number:" << jetMass->GetBin(i) << "Bin Content: " << jetMass->GetBinContent(i) << endl;
    JetMassBin = JetMassBin + jetMass->GetBinContent(i);
  }
  //cout << "Overflow: " << jetMass->GetBinContent(nbins+1) << endl;
  //cout << "Underflow: " << jetMass->GetBinContent(0) << endl;
  cout << "Events after 2 jet selection and 2 jet invariant mass cut: " << JetMassBin << endl;
  //cout << "integral: " << jetMass->Integral() << endl;
  jetEta->Write();
  H2JetPhi->Write();
  DijetRep->Write();
  PtMiss->Write();
  float nbinsM = PtMiss->GetNbinsX();
  for(int i=0; i<= nbinsM + 1; i++){
    PtMissBin = PtMissBin + PtMiss->GetBinContent(i);
  }
  cout << "Events After Cuts: " << PtMissBin << endl;
  //cout << "integral: " << PtMiss->Integral() << endl;
  DiMuRap->Write();
  P1MET0histo->Write();
  Jet1P->Write();
  JetMass1MET->Write();
  P2MET0histo->Write();
  Jet2P->Write();
  JetMass2MET->Write();
  JetMass85MET->Write();
  jjmmPtDiff->Write();
  DiffDimuPt->Write();
  DiffMF->Write();
  DiffMF2->Write();
  j1fitPt->Write();
  j2fitPt->Write();
  jjfitM->Write();
}
