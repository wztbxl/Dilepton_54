#include <cstdio> 
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "sys/types.h"
#include "dirent.h"
#include "math.h"
#include "string.h"

#ifndef __CINT__  
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TTimer.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TUnixSystem.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "miniDst.h"
#include "cuts.h"
#include "RefMfun.h"
#include "StRefMultCorr.h"

using namespace std;
#endif

void bookHistograms();
void writeHistograms(char* outFile);
void makeTags();
void makeRealPairs();
void makeMixPairs();
void copyCurrentToBuffer();
Bool_t Init();
Bool_t passEvent(miniDst* event);
Bool_t passTrack(miniDst* event, Int_t i);
Int_t  getCentralityBin9(Int_t cenBin16);
Double_t calReWeight(Double_t refMultCorr);
Double_t calCosTheta(TLorentzVector eVec,TLorentzVector eeVec);
Double_t reCalEventPlane(miniDst* event, Bool_t rejElectron = kFALSE);
Double_t phiVAngle(TLorentzVector e1, TLorentzVector e2, Int_t q1, Int_t q2);

TTimer   *timer;
TRandom3 *myRandom;

//variables 
Int_t dayIndex;
Int_t runIndex;
Int_t mCentrality;
map<Int_t,Int_t> mTotalDayId;
map<Int_t,Int_t> mTotalRunId;
map<Int_t,Int_t> mBadRunId_001;
map<Int_t,Int_t> mBadRunId_021;

Float_t bField;
Float_t reWeight;
Int_t iran = 0;

const Float_t pairPtCut = 0.1;

//StRefMultCorr* refMultCorrUtil;
//default categories for mixevent
//mCenBins=16; mVzBins=10; mEveBins=24; mMaxEventsInBuffer=100;
const Int_t mCenBins = 9; //16; //9;
const Int_t mVzBins = 10; //10; //6;
const Int_t mEveBins = 12; //24; //12;
const Int_t mMaxEventsInBuffer = 350; //100; //50;
const Int_t mMaxElectrons = 30;
const Float_t mPhiVCutMRange = 0.2;
Float_t current_EQx[mMaxElectrons],current_EQy[mMaxElectrons];
Int_t current_nE,current_nEPlus,current_nEMinus;
TLorentzVector current_ePlus[mMaxElectrons];
TLorentzVector current_eMinus[mMaxElectrons];
int current_ePlus_CellID[mMaxElectrons];
int current_eMinus_CellID[mMaxElectrons];
int current_ePlus_tags[mMaxElectrons];
int current_eMinus_tags[mMaxElectrons];
Int_t cenBufferPointer, vzBufferPointer, eveBufferPointer;
Int_t nEventsInBuffer[mCenBins][mVzBins][mEveBins];
Bool_t bufferFullFlag[mCenBins][mVzBins][mEveBins];
Int_t buffer_nEPlus[mCenBins][mVzBins][mEveBins][mMaxEventsInBuffer];
Int_t buffer_nEMinus[mCenBins][mVzBins][mEveBins][mMaxEventsInBuffer];
TLorentzVector buffer_ePlus[mCenBins][mVzBins][mEveBins][mMaxEventsInBuffer][mMaxElectrons];
TLorentzVector buffer_eMinus[mCenBins][mVzBins][mEveBins][mMaxEventsInBuffer][mMaxElectrons];
int buffer_ePlus_CellID[mCenBins][mVzBins][mEveBins][mMaxEventsInBuffer][mMaxElectrons];
int buffer_eMinus_CellID[mCenBins][mVzBins][mEveBins][mMaxEventsInBuffer][mMaxElectrons];

//***** constrain the bad dedx calibration geometry *****
TF1 *funPosHi;
TF1 *funPosLow;
TF1 *funNegHi;
TF1 *funNegLow;
Float_t par[4][4];
Float_t parErr[4][4];

//********* define function and histograms *********
TF1 *phiVcut;

//in Init function
TProfile2D *ShiftFactorcos[mArrayLength];
TProfile2D *ShiftFactorsin[mArrayLength];
TProfile2D *etapluszplusQx;
TProfile2D *etapluszminusQx;
TProfile2D *etaminuszplusQx;
TProfile2D *etaminuszminusQx;
TProfile2D *etapluszplusQy;
TProfile2D *etapluszminusQy;
TProfile2D *etaminuszplusQy;
TProfile2D *etaminuszminusQy;
//in passEvent function
TH1F *hCentrality9;
TH1F *hRefMult;
TH1F *hVertexZ;
TH1F *hVzDiff;
TH1F *hBField;
//eventPlane
TH1F *hRawEventPlane;
TH1F *hNewEventPlane;
TH1F *hFinalEventPlane;
TH1F *hReCenterEventPlane;
TH2F *hInclusiveEPhivsPt;
TH2F *hExclusiveEPhivsPt;
TH2F *hCut3EPhivsPt;
TH2F *hnEMinusvsEPlus;
//angleV 
TH2F *hULAngleVvsM;
TH2F *hLPosAngleVvsM;
TH2F *hLNegAngleVvsM;
TH2F *hMixULAngleVvsM;
TH2F *hMixLPosAngleVvsM;
TH2F *hMixLNegAngleVvsM;
//without phiV cut
// TH2F *hULMvsPtwophiV;
// TH2F *hLPosMvsPtwophiV;
// TH2F *hLNegMvsPtwophiV;
// TH2F *hMixULMvsPtwophiV;
// TH2F *hMixLPosMvsPtwophiV;
// TH2F *hMixLNegMvsPtwophiV;
//with phiV cut
TH2F *hULMvsPt;
TH2F *hLPosMvsPt;
TH2F *hLNegMvsPt;
TH2F *hMixULMvsPt;
TH2F *hMixLPosMvsPt;
TH2F *hMixLNegMvsPt;
TH2F *hElectronCenvsPt;
TH2F *hPositronCenvsPt;
//**********************
//add centrality dimension
TH3F *hULMvsPtCen;
TH3F *hLPosMvsPtCen;
TH3F *hLNegMvsPtCen;
TH3F *hULMvsPtCen_CutedbyPhiV;
TH3F *hLPosMvsPtCen_CutedbyPhiV;
TH3F *hLNegMvsPtCen_CutedbyPhiV;
TH3F *hMixULMvsPtCen;
TH3F *hMixLPosMvsPtCen;
TH3F *hMixLNegMvsPtCen;
TH3F *hMixLPosMvsPtCen_CutedbyPhiV;
TH3F *hMixLNegMvsPtCen_CutedbyPhiV;
TH3F *hULMvsPhiCen;
TH3F *hLPosMvsPhiCen;
TH3F *hLNegMvsPhiCen;

TH1D* hCellIDDiff;

// TH3F *hULCosThetavsMvsCen;
// TH3F *hLPosCosThetavsMvsCen;
// TH3F *hLNegCosThetavsMvsCen;
// TH3F *hMixULCosThetavsMvsCen;
// TH3F *hMixLPosCosThetavsMvsCen;
// TH3F *hMixLNegCosThetavsMvsCen;

//*** for QA purpose;
// TH3F *hULePtvsMeevsCen;
// TH3F *hLSePtvsMeevsCen;
// TH3F *hMixULePtvsMeevsCen;
// TH3F *hMixLSePtvsMeevsCen;

int main(int argc, char** argv)
{
	if(argc!=1&&argc!=3) return -1;

	TString inFile="test.list";
	char outFile[1024];
	sprintf(outFile,"test/test");
	if(argc==3){
		inFile = argv[1];
		sprintf(outFile,"%s",argv[2]);
	}

	//+---------------------------------+
	//| open files and add to the chain |
	//+---------------------------------+
	TChain *chain = new TChain("miniDst");

	Int_t ifile=0;
	char filename[512];
	ifstream *inputStream = new ifstream;
	inputStream->open(inFile.Data());
	if (!(inputStream)) {
		printf("can not open list file\n");
		return 0;
	}
	for(;inputStream->good();){
		inputStream->getline(filename,512);
		if(inputStream->good()) {
			TFile *ftmp = new TFile(filename);
			if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) {
				cout<<"something wrong"<<endl;
			} else {
				cout<<"read in "<<ifile<<"th file: "<<filename<<endl;
				chain->Add(filename);
				ifile++;
			}
			delete ftmp;
		}
	}
	delete inputStream;

	//intialization
	bookHistograms();

	if( !Init() ){
		cout<<"The initialization is failed !!!"<<endl;
		return 0;
	}else{
		timer = new TTimer();
		myRandom = new TRandom3();

		//+-------------------+
		//| initialize buffer |
		//+-------------------+
		memset(nEventsInBuffer,0,sizeof(nEventsInBuffer));
		memset(bufferFullFlag,0,sizeof(bufferFullFlag));
		memset(buffer_nEPlus,0,sizeof(buffer_nEPlus));
		memset(buffer_nEMinus,0,sizeof(buffer_nEMinus));
	}

	//+-------------+
	//| loop events |
	//+-------------+
	miniDst *event = new miniDst(chain);
	Int_t nEvts = chain->GetEntries();
	cout<<nEvts<<" events"<<endl;
	//refMultCorrUtil = new StRefMultCorr("refmult");
	for(int i=0;i<nEvts;i++){

		if(i%(nEvts/10)==0) cout << "begin " << i << "th entry...." << endl;
		event->GetEntry(i);
		Int_t runId = event->mRunId;

		map<Int_t,Int_t>::iterator iter = mTotalDayId.find((runId/1000)%1000);
		if(iter != mTotalDayId.end())
			dayIndex = iter->second;
		else{
			dayIndex = -1;
			cout<<"Can not find the dayId in the mTotalDayId list"<<endl;
			continue;
		}
    // cout << "begin " << i << "th entry...." << endl;

		if(dayIndex<0) continue;

		iter = mTotalRunId.find(runId);
		if(iter != mTotalRunId.end()){
			runIndex = iter->second;
		}
		else{
			cout<<"runNumber:"<<runId<<endl;
			cout<<"Can not find the runNumber in the mTotalRunId list"<<endl;
			continue;

		}

		if(i%1000==0){
			long long tmp = (long long)timer->GetAbsTime();
			UInt_t seed = tmp/myRandom->Rndm();
			myRandom->SetSeed(seed);
			//cout<<"random number:"<<myRandom->Uniform(-1,1)<<endl;
		}

		if(!passEvent(event)) continue; 

		current_nE=0;
		current_nEPlus=0;
		current_nEMinus=0;
		Int_t npTrks = event->mNTrks;
		for(int j=0;j<npTrks;j++) passTrack(event,j);
		hnEMinusvsEPlus->Fill(current_nEPlus,current_nEMinus);

		Double_t finalEventPlane = reCalEventPlane(event);
		if(finalEventPlane<0) continue;
		eveBufferPointer = (Int_t)(finalEventPlane/TMath::Pi()*mEveBins);
		////cout<<"eveBufferPointer:"<<eveBufferPointer<<endl;
		if(eveBufferPointer<0 || eveBufferPointer>=mEveBins) continue;

		for (int i = 0; i < mMaxElectrons; i++)
		{
			current_ePlus_tags[i] = 1;
			current_eMinus_tags[i] = 1;
		}

		makeTags();
		// cout << "after tags " << endl;
		makeRealPairs();
		// cout << "after real pairs " << endl;
		makeMixPairs();
		// cout << "after mixed pairs " << endl;
		copyCurrentToBuffer();
		// cout << "after copy ro buffer " << endl;
	}

	writeHistograms(outFile);
	delete chain;

	cout<<"end of program"<<endl;
	return 0;
}
//________________________________________________________________
Bool_t passEvent(miniDst* event)
{
	Int_t runId  = event->mRunId;
	Float_t vx = event->mVertexX;
	Float_t vy = event->mVertexY;
	Float_t vz = event->mVertexZ;
	Float_t vr = sqrt(vx*vx+vy*vy);
	Float_t vpdVz = event->mVpdVz;
	Float_t ZDCrate = event->mZDCRate;

	Float_t vzDiff = vz - vpdVz;
	Int_t refMult = event->mRefMult;
	int  nTrigs = event->mNTrigs;
	bool fireTrigger = kFALSE;
	bool RefMVzCorFlag = kFALSE;
	Bool_t is001Trigger = kFALSE;
	Bool_t is021Trigger = kFALSE;
	for(int i=0; i< nTrigs; i++){
		int trigId = event->mTrigId[i];
		if(trigId == mTrigId[0] || trigId == mTrigId[2]){
			fireTrigger = kTRUE;
		}
		if(trigId == mTrigId[0])RefMVzCorFlag = kTRUE, is001Trigger = kTRUE;
		if(trigId == mTrigId[2])is021Trigger = kTRUE;
	}
	if(!fireTrigger) return kFALSE;
	bField = event->mBField;
	mCentrality = event->mCentrality;

	map<Int_t, Int_t>::iterator iter_001 = mBadRunId_001.find(runId);
	if(iter_001 != mBadRunId_001.end() && is001Trigger){
		//cout<<"bad run, continue"<<endl;
		return kFALSE;
	}

	map<Int_t, Int_t>::iterator iter_021 = mBadRunId_021.find(runId);
	if(iter_021 != mBadRunId_021.end() && is021Trigger){
		//cout<<"bad run, continue"<<endl;
		return kFALSE;
	}

	// reWeight = 1.;
	Double_t RefMultCorr = refMult;
	mCentrality = event->mCentrality;
  mCentrality = mCentrality+1;
  cenBufferPointer = mCentrality;
  RefMultCorr = event->mGRefMultCorr;
  reWeight = event->mEvtWeight;

  // if(RefMVzCorFlag)RefMultCorr = GetRefMultCorr(refMult, vz);
	// reWeight = GetWeight(RefMultCorr);
	// mCentrality = GetCentrality(RefMultCorr);
	// cenBufferPointer = mCentrality-1;
	if (cenBufferPointer <0 || cenBufferPointer >8) return kFALSE;
	//cout << cenBufferPointer<<endl;
	
	
	// refMultCorrUtil->init(runId);
	// refMultCorrUtil->initEvent(refMult,vz,ZDCrate);
	// mCentrality = refMultCorrUtil->getCentralityBin9();
	// reWeight = refMultCorrUtil->getWeight();
	// double refMultCorr = refMultCorrUtil->getRefMultCorr();
	// cenBufferPointer = mCentrality;
	// if (cenBufferPointer <0 || cenBufferPointer >8) return kFALSE;

	//cout << cenBufferPointer<<endl;
	//if(refMult<300) cout<<"reWeight: "<<reWeight<<endl;

	if(TMath::Abs(vx)<1.e-5 && TMath::Abs(vy)<1.e-5 && TMath::Abs(vz)<1.e-5) return kFALSE;
	if(vr>=mVrCut) return kFALSE;
	if(TMath::Abs(vz)>=mVzCut) return kFALSE;//vz should also be in the range listed in the parameters file to do the refMult correction
	if(TMath::Abs(vzDiff)>=mVzDiffCut) return kFALSE;

	hBField->Fill(bField);
	hVertexZ->Fill(vz);
	hVzDiff->Fill(vzDiff);

	hRefMult->Fill(refMult,reWeight);

	Int_t centrality9 = mCentrality;
	hCentrality9->Fill(centrality9,reWeight);

	vzBufferPointer = (Int_t)((vz+mVzCut)/(2*mVzCut)*mVzBins);
	if(vzBufferPointer<0 || vzBufferPointer>=mVzBins) return kFALSE;

	return kTRUE;
}
//______________________________________________________________
Bool_t passTrack(miniDst* event, Int_t i)
{
	Int_t charge = event->mCharge[i];
	Int_t nHitsFit = event->mNHitsFit[i];
	Int_t nHitsDedx = event->mNHitsDedx[i];
	Int_t nHitsPoss = event->mNHitsPoss[i];
	Float_t nSigmaE = event->mNSigmaE[i];
	Float_t dca = event->mDca[i];
	Float_t pt = event->mPt[i];
	Float_t eta = event->mEta[i];
	Float_t phi = event->mPhi[i];
	Float_t beta2TOF = event->mBeta2TOF[i];
	Float_t ratio = 1.0*nHitsFit/nHitsPoss;
	int CellID = event->mTOFCellID[i];
	TVector3 mom;
	mom.SetPtEtaPhi(pt,eta,phi);
	Float_t p = mom.Mag();

	//if(charge<0)cout<<"charge="<<charge<<endl;
	if(pt<mTpcePtCut[0] || pt>mTpcePtCut[1]) return kFALSE;
	if(nHitsFit<20) return kFALSE;
	// if(nHitsFit<mTpceNHitsFitCut) return kFALSE;
	if(ratio<mTpceNHitsFitRatioCut) return kFALSE;
	// if(nHitsDedx<20) return kFALSE;
	if(nHitsDedx<mTpceNHitsDedxCut) return kFALSE;
	// if(dca>0.8) return kFALSE;
	if(dca>mTpceDcaCut) return kFALSE;
	if(TMath::Abs(eta)>mTpceEtaCut) return kFALSE;
	hInclusiveEPhivsPt->Fill(charge*pt,phi);
	if(beta2TOF<=0. || TMath::Abs(1.-1./beta2TOF)>mTpceBeta2TOFCut) return kFALSE;

	hExclusiveEPhivsPt->Fill(charge*pt,phi);
	Float_t mTpceNSigmaECutLow;
	if(p<.8){
		mTpceNSigmaECutLow = 3.0*p - 3.15; 
	}else{
		mTpceNSigmaECutLow = mTpceNSigmaECut[0];
	}
	if(nSigmaE<mTpceNSigmaECutLow+mNSigmaEShift || nSigmaE>mTpceNSigmaECut[1]+mNSigmaEShift) return kFALSE;
	hCut3EPhivsPt->Fill(charge*pt,phi);


	current_EQx[current_nE] = event->mEtaPlusQx;
	current_EQy[current_nE] = event->mEtaPlusQy;
	current_nE++;


	if(charge==1){
		current_ePlus[current_nEPlus].SetPtEtaPhiM(pt,eta,phi,Melectron);
		current_ePlus_CellID[current_nEPlus] = CellID;
		hPositronCenvsPt->Fill(pt,cenBufferPointer);
		current_nEPlus++;
	}
	else if(charge==-1){
		current_eMinus[current_nEMinus].SetPtEtaPhiM(pt,eta,phi,Melectron);
		current_eMinus_CellID[current_nEMinus] = CellID;
		hElectronCenvsPt->Fill(pt,cenBufferPointer);
		current_nEMinus++;
	}

	return kTRUE;
}

void makeTags()
{
	// cout << " nElectron = " << current_nEMinus << endl;
	// cout << " nPositron = " << current_nEPlus << endl;

	TLorentzVector pair(0,0,0,0);
	// e+e- and turn the electron under cuts tag to 0
	for (int i = 0; i < current_nEPlus; i++)
	{
		for ( int j = 0; j < current_nEMinus; j++)
		{
			pair = current_ePlus[i]+current_eMinus[j];
			// cout << "pair PseudoRapidity = " << pair.PseudoRapidity() << endl;
			if(TMath::Abs(pair.Rapidity())<=mPairYCut)
			{
				Double_t angleVcut = phiVcut->Eval(pair.M());
				Double_t angleV = phiVAngle(current_ePlus[i],current_eMinus[j],1,-1);
				// if( angleV < angleVcut && pair.M()<mPhiVCutMRange )
				// {
				// 	current_ePlus_tags[i] = 0;
				// 	current_eMinus_tags[j] = 0;
				// }
				if( pair.M() < 0.055 )
				{
					current_ePlus_tags[i] = 0;
					current_eMinus_tags[j] = 0;
				}
			}
			
		}
	}
	// for (int i = 0; i < current_nEPlus; i++)
	// {
	// 	for (int j = i+1; j < current_nEPlus; j++)
	// 	{
	// 		pair = current_ePlus[i]+current_ePlus[j];
	// 		// cout << "pair PseudoRapidity = " << pair.PseudoRapidity() << endl;
	// 		if(TMath::Abs(pair.Rapidity())<=mPairYCut)
	// 		{

	// 			Double_t angleVcut = phiVcut->Eval(pair.M());
	// 			Double_t angleV = phiVAngle(current_ePlus[i],current_ePlus[j],1,1);
	// 			// if( angleV < angleVcut && pair.M()<mPhiVCutMRange )
	// 			// {
	// 			// 	current_ePlus_tags[i] = 1;
	// 			// 	current_ePlus_tags[j] = 1;
	// 			// }
	// 			if( pair.M() < mMassCut )
	// 			{
	// 				current_ePlus_tags[i] = 0;
	// 				current_ePlus_tags[j] = 0;
	// 			}
	// 		}
	// 	}
	// }
	// for (int i = 0; i < current_nEMinus; i++)
	// {	
	// 	for (int j = i+1; j < current_nEMinus; j++)
	// 	{
	// 		pair = current_eMinus[i]+current_eMinus[j];
	// 		// cout << "pair PseudoRapidity = " << pair.PseudoRapidity() << endl;
	// 		if(TMath::Abs(pair.Rapidity())<=mPairYCut)
	// 		{
	// 			Double_t angleVcut = phiVcut->Eval(pair.M());
	// 			Double_t angleV = phiVAngle(current_eMinus[i],current_eMinus[j],-1,-1);
	// 			// if( angleV < angleVcut && pair.M()<mPhiVCutMRange )
	// 			// {
	// 			// 	current_eMinus_tags[i] = 1;
	// 			// 	current_eMinus_tags[j] = 1;
	// 			// }
	// 			if( pair.M() < mMassCut )
	// 			{
	// 				current_eMinus_tags[i] = 0;
	// 				current_eMinus_tags[j] = 0;
	// 			}
	// 		}
			
	// 	}
	// }


}

void makeRealPairs()
{
	//+--------------------------+
	//| current e+  + current e- |
	//+--------------------------+
	TLorentzVector pair(0,0,0,0);
	for(Int_t i=0;i<current_nEPlus;i++){
		if ( current_ePlus_tags[i] == 0) continue;
		for(Int_t j=0;j<current_nEMinus;j++){
			if ( current_eMinus_tags[j] == 0 ) continue;
			pair = current_ePlus[i]+current_eMinus[j];
			if(TMath::Abs(pair.Rapidity())<=mPairYCut){
				// hULMvsPtwophiV->Fill(pair.Pt(),pair.M(),reWeight);
				//hULMvsPtwophiV->Fill(pair.Pt(),pair.M());

				Double_t angleVcut = phiVcut->Eval(pair.M());
				Double_t angleV = phiVAngle(current_ePlus[i],current_eMinus[j],1,-1);
				hULAngleVvsM->Fill(pair.M(),angleV,reWeight);
				if( (angleV<angleVcut && pair.M()<mPhiVCutMRange) ) hULMvsPtCen_CutedbyPhiV->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
				if( (angleV>angleVcut && pair.M()<mPhiVCutMRange) || pair.M()>=mPhiVCutMRange ){
					hULMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
					//hULMvsPt->Fill(pair.Pt(),pair.M());
					hULMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
					hULMvsPhiCen->Fill(pair.Phi(),cenBufferPointer,pair.M(),reWeight);

					if(pair.Pt()<pairPtCut){
						Double_t costheta = calCosTheta(current_ePlus[i], pair);
						// hULMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
						//hULMvsPt->Fill(pair.Pt(),pair.M());
						// hULMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
						//hULCosThetavsMvsCen->Fill(cenBufferPointer, pair.M(), TMath::Abs(costheta), reWeight);

						// hULePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_ePlus[i].Pt(), reWeight);
						// hULePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_eMinus[j].Pt(), reWeight);
					}
				}
			}
		}//end of e- loop
	}//end of e+ loop

	//+--------------------------+
	//| current e+  + current e+ |
	//+--------------------------+
	for(Int_t i=0;i<current_nEPlus;i++){
		if ( current_ePlus_tags[i] == 0 ) continue;
		for(Int_t j=i+1;j<current_nEPlus;j++){
			if ( current_ePlus_tags[j] == 0 ) continue;
			pair = current_ePlus[i]+current_ePlus[j];
			if(TMath::Abs(pair.Rapidity())<=mPairYCut){
				// hLPosMvsPtwophiV->Fill(pair.Pt(),pair.M(),reWeight);
				//hLPosMvsPtwophiV->Fill(pair.Pt(),pair.M());

				// double TOF1x = cos(current_ePlus[i].Phi())*220;
				// double TOF1y = sin(current_ePlus[i].Phi())*220;
				// double TOF2x = cos(current_ePlus[j].Phi())*220;
				// double TOF2y = sin(current_ePlus[j].Phi())*220;
				// if ( sqrt( (TOF1x-TOF2x)*(TOF1x-TOF2x) + (TOF1y-TOF2y)*(TOF1y-TOF2y) ) < 6 ) continue;

				Double_t angleVcut = phiVcut->Eval(pair.M());
				Double_t angleV = phiVAngle(current_ePlus[i],current_ePlus[j],1,1);
				hLPosAngleVvsM->Fill(pair.M(),angleV,reWeight);
				if( (angleV<angleVcut && pair.M()<mPhiVCutMRange) ) hLPosMvsPtCen_CutedbyPhiV->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
				if( (angleV>angleVcut && pair.M()<mPhiVCutMRange) || pair.M()>=mPhiVCutMRange ){
				// if( (angleV>angleVcut && pair.M()<mPhiVCutMRange) || pair.M()>=mPhiVCutMRange ){
					hLPosMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
					//hLPosMvsPt->Fill(pair.Pt(),pair.M());
					hLPosMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
					hLPosMvsPhiCen->Fill(pair.Phi(),cenBufferPointer,pair.M(),reWeight);

					if(pair.Pt()<pairPtCut){
						Double_t costheta = calCosTheta(current_ePlus[i], pair);
						//hLPosCosThetavsMvsCen->Fill(cenBufferPointer, pair.M(), TMath::Abs(costheta), reWeight);

						// hLSePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_ePlus[i].Pt(), reWeight);
						// hLSePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_ePlus[j].Pt(), reWeight);
						// hLPosMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
						//hLPosMvsPt->Fill(pair.Pt(),pair.M());
						// hLPosMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
					}
				}
			}
		}//end of e+ loop
	}//end of e+ loop

	//+--------------------------+
	//| current e-  + current e- |
	//+--------------------------+
	for(Int_t i=0;i<current_nEMinus;i++){
		if ( current_eMinus_tags[i] == 0 ) continue;
		for(Int_t j=i+1;j<current_nEMinus;j++){
			if ( current_eMinus_tags[j] == 0 ) continue;
			pair = current_eMinus[i]+current_eMinus[j];
			if(TMath::Abs(pair.Rapidity())<=mPairYCut){
				// hLNegMvsPtwophiV->Fill(pair.Pt(),pair.M(),reWeight);
				//hLNegMvsPtwophiV->Fill(pair.Pt(),pair.M());

				// double TOF1x = cos(current_eMinus[i].Phi())*220;
				// double TOF1y = sin(current_eMinus[i].Phi())*220;
				// double TOF2x = cos(current_eMinus[j].Phi())*220;
				// double TOF2y = sin(current_eMinus[j].Phi())*220;
				// if ( sqrt( (TOF1x-TOF2x)*(TOF1x-TOF2x) + (TOF1y-TOF2y)*(TOF1y-TOF2y) ) < 6 ) continue;

				Double_t angleVcut = phiVcut->Eval(pair.M());
				Double_t angleV = phiVAngle(current_eMinus[i],current_eMinus[j],-1,-1);
				hLNegAngleVvsM->Fill(pair.M(),angleV,reWeight);
				if( (angleV<angleVcut && pair.M()<mPhiVCutMRange) ) hLNegMvsPtCen_CutedbyPhiV->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
				if( (angleV>angleVcut && pair.M()<mPhiVCutMRange) || pair.M()>=mPhiVCutMRange ){
					hLNegMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
					//hLNegMvsPt->Fill(pair.Pt(),pair.M());
					hLNegMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
					hLNegMvsPhiCen->Fill(pair.Phi(),cenBufferPointer,pair.M(),reWeight);

					if(pair.Pt()<pairPtCut){
						Double_t costheta = calCosTheta(current_eMinus[i], pair);
						//hLNegCosThetavsMvsCen->Fill(cenBufferPointer, pair.M(), TMath::Abs(costheta), reWeight);

						// hLSePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_eMinus[i].Pt(), reWeight);
						// hLSePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_eMinus[j].Pt(), reWeight);
						// hLNegMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
						//hLNegMvsPt->Fill(pair.Pt(),pair.M());
						// hLNegMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
					}
				}
			}
		}//end of e- loop
	}//end of e- loop
}
//_____________________________________________________________________________
void makeMixPairs()
{
	TLorentzVector pair(0,0,0,0);
	for(Int_t iBufferEvent=0;iBufferEvent<nEventsInBuffer[cenBufferPointer][vzBufferPointer][eveBufferPointer];iBufferEvent++){
		//+-------------------------+
		//| current e+  + buffer e- |
		//+-------------------------+
		for(Int_t i=0;i<current_nEPlus;i++){
			for(Int_t j=0;j<buffer_nEMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++){
				if ( current_ePlus_tags[i] == 0 ) continue;
				pair = current_ePlus[i] + buffer_eMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
				if ( abs( current_ePlus_CellID[i]-buffer_eMinus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] ) == 0 ) continue;
				// if ( abs( current_ePlus_CellID[i]-buffer_eMinus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] ) <= 1 ) continue;
				if(TMath::Abs(pair.Rapidity())<=mPairYCut){
					// hMixULMvsPtwophiV->Fill(pair.Pt(),pair.M(),reWeight);
					Double_t angleVcut = phiVcut->Eval(pair.M());
					Double_t angleV = phiVAngle(current_ePlus[i],buffer_eMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j],1,-1);
					hMixULAngleVvsM->Fill(pair.M(),angleV,reWeight);
					if( (angleV>angleVcut && pair.M()<mPhiVCutMRange) || pair.M()>=mPhiVCutMRange ){
						hMixULMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
						hMixULMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);

						if(pair.Pt()<pairPtCut){
							Double_t costheta = calCosTheta(current_ePlus[i], pair);
							//hMixULCosThetavsMvsCen->Fill(cenBufferPointer, pair.M(), TMath::Abs(costheta), reWeight);

							// hMixULePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_ePlus[i].Pt(), reWeight);
							// hMixULePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), buffer_eMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j].Pt(), reWeight);
							// hMixULMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
							// hMixULMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
						}
					}
				}
			}//end of buffer e- loop
		}//end of current e+ loop

		//+-------------------------+
		//| current e-  + buffer e+ |
		//+-------------------------+
		for(Int_t i=0;i<current_nEMinus;i++){
			for(Int_t j=0;j<buffer_nEPlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++){
				if ( current_eMinus_tags[i] == 0) continue;
				pair = current_eMinus[i] + buffer_ePlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];
				if ( abs( current_eMinus_CellID[i]-buffer_ePlus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] ) == 0 ) continue;
				// if ( abs( current_eMinus_CellID[i]-buffer_ePlus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] ) <= 1 ) continue;

				if(TMath::Abs(pair.Rapidity())<=mPairYCut){
					// hMixULMvsPtwophiV->Fill(pair.Pt(),pair.M(),reWeight);
					
					Double_t angleVcut = phiVcut->Eval(pair.M());
					Double_t angleV = phiVAngle(current_eMinus[i],buffer_ePlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j],-1,1);
					hMixULAngleVvsM->Fill(pair.M(),angleV,reWeight);
					if( (angleV>angleVcut && pair.M()<mPhiVCutMRange) || pair.M()>=mPhiVCutMRange ){
						hMixULMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
						hMixULMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);

						if(pair.Pt()<pairPtCut){
							Double_t costheta = calCosTheta(buffer_ePlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j], pair);
							//hMixULCosThetavsMvsCen->Fill(cenBufferPointer, pair.M(), TMath::Abs(costheta), reWeight);

							// hMixULePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_eMinus[i].Pt(), reWeight);
							// hMixULePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), buffer_ePlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j].Pt(), reWeight);
							// hMixULMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
							// hMixULMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
						}
					}
				}
			}//end of buffer e+ loop
		}//end of current e- loop

		//+-------------------------+
		//| current e+  + buffer e+ |
		//+-------------------------+
		for(Int_t i=0;i<current_nEPlus;i++){
			for(Int_t j=0;j<buffer_nEPlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++){
				if ( current_ePlus_tags[i] == 0 )  continue;
				pair = current_ePlus[i] + buffer_ePlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];

				if ( abs( current_ePlus_CellID[i]-buffer_ePlus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] ) == 0 ) continue;
				// if ( abs( current_ePlus_CellID[i]-buffer_ePlus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] ) <= 1 ) continue;
				hCellIDDiff->Fill(  current_ePlus_CellID[i]-buffer_ePlus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] );

				if(TMath::Abs(pair.Rapidity())<=mPairYCut){
					// hMixLPosMvsPtwophiV->Fill(pair.Pt(),pair.M(),reWeight);
					Double_t angleVcut = phiVcut->Eval(pair.M());
					Double_t angleV = phiVAngle(current_ePlus[i],buffer_ePlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j],1,1);
					hMixLPosAngleVvsM->Fill(pair.M(),angleV,reWeight);
					if( (angleV<angleVcut && pair.M()<mPhiVCutMRange) ) hMixLPosMvsPtCen_CutedbyPhiV->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
					if( (angleV>angleVcut && pair.M()<mPhiVCutMRange) || pair.M()>=mPhiVCutMRange ){
						hMixLPosMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
						hMixLPosMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);

						if(pair.Pt()<pairPtCut){
							Double_t costheta = calCosTheta(current_ePlus[i], pair);
							// hMixLPosCosThetavsMvsCen->Fill(cenBufferPointer, pair.M(), TMath::Abs(costheta), reWeight);

							// hMixLSePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_ePlus[i].Pt(), reWeight);
							// hMixLSePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), buffer_ePlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j].Pt(), reWeight);
							// hMixLPosMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
							// hMixLPosMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
						}

					}
				}
			}//end of buffer e+ loop
		}//endl of current e+ loop

		//+-------------------------+
		//| current e-  + buffer e- |
		//+-------------------------+
		for(Int_t i=0;i<current_nEMinus;i++){
			for(Int_t j=0;j<buffer_nEMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent];j++){
				if ( current_eMinus_tags[i] == 0 ) continue;
				pair = current_eMinus[i] + buffer_eMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j];

				if ( abs( current_eMinus_CellID[i]-buffer_eMinus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] ) <= 1 ) continue;
				hCellIDDiff->Fill(  current_eMinus_CellID[i]-buffer_eMinus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j] );

				if(TMath::Abs(pair.Rapidity())<=mPairYCut){
					// hMixLNegMvsPtwophiV->Fill(pair.Pt(),pair.M(),reWeight);
					Double_t angleVcut = phiVcut->Eval(pair.M());
					Double_t angleV = phiVAngle(current_eMinus[i],buffer_eMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j],-1,-1);
					hMixLNegAngleVvsM->Fill(pair.M(),angleV,reWeight);
					if( (angleV<angleVcut && pair.M()<mPhiVCutMRange) ) hMixLNegMvsPtCen_CutedbyPhiV->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
					if( (angleV>angleVcut && pair.M()<mPhiVCutMRange) || pair.M()>=mPhiVCutMRange ){
						hMixLNegMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
						hMixLNegMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);

						if(pair.Pt()<pairPtCut){
							Double_t costheta = calCosTheta(current_eMinus[i], pair);
							//hMixLNegCosThetavsMvsCen->Fill(cenBufferPointer, pair.M(), TMath::Abs(costheta), reWeight);

							// hMixLSePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), current_eMinus[i].Pt(), reWeight);
							// hMixLSePtvsMeevsCen->Fill(cenBufferPointer, pair.M(), buffer_eMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][iBufferEvent][j].Pt(), reWeight);
							// hMixLNegMvsPt->Fill(pair.Pt(),pair.M(),reWeight);
							// hMixLNegMvsPtCen->Fill(pair.Pt(),cenBufferPointer,pair.M(),reWeight);
						}
					}
				}
			}//end of buffer e- loop
		}//endl of current e- loop
	}
}
//_____________________________________________________________________________
void copyCurrentToBuffer()
{
	if(nEventsInBuffer[cenBufferPointer][vzBufferPointer][eveBufferPointer]>=mMaxEventsInBuffer) bufferFullFlag[cenBufferPointer][vzBufferPointer][eveBufferPointer] = kTRUE;
	Int_t eventPointer = -1;
	if(bufferFullFlag[cenBufferPointer][vzBufferPointer][eveBufferPointer]){
		eventPointer = (Int_t)myRandom->Uniform(0,mMaxEventsInBuffer-1.e-6);
	}else{
		eventPointer = nEventsInBuffer[cenBufferPointer][vzBufferPointer][eveBufferPointer];
	}

	buffer_nEPlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer] = current_nEPlus;
	buffer_nEPlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer] = current_nEPlus;
	int nTrks = 0;
	for(Int_t i=0;i<current_nEPlus; i++){
		if (current_ePlus_tags[i] == 0 ) continue;
		buffer_ePlus[cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_ePlus[i];
		buffer_ePlus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_ePlus_CellID[i];
		nTrks++;
	}

	buffer_nEMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer] = current_nEMinus;
	nTrks = 0;
	for(Int_t i=0;i<current_nEMinus;i++){
		if (current_eMinus_tags[i] == 0 ) continue;
		buffer_eMinus[cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_eMinus[i];
		buffer_eMinus_CellID[cenBufferPointer][vzBufferPointer][eveBufferPointer][eventPointer][i] = current_eMinus_CellID[i];
		nTrks++;
	}

	if(nEventsInBuffer[cenBufferPointer][vzBufferPointer][eveBufferPointer]<mMaxEventsInBuffer){
		nEventsInBuffer[cenBufferPointer][vzBufferPointer][eveBufferPointer]++;
	}
}
//_____________________________________________________________________________
Int_t getCentralityBin9(Int_t cenBin16)
{
	if(cenBin16<0 || cenBin16>15) return -1;
	else{
		if(cenBin16==15) return 8;
		else if(cenBin16==14) return 7;
		else return (Int_t)(0.5*cenBin16);
	}
}
//_____________________________________________________________________________
Double_t reCalEventPlane(miniDst* event, Bool_t rejElectron)
{
	Float_t vz = event->mVertexZ;
	Float_t mPlusQx = event->mEtaPlusQx;
	Float_t mPlusQy = event->mEtaPlusQy;
	Float_t mMinusQx = event->mEtaMinusQx;
	Float_t mMinusQy = event->mEtaMinusQy;
	Int_t mEtaPlusNTrks = event->mEtaPlusNTrks;
	Int_t mEtaMinusNTrks = event->mEtaMinusNTrks;
	Float_t Qx = mPlusQx + mMinusQx; 
	Float_t Qy = mMinusQy + mPlusQy;

	TVector2 mRawQ(Qx,Qy);
	Double_t rawEP = 0.5*mRawQ.Phi();
	if(rawEP<0.) rawEP += TMath::Pi();
	hRawEventPlane->Fill(rawEP);

	if(rejElectron){ //reject the contribution of electron
		for(Int_t i=0;i<current_nE;i++){
			Qx -= current_EQx[i];
			Qy -= current_EQy[i];
		}
	}
	Double_t eventPlane = -1;
	TVector2 Q(Qx,Qy);
	if((Q.Mod())!=0.){
		eventPlane = 0.5*Q.Phi();
		if(eventPlane<0.) eventPlane +=TMath::Pi();
	}
	hNewEventPlane->Fill(eventPlane);
	if(eventPlane<0.) return eventPlane;

	//********* get recenter number and recenter *********
	Double_t mReCenterQx, mReCenterQy;
	if(vz>0){
		mReCenterQx = Qx - mEtaPlusNTrks*etapluszplusQx->GetBinContent(runIndex+1, mCentrality+1) - mEtaMinusNTrks*etaminuszplusQx->GetBinContent(runIndex+1, mCentrality+1);
		mReCenterQy = Qy - mEtaPlusNTrks*etapluszplusQy->GetBinContent(runIndex+1, mCentrality+1) - mEtaMinusNTrks*etaminuszplusQy->GetBinContent(runIndex+1, mCentrality+1);
	}
	else{
		mReCenterQx = Qx - mEtaPlusNTrks*etapluszminusQx->GetBinContent(runIndex+1, mCentrality+1) - mEtaMinusNTrks*etaminuszminusQx->GetBinContent(runIndex+1, mCentrality+1);
		mReCenterQy = Qy - mEtaPlusNTrks*etapluszminusQy->GetBinContent(runIndex+1, mCentrality+1) - mEtaMinusNTrks*etaminuszminusQy->GetBinContent(runIndex+1, mCentrality+1);
	}
	Double_t recenterEP;
	TVector2 *mReCenterQ = new TVector2(mReCenterQx, mReCenterQy);
	if(mReCenterQ->Mod() > 0){
		recenterEP = 0.5*mReCenterQ->Phi();
		if(recenterEP<0.) recenterEP += TMath::Pi();
		hReCenterEventPlane->Fill(recenterEP);
	}

	//*********  get shift factor and add shift deltaPhi *********
	Float_t shiftCorrcos[mArrayLength];
	Float_t shiftCorrsin[mArrayLength];
	for(Int_t i=0;i<mArrayLength;i++){
		shiftCorrcos[i] = ShiftFactorcos[i]->GetBinContent(dayIndex+1,mCentrality+1);
		shiftCorrsin[i] = ShiftFactorsin[i]->GetBinContent(dayIndex+1,mCentrality+1);
	}
	Double_t deltaPhi=0;
	for(Int_t i=0;i<mArrayLength;i++){
		deltaPhi += 1./(i+1)*(-1.*shiftCorrsin[i]*cos(2.*(i+1)*recenterEP) + shiftCorrcos[i]*sin(2.*(i+1)*recenterEP));
	}
	if(deltaPhi<0.) deltaPhi += TMath::Pi();
	if(deltaPhi>=TMath::Pi()) deltaPhi -= TMath::Pi();
	recenterEP += deltaPhi;
	if(recenterEP<0.) recenterEP += TMath::Pi();
	if(recenterEP>=TMath::Pi()) recenterEP -= TMath::Pi();
	hFinalEventPlane->Fill(recenterEP);

	return recenterEP;
}
//____________________________________________________________
Double_t phiVAngle(TLorentzVector e1, TLorentzVector e2, Int_t q1, Int_t q2)
{
	Double_t pt1 = e1.Pt();
	Double_t eta1 = e1.Eta();
	Double_t phi1 = e1.Phi();

	Double_t pt2 = e2.Pt();
	Double_t eta2 = e2.Eta();
	Double_t phi2 = e2.Phi();

	TVector3 e1Mom,e2Mom;
	if(q1>0&&q2<0){
		e2Mom.SetPtEtaPhi(pt1,eta1,phi1);//e+
		e1Mom.SetPtEtaPhi(pt2,eta2,phi2);//e-
	}else if(q1<0&&q2>0){
		e2Mom.SetPtEtaPhi(pt2,eta2,phi2);//e+
		e1Mom.SetPtEtaPhi(pt1,eta1,phi1);//e-
	}else if(q1==q2&&TMath::Abs(q1)==1){
		Double_t ran = myRandom->Uniform(-1,1);
		if(ran>0){
			e2Mom.SetPtEtaPhi(pt1,eta1,phi1);
			e1Mom.SetPtEtaPhi(pt2,eta2,phi2);
		}
		else{
			e2Mom.SetPtEtaPhi(pt2,eta2,phi2);
			e1Mom.SetPtEtaPhi(pt1,eta1,phi1);
		}
	}else return -1;
	Double_t mN = 0.;
	if(bField<0.) mN = -1.;
	if(bField>0.) mN = 1.;

	TVector3 pu=e1Mom+e2Mom;
	TVector3 pv=e1Mom.Cross(e2Mom);
	TVector3 pw=pu.Cross(pv);
	TVector3 pnz(0.,0.,mN);
	TVector3 pwc=pu.Cross(pnz);

	Double_t angleV = pw.Angle(pwc);

	return angleV;
}
//____________________________________________________________
Double_t calCosTheta(TLorentzVector eVec,TLorentzVector eeVec)
{
	//eVec: positron TLorentzVector  eeVec: ee pair LorentzVector
	TLorentzVector positron(eVec); //positron
	TLorentzVector beam(0., 0., sqrt(pow(96.5,2)-pow(Mproton,2)), 96.5); // UU@193 GeV

	TVector3 dir = eeVec.BoostVector();

	positron.Boost(-1*dir);
	beam.Boost(-1*dir);

	Float_t theta = positron.Angle(beam.Vect());
	return TMath::Cos(theta);
}
//____________________________________________________________
void bookHistograms()
{
	char buf[500];
	for(int i=0;i<mArrayLength;i++){
		sprintf(buf,"ShiftFactorcos_%d",i);
		ShiftFactorcos[i] = new TProfile2D(buf,buf,mTotalDay,0,mTotalDay,mTotalCentrality,0,mTotalCentrality);
		sprintf(buf,"ShiftFactorsin_%d",i);
		ShiftFactorsin[i] = new TProfile2D(buf,buf,mTotalDay,0,mTotalDay,mTotalCentrality,0,mTotalCentrality);
	}

	hCentrality9 = new TH1F("hCentrality9","hCentrality9;Centrality;Counts",16,0,16);
	hRefMult = new TH1F("hRefMult","hRefMult;dN_{ch}/d#eta;Counts",1000,0,1000);
	hVertexZ = new TH1F("hVertexZ","hVertexZ;TPC VertexZ (cm);Counts",2000,-100,100);
	hVzDiff = new TH1F("hVzDiff","hVzDiff;Vz_{TPC} - Vz_{VPD} (cm);Counts",200,-10,10);
	hBField = new TH1F("hBField","hBField;Magnetic Filed (KiloGauss);Counts",400,-10,10);

	const Int_t    nPtBins   = 500;
	const Double_t ptLow     = 0;
	const Double_t ptHi      = 5;
	const Int_t    nMassBins = 800;
	const Double_t massLow   = 0;
	const Double_t massHi    = 4;

	//eventPlane
	hRawEventPlane = new TH1F("hRawEventPlane","hRawEventPlane;Reaction Plane (rad); Counts",300,0,TMath::Pi());
	hNewEventPlane = new TH1F("hNewEventPlane","hNewEventPlane;Reaction Plane (rad); Counts",300,0,TMath::Pi());
	hReCenterEventPlane = new TH1F("hReCenterEventPlane","hReCenterEventPlane;Reaction Plane (rad); Counts",300,0,TMath::Pi());
	hFinalEventPlane = new TH1F("hFinalEventPlane","hFinalEventPlane;Reaction Plane (rad); Counts",300,0,TMath::Pi());

	hInclusiveEPhivsPt = new TH2F("hInclusiveEPhivsPt","hInclusiveEPhivsPt;q*p_{T} (GeV/c); #phi",200,-10,10,600,-TMath::Pi(),TMath::Pi());
	hExclusiveEPhivsPt = new TH2F("hExclusiveEPhivsPt","hExclusiveEPhivsPt;q*p_{T} (GeV/c); #phi",200,-10,10,600,-TMath::Pi(),TMath::Pi());
	hCut3EPhivsPt = new TH2F("hCut3EPhivsPt","hCut3EPhivsPt;q*p_{T} (GeV/c); #phi",200,-10,10,600,-TMath::Pi(),TMath::Pi());

	hnEMinusvsEPlus = new TH2F("hnEMinusvsEPlus","hnEMinusvsEPlus;# e^{+};# e^{-}",30,0,30,30,0,30);

	//angleV 
	hULAngleVvsM = new TH2F("hULAngleVvsM","hULAngleVvsM;M_{ee} (GeV/c^{2});#phi_{V} (rad)",1000,0,1.,300,0,TMath::Pi());
	hLPosAngleVvsM = new TH2F("hLPosAngleVvsM","hLPosAngleVvsM;M_{ee} (GeV/c^{2});#phi_{V} (rad)",1000,0,1.,300,0,TMath::Pi());
	hLNegAngleVvsM = new TH2F("hLNegAngleVvsM","hLNegAngleVvsM;M_{ee} (GeV/c^{2});#phi_{V} (rad)",1000,0,1.,300,0,TMath::Pi());
	hMixULAngleVvsM = new TH2F("hMixULAngleVvsM","hMixULAngleVvsM;M_{ee} (GeV/c^{2});#phi_{V} (rad)",1000,0,1.,300,0,TMath::Pi());
	hMixLPosAngleVvsM = new TH2F("hMixLPosAngleVvsM","hMixLPosAngleVvsM;M_{ee} (GeV/c^{2});#phi_{V} (rad)",1000,0,1.,300,0,TMath::Pi());
	hMixLNegAngleVvsM = new TH2F("hMixLNegAngleVvsM","hMixLNegAngleVvsM;M_{ee} (GeV/c^{2});#phi_{V} (rad)",1000,0,1.,300,0,TMath::Pi());

	//without phiV cut
	// hULMvsPtwophiV = new TH2F("hULMvsPtwophiV","hULMvsPtwophiV;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	// hLPosMvsPtwophiV = new TH2F("hLPosMvsPtwophiV","hLPosMvsPtwophiV;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	// hLNegMvsPtwophiV = new TH2F("hLNegMvsPtwophiV","hLNegMvsPtwophiV;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	// hMixULMvsPtwophiV = new TH2F("hMixULMvsPtwophiV","hMixULMvsPtwophiV;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	// hMixLPosMvsPtwophiV = new TH2F("hMixLPosMvsPtwophiV","hMixLPosMvsPtwophiV;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	// hMixLNegMvsPtwophiV = new TH2F("hMixLNegMvsPtwophiV","hMixLNegMvsPtwophiV;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);

	//with phiV cut
	hULMvsPt = new TH2F("hULMvsPt","hULMvsPt;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	hLPosMvsPt = new TH2F("hLPosMvsPt","hLPosMvsPt;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	hLNegMvsPt = new TH2F("hLNegMvsPt","hLNegMvsPt;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	hMixULMvsPt = new TH2F("hMixULMvsPt","hMixULMvsPt;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	hMixLPosMvsPt = new TH2F("hMixLPosMvsPt","hMixLPosMvsPt;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	hMixLNegMvsPt = new TH2F("hMixLNegMvsPt","hMixLNegMvsPt;p_{T} (GeV/c);M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,nMassBins,massLow,massHi);
	hElectronCenvsPt = new TH2F("hElectronCenvsPt","hElectronCenvsPt;p_{T} (GeV/c);Centrality",nPtBins,ptLow,ptHi,16,0,16);
	hPositronCenvsPt = new TH2F("hPositronCenvsPt","hPositronCenvsPt;p_{T} (GeV/c);Centrality",nPtBins,ptLow,ptHi,16,0,16);

	//add centrality dimension
	hULMvsPtCen = new TH3F("hULMvsPtCen","hULMvsPtCen;p_{T} (GeV/c);Centrality;M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,16,0,16,nMassBins,massLow,massHi);
	hLPosMvsPtCen = new TH3F("hLPosMvsPtCen","hLPosMvsPtCen;p_{T} (GeV/c);Centrality;M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,16,0,16,nMassBins,massLow,massHi);
	hLNegMvsPtCen = new TH3F("hLNegMvsPtCen","hLNegMvsPtCen;p_{T} (GeV/c);Centrality;M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,16,0,16,nMassBins,massLow,massHi);
	hMixULMvsPtCen = new TH3F("hMixULMvsPtCen","hMixULMvsPtCen;p_{T} (GeV/c);Centrality;M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,16,0,16,nMassBins,massLow,massHi);
	hMixLPosMvsPtCen = new TH3F("hMixLPosMvsPtCen","hMixLPosMvsPtCen;p_{T} (GeV/c);Centrality;Centrality;M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,16,0,16,nMassBins,massLow,massHi);
	hMixLNegMvsPtCen = new TH3F("hMixLNegMvsPtCen","hMixLNegMvsPtCen;p_{T} (GeV/c);Centrality;M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,16,0,16,nMassBins,massLow,massHi);
	hMixLPosMvsPtCen_CutedbyPhiV = new TH3F("hMixLPosMvsPtCen_CutedbyPhiV","hMixLPosMvsPtCen_CutedbyPhiV;p_{T} (GeV/c);Centrality;Centrality;M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,16,0,16,nMassBins,massLow,massHi);
	hMixLNegMvsPtCen_CutedbyPhiV = new TH3F("hMixLNegMvsPtCen_CutedbyPhiV","hMixLNegMvsPtCen_CutedbyPhiV;p_{T} (GeV/c);Centrality;M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,16,0,16,nMassBins,massLow,massHi);
	hULMvsPhiCen = new TH3F("hULMvsPhiCen","hULMvsPhiCen;Phi;Centrality;M_{ee} (GeV/c^{2})",600,-3.14-0.1,3.14+0.1,16,0,16,1000,massLow,massHi);
	hLPosMvsPhiCen = new TH3F("hLPosMvsPhiCen","hLPosMvsPhiCen;Phi;Centrality;M_{ee} (GeV/c^{2})",600,-3.14-0.1,3.14+0.1,16,0,16,1000,massLow,massHi);
	hLNegMvsPhiCen = new TH3F("hLNegMvsPhiCen","hLNegMvsPhiCen;Phi;Centrality;M_{ee} (GeV/c^{2})",600,-3.14-0.1,3.14+0.1,16,0,16,1000,massLow,massHi);
	hULMvsPtCen_CutedbyPhiV = new TH3F("hULMvsPtCen_CutedbyPhiV","hULMvsPtCen;p_{T} (GeV/c^{2});Centrality;M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,16,0,16,nMassBins,massLow,massHi);
	hLPosMvsPtCen_CutedbyPhiV = new TH3F("hLPosMvsPtCen_CutedbyPhiV","hLPosMvsPtCen;p_{T} (GeV/c^{2});Centrality;M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,16,0,16,nMassBins,massLow,massHi);
	hLNegMvsPtCen_CutedbyPhiV = new TH3F("hLNegMvsPtCen_CutedbyPhiV","hLNegMvsPtCen;p_{T} (GeV/c^{2});Centrality;M_{ee} (GeV/c^{2})",nPtBins,ptLow,ptHi,16,0,16,nMassBins,massLow,massHi);
	hCellIDDiff = new TH1D("hCellIDDiff","hCellIDDiff;ID Diff ;counts;",32,-16,16);

	// hULMvsPtCen->Sumw2();
	// hLPosMvsPtCen->Sumw2();
	// hLNegMvsPtCen->Sumw2();
	// hMixULMvsPtCen->Sumw2();
	// hMixLPosMvsPtCen->Sumw2();
	// hMixLNegMvsPtCen->Sumw2();

	//low ee pair pT cos_{theta} distribution
	// hULCosThetavsMvsCen = new TH3F("hULCosThetavsMvsCen","hULCosThetavsMvsCen; Centrality; M_{ee} (GeV/c^{2}); cos_{#theta}",16,0,16,400,0,4,100,0,1);
	// hLPosCosThetavsMvsCen = new TH3F("hLPosCosThetavsMvsCen","hLPosCosThetavsMvsCen; Centrality; M_{ee} (GeV/c^{2}); cos_{#theta}",16,0,16,400,0,4,100,0,1);
	// hLNegCosThetavsMvsCen = new TH3F("hLNegCosThetavsMvsCen","hLNegCosThetavsMvsCen; Centrality; M_{ee} (GeV/c^{2}); cos_{#theta}",16,0,16,400,0,4,100,0,1);
	// hMixULCosThetavsMvsCen = new TH3F("hMixULCosThetavsMvsCen","hMixULCosThetavsMvsCen; Centrality; M_{ee} (GeV/c^{2}); cos_{#theta}",16,0,16,400,0,4,100,0,1);
	// hMixLPosCosThetavsMvsCen = new TH3F("hMixLPosCosThetavsMvsCen","hMixLPosCosThetavsMvsCen; Centrality; M_{ee} (GeV/c^{2}); cos_{#theta}",16,0,16,400,0,4,100,0,1);
	// hMixLNegCosThetavsMvsCen = new TH3F("hMixLNegCosThetavsMvsCen","hMixLNegCosThetavsMvsCen; Centrality; M_{ee} (GeV/c^{2}); cos_{#theta}",16,0,16,400,0,4,100,0,1);

	//low ee pair pT indivadual electron pT distribution
	// hULePtvsMeevsCen = new TH3F("hULePtvsMeevsCen", "hULePtvsMeevsCen; Centrality; M_{ee} (GeV/c^{2}); eletron p_{T} (GeV/c)",16,0,16,400,0,4,1000,0,10);
	// hLSePtvsMeevsCen = new TH3F("hLSePtvsMeevsCen", "hLSePtvsMeevsCen; Centrality; M_{ee} (GeV/c^{2}); eletron p_{T} (GeV/c)",16,0,16,400,0,4,1000,0,10);
	// hMixULePtvsMeevsCen = new TH3F("hMixULePtvsMeevsCen", "hMixULePtvsMeevsCen; Centrality; M_{ee} (GeV/c^{2}); eletron p_{T} (GeV/c)",16,0,16,400,0,4,1000,0,10);
	// hMixLSePtvsMeevsCen = new TH3F("hMixLSePtvsMeevsCen", "hMixLSePtvsMeevsCen; Centrality; M_{ee} (GeV/c^{2}); eletron p_{T} (GeV/c)",16,0,16,400,0,4,1000,0,10);
}
//=======================================================================================
void writeHistograms(char* outFile)
{
	char buf[1024];
	sprintf(buf,"%s.histo.root",outFile);
	cout<<"Writing histograms into "<<buf<<endl;
	TFile *mFile = new TFile(buf,"recreate");
	mFile->cd();

	//in passEvent function
	hCentrality9->Write();
	hRefMult->Write();
	hVertexZ->Write();
	hVzDiff->Write();
	hBField->Write();

	//eventPlane
	hRawEventPlane->Write();
	hNewEventPlane->Write();
	hFinalEventPlane->Write();
	hReCenterEventPlane->Write();
	hInclusiveEPhivsPt->Write();
	hExclusiveEPhivsPt->Write();
	hCut3EPhivsPt->Write();

	hnEMinusvsEPlus->Write();

	//angleV 
	hULAngleVvsM->Write();
	hLPosAngleVvsM->Write();
	hLNegAngleVvsM->Write();
	hMixULAngleVvsM->Write();
	hMixLPosAngleVvsM->Write();
	hMixLNegAngleVvsM->Write();

	//without phiV cut
	// hULMvsPtwophiV->Write();
	// hLPosMvsPtwophiV->Write();
	// hLNegMvsPtwophiV->Write();
	// hMixULMvsPtwophiV->Write();
	// hMixLPosMvsPtwophiV->Write();
	// hMixLNegMvsPtwophiV->Write();

	//with phiV cut
	hULMvsPt->Write();
	hLPosMvsPt->Write();
	hLNegMvsPt->Write();
	hMixULMvsPt->Write();
	hMixLPosMvsPt->Write();
	hMixLNegMvsPt->Write();

	//add centrality dimension
	hULMvsPtCen->Write();
	hLPosMvsPtCen->Write();
	hLNegMvsPtCen->Write();
	hMixULMvsPtCen->Write();
	hMixLPosMvsPtCen->Write();
	hMixLNegMvsPtCen->Write();
	hULMvsPhiCen->Write();
  	hLPosMvsPhiCen->Write();
  	hLNegMvsPhiCen->Write();
	hULMvsPtCen_CutedbyPhiV->Write();
	hLPosMvsPtCen_CutedbyPhiV->Write();
	hLNegMvsPtCen_CutedbyPhiV->Write();
	hMixLPosMvsPtCen_CutedbyPhiV->Write();
	hMixLNegMvsPtCen_CutedbyPhiV->Write();
	hElectronCenvsPt->Write();
	hPositronCenvsPt->Write();
	hCellIDDiff->Write();
	

	// hULCosThetavsMvsCen->Write();
	// hLPosCosThetavsMvsCen->Write();
	// hLNegCosThetavsMvsCen->Write();
	// hMixULCosThetavsMvsCen->Write();
	// hMixLPosCosThetavsMvsCen->Write();
	// hMixLNegCosThetavsMvsCen->Write();

	// hULePtvsMeevsCen->Write();
	// hLSePtvsMeevsCen->Write();
	// hMixULePtvsMeevsCen->Write();
	// hMixLSePtvsMeevsCen->Write();
}
//==============================================================================================
Bool_t Init()
{
	cout<<endl;

	ifstream indata;

	indata.open("/star/u/tc88qy/AuAu/run17/54GeV/QA/runList/output_all/GetRun/runList/mTotalDayList.dat");
	mTotalDayId.clear();
	if(indata.is_open()){
		cout<<"read in day number list and recode day number ...";
		Int_t oldId;
		Int_t newId=0;
		while(indata>>oldId){
			mTotalDayId[oldId] = newId;
			newId++;
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the day number list !!!"<<endl;
		return kFALSE;
	}
	indata.close();

	indata.open("/star/u/tc88qy/AuAu/run17/54GeV/QA/runList/output_all/GetRun/runList/mTotalRunList.dat");
	mTotalRunId.clear();
	if(indata.is_open()){
		cout<<"read in total run number list and recode run number ...";
		Int_t oldId;
		Int_t newId=0;
		while(indata>>oldId){
			mTotalRunId[oldId] = newId;
			newId++;
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the total run number list !!!"<<endl;
		return kFALSE;
	}
	indata.close();

	//read in bad run for 580001 and 580021
	ifstream indata_001;
	indata_001.open("/star/u/tc88qy/AuAu/run17/54GeV/EventPlane/reCenter/Badrun/mBadRunList_001.dat");
	mBadRunId_001.clear();
	if(indata_001.is_open()){
		cout<<"read in bad run list for trigger 580001 ...";
		Int_t oldId;
		Int_t newId=0;
		while(indata_001>>oldId){
			mBadRunId_001[oldId] = newId;
			newId++;
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the total run number list !!!"<<endl;
		return kFALSE;
	}
	indata_001.close();

	//read in bad run for 580001 and 580021
	ifstream indata_021;
	indata_021.open("/star/u/tc88qy/AuAu/run17/54GeV/EventPlane/reCenter/Badrun/mBadRunList_021.dat");
	mBadRunId_021.clear();
	if(indata_021.is_open()){
		cout<<"read in bad run list for trigger 580021 ...";
		Int_t oldId;
		Int_t newId=0;
		while(indata_021>>oldId){
			mBadRunId_021[oldId] = newId;
			newId++;
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the total run number list !!!"<<endl;
		return kFALSE;
	}
	indata_021.close();

	cout<<"bad run for trigger 580001"<<endl;
	for(map<Int_t,Int_t>::iterator iter=mBadRunId_001.begin();iter!=mBadRunId_001.end();iter++)
		cout<<iter->second<<" \t"<<iter->first<<endl;

	cout<<"bad run for trigger 580021"<<endl;
	for(map<Int_t,Int_t>::iterator iter=mBadRunId_021.begin();iter!=mBadRunId_021.end();iter++)
		cout<<iter->second<<" \t"<<iter->first<<endl;

	TFile *fReCenter = TFile::Open("/star/u/tc88qy/AuAu/run17/54GeV/EventPlane/reCenter/output_all/reCenter.root");
	if(fReCenter->IsOpen()){
		cout<<"read in re-center root file ...";
		etapluszplusQx   = (TProfile2D *)fReCenter->Get("etapluszplusQx");
		etapluszminusQx  = (TProfile2D *)fReCenter->Get("etapluszminusQx");
		etaminuszplusQx  = (TProfile2D *)fReCenter->Get("etaminuszplusQx");
		etaminuszminusQx = (TProfile2D *)fReCenter->Get("etaminuszminusQx");
		etapluszplusQy   = (TProfile2D *)fReCenter->Get("etapluszplusQy");
		etapluszminusQy  = (TProfile2D *)fReCenter->Get("etapluszminusQy");
		etaminuszplusQy  = (TProfile2D *)fReCenter->Get("etaminuszplusQy");
		etaminuszminusQy = (TProfile2D *)fReCenter->Get("etaminuszminusQy");
	}

	TFile *fShift = TFile::Open("/star/u/tc88qy/AuAu/run17/54GeV/EventPlane/shift/output_all/shift.histo.root");
	if(fShift->IsOpen()){
		cout<<"read in shiftfactor root file ...";
		for(int i=0;i<mArrayLength;i++){
			ShiftFactorcos[i] = (TProfile2D*)fShift->Get(Form("shiftfactorcos_%d",i));
			ShiftFactorsin[i] = (TProfile2D*)fShift->Get(Form("shiftfactorsin_%d",i));
		}
		cout<<" [OK]"<<endl;
	}


	phiVcut=new TF1("phiVcut","0.84326*exp((-49.4819)*x)+(-0.996609)*x+(0.19801)",0.,1.0); //jie's cut


	cout<<"Initialization DONE !!!"<<endl;
	cout<<endl;

	return kTRUE;
}
