#include <Higgsbb.hh>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/MCParticleImpl.h>
#include <EVENT/LCRelation.h>
#include <EVENT/Vertex.h>
#include <IMPL/LCRelationImpl.h>
#include <values.h>
#include <string>
#include <iostream>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <stdexcept>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TRandom.h>
#include <Rtypes.h>
#include <sstream>
#include <cmath>
#include <vector>
#include <TMath.h>
#include "TLorentzVector.h"
#include <UTIL/PIDHandler.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/Config.h"

using namespace std;

Higgsbb a_Higgsbb_instance;


Higgsbb::Higgsbb()
: Processor("Higgsbb"),
_output(0)
{
    _description = "Print MC Truth" ;
    
    _treeFileName="MCTruth.root";
    registerProcessorParameter( "TreeOutputFile" ,
                               "The name of the file to which the ROOT tree will be written" ,
                               _treeFileName ,
                               _treeFileName);
    
    _treeName="Tau";
    registerProcessorParameter( "TreeName" ,
                               "The name of the ROOT tree" ,
                               _treeName ,
                               _treeName);
    
    _overwrite=0;
    registerProcessorParameter( "OverwriteFile" ,
                               "If zero an already existing file will not be overwritten." ,
                               _overwrite ,
                               _overwrite);
    
}


void findCombine(int index, vector<int> &a, vector<int>& tmp, const int& n, vector<vector<int> >& res){
	if(tmp.size() == n/2){
		res.push_back(tmp);
		return;
	}
	for(int i=index; i<=n; i++){
		tmp.push_back(a[i-1]);
		findCombine(i+1, a, tmp, n, res);
		tmp.pop_back();
	}
	return;
}
vector<vector<int> > Combine(int n, vector<int> a){
	vector<vector<int> > res;
	if(n<1 || n%2!=0 || a.size()<1)
		return res;
	vector<int> tmp;
	findCombine(1, a, tmp, n, res);
	return res;
}
vector<int> DifferenceSet(vector<int> total, vector<int> original, vector<int> left){
    for(int i = 0; i<total.size(); i++)
    {
        int element = total[i];
        vector<int >::iterator it;
        it = find(original.begin(), original.end(), element);
        if(it == original.end()){
            left.push_back(element);
        }
    }
    return left;
}
std::vector<vector<int> > pair4jets(int numjets);

static bool sortEn(ReconstructedParticle* a1, ReconstructedParticle* a2){
    return a1->getEnergy() >= a2->getEnergy();
}

ReconstructedParticle* getMulti( std::vector<ReconstructedParticle* > vec, std::vector<double> &result );
void CalcuThrust(std::vector<TLorentzVector > UsedForThrust, std::vector<double> &result);
double CalcuSphericity(std::vector<TLorentzVector> UsedForSphericity, std::vector<double> &result);

/*
void CalcuXYVector(TVector3 V3, std::vector<TVector3> &result);

int ISBbar(int mc);
int ISB(int mc);
int ISCbar(int mc);
int ISC(int mc);
int ISS(int mc);
*/

void Higgsbb::init() {
    
    printParameters();
    TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));
    if (!tree_file->IsOpen()) {
        delete tree_file;
        tree_file=new TFile(_treeFileName.c_str(),"NEW");
    }
    
    _outputTree = new TTree(_treeName.c_str(),_treeName.c_str());
    _outputTree->SetAutoSave(32*1024*1024);  // autosave every 32MB
    _outputTree->Branch("EventNr",     &eventNr,       "EventNr/I");
    _outputTree->Branch("Num",         &Num,           "Num/I");
    _outputTree->Branch("num_lepton",  &num_lepton,    "num_lepton/I");
    _outputTree->Branch("num_quark",   &num_quark,     "num_quark/I");
    _outputTree->Branch("num_neutrino", &num_neutrino, "num_neutrino/I");
    _outputTree->Branch("ISR_En",      &ISR_En,        "ISR_En/D");
    _outputTree->Branch("count94",     &count94,       "count94/I");
    _outputTree->Branch("mass941",     &mass941,       "mass941/D");
    _outputTree->Branch("mass942",     &mass942,       "mass942/D");
    _outputTree->Branch("boson1Quark", &boson1Quark);
    _outputTree->Branch("boson2Quark", &boson2Quark);
    _outputTree->Branch("noHQPDG",     noHQPDG,        "noHQPDG[4]/I");
    _outputTree->Branch("HQPDG",       HQPDG,          "HQPDG[4]/I");
    _outputTree->Branch("B1Quark",     B1Quark,        "B1Quark[2]/I");
    _outputTree->Branch("B2Quark",     B2Quark,        "B2Quark[2]/I");
    _outputTree->Branch("OriLepPDG",   OriLepPDG,      "OriLepPDG[4]/I");
    _outputTree->Branch("HbbL",        &HbbL,           "HbbL[2]/F");

    _outputTree->Branch("HccL",        &HccL,           "HccL[2]/F");

    _outputTree->Branch("HooL",        HooL,           "HooL[2]/F");

    _outputTree->Branch("ZbbL",        &ZbbL,           "ZbbL[2]/F");

    _outputTree->Branch("ZccL",        &ZccL,           "ZccL[2]/F");

    _outputTree->Branch("ZooL",        ZooL,            "ZooL[2]/F");

    _outputTree->Branch("ISREn",       ISREn,          "ISREn[5]/F");
    _outputTree->Branch("ISRCosTheta", ISRCosTheta,    "ISRCosTheta[5]/F");
    _outputTree->Branch("HiggsJetsAngle",     &HiggsJetsAngle,     "HiggsJetsAngle/D");
    _outputTree->Branch("ZJetsAngle",         &ZJetsAngle,         "ZJetsAngle/D");
    _outputTree->Branch("ZHAngle",            &ZHAngle,            "ZHAngle/D");
    _outputTree->Branch("HCharge",            &HCharge,            "HCharge/I");
    _outputTree->Branch("ZCharge",            &ZCharge,            "ZCharge/I");

    _outputTree->Branch("B1Costheta",         &B1Costheta,         "B1Costheta/D");
    _outputTree->Branch("B2Costheta",         &B2Costheta,         "B2Costheta/D");

    _outputTree->Branch("numHiggsParticles",  &numHiggsParticles,  "numHiggsParticles/I");
    _outputTree->Branch("HiggsX",             HiggsX,              "HiggsX[4]/F");
    _outputTree->Branch("HiggsSumL",          HiggsSumL,           "HiggsSumL[4]/F");
    _outputTree->Branch("jetPx",              jetPx,           "jetPx[4]/F");
    _outputTree->Branch("jetPy",              jetPy,           "jetPy[4]/F");
    _outputTree->Branch("jetPz",              jetPz,           "jetPz[4]/F");
    _outputTree->Branch("jetEn",              jetEn,           "jetEn[4]/F");

    _outputTree->Branch("numChargeOfJets",    numChargeOfJets,     "numChargeOfJets[4]/I");
    _outputTree->Branch("maxHiggsX",          &maxHiggsX,          "maxHiggsX/F");
    _outputTree->Branch("miniAngleAmongJets", &miniAngleAmongJets, "miniAngleAmongJets/F");
    _outputTree->Branch("RminDif",            &RminDif,            "RminDif/D");
    _outputTree->Branch("RminDifZH",          &RminDifZH,            "RminDifZH/D");

    
    
    _outputTree->Branch("num_Vertex",   &num_Vertex,    "num_Vertex/I");
    _outputTree->Branch("num_jet",      &num_jet,       "num_jet/I");
    _outputTree->Branch("num_charge",   &num_charge,    "num_charge/I");
    _outputTree->Branch("Pmax",         &Pmax,          "Pmax/D");
    _outputTree->Branch("Pt",           &Pt,            "Pt/D");
    _outputTree->Branch("Pl",           &Pl,            "Pl/D");

    _outputTree->Branch("sumb",         &sumb,          "sumb/D");
    _outputTree->Branch("sum2b",        &sum2b,         "sum2b/D");
    _outputTree->Branch("sumc",         &sumc,          "sumc/D");
    _outputTree->Branch("sum2c",        &sum2c,         "sum2c/D");
    _outputTree->Branch("massB1b",      &massB1b,       "massB1b/D");
    _outputTree->Branch("massB2b",      &massB2b,       "massB2b/D");
    _outputTree->Branch("massRecoil1b", &massRecoil1b,  "massRecoil1b/D");
    _outputTree->Branch("massRecoil2b", &massRecoil2b,  "massRecoil2b/D");
    _outputTree->Branch("massB1c",      &massB1c,       "massB1c/D");
    _outputTree->Branch("massB2c",      &massB2c,       "massB2c/D");
    _outputTree->Branch("massRecoil1c", &massRecoil1c,  "massRecoil1c/D");
    _outputTree->Branch("massRecoil2c", &massRecoil2c,  "massRecoil2c/D");
    _outputTree->Branch("alpha",        &alpha,         "alpha/D");
    _outputTree->Branch("delta_R1",     &delta_R1,      "delta_R1/D");
    _outputTree->Branch("delta_R2",     &delta_R2,      "delta_R2/D");
    _outputTree->Branch("visEn",        &visEn,         "visEn/D");
    _outputTree->Branch("multiplicity", &multiplicity,  "multiplicity/I");
    _outputTree->Branch("LepNum",       &LepNum,        "LepNum/I");
    _outputTree->Branch("LeadLepEn",    &LeadLepEn,     "LeadLepEn/D");
    _outputTree->Branch("AvLepEn",      &AvLepEn,       "AvLepEn/D");
    _outputTree->Branch("hadNum",       &hadNum,        "hadNum/I");
    _outputTree->Branch("LeadhadEn",    &LeadhadEn,     "LeadhadEn/D");
    _outputTree->Branch("AvhadEn",      &AvhadEn,       "AvhadEn/D");
    _outputTree->Branch("gamaNum",      &gamaNum,       "gamaNum/I");
    _outputTree->Branch("LeadgamaEn",   &LeadgamaEn,    "LeadgamaEn/D");
    _outputTree->Branch("AvgamaEn",     &AvgamaEn,      "AvgamaEn/D");
    _outputTree->Branch("SubLeadgamaEn",&SubLeadgamaEn, "SubLeadgamaEn/D");
    _outputTree->Branch("LeadneutralEn",&LeadneutralEn, "LeadneutralEn/D");
    _outputTree->Branch("LeadElecEn",   &LeadElecEn,    "LeadElecEn/D");
    _outputTree->Branch("ElecNum",      &ElecNum,       "ElecNum/I");
    _outputTree->Branch("AvElecEn",     &AvElecEn,      "AvElecEn/D");
    _outputTree->Branch("MuonNum",      &MuonNum,       "MuonNum/I");
    _outputTree->Branch("LeadMuonEn",   &LeadMuonEn,    "LeadMuonEn/D");
    _outputTree->Branch("AvMuonEn",     &AvMuonEn,      "AvMuonEn/D");
    
    _outputTree->Branch("Thrust",        &Thrust,        "Thrust/D");
    _outputTree->Branch("MaxBroadening", &MaxBroadening, "MaxBroadening/D");
    _outputTree->Branch("HeavyMass",     &HeavyMass,     "HeavyMass/D");
    _outputTree->Branch("Sphericity",    &Sphericity,    "Sphericity/D");
    _outputTree->Branch("CPara",         &CPara,         "CPara/D");
    _outputTree->Branch("DPara",         &DPara,         "DPara/D");
    _outputTree->Branch("cosTthrust",    &cosTthrust,    "cosTthrust/D");
    _outputTree->Branch("Y01",           &Y01,           "Y01/D");
    _outputTree->Branch("Y12",           &Y12,           "Y12/D");
    _outputTree->Branch("Y23",           &Y23,           "Y23/D");
    _outputTree->Branch("Y34",           &Y34,           "Y34/D");
    _outputTree->Branch("Y45",           &Y45,           "Y45/D");
    _outputTree->Branch("boostThrust",   &boostThrust,   "boostThrust/D");
    _outputTree->Branch("boostH1En",     &boostH1En,     "boostH1En/D");
    _outputTree->Branch("boostH2En",     &boostH2En,     "boostH2En/D");
    _outputTree->Branch("boostH1Mass",   &boostH1Mass,   "boostH1Mass/D");
    _outputTree->Branch("boostH2Mass",   &boostH2Mass,   "boostH2Mass/D");
    
    _outputTree->Branch("LLRecoilMass",    &LLRecoilMass,    "LLRecoilMass/D");
    _outputTree->Branch("LLInvMass",       &LLInvMass,       "LLInvMass/D");
    _outputTree->Branch("JetsInvMass",     &JetsInvMass,     "JetsInvMass/D");
    _outputTree->Branch("JetsRecoilMass",  &JetsRecoilMass,  "JetsRecoilMass/D");
    _outputTree->Branch("haveMuplusminus", &haveMuplusminus, "haveMuplusminus/I");
    _outputTree->Branch("haveEplusminus",  &haveEplusminus,  "haveEplusminus/I");
    
    /*
    _outputTree->Branch("P1HBMass", &P1HBMass, "P1HBMass/D");
    _outputTree->Branch("P1LBMass", &P1LBMass, "P1LBMass/D");
    _outputTree->Branch("P1HF",     &P1HF,     "P1HF/I");
    _outputTree->Branch("P1LF",     &P1LF,     "P1LF/I");
    _outputTree->Branch("P1HLL",    &P1HLL,    "P1HLL/D");
    _outputTree->Branch("P1LLL",    &P1LLL,    "P1LLL/D");
    
    _outputTree->Branch("P2HBMass", &P2HBMass, "P2HBMass/D");
    _outputTree->Branch("P2LBMass", &P2LBMass, "P2LBMass/D");
    _outputTree->Branch("P2HF",     &P2HF,     "P2HF/I");
    _outputTree->Branch("P2LF",     &P2LF,     "P2LF/I");
    _outputTree->Branch("P2HLL",    &P2HLL,    "P2HLL/D");
    _outputTree->Branch("P2LLL",    &P2LLL,    "P2LLL/D");
*/
    _outputTree->Branch("pairMB1",   &pairMB1,    "pairMB1/D");
    _outputTree->Branch("pairMB2",   &pairMB2,    "pairMB2/D");
    

 
    Num = 0;
}

void Higgsbb::processEvent( LCEvent * evtP )
{
    
    if (evtP)
    {
        try{
            
            cout<<"Next Event *******************************************************************************************************"<<endl;
            eventNr = evtP->getEventNumber();
            cout<<"eventNr : "<<eventNr<<" Num : "<<Num<<endl;
            
            //ArborPFOs
            LCCollection* col_PFO = evtP->getCollection( "ArborPFOs" );
            int nPFO = col_PFO->getNumberOfElements();
            cout<<"nPFO : "<<nPFO<<endl;
            Pmax = 0; Pt = 0; Pl = 0;
            num_charge = 0;
            TLorentzVector TLPFO(0,0,0,0);
            std::vector<ReconstructedParticle*> vLepton;    vLepton.clear();
            std::vector<ReconstructedParticle*> vElec;      vElec.clear();
            std::vector<ReconstructedParticle*> vMuon;      vMuon.clear();
            std::vector<ReconstructedParticle*> vHadron;    vHadron.clear();
            std::vector<ReconstructedParticle*> vGamma;     vGamma.clear();
            std::vector<ReconstructedParticle*> vNeutral;   vNeutral.clear();
            std::vector<ReconstructedParticle*> vCHad;      vCHad.clear();
            std::vector<TLorentzVector>         vArborTL;   vArborTL.clear();
            std::vector<ReconstructedParticle*> muplusvec;  muplusvec.clear();
            std::vector<ReconstructedParticle*> muminusvec; muminusvec.clear();
            std::vector<ReconstructedParticle*> eplusvec;   eplusvec.clear();
            std::vector<ReconstructedParticle*> eminusvec;  eminusvec.clear();
            for(int i = 0; i<nPFO; i++){
                ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(i));
                TLorentzVector temp(pfo->getMomentum(), pfo->getEnergy());
                TVector3 TVtemp = temp.Vect();
                vArborTL.push_back(temp);
                TLPFO += temp;
                int type = abs(pfo->getType());
                if(type == 11 || type == 13){
                    vLepton.push_back(pfo);
                    if(type == 11){vElec.push_back(pfo);}
                    else if(type == 13){vMuon.push_back(pfo);}
                }
                else if(type == 22){vGamma.push_back(pfo); vNeutral.push_back(pfo);}
                else {
                    vHadron.push_back(pfo);
                    if(pfo->getCharge() == 0){vNeutral.push_back(pfo);}
                    else{ vCHad.push_back(pfo); }
                }
                if(pfo->getType() == 11){eminusvec.push_back(pfo);}
                else if(pfo->getType() == -11){eplusvec.push_back(pfo);}
                else if(pfo->getType() == 13){muminusvec.push_back(pfo);}
                else if(pfo->getType() == -13){muplusvec.push_back(pfo);}
                if(pfo->getCharge() != 0){
                    if(TVtemp.Mag() > Pmax){
                        Pmax = TVtemp.Mag();
                    }
                    num_charge += 1;
                }
            }
            cout<<"the total energy of ArborPFOs is : "<<TLPFO.E()<<endl;
            TVector3 TVPFO = TLPFO.Vect();
            Pt = TVPFO.Perp();
            Pl = pow(TVPFO.Mag2() - TVPFO.Perp2(), 0.5);
            TLorentzVector TLtest(0,0,0,0);
            TVector3 neg_tot_R = -1*TLPFO.BoostVector();
            std::vector<TLorentzVector>      vArborTLBoost; vArborTLBoost.clear();
            for(int i = 0; i<nPFO; i++){
                ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(col_PFO->getElementAt(i));
                TLorentzVector temp(pfo->getMomentum(), pfo->getEnergy());
                temp.Boost(neg_tot_R);
                TLtest += temp;
                vArborTLBoost.push_back(temp);
            }
//            cout<<"TLPFO.E() : "<<TLPFO.E()<<" TLtest.E() : "<<TLtest.E()<<endl;
            
            // the following code used to analysis LLH
//            cout<<"eminusvec.size() : "<<eminusvec.size()<<" eplusvec.size() : "<<eplusvec.size()<<" muplusvec.size() : "<<muplusvec.size()<<" muminusvec.size() : "<<muminusvec.size()<<endl;
            double muPlusMinusMass = 0, ePlusMinusMass = 0, muRecoilMass = 0, eRecoilMass = 0;
            LLRecoilMass = 0; LLInvMass = 0;
            TLorentzVector TLT240(0, 0, 0, 240);
            TLorentzVector TLTmu(0, 0, 0, 0);
            TLorentzVector TLTe(0,0,0,0);
            haveEplusminus = 0; haveMuplusminus = 0;
            if( eplusvec.size() != 0 && eminusvec.size() != 0 ){
                haveEplusminus = 1;
                sort(eplusvec.begin(), eplusvec.end(), sortEn);
                sort(eminusvec.begin(), eminusvec.end(), sortEn);
                ReconstructedParticle* eplus = eplusvec.at(0);
                ReconstructedParticle* eminus = eminusvec.at(0);
                TLorentzVector TLeplus(eplus->getMomentum(), eplus->getEnergy());
                TLorentzVector TLeminus(eminus->getMomentum(), eminus->getEnergy());
                TLorentzVector TLT = TLeplus + TLeminus;
                TLTe = TLT;
                TLorentzVector TLRecoil = TLT240 - TLT;
                ePlusMinusMass = TLT.M();
                eRecoilMass = TLRecoil.M();
            }
            if( muplusvec.size() != 0 && muminusvec.size() != 0 ){
                haveMuplusminus = 1;
                sort(muplusvec.begin(), muplusvec.end(), sortEn);
                sort(muminusvec.begin(), muminusvec.end(), sortEn);
                ReconstructedParticle* muplus = muplusvec.at(0);
                ReconstructedParticle* muminus = muminusvec.at(0);
                TLorentzVector TLmuplus(muplus->getMomentum(), muplus->getEnergy());
                TLorentzVector TLmuminus(muminus->getMomentum(), muminus->getEnergy());
                TLorentzVector TLT = TLmuplus + TLmuminus;
                TLTmu = TLT;
                TLorentzVector TLRecoil = TLT240 - TLT;
                muPlusMinusMass = TLT.M();
                muRecoilMass = TLRecoil.M();
            }
//            cout<<"have eplus and eminus : "<<haveEplusminus<<"  have muplus and muminus : "<<haveMuplusminus<<endl;
            TLorentzVector TLLL(0,0,0,0);
//            cout<<"ePlusMinusMass : "<<ePlusMinusMass<<" muPlusMinusMass : "<<muPlusMinusMass<<endl;
            if(abs(ePlusMinusMass - 91.2) <= abs(muPlusMinusMass - 91.2)){
                LLRecoilMass = eRecoilMass;
                LLInvMass = ePlusMinusMass;
                TLLL = TLTe;
            }
            else if(abs(ePlusMinusMass - 91.2) > abs(muPlusMinusMass - 91.2)){
                LLRecoilMass = muRecoilMass;
                LLInvMass = muPlusMinusMass;
                TLLL = TLTmu;
            }
                
                
            //the following code used to analysis the multiplicity information of final state particles
            visEn = 0; multiplicity = 0;
            visEn = TLPFO.E();
            multiplicity = nPFO;
            std::vector<double> leptonResult; leptonResult.clear();
            ReconstructedParticle* leadLepton = getMulti(  vLepton, leptonResult );
            LepNum     = leptonResult.at(0);
            LeadLepEn  = leptonResult.at(1);
            AvLepEn    = leptonResult.at(2);
            std::vector<double> ElecResult; ElecResult.clear();
            ReconstructedParticle* leadElec = getMulti( vElec, ElecResult );
            ElecNum    = ElecResult.at(0);
            LeadElecEn = ElecResult.at(1);
            AvElecEn   = ElecResult.at(2);
            std::vector<double> MuonResult; MuonResult.clear();
            ReconstructedParticle* leadMuon = getMulti( vMuon, MuonResult );
            MuonNum    = MuonResult.at(0);
            LeadMuonEn = MuonResult.at(1);
            AvMuonEn   = MuonResult.at(2);
            std::vector<double> hadResult; hadResult.clear();
            ReconstructedParticle* leadHadron = getMulti( vHadron, hadResult );
            hadNum     = hadResult.at(0);
            LeadhadEn  = hadResult.at(1);
            AvhadEn    = hadResult.at(2);
            std::vector<double> gammaResult; gammaResult.clear();
            ReconstructedParticle* leadGamma = getMulti( vGamma, gammaResult );
            gamaNum    = gammaResult.at(0);
            LeadgamaEn = gammaResult.at(1);
            AvgamaEn   = gammaResult.at(2);
            SubLeadgamaEn = gammaResult.at(3);
            std::vector<double> NeutralResult; NeutralResult.clear();
            ReconstructedParticle* leadNeutral = getMulti( vNeutral, NeutralResult );
            LeadneutralEn = NeutralResult.at(1);
            
            

            //the following code used to get the eventshape variables
            Thrust = 0; HeavyMass = 0; MaxBroadening = 0; Sphericity = 0; DPara = 0; CPara = 0;
            std::vector<double> shape1Result; shape1Result.clear();
            CalcuThrust(vArborTL, shape1Result);
            Thrust = shape1Result.at(0);
            double HemiMass1 = shape1Result.at(1);
            double HemiMass2 = shape1Result.at(2);
            HeavyMass = (HemiMass1 > HemiMass2) ? HemiMass1 : HemiMass2;
            double HemiBroadening1 = shape1Result.at(3);
            double HemiBroadening2 = shape1Result.at(4);
            MaxBroadening = (HemiBroadening1 > HemiBroadening2) ? HemiBroadening1 : HemiBroadening2;
            cosTthrust = shape1Result.at(5);
            
            std::vector<double> boostResult; boostResult.clear();
            CalcuThrust(vArborTLBoost, boostResult);
            boostThrust = boostResult.at(0); boostH1En = boostResult.at(6); boostH2En = boostResult.at(7);
            boostH1Mass = boostResult.at(8); boostH2Mass = boostResult.at(9);
            
            std::vector<double> shape2Result; shape2Result.clear();
            CalcuSphericity(vArborTL, shape2Result);
            Sphericity = shape2Result.at(0);
            DPara = shape2Result.at(1);
            CPara = shape2Result.at(2);
            
            

            

            //RecoJet
            std::map<double, TLorentzVector, greater<double> > bTLmap;       bTLmap.clear();
            std::map<double, TLorentzVector, greater<double> > cTLmap;       cTLmap.clear();
            std::map<double, TLorentzVector, greater<double> > oTLmap;       oTLmap.clear();
            std::map< int,   std::array<double, 3>           > mapProb;      mapProb.clear();

            Y01 = 999, Y12 = 999, Y23 = 999, Y34 = 999, Y45 = 999;
            LCCollection* col_Jet = evtP->getCollection( "RefinedJets" );
            num_jet = col_Jet->getNumberOfElements();
            cout<<"num_jet : "<<num_jet<<endl;
            JetsInvMass = 0;
            TLorentzVector TLJets(0,0,0,0);
            miniAngleAmongJets = 999;
            if(num_jet >= 1){
                if(num_jet >= 2){
                    for(int i = 0; i<(num_jet-1); i++){
                        for(int j = i+1; j<num_jet; j++){
                            ReconstructedParticle* jeti = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(i));
                            TLorentzVector TLtempi(jeti->getMomentum(), jeti->getEnergy());
                            TVector3 TVtempi = TLtempi.Vect();
                            ReconstructedParticle* jetj = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(j));
                            TLorentzVector TLtempj(jetj->getMomentum(), jetj->getEnergy());
                            TVector3 TVtempj = TLtempj.Vect();
                            if(TVtempi.Angle(TVtempj) < miniAngleAmongJets){miniAngleAmongJets = TVtempi.Angle(TVtempj);}
                        }
                    }
                }
                PIDHandler pidh(col_Jet);
                int algo   = 0;
                int ibtag  =-1;
                int ictag  =-1;
                int ibctag =-1;
                int icat   =-1;
                algo      = pidh.getAlgorithmID("lcfiplus");
                ibtag     = pidh.getParameterIndex (algo,  "BTag");
                ictag     = pidh.getParameterIndex (algo,  "CTag");
                ibctag    = pidh.getParameterIndex (algo,  "BCTag");
                icat      = pidh.getParameterIndex (algo,  "Category");
                for(int i = 0; i<num_jet; i++){
                    ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(i));
                    TLorentzVector TLtemp(jet->getMomentum(), jet->getEnergy());
                    TLJets += TLtemp;
                    const ParticleID &pid = pidh.getParticleID(jet, algo);
                    double btag  = pid.getParameters()[ibtag];
                    double ctag  = pid.getParameters()[ictag];
                    double bctag = pid.getParameters()[ibctag];
                    double cat   = pid.getParameters()[icat];
                    double otag = 1 - btag - ctag;
                    std::array<double, 3> myarray; myarray = { btag, ctag, otag };
                    mapProb.insert( std::make_pair(i, myarray) );
//                    mapTLFlavor.insert( std::make_pair(TLtemp, myarray) );
                    
                    cout<<"btag : "<<btag<<" ctag : "<<ctag<<" otag : "<<otag<<endl;
                    if( bTLmap.find( btag ) == bTLmap.end() ){ bTLmap[btag] = TLtemp;}
                    else{ btag = btag + 0.000001 * i; bTLmap[btag] = TLtemp;}
                    if( cTLmap.find( ctag ) == cTLmap.end() ){ cTLmap[ctag] = TLtemp; }
                    else{ ctag = ctag + 0.000001 * i; cTLmap[ctag] = TLtemp; }
                    if( oTLmap.find( otag ) == oTLmap.end() ){ oTLmap[otag] = TLtemp; }
                    else{ otag = otag + 0.000001 * i; oTLmap[otag] = TLtemp; }
                }
                algo = pidh.getAlgorithmID("yth");
                int iy01 = -1, iy12 = -1, iy23 = -1, iy34 = -1, iy45 = -1;
                iy01 = pidh.getParameterIndex (algo, "y01" );
                iy12 = pidh.getParameterIndex (algo, "y12" );
                iy23 = pidh.getParameterIndex (algo, "y23" );
                iy34 = pidh.getParameterIndex (algo, "y34" );
                iy45 = pidh.getParameterIndex (algo, "y45" );
                ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(0));
                const ParticleID &pid = pidh.getParticleID(jet, algo);
                Y01 = pid.getParameters()[iy01];
                Y12 = pid.getParameters()[iy12];
                Y23 = pid.getParameters()[iy23];
                Y34 = pid.getParameters()[iy34];
                Y45 = pid.getParameters()[iy45];
            }
            cout<<"Y01 : "<<Y01<<" Y12 : "<<Y12<<" Y23 : "<<Y23<<" Y34 : "<<Y34<<" Y45 : "<<Y45<<endl;
            JetsInvMass = TLJets.M();
            JetsRecoilMass = 0;
            JetsRecoilMass = (TLT240 - TLJets).M();
            cout<<"ArborPFOs, total energy : "<<(TLJets).E()<<endl;

            std::vector<double> vblikelihood; vblikelihood.clear();
            std::vector<TLorentzVector> vbTL; vbTL.clear();
            std::vector<double> vclikelihood; vclikelihood.clear();
            std::vector<TLorentzVector> vcTL; vcTL.clear();
            std::vector<double> volikelihood; volikelihood.clear();
            std::vector<TLorentzVector> voTL; voTL.clear();

            if(bTLmap.size() >= 1){
                for(auto it = bTLmap.begin(); it != bTLmap.end(); it++){
//                    cout<<"vblikelihood : "<<it->first<<endl;
                    vblikelihood.push_back(it->first);
                    vbTL.push_back(it->second);
                }
                for(auto it = cTLmap.begin(); it != cTLmap.end(); it++){
                    vclikelihood.push_back(it->first);
                    vcTL.push_back(it->second);
                }
                for(auto it = oTLmap.begin(); it != oTLmap.end(); it++){
                    volikelihood.push_back(it->first);
                    voTL.push_back(it->second);
                }
            }

            TLorentzVector TL2b(0,0,0,0); TLorentzVector TLOther2b(0,0,0,0);
            TLorentzVector TL2c(0,0,0,0); TLorentzVector TLOther2c(0,0,0,0);
            sumb = 0, sum2b = 0, massB1b = 0, massB2b = 0, massRecoil1b = 0, massRecoil2b = 0;
            sumc = 0, sum2c = 0, massB1c = 0, massB2c = 0, massRecoil1c = 0, massRecoil2c = 0;
            for(int i = 0; i<vbTL.size(); i++){
                sumb += vblikelihood.at(i);
                sumc += vclikelihood.at(i);
                if(i<2){
                    sum2b += vblikelihood.at(i);
                    TL2b += vbTL.at(i);
                    sum2c += vclikelihood.at(i);
                    TL2c += vcTL.at(i);
                }
                if(i>=2){TLOther2b += vbTL.at(i); TLOther2c += vcTL.at(i);}
            }
            cout<<"sumb : "<<sumb<<endl;
            massB1b = TL2b.M(); massB2b = TLOther2b.M(); massRecoil1b = (TLT240 - TL2b).M(); massRecoil2b = (TLT240 - TLOther2b).M();
            massB1c = TL2c.M(); massB2c = TLOther2c.M(); massRecoil1c = (TLT240 - TL2c).M(); massRecoil2c = (TLT240 - TLOther2c).M();
            
            //vertex analysis
            LCCollection* col_BuildUpVertexRP = evtP->getCollection( "BuildUpVertex_RP" );
            num_Vertex = 0;
            num_Vertex = col_BuildUpVertexRP->getNumberOfElements();

            std::vector<MCParticle* > allQuark;  allQuark.clear();
            std::vector<MCParticle* > HQuark;    HQuark.clear();
            std::vector<MCParticle* > noHQuark;  noHQuark.clear();
            std::vector<MCParticle* > OriLepton; OriLepton.clear();
            std::vector<MCParticle* > HLepton;   HLepton.clear();


            TLorentzVector TLISR(0,0,0,0);
            count94 = 0;
            mass941 = 0; mass942 = 0;
            boson1Quark.clear(); boson2Quark.clear();
            TLorentzVector TL941(0,0,0,0); TLorentzVector TL942(0,0,0,0);
            std::map<double, float, greater<double> > ISRmap; ISRmap.clear();
            LCCollection* col_MCP = evtP->getCollection( "MCParticle" );
            int n_MCP = col_MCP->getNumberOfElements();
            num_lepton = 0, num_quark = 0, num_neutrino = 0;
            B1Quark[0] = 0; B1Quark[1] = 0; B2Quark[0] = 0; B2Quark[1] = 0;
            for(int i = 0; i<n_MCP; i++)
            {
                MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(i));
                int NParents      = a_MCP->getParents().size();
                int NDaughters    = a_MCP->getDaughters().size();
                int PDG           = a_MCP->getPDG();
                TLorentzVector temp(a_MCP->getMomentum(), a_MCP->getEnergy());
                TVector3 tempV3 = temp.Vect();
                
                if(NParents == 0 && (abs(PDG) == 1 || abs(PDG) == 2 || abs(PDG) == 3 || abs(PDG) == 4 || abs(PDG) == 5 || abs(PDG) == 6 || abs(PDG) == 21   ) )
                {
                    allQuark.push_back(a_MCP);
                    noHQuark.push_back(a_MCP);
                    num_quark += 1;
                }
                if(NParents == 0 && (abs(PDG) == 11 || abs(PDG) == 12 || abs(PDG) == 13 || abs(PDG) == 14 || abs(PDG) == 15 || abs(PDG) == 16) ){
                    num_lepton += 1;
                    OriLepton.push_back( a_MCP );
                    if( abs(PDG) == 12 || abs(PDG) == 14 || abs(PDG) == 16 ){num_neutrino += 1;}
                }
                
                if(PDG == 25 && NDaughters == 2){
                    MCParticle* HDau1 = a_MCP->getDaughters()[0];
                    MCParticle* HDau2 = a_MCP->getDaughters()[1];
                    cout<<"HDau1->getPDG() : "<<HDau1->getPDG()<<" HDau2->getPDG() : "<<HDau2->getPDG()<<endl;
                    if(abs(HDau1->getPDG()) == 1 || abs(HDau1->getPDG()) == 2 || abs(HDau1->getPDG()) == 3 || abs(HDau1->getPDG()) == 4 || abs(HDau1->getPDG()) == 5 || abs(HDau1->getPDG()) == 6 || abs(HDau1->getPDG()) == 21){
                        allQuark.push_back(HDau1);
                        allQuark.push_back(HDau2);
                        HQuark.push_back(HDau1);
                        HQuark.push_back(HDau2);
//                        cout<<"HDau1->getPDG() : "<<HDau1->getPDG()<<" HDau2->getPDG() : "<<HDau2->getPDG()<<endl;
                        num_quark += 2;
                    }
                    if(abs(HDau1->getPDG()) == 11 || abs(HDau1->getPDG()) == 12 || abs(HDau1->getPDG()) == 13 || abs(HDau1->getPDG()) == 14 || abs(HDau1->getPDG()) == 15 || abs(HDau1->getPDG()) == 16 ){
                        num_lepton += 2;
                    }
                    if(abs(HDau1->getPDG()) == 23 || abs(HDau1->getPDG()) == 24){
                        if(HDau1->getDaughters().size() >= 2){
                            MCParticle* wzDau = HDau1->getDaughters()[0];
                            int wzDauPDG = abs(wzDau->getPDG());
                            if(abs(wzDauPDG) == 11 || abs(wzDauPDG) == 12 || abs(wzDauPDG) == 13 || abs(wzDauPDG) == 14 || abs(wzDauPDG) == 15 || abs(wzDauPDG) == 16 ){
                                num_lepton += 2;
                            }
                            if(abs(wzDauPDG) == 1 || abs(wzDauPDG) == 2 || abs(wzDauPDG) == 3 || abs(wzDauPDG) == 4 || abs(wzDauPDG) == 5 || abs(wzDauPDG) == 6 || wzDauPDG == 21 ){
                                num_quark += 2;
                                HQuark.push_back(HDau1->getDaughters()[0]); HQuark.push_back(HDau1->getDaughters()[1]);
                            }
                        }
                        if(HDau2->getDaughters().size() >= 2){
                            MCParticle* wzDau2 = HDau2->getDaughters()[0];
                            int wzDauPDG2 = abs(wzDau2->getPDG());
                            if(abs(wzDauPDG2) == 11 || abs(wzDauPDG2) == 12 || abs(wzDauPDG2) == 13 || abs(wzDauPDG2) == 14 || abs(wzDauPDG2) == 15 || abs(wzDauPDG2) == 16 ){
                                num_lepton += 2;
                            }
                            if(abs(wzDauPDG2) == 1 || abs(wzDauPDG2) == 2 || abs(wzDauPDG2) == 3 || abs(wzDauPDG2) == 4 || abs(wzDauPDG2) == 5 || abs(wzDauPDG2) == 6 || wzDauPDG2 == 21 ){
                                num_quark += 2;
                                HQuark.push_back(HDau2->getDaughters()[0]); HQuark.push_back(HDau2->getDaughters()[1]);
                            }
                        }
                    }
                }
                
                
                if(PDG == 94){
                    count94 += 1;
                    if(count94 == 1){
                        mass941 = a_MCP->getMass();
                        TL941 = temp;
                        MCParticle* par941 = a_MCP->getParents()[0];
                        MCParticle* par942 = a_MCP->getParents()[1];
                        boson1Quark.push_back( abs(par941->getPDG()) );
                        boson1Quark.push_back( abs(par942->getPDG()) );
                        B1Quark[0] = abs(par941->getPDG());
                        B1Quark[1] = abs(par942->getPDG());
                    }
                    if(count94 == 2){
                        mass942 = a_MCP->getMass();
                        TL942 = temp;
                        MCParticle* par941 = a_MCP->getParents()[0];
                        MCParticle* par942 = a_MCP->getParents()[1];
                        boson2Quark.push_back( abs(par941->getPDG()) );
                        boson2Quark.push_back( abs(par942->getPDG()) );
                        B2Quark[0] = abs(par941->getPDG());
                        B2Quark[1] = abs(par942->getPDG());
                    }
                }
                
                if(NParents == 0 && PDG == 22)
                {
                    TLISR += temp;
                    ISRmap[a_MCP->getEnergy()] = abs(tempV3.CosTheta());
                }
                
            }

            
            for(int i = 0; i<5; i++){
                ISREn[i] = 999;
                ISRCosTheta[i] = 999;
            }
            std::vector<double> vISREn;     vISREn.clear();
            std::vector<float>  vISRCTheta; vISRCTheta.clear();
            int count = 0;
            for(auto it = ISRmap.begin(); it != ISRmap.end(); it++){
                ISREn[count] = it->first;
                ISRCosTheta[count] = it->second;
                count += 1;
            }
            ISR_En = 0; ISR_En = TLISR.E();
            
            for(int i = 0; i<4; i++){
                HQPDG[i] = 0;
                noHQPDG[i] = 0;
                OriLepPDG[i] = 0;
            }
            
            if(HQuark.size() != 0){
                for(int i = 0; i<HQuark.size(); i++){
                    MCParticle* a_MCP = HQuark.at(i);
                    HQPDG[i] = abs(a_MCP->getPDG());
                }
            }
            if(noHQuark.size() != 0){
                for(int i = 0; i<noHQuark.size(); i++){
                    MCParticle* a_MCP = noHQuark.at(i);
                    noHQPDG[i] = abs(a_MCP->getPDG());
                }
            }
            if(OriLepton.size() != 0){
                for(int i = 0; i<OriLepton.size(); i++){
                    MCParticle* a_MCP = OriLepton.at(i);
                    OriLepPDG[i] = abs(a_MCP->getPDG());
                }
            }

            

            for(int i = 0; i<2; i++){

                HbbL[i] = 999;
                HccL[i] = 999;
                ZbbL[i] = 999;
                ZccL[i] = 999;
                
                HooL[i] = 999;
                HbcL[i] = 999;


                ZooL[i] = 999;
                ZbcL[i] = 999;

            }
            double wmass = 80.4, zmass = 91.2, hmass = 125;
            std::vector<vector<int> > vect1 = pair4jets(4);
        
            //pair 4 jets according to invariant mass
            pairMB1 = 0; pairMB2 = 0; HiggsJetsAngle = 0; numHiggsParticles = 999; ZJetsAngle = 0; ZHAngle = 0;
            std::vector<double> HiggsbbL; HiggsbbL.clear();
            std::vector<double> HiggsccL; HiggsccL.clear();
            std::vector<double> HiggsooL; HiggsooL.clear();
            std::vector<double> HiggsbcL; HiggsbcL.clear();

            for(int i = 0; i<4; i++){
                HiggsX[i] = -1;
                HiggsSumL[i] = 0;
                numChargeOfJets[i] = 0;
                jetPx[i] = 0; jetPy[i] = 0; jetPz[i] = 0; jetEn[i] = 0;
            }
            maxHiggsX = -1;
            int checkJetsPs = 0;
            if(num_jet == 2){
                PIDHandler pidh(col_Jet);
                int algo   = 0;
                int ibtag  =-1;
                int ictag  =-1;
                int ibctag =-1;
                int icat   =-1;
                algo      = pidh.getAlgorithmID("lcfiplus");
                ibtag     = pidh.getParameterIndex (algo,  "BTag");
                ictag     = pidh.getParameterIndex (algo,  "CTag");
                ibctag    = pidh.getParameterIndex (algo,  "BCTag");
                icat      = pidh.getParameterIndex (algo,  "Category");

                ReconstructedParticle* jet1 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(0));
                jetPx[0] = jet1->getMomentum()[0]; jetPy[0] = jet1->getMomentum()[1]; jetPz[0] = jet1->getMomentum()[2]; jetEn[0] = jet1->getEnergy();
                const ParticleID &pid1 = pidh.getParticleID(jet1, algo);
                double btag1  = pid1.getParameters()[ibtag];
                double ctag1  = pid1.getParameters()[ictag];
                double bctag1 = pid1.getParameters()[ibctag];
                double cat1   = pid1.getParameters()[icat];
                double otag1  = 1 - btag1 - ctag1;
                TLorentzVector TLtemp1(jet1->getMomentum(), jet1->getEnergy());
                TVector3 TVtemp1 = TLtemp1.Vect();
                ReconstructedParticle* jet2 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(1));
                jetPx[1] = jet2->getMomentum()[0]; jetPy[1] = jet2->getMomentum()[1]; jetPz[1] = jet2->getMomentum()[2]; jetEn[1] = jet2->getEnergy();
                const ParticleID &pid2 = pidh.getParticleID(jet2, algo);
                double btag2  = pid2.getParameters()[ibtag];
                double ctag2  = pid2.getParameters()[ictag];
                double bctag2 = pid2.getParameters()[ibctag];
                double cat2   = pid2.getParameters()[icat];
                double otag2  = 1 - btag2 - ctag2;
                HbbL[0] = btag1; HbbL[1] = btag2;
                HccL[0] = ctag1; HccL[1] = ctag2;
                HooL[0] = otag1; HooL[1] = otag2;
                HbcL[0] = bctag1; HbcL[1] = bctag2;
                
                if(btag1 + ctag1 > 1.1){
                    cout<<"***********************************************************************************************************"<<endl;
                    cout<<"***********************************************************************************************************"<<endl;
                    cout<<"***********************************************************************************************************"<<endl;

                    cout<<"btag1 + ctag1 : "<<btag1 + ctag1<<endl;}
                if(btag2 + ctag2 > 1.1){
                    cout<<"***********************************************************************************************************"<<endl;
                    cout<<"***********************************************************************************************************"<<endl;
                    cout<<"***********************************************************************************************************"<<endl;
                    cout<<"btag2 + ctag2 : "<<btag2 + ctag2<<endl;}
                
//                double Xb = HbbL[0]*HbbL[1]/(HbbL[0]*HbbL[1] + (1-HbbL[0])*(1-HbbL[1]));
//                double Xc = HccL[0]*HccL[1]/(HccL[0]*HccL[1] + (1-HccL[0])*(1-HccL[1]));
//                double Xo = HooL[0]*HooL[1]/(HooL[0]*HooL[1] + (1-HooL[0])*(1-HooL[1]));
//                HiggsX[0] = Xb; HiggsX[1] = Xc; HiggsX[2] = Xo;
//                cout<<HbbL0 + HccL[0] + HooL[0]<<endl;
//                cout<<HbbL1 + HccL[1] + HooL[1]<<endl;
            }
            cout<<"checkJetsPs : "<<checkJetsPs<<" nPFO : "<<nPFO<<endl;
                                                                                       
                                                                                       
            HCharge = 0; ZCharge = 0;                                                                                      
            B1Costheta = 999; B2Costheta = 999;                                                       
                                                                                       
            if(num_jet == 4){
                PIDHandler pidh(col_Jet);
                int algo   = 0;
                int ibtag  =-1;
                int ictag  =-1;
                int ibctag =-1;
                int icat   =-1;
                algo      = pidh.getAlgorithmID("lcfiplus");
                ibtag     = pidh.getParameterIndex (algo,  "BTag");
                ictag     = pidh.getParameterIndex (algo,  "CTag");
                ibctag    = pidh.getParameterIndex (algo,  "BCTag");
                icat      = pidh.getParameterIndex (algo,  "Category");
                
                int slcti3 = 999;
                TLorentzVector TLRecoB1(0,0,0,0);
                TLorentzVector TLRecoB2(0,0,0,0);
                RminDif = 9999.0;
                RminDifZH = 9999.0;
                for(int i = 0; i<(vect1.size() - 1); i = i + 2){
                    TLorentzVector TL12(0,0,0,0);
                    TLorentzVector TL34(0,0,0,0);
                    for(int j = 0; j<vect1[i].size(); j++){
                        ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(vect1[i][j]));
                        TLorentzVector TLtemp(jet->getMomentum(), jet->getEnergy());
                        TL12 += TLtemp;
                    }
                    for(int j = 0; j<vect1[i+1].size(); j++){
                        ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(vect1[i+1][j]));
                        TLorentzVector TLtemp(jet->getMomentum(), jet->getEnergy());
                        TL34 += TLtemp;
                    }
                    double compareW = (TL12.M() - wmass)*(TL12.M() - wmass)/14.44 + (TL34.M() - wmass)*(TL34.M() - wmass)/14.44;
                    double compareZ = (TL12.M() - zmass)*(TL12.M() - zmass)/19.36 + (TL34.M() - zmass)*(TL34.M() - zmass)/19.36;
                    double compareZH= (TL12.M() - zmass)*(TL12.M() - zmass)/19.36 + (TL34.M() - hmass)*(TL34.M() - hmass)/25;
                    double compareHZ= (TL12.M() - hmass)*(TL12.M() - hmass)/25    + (TL34.M() - zmass)*(TL34.M() - zmass)/19.36;
                    if(compareW  < RminDif){RminDif = compareW;  TLRecoB1 = TL12; TLRecoB2 = TL34; slcti3 = i;}
                    if(compareZ  < RminDif){RminDif = compareZ;  TLRecoB1 = TL12; TLRecoB2 = TL34; slcti3 = i;}
                    if(compareZH < RminDif){RminDif = compareZH; TLRecoB1 = TL12; TLRecoB2 = TL34; slcti3 = i;}
                    if(compareHZ < RminDif){RminDif = compareHZ; TLRecoB1 = TL12; TLRecoB2 = TL34; slcti3 = i;}
                    cout<<"RminDif : "<<RminDif<<" TLRecoB1.M() : "<<TLRecoB1.M()<<" TLRecoB2.M() : "<<TLRecoB2.M()<<endl;
                }
                if(TLRecoB1.M() > TLRecoB2.M()){
                    RminDifZH = (TLRecoB2.M() - zmass)*(TLRecoB2.M() - zmass)/19.36 + (TLRecoB1.M() - hmass)*(TLRecoB1.M() - hmass)/25;
                    int jet1Charge = 0, jet2Charge = 0, jet3Charge = 0, jet4Charge = 0;
                    pairMB1 = TLRecoB1.M();   pairMB2 = TLRecoB2.M();
                    ReconstructedParticle* jet1 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(vect1[slcti3][0]));
                    jetPx[0] = jet1->getMomentum()[0]; jetPy[0] = jet1->getMomentum()[1]; jetPz[0] = jet1->getMomentum()[2]; jetEn[0] = jet1->getEnergy();
                    const ParticleID &pid1 = pidh.getParticleID(jet1, algo);
                    double btag1  = pid1.getParameters()[ibtag];
                    double ctag1  = pid1.getParameters()[ictag];
                    double bctag1 = pid1.getParameters()[ibctag];
                    double cat1   = pid1.getParameters()[icat];
                    double otag1  = 1 - btag1 - ctag1;
                    TLorentzVector TLtemp1(jet1->getMomentum(), jet1->getEnergy());
                    TVector3 TVtemp1 = TLtemp1.Vect();
                    ReconstructedParticle* jet2 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(vect1[slcti3][1]));
                    jetPx[1] = jet2->getMomentum()[0]; jetPy[1] = jet2->getMomentum()[1]; jetPz[1] = jet2->getMomentum()[2]; jetEn[1] = jet2->getEnergy();
                    const ParticleID &pid2 = pidh.getParticleID(jet2, algo);
                    double btag2  = pid2.getParameters()[ibtag];
                    double ctag2  = pid2.getParameters()[ictag];
                    double bctag2 = pid2.getParameters()[ibctag];
                    double cat2   = pid2.getParameters()[icat];
                    double otag2  = 1 - btag2 - ctag2;
                    HbbL[0] = btag1; HbbL[1] = btag2;
                    HccL[0] = ctag1; HccL[1] = ctag2;
                    HooL[0] = otag1; HooL[1] = otag2;
                    HbcL[0] = bctag1; HbcL[1] = bctag2;

                    if(btag1 + ctag1 > 1.1){
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;

                        cout<<"btag1 + ctag1 : "<<btag1 + ctag1<<endl;}
                    if(btag2 + ctag2 > 1.1){
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"btag2 + ctag2 : "<<btag2 + ctag2<<endl;}
                    
                    
                    if(HbbL[0] + HccL[0] > 1.1){
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;

                        cout<<"HbbL[0] + HccL[0] : "<<HbbL[0] + HccL[0]<<endl;}
                    if(HbbL[1] + HccL[1] > 1.1){
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"HbbL[1] + HccL[1] : "<<HbbL[1] + HccL[1]<<endl;}
                    
                    
                    TLorentzVector TLtemp2(jet2->getMomentum(), jet2->getEnergy());
                    TVector3 TVtemp2 = TLtemp2.Vect();
                    HiggsJetsAngle = TVtemp1.Angle( TVtemp2 );
                    ReconstructedParticleVec components1 = jet1->getParticles();
                    ReconstructedParticleVec components2 = jet2->getParticles();
                    for(int i = 0; i<components1.size(); i++){
                        ReconstructedParticle* compi = components1.at(i);
                        if(compi->getCharge() != 0){jet1Charge += 1; HCharge += compi->getCharge(); }
                    }
                    for(int i = 0; i<components2.size(); i++){
                        ReconstructedParticle* compi = components2.at(i);
                        if(compi->getCharge() != 0){jet2Charge += 1; HCharge += compi->getCharge(); }
                    }
                    numHiggsParticles = components1.size() + components2.size();
                    double Xb = btag1*btag2/(btag1*btag2 + (1-btag1)*(1-btag2));
                    double Xc = ctag1*ctag2/(ctag1*ctag2 + (1-ctag1)*(1-ctag2));
                    double Xo = otag1*otag2/(otag1*otag2 + (1-otag1)*(1-otag2));
                    double Xbc = bctag1*bctag2/(bctag1*bctag2 + (1-bctag1)*(1-bctag2));
                    HiggsbbL.push_back(btag1); HiggsbbL.push_back(btag2); HiggsccL.push_back(ctag1); HiggsccL.push_back(ctag2); HiggsooL.push_back(otag1); HiggsooL.push_back(otag2); HiggsbcL.push_back(bctag1); HiggsbcL.push_back(bctag2);
                    HiggsX[0] = Xb; HiggsX[1] = Xc; HiggsX[2] = Xo; HiggsX[3] = Xbc;
                    HiggsSumL[0] = btag1 + btag2; HiggsSumL[1] = ctag1 + ctag2; HiggsSumL[2] = otag1 + otag2; HiggsSumL[3] = bctag1 + bctag2;
                    ReconstructedParticle* jet3 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(vect1[slcti3+1][0]));
                    const ParticleID &pid3 = pidh.getParticleID(jet3, algo);
                    double btag3  = pid3.getParameters()[ibtag];
                    double ctag3  = pid3.getParameters()[ictag];
                    double bctag3 = pid3.getParameters()[ibctag];
                    double cat3   = pid3.getParameters()[icat];
                    double otag3  = 1 - btag3 - ctag3;
                    jetPx[2] = jet3->getMomentum()[0]; jetPy[2] = jet3->getMomentum()[1]; jetPz[2] = jet3->getMomentum()[2]; jetEn[2] = jet3->getEnergy();
                    TLorentzVector TLtemp3(jet3->getMomentum(), jet3->getEnergy());
                    ReconstructedParticle* jet4 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(vect1[slcti3+1][1]));
                    const ParticleID &pid4 = pidh.getParticleID(jet4, algo);
                    double btag4  = pid4.getParameters()[ibtag];
                    double ctag4  = pid4.getParameters()[ictag];
                    double bctag4 = pid4.getParameters()[ibctag];
                    double cat4   = pid4.getParameters()[icat];
                    double otag4  = 1 - btag4 - ctag4;
                    ZbbL[0] = btag3; ZbbL[1] = btag4;
                    ZccL[0] = ctag3; ZccL[1] = ctag4;
                    ZooL[0] = otag3; ZooL[1] = otag4;
                    ZbcL[0] = bctag3; ZbcL[1] = bctag4;
                    
                    
                    if(ZbbL[0] + ZccL[0] > 1.1){
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;

                        cout<<"ZbbL[0] + ZccL[0] : "<<ZbbL[0] + ZccL[0]<<endl;}
                    if(ZbbL[1] + ZccL[1] > 1.1){
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"ZbbL[1] + ZccL[1] : "<<ZbbL[1] + ZccL[1]<<endl;}
                    
                    
                    jetPx[3] = jet4->getMomentum()[0]; jetPy[3] = jet4->getMomentum()[1]; jetPz[3] = jet4->getMomentum()[2]; jetEn[3] = jet4->getEnergy();
                    TLorentzVector TLtemp4(jet4->getMomentum(), jet4->getEnergy());
                    ZJetsAngle = (TLtemp3.Vect()).Angle( TLtemp4.Vect() );
                    ReconstructedParticleVec components3 = jet3->getParticles();
                    ReconstructedParticleVec components4 = jet4->getParticles();
                    for(int i = 0; i<components3.size(); i++){
                        ReconstructedParticle* compi = components3.at(i);
                        if(compi->getCharge() != 0){jet3Charge += 1; ZCharge += compi->getCharge();}
                    }
                    for(int i = 0; i<components4.size(); i++){
                        ReconstructedParticle* compi = components4.at(i);
                        if(compi->getCharge() != 0){jet4Charge += 1; ZCharge += compi->getCharge();}
                    }
                    numChargeOfJets[0] = jet1Charge; numChargeOfJets[1] = jet2Charge; numChargeOfJets[2] = jet3Charge; numChargeOfJets[3] = jet4Charge;
                }
                else{
                    RminDifZH = (TLRecoB1.M() - zmass)*(TLRecoB1.M() - zmass)/19.36 + (TLRecoB2.M() - hmass)*(TLRecoB2.M() - hmass)/25;
                    int jet1Charge = 0, jet2Charge = 0, jet3Charge = 0, jet4Charge = 0;
                    pairMB2 = TLRecoB1.M();   pairMB1 = TLRecoB2.M();
                    ReconstructedParticle* jet1 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(vect1[slcti3+1][0]));
                    jetPx[0] = jet1->getMomentum()[0]; jetPy[0] = jet1->getMomentum()[1]; jetPz[0] = jet1->getMomentum()[2]; jetEn[0] = jet1->getEnergy();
                    const ParticleID &pid1 = pidh.getParticleID(jet1, algo);
                    double btag1  = pid1.getParameters()[ibtag];
                    double ctag1  = pid1.getParameters()[ictag];
                    double bctag1 = pid1.getParameters()[ibctag];
                    double cat1   = pid1.getParameters()[icat];
                    double otag1 = 1 - btag1 - ctag1;
                    TLorentzVector TLtemp1(jet1->getMomentum(), jet1->getEnergy());
                    TVector3 TVtemp1 = TLtemp1.Vect();
                    ReconstructedParticle* jet2 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(vect1[slcti3+1][1]));
                    jetPx[1] = jet2->getMomentum()[0]; jetPy[1] = jet2->getMomentum()[1]; jetPz[1] = jet2->getMomentum()[2]; jetEn[1] = jet2->getEnergy();
                    const ParticleID &pid2 = pidh.getParticleID(jet2, algo);
                    double btag2  = pid2.getParameters()[ibtag];
                    double ctag2  = pid2.getParameters()[ictag];
                    double bctag2 = pid2.getParameters()[ibctag];
                    double cat2   = pid2.getParameters()[icat];
                    double otag2 = 1 - btag2 - ctag2;
                    HbbL[0] = btag1; HbbL[1] = btag2;
                    HccL[0] = ctag1; HccL[1] = ctag2;
                    HooL[0] = otag1; HooL[1] = otag2;
                    HbcL[0] = bctag1; HbcL[1] = bctag2;
                    
                    if(HbbL[0] + HccL[0] > 1.1){
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;

                        cout<<"HbbL[0] + HccL[0] : "<<HbbL[0] + HccL[0]<<endl;}
                    if(HbbL[1] + HccL[1] > 1.1){
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"HbbL[1] + HccL[1] : "<<HbbL[1] + HccL[1]<<endl;}
                    
                    TLorentzVector TLtemp2(jet2->getMomentum(), jet2->getEnergy());
                    TVector3 TVtemp2 = TLtemp2.Vect();
                    HiggsJetsAngle = TVtemp1.Angle( TVtemp2 );
                    ReconstructedParticleVec components1 = jet1->getParticles();
                    ReconstructedParticleVec components2 = jet2->getParticles();
                    numHiggsParticles = components1.size() + components2.size();
                    for(int i = 0; i<components1.size(); i++){
                        ReconstructedParticle* compi = components1.at(i);
                        if(compi->getCharge() != 0){jet1Charge += 1; HCharge += compi->getCharge();}
                    }
                    for(int i = 0; i<components2.size(); i++){
                        ReconstructedParticle* compi = components2.at(i);
                        if(compi->getCharge() != 0){jet2Charge += 1; HCharge += compi->getCharge();}
                    }
                    double Xb = btag1*btag2/(btag1*btag2 + (1-btag1)*(1-btag2));
                    double Xc = ctag1*ctag2/(ctag1*ctag2 + (1-ctag1)*(1-ctag2));
                    double Xo = otag1*otag2/(otag1*otag2 + (1-otag1)*(1-otag2));
                    double Xbc = bctag1*bctag2/(bctag1*bctag2 + (1-bctag1)*(1-bctag2));
                    HiggsbbL.push_back(btag1); HiggsbbL.push_back(btag2); HiggsccL.push_back(ctag1); HiggsccL.push_back(ctag2); HiggsooL.push_back(otag1); HiggsooL.push_back(otag2); HiggsbcL.push_back(bctag1); HiggsbcL.push_back(bctag2);
                    HiggsX[0] = Xb; HiggsX[1] = Xc; HiggsX[2] = Xo; HiggsX[3] = Xbc;
                    HiggsSumL[0] = btag1 + btag2; HiggsSumL[1] = ctag1 + ctag2; HiggsSumL[2] = otag1 + otag2; HiggsSumL[3] = bctag1 + bctag2;
                    ReconstructedParticle* jet3 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(vect1[slcti3][0]));
                    const ParticleID &pid3 = pidh.getParticleID(jet3, algo);
                    double btag3  = pid3.getParameters()[ibtag];
                    double ctag3  = pid3.getParameters()[ictag];
                    double bctag3 = pid3.getParameters()[ibctag];
                    double cat3   = pid3.getParameters()[icat];
                    double otag3  = 1 - btag3 - ctag3;
                    jetPx[2] = jet3->getMomentum()[0]; jetPy[2] = jet3->getMomentum()[1]; jetPz[2] = jet3->getMomentum()[2]; jetEn[2] = jet3->getEnergy();
                    TLorentzVector TLtemp3(jet3->getMomentum(), jet3->getEnergy());
                    ReconstructedParticle* jet4 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt(vect1[slcti3][1]));
                    const ParticleID &pid4 = pidh.getParticleID(jet4, algo);
                    double btag4  = pid4.getParameters()[ibtag];
                    double ctag4  = pid4.getParameters()[ictag];
                    double bctag4 = pid4.getParameters()[ibctag];
                    double cat4   = pid4.getParameters()[icat];
                    double otag4  = 1 - btag4 - ctag4;
                    ZbbL[0] = btag3; ZbbL[1] = btag4;
                    ZccL[0] = ctag3; ZccL[1] = ctag4;
                    ZooL[0] = otag3; ZooL[1] = otag4;
                    ZbcL[0] = bctag3; ZbcL[1] = bctag4;
                    
                    if(ZbbL[0] + ZccL[0] > 1.1){
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;

                        cout<<"ZbbL[0] + ZccL[0] : "<<ZbbL[0] + ZccL[0]<<endl;}
                    if(ZbbL[1] + ZccL[1] > 1.1){
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"***********************************************************************************************************"<<endl;
                        cout<<"ZbbL[1] + ZccL[1] : "<<ZbbL[1] + ZccL[1]<<endl;}
                    
                    jetPx[3] = jet4->getMomentum()[0]; jetPy[3] = jet4->getMomentum()[1]; jetPz[3] = jet4->getMomentum()[2]; jetEn[3] = jet4->getEnergy();
                    TLorentzVector TLtemp4(jet4->getMomentum(), jet4->getEnergy());
                    ZJetsAngle = (TLtemp3.Vect()).Angle( TLtemp4.Vect() );
                    ReconstructedParticleVec components3 = jet3->getParticles();
                    ReconstructedParticleVec components4 = jet4->getParticles();
                    for(int i = 0; i<components3.size(); i++){
                        ReconstructedParticle* compi = components3.at(i);
                        if(compi->getCharge() != 0){jet3Charge += 1; ZCharge += compi->getCharge();}
                    }
                    for(int i = 0; i<components4.size(); i++){
                        ReconstructedParticle* compi = components4.at(i);
                        if(compi->getCharge() != 0){jet4Charge += 1; ZCharge += compi->getCharge();}
                    }
                    numChargeOfJets[0] = jet1Charge; numChargeOfJets[1] = jet2Charge; numChargeOfJets[2] = jet3Charge; numChargeOfJets[3] = jet4Charge;
                }
                for(int i = 0; i<4; i++){
                    if(HiggsX[i] > maxHiggsX){maxHiggsX = HiggsX[i];}
                }
//                sort(HiggsbbL.begin(), HiggsbbL.end());
//                sort(HiggsccL.begin(), HiggsccL.end());
//                sort(HiggsooL.begin(), HiggsooL.end());
//                sort(HiggsbcL.begin(), HiggsbcL.end());
//                cout<<"HiggsbbL.at(0) : "<<HiggsbbL.at(0)<<" : "<<HiggsbbL.at(1)<<endl;
//                cout<<"HiggsccL.at(0) : "<<HiggsccL.at(0)<<" : "<<HiggsccL.at(1)<<endl;
//                cout<<"HiggsooL.at(0) : "<<HiggsooL.at(0)<<" : "<<HiggsooL.at(1)<<endl;
//                HbbL[0] = HiggsbbL.at(0); HbbL[1] = HiggsbbL.at(1);
//                HccL[0] = HiggsccL.at(0); HccL[1] = HiggsccL.at(1);
//                HooL[0] = HiggsooL.at(0); HooL[1] = HiggsooL.at(1);
//                HbcL[0] = HiggsbcL.at(0); HbcL[1] = HiggsbcL.at(1);

                cout<<"HiggsJetsAngle : "<<HiggsJetsAngle<<" numHiggsParticles : "<<numHiggsParticles<<endl;
                cout<<"pairMB1 : "<<pairMB1<<" pairMB2 : "<<pairMB2<<endl;
                TVector3 B1TV = TL941.Vect();          TVector3 B2TV = TL942.Vect();
                TVector3 R1TV = TLRecoB1.Vect();       TVector3 R2TV = TLRecoB2.Vect();
                ZHAngle = R1TV.Angle(R2TV);

		B1Costheta = R1TV.CosTheta();
		B2Costheta = R2TV.CosTheta();

                float compare1 = B1TV.Angle(R1TV) + B2TV.Angle(R2TV);
                float compare2 = B1TV.Angle(R2TV) + B2TV.Angle(R1TV);
                cout<<"compare1 : "<<compare1<<" compare2 : "<<compare2<<endl;
                cout<<"ZHAngle : "<<ZHAngle<<endl;
                delta_R1 = 999, delta_R2 = 999;
                if(compare1 <= compare2){
                    delta_R1 = B1TV.Angle(R1TV);
                    delta_R2 = B2TV.Angle(R2TV);
                }
                else {
                    delta_R1 = B1TV.Angle(R2TV);
                    delta_R2 = B2TV.Angle(R1TV);
                }
                cout<<"delta_R1 : "<<delta_R1<<" delta_R2 : "<<delta_R2<<endl;
                alpha = log10(delta_R1 * delta_R2);
            }

            

           
              /*
            double miniDif = 0, miniDif2 = 0;
            int slcti = 999, slctj = 999, slctm = 999, slcti2 = 999, slctj2 = 999, slctm2 = 999;
            double sumLLj = 0, sumLLm = 0, Xj = 0, Xm = 0;
            P1HBMass = 999; P1LBMass = 999; P1HF = 999; P1LF = 999; P1HLL = 999; P1LLL = 999;
            P2HBMass = 999; P2LBMass = 999; P2HF = 999; P2LF = 999; P2HLL = 999; P2LLL = 999;
            if(num_jet == 4){
                for(int i = 0; i<(vect1.size() - 1); i = i+2){
                    std::array<double, 3> prob0 = mapProb[ vect1[i][0] ];
                    std::array<double, 3> prob1 = mapProb[ vect1[i][1] ];
                    std::array<double, 3> prob2 = mapProb[ vect1[i+1][0] ];
                    std::array<double, 3> prob3 = mapProb[ vect1[i+1][1] ];
                    for(int j = 0; j<3; j++){
                        for(int m = 0; m<3; m++){
//                            double dif = abs( 1.0/(prob0[j] + prob1[j]) + 1.0/(prob2[m] + prob3[m]) - 1 );
                            double dif = prob0[j] + prob1[j] + prob2[m] + prob3[m];
                            double dif2 = prob0[j]*prob1[j]/(prob0[j]*prob1[j] + (1-prob0[j])*(1-prob1[j])) + prob2[j]*prob3[j]/(prob2[j]*prob3[j] + (1-prob2[j])*(1-prob3[j]));
//                            if(dif < miniDif){miniDif = dif; slcti= i; slctj = j; slctm = m; sumLLj = prob0[j] + prob1[j]; sumLLm = prob2[m] + prob3[m];}
                            if(dif > miniDif){miniDif = dif; slcti= i; slctj = j; slctm = m; sumLLj = prob0[j] + prob1[j]; sumLLm = prob2[m] + prob3[m];}
//                            cout<<"dif : "<<dif<<" miniDif : "<<miniDif<<endl;
                            if(dif2 > miniDif2){miniDif2 = dif2; slcti2 = i; slctj2 = j; slctm2 = m;
                                Xj = prob0[j]*prob1[j]/(prob0[j]*prob1[j] + (1-prob0[j])*(1-prob1[j]));
                                Xm = prob2[j]*prob3[j]/(prob2[j]*prob3[j] + (1-prob2[j])*(1-prob3[j]));
                            }
//                            cout<<"dif2 : "<<dif2<<" miniDif2 : "<<miniDif2<<endl;
                        }
                    }
                }
//                cout<<" miniDif : "<<miniDif<<endl;
//                cout<<" miniDif2 : "<<miniDif2<<endl;

                ReconstructedParticle* jet0 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt( vect1[slcti][0] ));
                TLorentzVector TLtemp0(jet0->getMomentum(), jet0->getEnergy());
                ReconstructedParticle* jet1 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt( vect1[slcti][1] ));
                TLorentzVector TLtemp1(jet1->getMomentum(), jet1->getEnergy());
                double newMass1 = (TLtemp0 + TLtemp1).M();
                ReconstructedParticle* jet2 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt( vect1[slcti+1][0] ));
                TLorentzVector TLtemp2(jet2->getMomentum(), jet2->getEnergy());
                ReconstructedParticle* jet3 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt( vect1[slcti+1][1] ));
                TLorentzVector TLtemp3(jet3->getMomentum(), jet3->getEnergy());
                double newMass2 = (TLtemp2 + TLtemp3).M();
                if(newMass1 >= newMass2){ P1HBMass = newMass1; P1LBMass = newMass2; P1HF = slctj; P1LF = slctm; P1HLL = sumLLj; P1LLL = sumLLm;}
                else{P1HBMass = newMass2; P1LBMass = newMass1; P1HF = slctm; P1LF = slctj; P1HLL = sumLLm; P1LLL = sumLLj;}
                cout<<"P1HBMass : "<<P1HBMass<<" P1LBMass : "<<P1LBMass<<" P1HF : "<<P1HF<<" P1LF : "<<P1LF<<" P1HLL : "<<P1HLL<<" P1LLL : "<<P1LLL<<endl;
                
                ReconstructedParticle* jet02 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt( vect1[slcti2][0] ));
                TLorentzVector TLtemp02(jet02->getMomentum(), jet02->getEnergy());
                ReconstructedParticle* jet12 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt( vect1[slcti2][1] ));
                TLorentzVector TLtemp12(jet12->getMomentum(), jet12->getEnergy());
                double newMass12 = (TLtemp02 + TLtemp12).M();
                ReconstructedParticle* jet22 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt( vect1[slcti2+1][0] ));
                TLorentzVector TLtemp22(jet22->getMomentum(), jet22->getEnergy());
                ReconstructedParticle* jet32 = dynamic_cast<ReconstructedParticle*>(col_Jet->getElementAt( vect1[slcti2+1][1] ));
                TLorentzVector TLtemp32(jet32->getMomentum(), jet32->getEnergy());
                double newMass22 = (TLtemp22 + TLtemp32).M();
                if(newMass12 >= newMass22){ P2HBMass = newMass12; P2LBMass = newMass22; P2HF = slctj2; P2LF = slctm2; P2HLL = Xj; P2LLL = Xm;}
                else{P2HBMass = newMass22; P2LBMass = newMass12; P2HF = slctm2; P2LF = slctj2; P2HLL = Xm; P2LLL = Xj;}
                cout<<"P2HBMass : "<<P2HBMass<<" P2LBMass : "<<P2LBMass<<" P2HF : "<<P2HF<<" P2LF : "<<P2LF<<" P2HLL : "<<P2HLL<<" P2LLL : "<<P2LLL<<endl;
            }
*/
            
            
            

        }catch (lcio::DataNotAvailableException err) {  }

    }
    
    _outputTree->Fill();
    Num ++;
}



void Higgsbb::end()
{
    
    if (_outputTree) {
        
        TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
        //tree_file->cd();
        tree_file->Write();
        delete tree_file;
        //tree_file->Close();
    }
    
}


ReconstructedParticle* getMulti( std::vector<ReconstructedParticle* > vec, std::vector<double> &result ){
    ReconstructedParticle* leadEnPFO;
    result.clear();
    int num = 0;
    num = vec.size();
    double totalEn = 0;
    double LeadEn = 0;
    std::vector<double> energySort; energySort.clear();
    if(vec.size() != 0){
        for(int i = 0; i<vec.size(); i++){
            ReconstructedParticle* pfo = vec.at(i);
            totalEn += pfo->getEnergy();
            energySort.push_back(pfo->getEnergy());
            if(pfo->getEnergy() > LeadEn){LeadEn = pfo->getEnergy(); leadEnPFO = pfo;}
        }
    }
    double avEn = 0;
    if(num != 0){avEn = totalEn / num;}
    result.push_back( num );
    result.push_back( LeadEn );
    result.push_back( avEn );
    cout<<"LeadEn : "<<LeadEn<<endl;
    double subleadEn = 0;
    sort(energySort.begin(), energySort.end());
    if(energySort.size() >= 2){
        int count = energySort.size();
        cout<<""<<energySort.at( count - 1 )<<" "<<energySort.at( count - 2 )<<endl;
        subleadEn = energySort.at( count - 2 );
    }
    result.push_back(subleadEn);
    return leadEnPFO;
}

void CalcuThrust(std::vector<TLorentzVector > UsedForThrust, std::vector<double> &result){
    result.clear();
    double T = 0;
    //the following code used to find the thrust
    double thetaMin = TMath::Pi(), phiMin = 2*TMath::Pi();
    double thetaMin2 = 0, phiMin2 = 0;
    double thetaL = 0, phiL = 0, thetaR = TMath::Pi(), phiR = 2*TMath::Pi();
    int iter = 0;
    double Told = 0;
    double Tnew = 0;
    double cut = 1;
    double thetaRange = 0, phiRange = 0;
    do{
        iter += 1;
//        cout<<"iter : "<<iter<<endl;
        if(iter == 1){
            thetaRange = thetaR - thetaL, phiRange = phiR - phiL;
        }
        else if(iter != 1){
            
            thetaRange = 0.1*(thetaR - thetaL);
            phiRange = 0.1*(phiR - phiL);
            
            thetaL =  thetaMin - thetaRange;
            thetaR = thetaMin + thetaRange;
            phiL = phiMin - phiRange;
            phiR = phiMin + phiRange;
            thetaRange = thetaR - thetaL, phiRange = phiR - phiL;
            
//            cout<<"thetaL : "<<thetaL<<" thetaR : "<<thetaR<<endl;
//            cout<<"phiL : "<<phiL<<" phiR : "<<phiR<<endl;
        }
//        cout<<"thetaRange : "<<thetaRange<<" phiRange : "<<phiRange<<endl;
        for(double theta = thetaL; theta <= thetaR; theta += 0.1*thetaRange){   //in this round, find the max T
            for(double phi = phiL; phi <= phiR; phi += 0.1*phiRange){
                
                double x = sin(theta)*cos(phi);
                double y = sin(theta)*sin(phi);
                double z = cos(theta);
                
                double denominator = 0;
                double numerator = 0;
                for(int i = 0; i<UsedForThrust.size(); i++){
                    TLorentzVector TLtemp = UsedForThrust.at(i);
//                    TLorentzVector TLtemp(temp->getMomentum(), temp->getEnergy());
                    TVector3 TVtemp = TLtemp.Vect();
                    denominator += TVtemp.Mag();
                    numerator += abs(x*TVtemp(0) + y*TVtemp(1) + z*TVtemp(2));
                }
                double Ttemp = numerator/denominator;
                if(Ttemp > T){
                    thetaMin = theta;   phiMin = phi; T = Ttemp;
 //                   cout<<"*************"<<endl;
 //                   cout<<"T : "<<T<<"thetaMin : phiMin "<<thetaMin<<" : "<<phiMin<<endl;
 //                   cout<<"*************"<<endl;
                }
            }
        }
        if(iter == 1){Told = T; Tnew = T;}
        else if(T >= Tnew && iter != 1){
            Told = Tnew; Tnew = T; cut = (Tnew - Told)/Tnew;
        }
//        cout<<"cut : "<<cut<<endl;
    }
    while(cut >= 0.2);
    
//    result[3] = {T, phiMin, thetaMin};
    result.push_back(T);

    TVector3 tempThrust(0,0,0);
    tempThrust.SetXYZ(sin(thetaMin)*cos(phiMin), sin(thetaMin)*sin(phiMin), cos(thetaMin));

    //the following code used to get Hemisphere masses
    std::vector<TLorentzVector > hemisphere1;
    std::vector<TLorentzVector > hemisphere2;
    hemisphere1.clear(); hemisphere2.clear();
    double visEn = 0;
    double JetBroadeningDenominator = 0;
    TVector3 TVtotal(0,0,0);
    for(int i = 0; i<UsedForThrust.size(); i++){
        TLorentzVector TLtemp = UsedForThrust.at(i);
//        TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
        TVector3 TVtemp = TLtemp.Vect();
        TVtotal += TVtemp;
        if(TVtemp.Angle(tempThrust) > 0.5*TMath::Pi()){hemisphere1.push_back(TLtemp);}
        else {hemisphere2.push_back(TLtemp);}
        visEn += TLtemp.E();
        JetBroadeningDenominator += TVtemp.Mag();
    }
    double hemi1En = 0, hemi1Mass = 0, hemi2En = 0, hemi2Mass = 0;
    TLorentzVector TLsphere1(0,0,0,0);
    TLorentzVector TLsphere2(0,0,0,0);
    for(int i = 0; i<hemisphere1.size(); i++){
        TLsphere1 += hemisphere1.at(i);
    }
    for(int i = 0; i<hemisphere2.size(); i++){
        TLsphere2 += hemisphere2.at(i);
    }
    hemi1En = TLsphere1.E();   hemi1Mass = TLsphere1.M();
    hemi2En = TLsphere2.E();   hemi2Mass = TLsphere2.M();
    cout<<"hemi1En : "<<hemi1En<<" hemi2En : "<<hemi2En<<endl;
    cout<<"hemi1Mass : "<<hemi1Mass<<" hemi2Mass : "<<hemi2Mass<<endl;
//    cout<<"the number of particles in two hemispheres is "<<hemisphere1.size()+hemisphere2.size()<<endl;
    
    
    double JetBroadeningNumerator1 = 0, JetBroadeningNumerator2 =0;
    TLorentzVector TLHemi1(0,0,0,0);
    TLorentzVector TLHemi2(0,0,0,0);
    for(int i = 0; i<hemisphere1.size(); i++){
        TLorentzVector TLtemp = hemisphere1.at(i);
//        TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
        TVector3 TVtemp = TLtemp.Vect();
        TLHemi1 += TLtemp;
        JetBroadeningNumerator1 += abs(TVtemp.Mag() * tempThrust.Mag() * sin(TVtemp.Angle(tempThrust)));
        
    }
    for(int i = 0; i<hemisphere2.size(); i++){
        TLorentzVector TLtemp = hemisphere2.at(i);
//        TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
        TVector3 TVtemp = TLtemp.Vect();
        TLHemi2 += TLtemp;
        JetBroadeningNumerator2 += abs(TVtemp.Mag() * tempThrust.Mag() * sin(TVtemp.Angle(tempThrust)));
    }
    double hemiMass1 = 0, hemiMass2 = 0, hemiBroadening1 = 0, hemiBroadening2 = 0;
    if(visEn != 0){
        hemiMass1 = (TLHemi1.M())*(TLHemi1.M())/(visEn*visEn);
        hemiMass2 = (TLHemi2.M())*(TLHemi2.M())/(visEn*visEn);
    }
    if(JetBroadeningDenominator != 0){
        hemiBroadening1 = JetBroadeningNumerator1/(2*JetBroadeningDenominator);
        hemiBroadening2 = JetBroadeningNumerator2/(2*JetBroadeningDenominator);
    }

//    cout<<"hemiMass1 : "<<hemiMass1<<" hemiMass2 : "<<hemiMass2<<endl;
//    cout<<"hemiBroadening1 : "<<hemiBroadening1<<" hemiBroadening2 : "<<hemiBroadening2<<endl;
    result.push_back(hemiMass1);
    result.push_back(hemiMass2);
    result.push_back(hemiBroadening1);
    result.push_back(hemiBroadening2);
    
    
    
    //the following code used to get rapidity
    TVector3 TVParaThrust = tempThrust.Orthogonal();
    double transverseM = TVtotal.Perp(TVParaThrust);
    double rapidity = 0.5*log((visEn + transverseM)/(visEn - transverseM));
//    result.push_back(rapidity);
    result.push_back( tempThrust.CosTheta() );
    result.push_back( hemi1En ); result.push_back(hemi2En); result.push_back(hemi1Mass); result.push_back(hemi2Mass);
    
    cout<<"visEn : "<<visEn<<endl;
    
//    return 0;
}

double CalcuSphericity(std::vector<TLorentzVector> UsedForSphericity, std::vector<double> &result){
    result.clear();
    //the following code used to calculate sphericity
    double CPara = 0;
    double p2Min = 1e-20;
    double denom = 0.;
    double tt[4][4];
    double eVal1 = 0, eVal2 = 0, eVal3 = 0;
    for(int j = 1; j<4; ++j){
        for(int k = j; k<4; ++k){
            tt[j][k] = 0.;
        }
    }
    double newCPara = 0;
    double newdenom = 0;
    double theta11 = 0, theta22 = 0, theta33 = 0, theta12 = 0, theta23 = 0, theta13 = 0;
    for(int i = 0; i<UsedForSphericity.size(); i++){
        TLorentzVector TLtemp = UsedForSphericity.at(i);
//        TLorentzVector TLtemp(a_MCP->getMomentum(), a_MCP->getEnergy());
        TVector3 TVtemp = TLtemp.Vect();
        double pNow[4];
        pNow[1] = TLtemp.Px();
        pNow[2] = TLtemp.Py();
        pNow[3] = TLtemp.Pz();
        double p2Now = pow(pNow[1],2) + pow(pNow[2],2) + pow(pNow[3],2);
        double pWeight = 1./sqrt(max(p2Now, p2Min));
        for(int j = 1; j<4; ++j){
            for(int k=j; k<4; ++k){
                tt[j][k] += pWeight * pNow[j] * pNow[k];
                denom += pWeight * p2Now;
            }
        }
        newdenom += sqrt(max(p2Now, p2Min));
        theta11 += pNow[1]*pNow[1]/sqrt(max(p2Now, p2Min));
        theta22 += pNow[2]*pNow[2]/sqrt(max(p2Now, p2Min));
        theta33 += pNow[3]*pNow[3]/sqrt(max(p2Now, p2Min));
        theta12 += pNow[1]*pNow[2]/sqrt(max(p2Now, p2Min));
        theta13 += pNow[1]*pNow[3]/sqrt(max(p2Now, p2Min));
        theta23 += pNow[2]*pNow[3]/sqrt(max(p2Now, p2Min));
    }
    //calculate C parameter according to the definition from EERAD3
    if(newdenom != 0){
        newCPara = 3*(theta11*theta22 + theta22*theta33 + theta33*theta11 - theta12*theta12 - theta23*theta23 - theta13*theta13)/pow(newdenom, 2);
    }
    
    
    //Normalize tensor to trace = 1.
    for(int j = 1; j < 4; ++j){
        for(int k = j; k < 4; ++k){
            tt[j][k] /= denom;
        }
    }
    //find eigenvalues to matrix (third degree equation)
    double qCoef = ( tt[1][1] * tt[2][2] + tt[1][1] * tt[3][3]
                    + tt[2][2] * tt[3][3] - pow(tt[1][2],2) - pow(tt[1][3],2)
                    - pow(tt[2][3],2) ) / 3. - 1./9.;
    double qCoefRt = sqrt( -qCoef);
    double rCoef = -0.5 * ( qCoef + 1./9. + tt[1][1] * pow(tt[2][3],2)
                           + tt[2][2] * pow(tt[1][3],2) + tt[3][3] * pow(tt[1][2],2)
                           - tt[1][1] * tt[2][2] * tt[3][3] )
    + tt[1][2] * tt[1][3] * tt[2][3] + 1./27.;
    double pTemp = max( min( rCoef / pow(qCoefRt,3), 1.), -1.);
    double pCoef = cos( acos(pTemp) / 3.);
    double pCoefRt = sqrt( 3. * (1. - pow(pCoef,2)) );
    eVal1 = 1./3. + qCoefRt * max( 2. * pCoef,  pCoefRt - pCoef);
    eVal3 = 1./3. + qCoefRt * min( 2. * pCoef, -pCoefRt - pCoef);
    eVal2 = 1. - eVal1 - eVal3;
    CPara = 3*(eVal1*eVal2 + eVal2*eVal3 + eVal3*eVal1);
    
//    cout<<"C parameter is "<<CPara<<endl;
    result.push_back(CPara);
    double DPara = 27*eVal1*eVal2*eVal3;
    result.push_back(DPara);
    result.push_back(newCPara);
//    return 0;
//    return CPara;
}




/*
int ISBbar(int mc){
    if(mc == 511 || mc == 521 || mc == 10511 || mc == 10521 || mc == 513 || mc == 523 || mc == 10513 || mc == 10523 || mc == 20513 || mc == 20523 || mc == 515 || mc == 525 || mc == 531 || mc == 10531 || mc == 533 || mc == 10533 || mc == 20533 || mc == 535 || mc == 541 || mc == 10541 || mc == 543 || mc == 10543 || mc == 20543 || mc == 545 || mc == -5122 || mc == -5112 || mc == -5212 || mc == -5222 || mc == -5114 || mc == -5214 || mc == -5224 || mc == -5132 || mc == -5232 || mc == -5312 || mc == -5322 || mc == -5314 || mc == -5324 || mc == -5332 || mc == -5334 || mc == -5142 || mc == -5242 || mc == -5412 || mc == -5422 || mc == -5414 || mc == -5424 || mc == -5342 || mc == -5432 || mc == -5434 || mc == -5442 || mc == -5444 || mc == -5512 || mc == -5522 || mc == -5514 || mc == -5524 || mc == -5532 || mc == -5534 || mc == -5542 || mc == -5544 || mc == -5544 ) {return 1;}
    else {return 0;}
}

int ISB(int mc){
    if(mc == -511 || mc == -521 || mc == -10511 || mc == -10521 || mc == -513 || mc == -523 || mc == -10513 || mc == -10523 || mc == -20513 || mc == -20523 || mc == -515 || mc == -525 || mc == -531 || mc == -10531 || mc == -533 || mc == -10533 || mc == -20533 || mc == -535 || mc == -541 || mc == -10541 || mc == -543 || mc == -10543 || mc == -20543 || mc == -545 || mc == 5122 || mc == 5112 || mc == 5212 || mc == 5222 || mc == 5114 || mc == 5214 || mc == 5224 || mc == 5132 || mc == 5232 || mc == 5312 || mc == 5322 || mc == 5314 || mc == 5324 || mc == 5332 || mc == 5334 || mc == 5142 || mc == 5242 || mc == 5412 || mc == 5422 || mc == 5414 || mc == 5424 || mc == 5342 || mc == 5432 || mc == 5434 || mc == 5442 || mc == 5444 || mc == 5512 || mc == 5522 || mc == 5514 || mc == 5524 || mc == 5532 || mc == 5534 || mc == 5542 || mc == 5544 || mc == 5544 ) {return 1;}
    else {return 0;}

}

int ISBLong(int mc){
    if(mc == 511 || mc == 521 || mc == 531 || mc == 541 || mc == 5112 || mc == 5122 || mc == 5132 || mc == 5232 || mc == 5332)
    {return 1;}
    else {return 0;}

}

int ISC(int mc){
    if(mc == 411 || mc == 421 || mc == 10411 || mc == 10421 || mc == 413 || mc == 423 || mc == 10413 || mc == 10423 || mc == 20413 || mc == 20423 || mc == 415 || mc == 425 || mc == 431 || mc == 10431 || mc == 433 || mc == 10433 || mc == 20433 || mc == 435      || mc == 4122 || mc == 4222 || mc == 4212 || mc == 4112 || mc == 4224 || mc == 4214 || mc == 4114 || mc == 4232 || mc == 4132 || mc == 4322 || mc == 4312 || mc == 4324 || mc == 4314 || mc == 4332 || mc == 4334 || mc == 4412 || mc == 4422 || mc == 4414 || mc == 4424 || mc == 4432 || mc == 4434 || mc == 4444) {return 1;}
    else {return 0;}
}

int ISCbar(int mc){
    if(mc == -411 || mc == -421 || mc == -10411 || mc == -10421 || mc == -413 || mc == -423 || mc == -10413 || mc == -10423 || mc == -20413 || mc == -20423 || mc == -415 || mc == -425 || mc == -431 || mc == -10431 || mc == -433 || mc == -10433 || mc == -20433 || mc == -435      || mc == -4122 || mc == -4222 || mc == -4212 || mc == -4112 || mc == -4224 || mc == -4214 || mc == -4114 || mc == -4232 || mc == -4132 || mc == -4322 || mc == -4312 || mc == -4324 || mc == -4314 || mc == -4332 || mc == -4334 || mc == -4412 || mc == -4422 || mc == -4414 || mc == -4424 || mc == -4432 || mc == -4434 || mc == -4444) {return 1;}
    else {return 0;}
}

int ISS(int mc){
    if(mc == 130 || mc == 310 || mc == 311 || mc == 321 || mc == 10311 || mc == 10321 || mc == 100311 || mc == 100321 || mc == 200311 || mc == 200321 || mc == 9000311 || mc == 9000321 || mc == 313 || mc == 323 || mc == 10313 || mc == 10323 || mc == 20313 || mc == 20323 || mc == 100313 || mc == 100323 || mc == 9000313 || mc == 9000323 || mc == 30313 || mc == 30323 || mc == 315 || mc == 325 || mc == 9000315 || mc == 9000325 || mc == 10315 || mc == 10325 || mc == 20315 || mc == 20325 || mc == 100315 || mc == 100325 || mc == 9010315 || mc == 9010325 || mc == 317 || mc == 327 || mc == 9010317 || mc == 9010327 || mc == 319 || mc == 329 || mc == 9000319 || mc == 9000329       || mc == 3122 || mc == 3222 || mc == 3212 || mc == 3112 || mc == 3224 || mc == 3214 || mc == 3114 || mc == 3322 || mc == 3312 || mc == 3324 || mc == 3314 || mc == 3334) {return 1;}
    else {return 0;}
}


void CalcuXYVector(TVector3 V3, std::vector<TVector3> &result){
    result.clear();
    TVector3 XV3( 1, 1, (-V3(0) - V3(1))/V3(2) );
    TVector3 YV3( V3(1)*XV3(2) - V3(2)*XV3(1), V3(2)*XV3(0) - V3(0)*XV3(2), V3(0)*XV3(1) - V3(1)*XV3(0));
    TVector3 XV3Unit = XV3.Unit();
    TVector3 YV3Unit = YV3.Unit();
    result.push_back(XV3Unit);
    result.push_back(YV3Unit);
    
}
*/

std::vector<vector<int> > pair4jets(int numjets){
    vector<int > a(numjets);
    for(int i = 0; i<numjets; i++){
        a[i] = i;
    }
    
    vector<vector<int> > vect1;
    for(int i = 0; i<(numjets-1); i++)
    {
        for(int j = i+1; j<numjets; j++)
        {
                vector<int> temp(2);
                vector<int> result;
                
                temp[0] = a[i];
                temp[1] = a[j];
                vect1.push_back(temp);
                result = DifferenceSet(a, temp, result);
                vect1.push_back(result);
        }
    }
    return vect1;
}
