#ifndef _Higgsbb_hh_
#define _Higgsbb_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>

class TTree;

class Higgsbb  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new Higgsbb ; }

		Higgsbb();

		~Higgsbb() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;
		std::string _colAdcVals;
//		TFile *tree_file;

		int _overwrite;
		TTree *_outputTree;

		unsigned int eventNr;
		int Num;
        double ISR_En;
        int    num_Vertex;
        int num_jet, num_lepton, num_quark, num_charge, numHiggsParticles, num_neutrino;
        int count94;
        int HCharge, ZCharge;

	double B1Costheta, B2Costheta;

        double mass941, mass942;
        double sumb, sum2b, massB1b, massB2b, massRecoil1b, massRecoil2b, alpha;
    double delta_R1, delta_R2;
        double sumc, sum2c, massB1c, massB2c, massRecoil1c, massRecoil2c;
        int HQPDG[4], noHQPDG[4], numChargeOfJets[4], OriLepPDG[4], B1Quark[2], B2Quark[2];
        float ISRCosTheta[5], ISREn[5], HiggsX[4], maxHiggsX, miniAngleAmongJets, HiggsSumL[4], jetPx[4], jetPy[4], jetPz[4], jetEn[4];
    float HbbL[2], HccL[2], ZbbL[2], ZccL[2];
    float HbcL[2], HooL[2], ZbcL[2], ZooL[2];
    double jetsAngleOfHiggsbb, jetsAngleOfHiggscc, jetsAngleOfHiggsoo, HiggsJetsAngle, ZJetsAngle, ZHAngle;
    
    double visEn, LeadLepEn, AvLepEn, LeadhadEn, AvhadEn, LeadgamaEn, AvgamaEn, SubLeadgamaEn;
    double LeadElecEn, AvElecEn, LeadMuonEn, AvMuonEn;
    int    ElecNum, MuonNum;
    double LeadneutralEn;
    int    multiplicity, LepNum, hadNum, gamaNum;
    double Thrust, MaxBroadening, HeavyMass, Sphericity, CPara, DPara, cosTthrust;
    double boostThrust, boostH1En, boostH1Mass, boostH2En, boostH2Mass;
    double Y01, Y12, Y23, Y34, Y45;
    double Pmax, Pt, Pl;
    double LLRecoilMass, LLInvMass, JetsInvMass, JetsRecoilMass;
    int    haveEplusminus, haveMuplusminus;
    
    double P1HBMass, P1LBMass, P1HLL, P1LLL, P2HBMass, P2LBMass, P2HLL, P2LLL;
    int    P1HF, P1LF, P2HF, P2LF;
    double pairMB1, pairMB2, RminDif, RminDifZH;

    
    

    std::vector<Int_t> boson1Quark;
    std::vector<Int_t> boson2Quark;

    
		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


