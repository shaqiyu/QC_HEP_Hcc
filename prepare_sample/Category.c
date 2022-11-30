#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <utility>
#include "TFile.h"
#include "TH2.h"
#include "TString.h"
#include "data.h"
#include "TLorentzVector.h"
#include <iostream>
#include <iomanip>
#include <ctime>
#include <TTreeReaderArray.h>

using namespace TMath;
using namespace std;
using namespace RooFit;

void Cate(TChain &chainname, TString outputname)
{
    if (chainname.GetNtrees()==0) 
    {
        cout<<"No files yet!"<<endl;
            return;
                                
    }
    TTreeReader Reader(&chainname);
    cout << "now working on: " << outputname << endl;

    TTreeReaderValue<int> EventNr(Reader, "EventNr");
    TTreeReaderValue<int> Num(Reader, "Num");
    TTreeReaderValue<int> num_lepton(Reader, "num_lepton");
    TTreeReaderValue<int> num_quark(Reader, "num_quark");
    TTreeReaderValue<int> num_neutrino(Reader, "num_neutrino");
    TTreeReaderValue<int> num_Vertex(Reader, "num_Vertex");
    TTreeReaderValue<int> num_jet(Reader, "num_jet");
    TTreeReaderValue<int> num_charge(Reader, "num_charge");

    TTreeReaderValue<double> ISR_En(Reader, "ISR_En");
    TTreeReaderValue<int> count94(Reader, "count94");
    TTreeReaderValue<double> mass941(Reader, "mass941");
    TTreeReaderValue<double> mass942(Reader, "mass942");

    TTreeReaderValue<vector<int>> boson1Quark(Reader, "boson1Quark");
    TTreeReaderValue<vector<int>> boson2Quark(Reader, "boson2Quark");

    TTreeReaderArray<int> noHQPDG(Reader, "noHQPDG");
    TTreeReaderArray<int> HQPDG(Reader, "HQPDG");
    TTreeReaderArray<int> B1Quark(Reader, "B1Quark");
    TTreeReaderArray<int> B2Quark(Reader, "B2Quark");
    TTreeReaderArray<int> OriLepPDG(Reader, "OriLepPDG");

    TTreeReaderArray<float> HbbL(Reader, "HbbL");
    TTreeReaderArray<float> HccL(Reader, "HccL");
    TTreeReaderArray<float> HooL(Reader, "HooL");
   // TTreeReaderArray<float> HbcL(Reader, "HbcL");
    TTreeReaderArray<float> ZbbL(Reader, "ZbbL");
    TTreeReaderArray<float> ZccL(Reader, "ZccL");
    TTreeReaderArray<float> ZooL(Reader, "ZooL");
   // TTreeReaderArray<float> ZbcL(Reader, "ZbcL");

    TTreeReaderArray<float> ISREn(Reader, "ISREn");
    TTreeReaderArray<float> ISRCosTheta(Reader, "ISRCosTheta");

    TTreeReaderValue<double> HiggsJetsAngle(Reader, "HiggsJetsAngle");
    TTreeReaderValue<double> ZJetsAngle(Reader, "ZJetsAngle");
    TTreeReaderValue<double> ZHAngle(Reader, "ZHAngle");
    TTreeReaderValue<int> numHiggsParticles(Reader, "numHiggsParticles");

    TTreeReaderArray<float> HiggsX(Reader, "HiggsX");
    TTreeReaderArray<float> HiggsSumL(Reader, "HiggsSumL");
    TTreeReaderArray<float> jetPx(Reader, "jetPx");
    TTreeReaderArray<float> jetPy(Reader, "jetPy");
    TTreeReaderArray<float> jetPz(Reader, "jetPz");
    TTreeReaderArray<float> jetEn(Reader, "jetEn");
    TTreeReaderArray<int> numChargeOfJets(Reader, "numChargeOfJets");
    
    TTreeReaderValue<float> maxHiggsX(Reader, "maxHiggsX");
    TTreeReaderValue<float> miniAngleAmongJets(Reader, "miniAngleAmongJets");

    TTreeReaderValue<double> RminDif(Reader, "RminDif");
    TTreeReaderValue<double> RminDifZH(Reader, "RminDifZH");

    TTreeReaderValue<double> Pmax(Reader, "Pmax");
    TTreeReaderValue<double> Pt(Reader, "Pt");
    TTreeReaderValue<double> Pl(Reader, "Pl");

    TTreeReaderValue<double> sumb(Reader, "sumb");
    TTreeReaderValue<double> sum2b(Reader, "sum2b");
    TTreeReaderValue<double> sumc(Reader, "sumc");
    TTreeReaderValue<double> sum2c(Reader, "sum2c");
    TTreeReaderValue<double> massB1b(Reader, "sum2c");
    TTreeReaderValue<double> massB2b(Reader, "massB2b");
    TTreeReaderValue<double> massRecoil1b(Reader, "massRecoil1b");
    TTreeReaderValue<double> massRecoil2b(Reader, "massRecoil2b");
    TTreeReaderValue<double> massB1c(Reader, "massB1c");
    TTreeReaderValue<double> massB2c(Reader, "massB2c");
    TTreeReaderValue<double> massRecoil1c(Reader, "massRecoil1c");
    TTreeReaderValue<double> massRecoil2c(Reader, "massRecoil2c");
    
    TTreeReaderValue<double> alpha(Reader, "alpha");
    TTreeReaderValue<double> delta_R1(Reader, "delta_R1");
    TTreeReaderValue<double> delta_R2(Reader, "delta_R2");
    TTreeReaderValue<double> visEn(Reader, "visEn");

    TTreeReaderValue<int> multiplicity(Reader, "multiplicity");
    TTreeReaderValue<int> LepNum(Reader, "LepNum");  //???
    TTreeReaderValue<int> hadNum(Reader, "hadNum");  
    TTreeReaderValue<int> gamaNum(Reader, "gamaNum");  
    TTreeReaderValue<int> ElecNum(Reader, "ElecNum");
    TTreeReaderValue<int> MuonNum(Reader, "MuonNum");

    TTreeReaderValue<double> LeadLepEn(Reader, "LeadLepEn");
    TTreeReaderValue<double> AvLepEn(Reader, "AvLepEn");
    TTreeReaderValue<double> AvhadEn(Reader, "AvhadEn");
    TTreeReaderValue<double> LeadhadEn(Reader, "LeadhadEn");
    TTreeReaderValue<double> LeadgamaEn(Reader, "LeadgamaEn");
    TTreeReaderValue<double> SubLeadgamaEn(Reader, "SubLeadgamaEn");
    TTreeReaderValue<double> AvgamaEn(Reader, "AvgamaEn");
    TTreeReaderValue<double> LeadneutralEn(Reader, "LeadneutralEn");
    TTreeReaderValue<double> LeadElecEn(Reader, "LeadElecEn");
    TTreeReaderValue<double> AvElecEn(Reader, "AvElecEn");
    TTreeReaderValue<double> LeadMuonEn(Reader, "LeadMuonEn");
    TTreeReaderValue<double> AvMuonEn(Reader, "AvMuonEn");

    TTreeReaderValue<double> Thrust(Reader, "Thrust");
    TTreeReaderValue<double> MaxBroadening(Reader, "MaxBroadening");
    TTreeReaderValue<double> HeavyMass(Reader, "HeavyMass");
    TTreeReaderValue<double> Sphericity(Reader, "Sphericity");
    TTreeReaderValue<double> CPara(Reader, "CPara");
    TTreeReaderValue<double> DPara(Reader, "DPara");
    TTreeReaderValue<double> cosTthrust(Reader, "cosTthrust");
    TTreeReaderValue<double> boostThrust(Reader, "boostThrust");
    TTreeReaderValue<double> boostH1En(Reader, "boostH1En");
    TTreeReaderValue<double> boostH2En(Reader, "boostH2En");
    TTreeReaderValue<double> boostH1Mass(Reader, "boostH1Mass");
    TTreeReaderValue<double> boostH2Mass(Reader, "boostH2Mass");

    TTreeReaderValue<double> Y01(Reader, "Y01");
    TTreeReaderValue<double> Y12(Reader, "Y12");
    TTreeReaderValue<double> Y23(Reader, "Y23");
    TTreeReaderValue<double> Y34(Reader, "Y34");
    TTreeReaderValue<double> Y45(Reader, "Y45");
    
    TTreeReaderValue<double> LLRecoilMass(Reader, "LLRecoilMass");
    TTreeReaderValue<double> LLInvMass(Reader, "LLInvMass");
    TTreeReaderValue<double> JetsInvMass(Reader, "JetsInvMass");
    TTreeReaderValue<double> JetsRecoilMass(Reader, "JetsRecoilMass");

    TTreeReaderValue<int> haveMuplusminus(Reader, "haveMuplusminus");
    TTreeReaderValue<int> haveEplusminus(Reader, "haveEplusminus");
    
    TTreeReaderValue<double> pairMB1(Reader, "pairMB1");
    TTreeReaderValue<double> pairMB2(Reader, "pairMB2");
  
    TTreeReaderValue<float> weight(Reader, "weight");
    //TTreeReaderValue<int> filename0(Reader, "filename");

//output
    TString outputname0 = outputname + ".root";
    TFile f_output(outputname0, "RECREATE");
    TTree output("output", "output");

    
    int Types, DTypes;
    int O_EventNr, O_Num, O_num_lepton, O_num_quark, O_num_neutrino, O_num_Vertex, O_num_jet, O_num_charge;
    double O_ISR_En, O_mass941, O_mass942;
    int O_count94;
    vector<int> O_boson1Quark, O_boson2Quark;
    double O_HiggsJetsAngle, O_ZJetsAngle, O_ZHAngle;
    int O_numHiggsParticles;
    float O_maxHiggsX, O_miniAngleAmongJets;
    double O_RminDif, O_RminDifZH;
    double O_Pmax, O_Pt, O_Pl, O_alpha, O_delta_R1, O_delta_R2, O_visEn;
    double O_sumb, O_sum2b, O_sumc, O_sum2c, O_massB1b, O_massB2b, O_massRecoil1b, O_massRecoil2b, O_massB1c, O_massB2c, O_massRecoil1c, O_massRecoil2c;
    
    int O_multiplicity, O_LepNum, O_hadNum, O_gamaNum, O_ElecNum, O_MuonNum;
    double O_LeadLepEn, O_AvLepEn, O_AvhadEn, O_LeadhadEn, O_LeadgamaEn, O_SubLeadgamaEn, O_AvgamaEn, O_LeadneutralEn, O_LeadElecEn, O_AvElecEn, O_LeadMuonEn, O_AvMuonEn;
    double O_Thrust, O_MaxBroadening, O_HeavyMass, O_Sphericity, O_CPara, O_DPara, O_cosTthrust, O_boostThrust, O_boostH1En, O_boostH2En, O_boostH1Mass, O_boostH2Mass; 
    double O_Y01, O_Y12, O_Y23, O_Y34, O_Y45;
    double O_LLRecoilMass, O_LLInvMass, O_JetsInvMass, O_JetsRecoilMass;
    int O_haveMuplusminus, O_haveEplusminus;
    double O_pairMB1, O_pairMB2;
    float O_weight;
    long long eventID = 0;
    long long eventID_2 = 0;
    bool pass_sigcut;
    bool pass_sigcutbb;
    bool pass_sigcutcc;
    bool pass_sigcutgg;
    int O_noHQPDG[4], O_HQPDG[4], O_B1Quark[2], O_B2Quark[2], O_OriLepPDG[4];
    float O_HbbL[2], O_HccL[2], O_HooL[2], O_HbcL[2], O_ZbbL[2], O_ZccL[2], O_ZooL[2], O_ZbcL[2], O_ISREn[5], O_ISRCosTheta[5];
    float O_HiggsX[4], O_HiggsSumL[4], O_jetPx[4], O_jetPy[4], O_jetPz[4], O_jetEn[4], O_numChargeOfJets[4];

    output.Branch("EventNr",            &O_EventNr);
    output.Branch("Num",                &O_Num);
    output.Branch("num_lepton",         &O_num_lepton);
    output.Branch("num_quark",          &O_num_quark);
    output.Branch("num_neutrino",       &O_num_neutrino);
    output.Branch("num_Vertex",         &O_num_Vertex);
    output.Branch("num_jet",            &O_num_jet);
    output.Branch("num_charge",         &O_num_charge);

    output.Branch("ISR_En",             &O_ISR_En);
    output.Branch("count94",            &O_count94);
    output.Branch("mass941",            &O_mass941);
    output.Branch("mass942",            &O_mass942);

    output.Branch("boson1Quark",        &O_boson1Quark);
    output.Branch("boson2Quark",        &O_boson2Quark);

    output.Branch("HiggsJetsAngle",     &O_HiggsJetsAngle);
    output.Branch("ZJetsAngle",         &O_ZJetsAngle);
    output.Branch("ZHAngle",            &O_ZHAngle);
    output.Branch("numHiggsParticles",  &O_numHiggsParticles);
    
    output.Branch("maxHiggsX",          &O_maxHiggsX);
    output.Branch("miniAngleAmongJets", &O_miniAngleAmongJets);

    output.Branch("RminDif",            &O_RminDif);
    output.Branch("RminDifZH",          &O_RminDifZH);

    output.Branch("Pmax",               &O_Pmax);
    output.Branch("Pt",                 &O_Pt);
    output.Branch("Pl",                 &O_Pl);
   
    output.Branch("sumb",               &O_sumb);
    output.Branch("sum2b",              &O_sum2b);
    output.Branch("sumc",               &O_sumc);
    output.Branch("sum2c",              &O_sum2c);
    output.Branch("massB1b",            &O_sum2c);
    output.Branch("massB2b",            &O_massB2b);
    output.Branch("massRecoil1b",       &O_massRecoil1b);
    output.Branch("massRecoil2b",       &O_massRecoil2b);
    output.Branch("massB1c",            &O_massB1c);
    output.Branch("massB2c",            &O_massB2c);
    output.Branch("massRecoil1c",       &O_massRecoil1c);
    output.Branch("massRecoil2c",       &O_massRecoil2c);

    output.Branch("alpha",              &O_alpha);
    output.Branch("delta_R1",           &O_delta_R1);
    output.Branch("delta_R2",           &O_delta_R2);
    output.Branch("visEn",              &O_visEn);

    output.Branch("multiplicity",       &O_multiplicity);
    output.Branch("LepNum",             &O_LepNum);  //???
    output.Branch("hadNum",             &O_hadNum);  
    output.Branch("gamaNum",            &O_gamaNum);  
    output.Branch("ElecNum",            &O_ElecNum);
    output.Branch("MuonNum",            &O_MuonNum);

    output.Branch("LeadLepEn",          &O_LeadLepEn);
    output.Branch("AvLepEn",            &O_AvLepEn);
    output.Branch("AvhadEn",            &O_AvhadEn);
    output.Branch("LeadhadEn",          &O_LeadhadEn);
    output.Branch("LeadgamaEn",         &O_LeadgamaEn);
    output.Branch("SubLeadgamaEn",      &O_SubLeadgamaEn);
    output.Branch("AvgamaEn",           &O_AvgamaEn);
    output.Branch("LeadneutralEn",      &O_LeadneutralEn);
    output.Branch("LeadElecEn",         &O_LeadElecEn);
    output.Branch("AvElecEn",           &O_AvElecEn);
    output.Branch("LeadMuonEn",         &O_LeadMuonEn);
    output.Branch("AvMuonEn",           &O_AvMuonEn);

    output.Branch("Thrust",             &O_Thrust);
    output.Branch("MaxBroadening",      &O_MaxBroadening);
    output.Branch("HeavyMass",          &O_HeavyMass);
    output.Branch("Sphericity",         &O_Sphericity);
    output.Branch("CPara",              &O_CPara);
    output.Branch("DPara",              &O_DPara);
    output.Branch("cosTthrust",         &O_cosTthrust);
    output.Branch("boostThrust",        &O_boostThrust);
    output.Branch("boostH1En",          &O_boostH1En);
    output.Branch("boostH2En",          &O_boostH2En);
    output.Branch("boostH1Mass",        &O_boostH1Mass);
    output.Branch("boostH2Mass",        &O_boostH2Mass);

    output.Branch("Y01",                &O_Y01);
    output.Branch("Y12",                &O_Y12);
    output.Branch("Y23",                &O_Y23);
    output.Branch("Y34",                &O_Y34);
    output.Branch("Y45",                &O_Y45);
    
    output.Branch("LLRecoilMass",       &O_LLRecoilMass);
    output.Branch("LLInvMass",          &O_LLInvMass);
    output.Branch("JetsInvMass",        &O_JetsInvMass);
    output.Branch("JetsRecoilMass",     &O_JetsRecoilMass);

    output.Branch("haveMuplusminus",    &O_haveMuplusminus);
    output.Branch("haveEplusminus",     &O_haveEplusminus);
    
    output.Branch("pairMB1",            &O_pairMB1);
    output.Branch("pairMB2",            &O_pairMB2);
  
    output.Branch("weight",             &O_weight);

	output.Branch("pass_sigcut",        &pass_sigcut);
	output.Branch("pass_sigcutbb",      &pass_sigcutbb);
	output.Branch("pass_sigcutcc",      &pass_sigcutcc);
	output.Branch("pass_sigcutgg",      &pass_sigcutgg);
    output.Branch("Types",              &Types);
    output.Branch("DTypes",             &DTypes);

  
    Long64_t ievt = Reader.GetEntries(true);
    cout << ievt << endl;

    Long64_t kk = 0;

    long double count[50];
    for (int i = 0; i < 50; i++)
        count[i] = 0;

    TString filename = "";
    while (Reader.Next())
    {
        kk++;
    if (filename != TString(chainname.GetFile()->GetName()))
        {
            filename = TString(chainname.GetFile()->GetName());
            std::cout << "---Now running with " << filename << endl;
            TTree *temptree = (TTree *)(chainname.GetFile()->Get("Tau"));
            std::cout << "Entries: " << temptree->GetEntries() << endl;
        }
    if (kk % 100000 == 0)
        std::cout << "--- ... Processing event: " << kk << std::endl;

//cutflow
//clear
    Types = -1;
    DTypes = -1;
    O_weight = -1;
    O_EventNr            = -1;
    O_Num                = -1;
    O_num_lepton         = -1;
    O_num_quark          = -1;
    O_num_neutrino       = -1;
    O_num_Vertex         = -1;
    O_num_jet            = -1;
    O_num_charge         = -1;
    O_ISR_En             = -1;
    O_count94            = -1;
    O_mass941            = -1;
    O_mass942            = -1;
    O_boson1Quark.clear();
    O_boson2Quark.clear();
    O_HiggsJetsAngle     = -1;
    O_ZJetsAngle         = -1;
    O_ZHAngle            = -1;
    O_numHiggsParticles  = -1;
    O_maxHiggsX          = -1;
    O_miniAngleAmongJets = -1;
    O_RminDif            = -1;
    O_RminDifZH          = -1;
    O_Pmax               = -1;
    O_Pt                 = -1;
    O_Pl                 = -1;
    O_sumb               = -1;
    O_sum2b              = -1;
    O_sumc               = -1;
    O_sum2c              = -1;
    O_sum2c              = -1;
    O_massB2b            = -1;
    O_massRecoil1b       = -1;
    O_massRecoil2b       = -1;
    O_massB1c            = -1;
    O_massB2c            = -1;
    O_massRecoil1c       = -1;
    O_massRecoil2c       = -1;
    O_alpha              = -1;
    O_delta_R1           = -1;
    O_delta_R2           = -1;
    O_visEn              = -1;
    O_multiplicity       = -1;
    O_LepNum             = -1;
    O_hadNum             = -1;
    O_gamaNum            = -1;
    O_ElecNum            = -1;
    O_MuonNum            = -1;
    O_LeadLepEn          = -1;
    O_AvLepEn            = -1;
    O_AvhadEn            = -1;
    O_LeadhadEn          = -1;
    O_LeadgamaEn         = -1;
    O_SubLeadgamaEn      = -1;
    O_AvgamaEn           = -1;
    O_LeadneutralEn      = -1;
    O_LeadElecEn         = -1;
    O_AvElecEn           = -1;
    O_LeadMuonEn         = -1;
    O_AvMuonEn           = -1;
    O_Thrust             = -1;
    O_MaxBroadening      = -1;
    O_HeavyMass          = -1;
    O_Sphericity         = -1;
    O_CPara              = -1;
    O_DPara              = -1;
    O_cosTthrust         = -1;
    O_boostThrust        = -1;
    O_boostH1En          = -1;
    O_boostH2En          = -1;
    O_boostH1Mass        = -1;
    O_boostH2Mass        = -1;
    O_Y01                = -1;
    O_Y12                = -1;
    O_Y23                = -1;
    O_Y34                = -1;
    O_Y45                = -1;
    O_LLRecoilMass       = -1;
    O_LLInvMass          = -1;
    O_JetsInvMass        = -1;
    O_JetsRecoilMass     = -1;
    O_haveMuplusminus    = -1;
    O_haveEplusminus     = -1;
    O_pairMB1            = -1;
    O_pairMB2            = -1;

//give value
    O_EventNr                = *EventNr;                        
    O_Num                    = *Num;
    O_num_lepton             = *num_lepton;
    O_num_quark              = *num_quark;
    O_num_neutrino           = *num_neutrino;
    O_num_Vertex             = *num_Vertex;
    O_num_jet                = *num_jet;
    O_num_charge             = *num_charge;
    O_ISR_En                 = *ISR_En;
    O_count94                = *count94;
    O_mass941                = *mass941;
    O_mass942                = *mass942;
    O_boson1Quark            = *boson1Quark;
    O_boson2Quark            = *boson2Quark;
    O_HiggsJetsAngle         = *HiggsJetsAngle;
    O_ZJetsAngle             = *ZJetsAngle;
    O_ZHAngle                = *ZHAngle;
    O_numHiggsParticles      = *numHiggsParticles;
    O_maxHiggsX              = *maxHiggsX;
    O_miniAngleAmongJets     = *miniAngleAmongJets;
    O_RminDif                = *RminDif;
    O_RminDifZH              = *RminDifZH;
    O_Pmax                   = *Pmax;
    O_Pt                     = *Pt;
    O_Pl                     = *Pl;
    O_sumb                   = *sumb;
    O_sum2b                  = *sum2b;
    O_sumc                   = *sumc;
    O_sum2c                  = *sum2c;
    O_sum2c                  = *sum2c;
    O_massB2b                = *massB2b;
    O_massRecoil1b           = *massRecoil1b;
    O_massRecoil2b           = *massRecoil2b;
    O_massB1c                = *massB1c;
    O_massB2c                = *massB2c;
    O_massRecoil1c           = *massRecoil1c;
    O_massRecoil2c           = *massRecoil2c;
    O_alpha                  = *alpha;
    O_delta_R1               = *delta_R1;
    O_delta_R2               = *delta_R2;
    O_visEn                  = *visEn;
    O_multiplicity           = *multiplicity;
    O_LepNum                 = *LepNum;  //???
    O_hadNum                 = *hadNum;  
    O_gamaNum                = *gamaNum;  
    O_ElecNum                = *ElecNum;
    O_MuonNum                = *MuonNum;
    O_LeadLepEn              = *LeadLepEn;
    O_AvLepEn                = *AvLepEn;
    O_AvhadEn                = *AvhadEn;
    O_LeadhadEn              = *LeadhadEn;
    O_LeadgamaEn             = *LeadgamaEn;
    O_SubLeadgamaEn          = *SubLeadgamaEn;
    O_AvgamaEn               = *AvgamaEn;
    O_LeadneutralEn          = *LeadneutralEn;
    O_LeadElecEn             = *LeadElecEn;
    O_AvElecEn               = *AvElecEn;
    O_LeadMuonEn             = *LeadMuonEn;
    O_AvMuonEn               = *AvMuonEn;
    O_Thrust                 = *Thrust;
    O_MaxBroadening          = *MaxBroadening;
    O_HeavyMass              = *HeavyMass;
    O_Sphericity             = *Sphericity;
    O_CPara                  = *CPara;
    O_DPara                  = *DPara;
    O_cosTthrust             = *cosTthrust;
    O_boostThrust            = *boostThrust;
    O_boostH1En              = *boostH1En;
    O_boostH2En              = *boostH2En;
    O_boostH1Mass            = *boostH1Mass;
    O_boostH2Mass            = *boostH2Mass;
    O_Y01                    = *Y01;
    O_Y12                    = *Y12;
    O_Y23                    = *Y23;
    O_Y34                    = *Y34;
    O_Y45                    = *Y45;
    O_LLRecoilMass           = *LLRecoilMass;
    O_LLInvMass              = *LLInvMass;
    O_JetsInvMass            = *JetsInvMass;
    O_JetsRecoilMass         = *JetsRecoilMass;
    O_haveMuplusminus        = *haveMuplusminus;
    O_haveEplusminus         = *haveEplusminus;
    O_pairMB1                = *pairMB1;
    O_pairMB2                = *pairMB2;

    O_weight = *weight;

    pass_sigcut = false;
    pass_sigcutbb = false;
    pass_sigcutcc = false;
    pass_sigcutgg = false;

//cut

    if (O_num_lepton == 2 && O_num_quark == 2 && O_num_neutrino == 2 && HQPDG[0] != 0 && HQPDG[1] != 0)    pass_sigcut=true;
    if (O_num_lepton == 2 && O_num_quark == 2 && O_num_neutrino == 2 && HQPDG[0] == 5 && HQPDG[1] == 5)    pass_sigcutbb=true;
    if (O_num_lepton == 2 && O_num_quark == 2 && O_num_neutrino == 2 && HQPDG[0] == 4 && HQPDG[1] == 4)    pass_sigcutcc=true;
    if (O_num_lepton == 2 && O_num_quark == 2 && O_num_neutrino == 2 && HQPDG[0] == 21 && HQPDG[1] == 21)  pass_sigcutgg=true;

    //Only for ZH
   // if (O_num_lepton == 2 && O_num_quark == 2 && O_num_neutrino == 2 && HQPDG[0] != 0 && HQPDG[1] != 0) continue;

    count[0] += O_weight;
    if(!(O_JetsRecoilMass >= 74 && O_JetsRecoilMass <= 131)) continue;
    count[1] += O_weight;
    if(!(O_visEn >= 109 && O_visEn <= 143)) continue;
    count[2] += O_weight;
    if(!(O_LeadLepEn <= 42)) continue;
    count[3] += O_weight;
    if(!(O_multiplicity >= 40 && O_multiplicity <= 130)) continue;
    count[4] += O_weight;
    if(!(O_LeadneutralEn <= 41)) continue;
    count[5] += O_weight;
    if(!(O_Pt >= 20 && O_Pt <= 60)) continue;
    count[6] += O_weight;
    if(!(O_Pl <= 50)) continue;
    count[7] += O_weight;
    if(!(-log(O_Y23) >= 3.375)) continue;
    count[8] += O_weight;
    if(!(O_JetsInvMass >= 116 && O_JetsInvMass <= 134)) continue;
    count[9] += O_weight;

    if(pass_sigcut) count[10] += O_weight;
    if(pass_sigcutbb) count[11] += O_weight;
    if(pass_sigcutcc) count[12] += O_weight;
    if(pass_sigcutgg) count[13] += O_weight;

//sample types
    if (filename.Contains("nnh_X"))
        {DTypes = 0;}
    else if (filename.Contains("qqh_X"))
        {DTypes = 1;}
    else if (filename.Contains("bhabha.root"))
        {Types = 1;}
    else if (filename.Contains("e2e2.root"))
        {Types = 1;}
    else if (filename.Contains("e3e3.root"))
        {Types = 1;}
    else if (filename.Contains("nn.root"))
        {Types = 1;}
    else if (filename.Contains("qq.root"))
        {Types = 1;}
    else if (filename.Contains("sw"))
        {Types = 2;}
    else if (filename.Contains("sze"))
        {Types = 3;
        DTypes = 3;}
    else if (filename.Contains("sznu"))
        {Types = 3;
        DTypes = 4;}
    else if (filename.Contains("ww"))
        {Types = 4;}
    else if (filename.Contains("zz"))
        {Types = 5;}
    else if (filename.Contains("h_X"))
        {Types = 6;}
    else 
        {Types = 7;}
    
    eventID++;
        output.Fill();
    }
    cout<<endl;
    for(int i=0;i<14;i++){
        cout<<count[i]<<endl;
    }

    output.AutoSave("Overwrite");
    f_output.Write();

    output.SetBranchAddress("pass_sigcut", &pass_sigcut);
    output.SetBranchAddress("pass_sigcutbb", &pass_sigcutbb);
    output.SetBranchAddress("pass_sigcutcc", &pass_sigcutcc);
    output.SetBranchAddress("pass_sigcutgg", &pass_sigcutgg);
    output.SetBranchAddress("Types", &Types);
    ULong64_t NEntries=output.GetEntries();
    if (DTypes == 0 )
    {
    TString outputname00=outputname+"_sig.root";
    TFile f_output00(outputname00, "RECREATE");
    TTree *newTree00=output.CloneTree(0);

    for (ULong64_t i=0; i<NEntries; i++)
    {
        output.GetEntry(i);
        if (pass_sigcut) 
        {
            eventID_2++;
            newTree00->Fill();
        }
    }
    newTree00->AutoSave("Overwrite");
    f_output00.Close();
    cout<<"Done sigcut.root"<<endl;
    
    TString outputname01=outputname+"_bb.root";
    TFile f_output01(outputname01, "RECREATE");
    TTree *newTree01=output.CloneTree(0);

    for (ULong64_t i=0; i<NEntries; i++)
    {
        output.GetEntry(i);
        if (pass_sigcutbb) 
        {
            eventID_2++;
            newTree01->Fill();
        }
    }
    newTree01->AutoSave("Overwrite");
    f_output01.Close();
    cout<<"Done cut_bb.root"<<endl;

    TString outputname02=outputname+"_cc.root";
    TFile f_output02(outputname02, "RECREATE");
    TTree *newTree02=output.CloneTree(0);

    for (ULong64_t i=0; i<NEntries; i++)
    {
        output.GetEntry(i);
        if (pass_sigcutcc) 
        {
            eventID_2++;
            newTree02->Fill();
        }
    }
    newTree02->AutoSave("Overwrite");
    f_output02.Close();
    cout<<"Done cut_cc.root"<<endl;

    TString outputname03=outputname+"_gg.root";
    TFile f_output03(outputname03, "RECREATE");
    TTree *newTree03=output.CloneTree(0);

    for (ULong64_t i=0; i<NEntries; i++)
    {
        output.GetEntry(i);
        if (pass_sigcutgg) 
        {
            eventID_2++;
            newTree03->Fill();
        }
    } 
 
    newTree03->AutoSave("Overwrite");
    f_output03.Close();
    cout<<"Done cut_gg.root"<<endl;
    }

    f_output.Close();
    
}
void MCCategory()
{
    preparedata();

   // Cate(vvH,    "vvH");
    //Cate(ff,    "ff");
    //Cate(SW,    "SW");
   // Cate(SZ,    "SZ");
   // Cate(WW,    "WW");
   // Cate(ZZ,    "ZZ");
    Cate(ZH,    "ZH");
   // Cate(Mix,    "Mix");
    
}  
