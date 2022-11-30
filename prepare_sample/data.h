#include "TChain.h"
TChain vvH("Tau","Tau"); 
TChain qqH("Tau","Tau"); 

TChain ff("Tau","Tau"); 
TChain SW("Tau","Tau");
TChain SZ("Tau","Tau");
TChain WW("Tau","Tau");
TChain ZZ("Tau","Tau");
TChain ZH("Tau","Tau");
TChain Mix("Tau","Tau");


void preparedata()
{
   vvH.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/nnh_X.root"); 
   qqH.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/qqh_X.root");

   ff.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/bhabha.root");
   ff.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/e2e2.root");
   ff.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/e3e3.root");
   ff.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/nn.root");
   ff.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/qq.root");

   SW.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/sw_l0mu.root");
   SW.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/sw_l0tau.root");
   SW.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/sw_sl0qq.root");

   SZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/sze_l0e.root");
   SZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/sze_l0mu.root");
   SZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/sze_l0nunu.root");
   SZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/sze_l0tau.root");
   SZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/sze_sl0dd.root");
   SZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/sze_sl0uu.root");
   SZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/sznu_l0mumu.root");
   SZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/sznu_l0tautau.root");
   SZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/sznu_sl0nu_down.root");
   SZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/sznu_sl0nu_up.root");

   WW.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/ww_h0ccbs.root");
   WW.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/ww_h0ccds.root");
   WW.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/ww_h0cuxx.root");
   WW.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/ww_h0uubd.root");
   WW.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/ww_h0uusd.root");
   WW.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/ww_l0ll.root");
   WW.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/ww_sl0muq.root");
   WW.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/ww_sl0tauq.root");

   ZZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zz_h0cc_nots.root");
   ZZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zz_h0dtdt.root");
   ZZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zz_h0utut.root");
   ZZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zz_h0uu_notd.root");
   ZZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zz_l04mu.root");
   ZZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zz_l04tau.root");
   ZZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zz_l0mumu.root");
   ZZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zz_l0taumu.root");
   ZZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zz_l0tautau.root");
   ZZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zz_sl0mu_down.root");
   ZZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zz_sl0mu_up.root");
   ZZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zz_sl0nu_down.root");
   ZZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zz_sl0nu_up.root");
   ZZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zz_sl0tau_down.root");
   ZZ.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zz_sl0tau_up.root");

   ZH.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/nnh_X.root");
   ZH.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/e1e1h_X.root");
   ZH.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/e2e2h_X.root");
   ZH.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/e3e3h_X.root");
   ZH.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/qqh_X.root");

   Mix.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zzorww_l0mumu.root");
   Mix.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zzorww_l0tautau.root");
   Mix.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zzorww_h0cscs.root");
   Mix.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/zzorww_h0udud.root");
   Mix.Add("/cefs/higgs/shaqy/Quantum/sample_qq/Higgs_cc/vvH/sample/szeorsw_l0l.root");

}
