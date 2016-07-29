//----------------------------------------------------------------------------
// alias_ra2b - Sets aliases for the #$@%-named branches in TreeMaker
//----------------------------------------------------------------------------

#ifndef ROOT_VERSION
#include "alias_ra2b.hpp"
#endif

void setAliasRa2b(TTree *tree){
  tree->SetAlias("ht"	   ,"HT");
  tree->SetAlias("mht"	   ,"MHT");
  tree->SetAlias("met"	   ,"MET");
  tree->SetAlias("njets"   ,"NJets");
  tree->SetAlias("nbm"	   ,"BTags");
  tree->SetAlias("nels"	   ,"@Electrons.size()");
  tree->SetAlias("nmus"	   ,"@Muons.size()");
  tree->SetAlias("nleps"   ,"@Electrons.size()+@Muons.size()");
  tree->SetAlias("dphi1"   ,"DeltaPhi1");
  tree->SetAlias("dphi2"   ,"DeltaPhi2");
  tree->SetAlias("dphi3"   ,"DeltaPhi3");
  tree->SetAlias("dphi4"   ,"DeltaPhi4");
  tree->SetAlias("low_dphi" ,"(DeltaPhi1<0.5||DeltaPhi2<0.5||DeltaPhi3<0.3||DeltaPhi4<0.3)");
  tree->SetAlias("ntrks"   ,"isoMuonTracks+isoElectronTracks+isoPionTracks");
  tree->SetAlias("mc"	   ,"GenParticles");
  tree->SetAlias("mc_id"   ,"GenParticles_PdgId");
  tree->SetAlias("mc_mom"  ,"GenParticles_ParentId");
  tree->SetAlias("weight"  ,"1000*Weight");
  tree->SetAlias("trig"	   ,"TriggerPass");
  tree->SetAlias("mus_pt"  ,"Muons.Pt()");
  tree->SetAlias("els_pt"  ,"Electrons.Pt()");
  tree->SetAlias("pf_calo" ,"PFCaloMETRatio");
  tree->SetAlias("pass_nopf"   ,"globalTightHalo2016Filter==1&&HBHENoiseFilter==1&&HBHEIsoNoiseFilter==1&&eeBadScFilter==1&&EcalDeadCellTriggerPrimitiveFilter==1&&BadChargedCandidateFilter&&BadPFMuonFilter&&NVtx>0&&JetID");
  tree->SetAlias("pass"	   ,"PFCaloMETRatio<5&&globalTightHalo2016Filter==1&&HBHENoiseFilter==1&&HBHEIsoNoiseFilter==1&&eeBadScFilter==1&&EcalDeadCellTriggerPrimitiveFilter==1&&BadChargedCandidateFilter&&BadPFMuonFilter&&NVtx>0&&JetID");
}
