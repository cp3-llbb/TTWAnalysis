#pragma once

#include <utility>
#include <vector>
#include <limits>

#include <Math/PtEtaPhiE4D.h>
#include <Math/LorentzVector.h>
#include <Math/VectorUtil.h>

#include <cp3_llbb/TTWAnalysis/interface/Indices.h>

// Needed because of gcc bug when using typedef and std::map
#define myLorentzVector ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>>

namespace TTWAnalysis {

  // Forward declaration to use this here
  float DeltaEta(const myLorentzVector &v1, const myLorentzVector &v2);

  struct BaseObject {
    BaseObject(myLorentzVector p4): p4(p4) {}
    BaseObject() {}
    virtual ~BaseObject() {}

    myLorentzVector p4;
  };

  struct GenParticle: BaseObject {
    GenParticle() {}
    GenParticle(myLorentzVector p4, int16_t pdg_id, int32_t pruned_idx): BaseObject(p4), pdg_id(pdg_id), pruned_idx(pruned_idx) {}
    int16_t pdg_id;
    // This will be used to check the decay chains of the particles, and not stored in the tree
    int32_t pruned_idx;
  };

  /*
  struct Lepton: BaseObject {
    Lepton():
      ID(LepID::Count, false),
      iso(LepIso::Count, false)
    {}
    Lepton(myLorentzVector p4, uint16_t idx, int16_t charge, bool isEl, bool isMu, bool isVeto = false, bool isLoose = false, bool isMedium = false, bool isTight = false, float isoValue = 0, bool isoLoose = false, bool isoTight = false):
      BaseObject(p4),
      idx(idx),
      charge(charge),
      isoValue(isoValue),
      isEl(isEl),
      isMu(isMu),
      ID(LepID::Count, false),
      iso(LepIso::Count, false)
      {
        if(isEl)
          ID[LepID::V] = isVeto;
        else
          ID[LepID::V] = isLoose; // for muons, re-use Loose as Veto ID
        ID[LepID::L] = isLoose;
        ID[LepID::M] = isMedium;
        ID[LepID::T] = isTight;

        if(isMu){
          iso[LepIso::L] = isoLoose;
          iso[LepIso::T] = isoTight;
        }else{
          iso[LepIso::L] = true;
        }
      }

    uint16_t idx; // stores index to electron/muon arrays
    int16_t charge;
    float isoValue;
    int16_t hlt_idx = -1; // Index to the matched HLT object. -1 if no match
    bool isEl;
    bool isMu;
    std::vector<bool> ID; // lepton ID: veto-loose-medium-tight
    std::vector<bool> iso; // lepton Iso: loose-tight (only for muons -> electrons only have loose)

    float hlt_DR_matched_object;
    float hlt_DPt_matched_object;

    bool hlt_already_tried_matching = false; // Internal flag; if true, it means this lepton has already been matched to an online object, even if no match has been found.

    int8_t pdg_id() const {
        int8_t id = (isEl) ? 11 : 13;
        return charge * id;
    }
  };

  struct DiLepton: BaseObject {
    DiLepton():
      ID(LepID::Count*LepID::Count, false),
      iso(LepIso::Count*LepIso::Count, false)
      {}

    std::pair<uint16_t, uint16_t> idxs; // stores indices to electron/muon arrays
    std::pair<uint16_t, uint16_t> lidxs; // stores indices to Lepton array
    std::pair<int16_t, int16_t> hlt_idxs; // Stores indices of matched online objects
    bool isElEl, isElMu, isMuEl, isMuMu;
    bool isOS; // opposite sign
    bool isSF; // same flavour
    std::vector<bool> ID; // combination of two lepton IDs
    std::vector<bool> iso; // combination of two lepton isolations
    float DR;
    float DEta;
    float DPhi;
  };

  struct Jet: BaseObject {
    Jet():
      ID(JetID::Count, false),
      minDRjl_lepIDIso(LepID::Count*LepIso::Count, std::numeric_limits<float>::max()),
      BWP(BWP::Count, false)
      {}

    uint16_t idx; // index to jet array
    std::vector<bool> ID;
    std::vector<float> minDRjl_lepIDIso; // defined for each combination of a lepton ID and isolation
    float CSVv2;
    std::vector<bool> BWP;
  };

  struct DiJet: BaseObject {
    DiJet():
      minDRjl_lepIDIso(LepID::Count*LepIso::Count, std::numeric_limits<float>::max()),
      BWP(BWP::Count*BWP::Count, false)
      {}

    std::pair<uint16_t, uint16_t> idxs; // stores indices to jets array
    std::pair<uint16_t, uint16_t> jidxs; // stores indices to TTWAnalysis::Jet array
    std::vector<float> minDRjl_lepIDIso; // defined for each combination of a lepton ID and isolation
    std::vector<bool> BWP; // combination of two b-tagging working points
    float DR;
    float DEta;
    float DPhi;
  };

  struct DiLepDiJet: BaseObject {
    DiLepDiJet():
      diLepton(nullptr),
      diJet(nullptr)
      {}
    DiLepDiJet(const DiLepton& _diLepton, const int diLepIdx, const DiJet& _diJet, const int diJetIdx):
      BaseObject(_diLepton.p4 + _diJet.p4),
      diLepton(&_diLepton),
      diLepIdx(diLepIdx),
      diJet(&_diJet),
      diJetIdx(diJetIdx),
      DR_ll_jj( ROOT::Math::VectorUtil::DeltaR(_diLepton.p4, _diJet.p4) ),
      DEta_ll_jj( DeltaEta(_diLepton.p4, _diJet.p4) ),
      DPhi_ll_jj( ROOT::Math::VectorUtil::DeltaPhi(_diLepton.p4, _diJet.p4) )
      {}

    const DiLepton* diLepton;
    uint16_t diLepIdx;
    const DiJet* diJet;
    uint16_t diJetIdx;

    float DR_ll_jj, DEta_ll_jj, DPhi_ll_jj;

    float minDRjl, maxDRjl;
    float minDEtajl, maxDEtajl;
    float minDPhijl, maxDPhijl;
  };

  struct DiLepDiJetMet: DiLepDiJet {
    DiLepDiJetMet() {}
    DiLepDiJetMet(const DiLepDiJet& diLepDiJet, uint16_t diLepDiJetIdx, const myLorentzVector& MetP4, bool hasNoHFMet = false):
      DiLepDiJet(*diLepDiJet.diLepton, diLepDiJet.diLepIdx, *diLepDiJet.diJet, diLepDiJet.diJetIdx),
      diLepDiJetIdx(diLepDiJetIdx),
      hasNoHFMet(hasNoHFMet)
    {
      DiLepDiJet::minDRjl = diLepDiJet.minDRjl;
      DiLepDiJet::maxDRjl = diLepDiJet.maxDRjl;
      DiLepDiJet::minDEtajl = diLepDiJet.minDEtajl;
      DiLepDiJet::maxDEtajl = diLepDiJet.maxDEtajl;
      DiLepDiJet::minDPhijl = diLepDiJet.minDPhijl;
      DiLepDiJet::maxDPhijl = diLepDiJet.maxDPhijl;

      p4 += MetP4;

      DR_ll_Met = ROOT::Math::VectorUtil::DeltaR(diLepton->p4, MetP4);
      DR_jj_Met = ROOT::Math::VectorUtil::DeltaR(diJet->p4, MetP4);

      DEta_ll_Met = DeltaEta(diLepton->p4, MetP4);
      DEta_jj_Met = DeltaEta(diJet->p4, MetP4);

      DPhi_ll_Met = ROOT::Math::VectorUtil::DeltaPhi(diLepton->p4, MetP4);
      DPhi_jj_Met = ROOT::Math::VectorUtil::DeltaPhi(diJet->p4, MetP4);

      DR_lljj_Met = ROOT::Math::VectorUtil::DeltaR(diLepton->p4 + diJet->p4, MetP4);
      DEta_lljj_Met = DeltaEta(diLepton->p4 + diJet->p4, MetP4);
      DPhi_lljj_Met = ROOT::Math::VectorUtil::DeltaPhi(diLepton->p4 + diJet->p4, MetP4);
    }

    uint16_t diLepDiJetIdx;
    bool hasNoHFMet;

    float DR_ll_Met, DR_jj_Met;
    float DEta_ll_Met, DEta_jj_Met;
    float DPhi_ll_Met, DPhi_jj_Met;

    float DR_lljj_Met;
    float DEta_lljj_Met;
    float DPhi_lljj_Met;

    float minDR_l_Met, minDR_j_Met;
    float maxDR_l_Met, maxDR_j_Met;
    float minDEta_l_Met, minDEta_j_Met;
    float maxDEta_l_Met, maxDEta_j_Met;
    float minDPhi_l_Met, minDPhi_j_Met;
    float maxDPhi_l_Met, maxDPhi_j_Met;
  };
  */

}

