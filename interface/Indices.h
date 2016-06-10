#pragma once

#include <map>
#include <array>
#include <string>

namespace TTWAnalysis {

  // Lepton IDs
  namespace LepID {
    enum LepID{ V, L, M, T, Count };
    // Ugly way to allow iterating over all items in the enumeration ( for(const LepID::LepID& id: LepID::it) )
    const std::array<LepID, Count> it = {{ V, L, M, T }};
    // Is useful in categories to construct cut strings out of each working point
    const std::map<LepID, std::string> map = { {V, "V"}, {L, "L"}, {M, "M"}, {T, "T"} };
  }

  // Lepton Isolation
  namespace LepIso {
    enum LepIso{ L, T, Count };
    const std::array<LepIso, Count> it = {{ L, T }};
    const std::map<LepIso, std::string> map = { {L, "L"}, {T, "T"} };
  }

  // Combination of Lepton ID + Lepton Isolation for a single lepton
  uint16_t LepIDIso(const LepID::LepID& id, const LepIso::LepIso& iso);
  std::string LepIDIsoStr(const LepID::LepID& id, const LepIso::LepIso& iso);

  // Combination of Lepton ID for a DiLepton object
  uint16_t LepLepID(const LepID::LepID& id1, const LepID::LepID& id2);
  std::string LepLepIDStr(const LepID::LepID& id1, const LepID::LepID& id2);

  // Combination of Lepton Isolation for a DiLepton object
  uint16_t LepLepIso(const LepIso::LepIso& iso1, const LepIso::LepIso& iso2);
  std::string LepLepIsoStr(const LepIso::LepIso& iso1, const LepIso::LepIso& iso2);

  // Combination of Lepton ID + Lepton Isolation for a DiLepton object
  uint16_t LepLepIDIso(const LepID::LepID& id1, const LepIso::LepIso& iso1, const LepID::LepID& id2, const LepIso::LepIso& iso2);
  std::string LepLepIDIsoStr(const LepID::LepID& id1, const LepIso::LepIso& iso1, const LepID::LepID& id2, const LepIso::LepIso& iso2);

  // Jet ID
  namespace JetID {
    enum JetID{ L, T, TLV, Count };
    const std::array<JetID, Count> it = {{ L, T, TLV }};
    const std::map<JetID, std::string> map = { {L, "L"}, {T, "T"}, {TLV, "TLV"} };
  }

  // Combination of Jet IDs for two jets (NOTE: NOT USED FOR NOW)
  uint16_t JetJetID(const JetID::JetID& id1, const JetID::JetID& id2);
  std::string JetJetIDStr(const JetID::JetID& id1, const JetID::JetID& id2);

  // B-tagging working points
  namespace BWP {
    enum BWP{ L, M, T, Count };
    const std::array<BWP, Count> it = {{ L, M, T }};
    const std::map<BWP, std::string> map = { {L, "L"}, {M, "M"}, {T, "T"} };
  }

  // Combination of Jet ID and B-tagging working point (NOTE: NOT USED FOR NOW)
  uint16_t JetIDBWP(const JetID::JetID& id, const BWP::BWP& wp);
  std::string JetIDBWPStr(const JetID::JetID& id, const BWP::BWP& wp);

  // Combination of Lepton ID + Lepton Isolation (one lepton) and B-tagging working point for one jet
  uint16_t LepIDIsoJetBWP(const LepID::LepID& id, const LepIso::LepIso& iso, const BWP::BWP& wp);
  std::string LepIDIsoJetBWPStr(const LepID::LepID& id, const LepIso::LepIso& iso, const BWP::BWP& wp);

  // Combination of B-tagging working points for two jets
  uint16_t JetJetBWP(const BWP::BWP& wp1, const BWP::BWP& wp2);
  std::string JetJetBWPStr(const BWP::BWP& wp1, const BWP::BWP& wp2);

  // Combination of Jet ID and B-tagging working points for two jets (NOTE: NOT USED FOR NOW)
  uint16_t JetJetIDBWP(const JetID::JetID& id1, const BWP::BWP& wp1, const JetID::JetID& id2, const BWP::BWP& wp2);
  std::string JetJetIDBWPStr(const JetID::JetID& id1, const BWP::BWP wp1, const JetID::JetID& id2, const BWP::BWP wp2);

  // Combination of Lepton ID + Lepton Isolation (one lepton) and B-tagging working points for two jets
  uint16_t LepIDIsoJetJetBWP(const LepID::LepID& id, const LepIso::LepIso& iso, const BWP::BWP& wp1, const BWP::BWP& wp2);
  std::string LepIDIsoJetJetBWPStr(const LepID::LepID& id, const LepIso::LepIso& iso, const BWP::BWP& wp1, const BWP::BWP& wp2);

  // Combination of Lepton ID, Lepton Isolation, and B-tagging working points for a two-lepton-two-b-jets object
  uint16_t LepLepIDIsoJetJetBWP(const LepID::LepID& id1, const LepIso::LepIso& iso1, const LepID::LepID& id2, const LepIso::LepIso& iso2, const BWP::BWP& wp1, const BWP::BWP& wp2);
  std::string LepLepIDIsoJetJetBWPStr(const LepID::LepID& id1, const LepIso::LepIso& iso1, const LepID::LepID& id2, const LepIso::LepIso& iso2, const BWP::BWP& wp1, const BWP::BWP& wp2);


  enum TTDecayType {
    UnknownTT = -1,
    NotTT = 0,
    Hadronic,
    Semileptonic_e,
    Semileptonic_mu,
    Dileptonic_mumu,
    Dileptonic_ee,
    Dileptonic_mue,

    // With tau
    Semileptonic_tau,
    Dileptonic_tautau,
    Dileptonic_mutau,
    Dileptonic_etau
  };

  enum WDecayType {
    UnknownW = -1,
    NotW = 0,
    HadronicW,
    LeptonicW_e,
    LeptonicW_mu,
    // With tau
    LeptonicW_tau,
  };

  enum TTWDecayType {
    UnknownTTW = -1,
    NotTTW = 0,
    HadronicTTW,
    SingleLeptonW,
    SingleLeptonTT,
    DiLeptonTT,
    DiLeptonWOST,
    DiLeptonWSST,
    TriLepton
  };
  enum TTWLeptonCombo {
    ZeroLepton=0,
    Single_e,
    Single_mu,
    Di_ee,
    Di_emu,
    Di_mumu,
    Tri_eee,
    Tri_eemu,
    Tri_emumu,
    Tri_mumumu,
    // with tau
    Single_tau,
    Di_etau,
    Di_mutau,
    Di_tautau,
    Tri_eetau,
    Tri_emutau,
    Tri_etautau,
    Tri_mumutau,
    Tri_mutautau,
    Tri_tautautau
  };
}
