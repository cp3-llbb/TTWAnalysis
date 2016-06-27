#include <boost/container/flat_map.hpp>

#include "cp3_llbb/Framework/interface/Analyzer.h"
#include "cp3_llbb/Framework/interface/AnalyzersManager.h"
#include "cp3_llbb/TTWAnalysis/interface/NewTypes.h"
#include "cp3_llbb/TTWAnalysis/interface/Tools.h"
#include "cp3_llbb/TTWAnalysis/interface/GenStatusFlags.h"

#define TT_GEN_DEBUG (false)

#if TT_GEN_DEBUG
#define FILL_GEN_COLL( X ) \
    if (flags.isLastCopy()) { \
        genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
        gen_##X = gen_index; \
        gen_index++; \
        std::cout << "Assigning gen_" #X " = " << i << " (" << pdg_id << ")" << std::endl; \
    } \
    if (flags.isFirstCopy()) { \
        genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
        gen_##X##_beforeFSR = gen_index; \
        gen_index++; \
        std::cout << "Assigning gen_" #X "_beforeFSR = " << i << " (" << pdg_id << ")" << std::endl; \
    }
#else
#define FILL_GEN_COLL( X ) \
    if (flags.isLastCopy()) { \
        genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
        gen_##X = gen_index; \
        gen_index++; \
    } \
    if (flags.isFirstCopy()) { \
        genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
        gen_##X##_beforeFSR = gen_index; \
        gen_index++; \
    }
#endif

// Assign index to X if it's empty, or Y if not
#if TT_GEN_DEBUG
#define FILL_GEN_COLL2(X, Y, ERROR) \
    if (flags.isLastCopy()) { \
        if (gen_##X == -1){ \
            genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
            gen_##X = gen_index; \
            gen_index++; \
            std::cout << "Assigning gen_" #X " = " << i << " (" << pdg_id << ")" << std::endl; \
        } else if (gen_##Y == -1) { \
            genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
            gen_##Y = gen_index; \
            gen_index++; \
            std::cout << "Assigning gen_" #Y " = " << i << " (" << pdg_id << ")" << std::endl; \
        } else \
            std::cout << ERROR << std::endl; \
    } \
    if (flags.isFirstCopy()) { \
        if (gen_##X##_beforeFSR == -1) { \
            genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
            gen_##X##_beforeFSR = gen_index; \
            gen_index++; \
            std::cout << "Assigning gen_" #X "_beforeFSR = " << i << " (" << pdg_id << ")" << std::endl; \
        } else if (gen_##Y##_beforeFSR == -1) { \
            genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
            gen_##Y##_beforeFSR = gen_index; \
            gen_index++; \
            std::cout << "Assigning gen_" #Y "_beforeFSR = " << i << " (" << pdg_id << ")" << std::endl; \
        } else \
            std::cout << ERROR << std::endl; \
    }
#else
#define FILL_GEN_COLL2(X, Y, ERROR) \
    if (flags.isLastCopy()) { \
        if (gen_##X == -1){ \
            genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
            gen_##X = gen_index; \
            gen_index++; \
        } else if (gen_##Y == -1) { \
            genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
            gen_##Y = gen_index; \
            gen_index++; \
        } else \
            std::cout << ERROR << std::endl; \
    } \
    if (flags.isFirstCopy()) { \
        if (gen_##X##_beforeFSR == -1) { \
            genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
            gen_##X##_beforeFSR = gen_index; \
            gen_index++; \
        } else if (gen_##Y##_beforeFSR == -1) { \
            genParticles.push_back( GenParticle(gen_particles.pruned_p4[i], pdg_id, i) ); \
            gen_##Y##_beforeFSR = gen_index; \
            gen_index++; \
        } else \
            std::cout << ERROR << std::endl; \
    }
#endif

/**
 * Generator-level truth information about the ttbar system
 * and matching of the jets and leptons to top decay products
 */
class TTWTruthAnalyzer : public Framework::Analyzer {
  public:
    TTWTruthAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config)
      : Analyzer(name, tree_, config)
      {
        m_ttWName = config.getParameter<std::string>("TTWAnalyzer");
        m_jetWPs = config.getParameter<std::vector<std::string>>("JetWP");
        std::sort(std::begin(m_jetWPs), std::end(m_jetWPs));
        for ( const auto& iWP : m_jetWPs ) {
          gen_b_deltaR              [iWP] = &(tree["gen_b_deltaR_"              +iWP].write<std::vector<float>>());
          gen_bbar_deltaR           [iWP] = &(tree["gen_bbar_deltaR_"           +iWP].write<std::vector<float>>());
          gen_b_beforeFSR_deltaR    [iWP] = &(tree["gen_b_beforeFSR_deltaR_"    +iWP].write<std::vector<float>>());
          gen_bbar_beforeFSR_deltaR [iWP] = &(tree["gen_bbar_beforeFSR_deltaR_" +iWP].write<std::vector<float>>());

          gen_matched_b             [iWP] = &(tree["gen_matched_b_"             +iWP].write<int16_t>());
          gen_matched_bbar          [iWP] = &(tree["gen_matched_bbar_"          +iWP].write<int16_t>());
          gen_matched_b_beforeFSR   [iWP] = &(tree["gen_matched_b_beforeFSR_"   +iWP].write<int16_t>());
          gen_matched_bbar_beforeFSR[iWP] = &(tree["gen_matched_bbar_beforeFSR_"+iWP].write<int16_t>());
        }
      }

    virtual ~TTWTruthAnalyzer();

    virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const AnalyzersManager&, const CategoryManager&) override;
    virtual void registerCategories(CategoryManager& manager, const edm::ParameterSet&) override;

  private:
    std::string m_ttWName;
    std::vector<std::string> m_jetWPs;

    // Gen matching. All indexes are from the `genParticles` collection
    BRANCH(genParticles, std::vector<TTWAnalysis::GenParticle>);
    BRANCH(gen_t, int16_t); // Index of the top quark
    BRANCH(gen_t_beforeFSR, int16_t); // Index of the top quark, before any FSR
    BRANCH(gen_tbar, int16_t); // Index of the anti-top quark
    BRANCH(gen_tbar_beforeFSR, int16_t); // Index of the anti-top quark, before any FSR
    BRANCH(gen_t_tbar_deltaR, float); // DeltaR between the top and the anti-top quark
    BRANCH(gen_t_tbar_deltaEta, float); // DeltaEta between the top and the anti-top quark
    BRANCH(gen_t_tbar_deltaPhi, float); // DeltaPhi between the top and the anti-top quark

    BRANCH(gen_W, int16_t); // Index of the associated W-boson
    BRANCH(gen_W_beforeFSR, int16_t); // Index of the associated W-boson
    BRANCH(gen_W_t, int16_t); // Index of the W-boson from the top decay chain
    BRANCH(gen_W_t_beforeFSR, int16_t); // Index of the W-boson from the top decay chain
    BRANCH(gen_W_tbar, int16_t); // Index of the W-boson from the anti-top decay chain
    BRANCH(gen_W_tbar_beforeFSR, int16_t); // Index of the W-boson from the anti-top decay chain

    BRANCH(gen_b, int16_t); // Index of the b quark coming from the top decay
    BRANCH(gen_b_beforeFSR, int16_t); // Index of the b quark coming from the top decay, before any FSR
    BRANCH(gen_bbar, int16_t); // Index of the anti-b quark coming from the anti-top decay
    BRANCH(gen_bbar_beforeFSR, int16_t); // Index of the anti-b quark coming from the anti-top decay, before any FSR
    BRANCH(gen_b_bbar_deltaR, float); // DeltaR between the b and the anti-b quark

    BRANCH(gen_jet1_t, int16_t); // Index of the first jet from the top decay chain
    BRANCH(gen_jet1_t_beforeFSR, int16_t); // Index of the first jet from the top decay chain, before any FSR
    BRANCH(gen_jet2_t, int16_t); // Index of the second jet from the top decay chain
    BRANCH(gen_jet2_t_beforeFSR, int16_t); // Index of the second jet from the top decay chain, before any FSR

    BRANCH(gen_jet1_tbar, int16_t); // Index of the first jet from the anti-top decay chain
    BRANCH(gen_jet1_tbar_beforeFSR, int16_t); // Index of the first jet from the anti-top decay chain, before any FSR
    BRANCH(gen_jet2_tbar, int16_t); // Index of the second jet from the anti-top decay chain
    BRANCH(gen_jet2_tbar_beforeFSR, int16_t); // Index of the second jet from the anti-top decay chain, before any FSR

    BRANCH(gen_jet1_W, int16_t); // Index of the first jet from the associated W decay chain
    BRANCH(gen_jet1_W_beforeFSR, int16_t); // Index of the first jet from the associated W decay chain, before any FSR
    BRANCH(gen_jet2_W, int16_t); // Index of the second jet from the associated W decay chain
    BRANCH(gen_jet2_W_beforeFSR, int16_t); // Index of the second jet from the associated W decay chain, before any FSR

    BRANCH(gen_lepton_t, int16_t); // Index of the lepton from the top decay chain
    BRANCH(gen_lepton_t_beforeFSR, int16_t); // Index of the lepton from the top decay chain, before any FSR
    BRANCH(gen_neutrino_t, int16_t); // Index of the neutrino from the top decay chain
    BRANCH(gen_neutrino_t_beforeFSR, int16_t); // Index of the neutrino from the top decay chain, before any FSR

    BRANCH(gen_lepton_tbar, int16_t); // Index of the lepton from the anti-top decay chain
    BRANCH(gen_lepton_tbar_beforeFSR, int16_t); // Index of the lepton from the anti-top decay chain, before any FSR
    BRANCH(gen_neutrino_tbar, int16_t); // Index of the neutrino from the anti-top decay chain
    BRANCH(gen_neutrino_tbar_beforeFSR, int16_t); // Index of the neutrino from the anti-top decay chain, before any FSR

    BRANCH(gen_lepton_W, int16_t); // Index of the lepton from the associated W decay chain
    BRANCH(gen_lepton_W_beforeFSR, int16_t); // Index of the lepton from the associated W decay chain, before any FSR
    BRANCH(gen_neutrino_W, int16_t); // Index of the neutrino from the associated W decay chain
    BRANCH(gen_neutrino_W_beforeFSR, int16_t); // Index of the neutrino from the associated W decay chain, before any FSR

    BRANCH(gen_ttbar_decay_type, char); // Type of ttbar decay. Can take any values from TTDecayType enum
    BRANCH(gen_W_decay_type, char); // Type of W decay. Can take any values from WDecayType enum
    BRANCH(gen_ttbarW_decay_type, char); // Type of ttbar-W decay. Can take any values from TTWDecayType enum
    BRANCH(gen_ttbarW_flavour_type, char); // Type of ttbar-W decay. Can take any values from TTWLeptonCombo enum

    BRANCH(gen_ttbar_beforeFSR_p4, reco::Candidate::LorentzVector);
    BRANCH(gen_ttbar_p4, reco::Candidate::LorentzVector);
    BRANCH(gen_W_p4, reco::Candidate::LorentzVector);

    BRANCH(gen_forwardbackward_top, short);
    BRANCH(gen_forwardbackward_beauty, short);
    BRANCH(gen_forwardbackward_lepton, short);

    // Matching for the dileptonic case

    BRANCH(gen_b_lepton_t_deltaR, float); // DeltaR between the b quark and the lepton coming from the top decay chain
    BRANCH(gen_bbar_lepton_tbar_deltaR, float); // DeltaR between the b quark and the lepton coming from the top decay chain

    // These two vectors are indexed wrt LepLepId enum
    boost::container::flat_map<std::string,std::vector<float>*> gen_b_deltaR; // DeltaR between the gen b coming from the top decay and each selected jets. Indexed as `selectedJets_tightID_DRcut` array
    boost::container::flat_map<std::string,std::vector<float>*> gen_bbar_deltaR; // DeltaR between the gen bbar coming from the anti-top decay chain and each selected jets. Indexed as `selectedJets_tightID_DRcut` array

    // These two vectors are indexed wrt LepLepId enum
    boost::container::flat_map<std::string,std::vector<float>*> gen_b_beforeFSR_deltaR; // DeltaR between the gen b coming from the top decay and each selected jets. Indexed as `selectedJets_tightID_DRcut` array
    boost::container::flat_map<std::string,std::vector<float>*> gen_bbar_beforeFSR_deltaR; // DeltaR between the gen bbar coming from the anti-top decay chain and each selected jets. Indexed as `selectedJets_tightID_DRcut` array

    BRANCH(gen_lepton_W_deltaR, std::vector<float>); // DeltaR between the gen lepton coming from the anti-top decay chain and each selected lepton. Indexed as `leptons` array
    BRANCH(gen_lepton_t_deltaR, std::vector<float>); // DeltaR between the gen lepton coming from the top decay chain and each selected lepton. Indexed as `leptons` array
    BRANCH(gen_lepton_tbar_deltaR, std::vector<float>); // DeltaR between the gen lepton coming from the anti-top decay chain and each selected lepton. Indexed as `leptons` array

    // These two vectors are indexed wrt LepLepId enum
    boost::container::flat_map<std::string,int16_t*> gen_matched_b; // Index inside the `selectedJets_tightID_DRcut` collection of the jet with the smallest deltaR with the gen b coming from the top decay
    boost::container::flat_map<std::string,int16_t*> gen_matched_bbar; // Index inside the `selectedJets_tightID_DRcut` collection of the jet with the smallest deltaR with the gen bbar coming from the anti-top decay

    // These two vectors are indexed wrt LepLepId enum
    boost::container::flat_map<std::string,int16_t*> gen_matched_b_beforeFSR; // Index inside the `selectedJets_tightID_DRcut` collection of the jet with the smallest deltaR with the gen b coming from the top decay
    boost::container::flat_map<std::string,int16_t*> gen_matched_bbar_beforeFSR; // Index inside the `selectedJets_tightID_DRcut` collection of the jet with the smallest deltaR with the gen bbar coming from the anti-top decay

    BRANCH(gen_matched_lepton_t, int16_t); // Index inside the `leptons` collection of the lepton with the smallest deltaR with the gen lepton coming from the top decay chain
    BRANCH(gen_matched_lepton_tbar, int16_t); // Index inside the `leptons` collection of the lepton with the smallest deltaR with the gen lepton coming from the anti-top decay chain
    BRANCH(gen_matched_lepton_W, int16_t); // Index inside the `leptons` collection of the lepton with the smallest deltaR with the gen lepton coming from the associated W decay chain
};

////////////////////////////////////////////////////////////////////////////////
// Implementation of TTWTruthAnalyzer                                          //
////////////////////////////////////////////////////////////////////////////////

#include "cp3_llbb/Framework/interface/GenParticlesProducer.h"
#include "cp3_llbb/TTWAnalysis/interface/CandidatesProducer.h"
#include "cp3_llbb/TTWAnalysis/interface/TTWAnalyzer.h"

TTWTruthAnalyzer::~TTWTruthAnalyzer() {}

void TTWTruthAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup, const ProducersManager& producers, const AnalyzersManager& analyzers, const CategoryManager& categories)
{
  using namespace TTWAnalysis;
  // To access VectorUtil::DeltaR() more easily
  using namespace ROOT::Math;
  ///////////////////////////
  //       GEN INFO        //
  ///////////////////////////

  #ifdef _TT_DEBUG_
  std::cout << "Generator" << std::endl;
  #endif

  if (event.isRealData())
    return;

  const GenParticlesProducer& gen_particles = producers.get<GenParticlesProducer>("gen_particles");
  const TTWAnalyzer& ttW_ana = analyzers.get<TTWAnalyzer>(m_ttWName);
  auto leptons = ttW_ana.getObjList<Lepton>();
  auto jets = ttW_ana.getPtrList<pat::Jet>();

  // 'Pruned' particles are from the hard process
  // 'Packed' particles are stable particles

  std::function<bool(size_t, size_t)> pruned_decays_from = [&pruned_decays_from, &gen_particles](size_t particle_index, size_t mother_index) -> bool {
    // Iterator over all pruned particles to find if the particle `particle_index` has `mother_index` in its decay history
    if (gen_particles.pruned_mothers_index[particle_index].empty())
      return false;

    size_t index = gen_particles.pruned_mothers_index[particle_index][0];

    if (index == mother_index) {
      return true;
    }

    if (pruned_decays_from(index, mother_index))
      return true;

    return false;
  };

#if TT_GEN_DEBUG
  std::function<void(size_t)> print_mother_chain = [&gen_particles, &print_mother_chain](size_t p) {

    if (gen_particles.pruned_mothers_index[p].empty()) {
      std::cout << std::endl;
      return;
    }

    size_t index = gen_particles.pruned_mothers_index[p][0];
      std::cout << " <- #" << index << "(" << gen_particles.pruned_pdg_id[index] << ")";
      print_mother_chain(index);
  };

  for ( std::size_t i{0}; gen_particles.pruned_pdg_id.size() !=i; ++i ) {
    std::cout << "#" << i << "  " << gen_particles.pruned_pdg_id[i] << " (";
    GenStatusFlags flags{gen_particles.pruned_status_flags[i]};
    if ( flags.fromHardProcess() ) { std::cout << "fromHard, "; }
    if ( flags.fromHardProcessBeforeFSR() ) { std::cout << "fromHardBeforeFSR, "; }
    if ( flags.isHardProcess() ) { std::cout << "isHard, "; }
    if ( flags.isFirstCopy() ) { std::cout << "first, "; }
    if ( flags.isLastCopy() ) { std::cout << "last, "; }
    if ( flags.isLastCopyBeforeFSR() ) { std::cout << "lastBeforeFSR, "; }
    std::cout << ") ";
    print_mother_chain(i);
  }
#endif

  // We need to initialize everything to -1, since 0 is a valid entry in the tt_genParticles array
  gen_t = -1;
  gen_t_beforeFSR = -1;
  gen_tbar = -1;
  gen_tbar_beforeFSR = -1;
  gen_W = -1;
  gen_W_beforeFSR = -1;
  gen_W_t = -1;
  gen_W_t_beforeFSR = -1;
  gen_W_tbar = -1;
  gen_W_tbar_beforeFSR = -1;

  gen_b = -1;
  gen_b_beforeFSR = -1;
  gen_bbar = -1;
  gen_bbar_beforeFSR = -1;

  gen_jet1_t = -1;
  gen_jet1_t_beforeFSR = -1;
  gen_jet1_tbar = -1;
  gen_jet1_tbar_beforeFSR = -1;
  gen_jet1_W = -1;
  gen_jet1_W_beforeFSR = -1;
  gen_jet2_t = -1;
  gen_jet2_t_beforeFSR = -1;
  gen_jet2_tbar = -1;
  gen_jet2_tbar_beforeFSR = -1;
  gen_jet2_W = -1;
  gen_jet2_W_beforeFSR = -1;

  gen_lepton_t = -1;
  gen_lepton_t_beforeFSR = -1;
  gen_neutrino_t = -1;
  gen_neutrino_t_beforeFSR = -1;
  gen_lepton_tbar = -1;
  gen_lepton_tbar_beforeFSR = -1;
  gen_neutrino_tbar = -1;
  gen_neutrino_tbar_beforeFSR = -1;
  gen_lepton_W = -1;
  gen_lepton_W_beforeFSR = -1;
  gen_neutrino_W = -1;
  gen_neutrino_W_beforeFSR = -1;

  gen_matched_lepton_t = -1;
  gen_matched_lepton_tbar = -1;
  gen_matched_lepton_W = -1;

  size_t gen_index(0);
  for (size_t i = 0; i < gen_particles.pruned_pdg_id.size(); i++) {

    int16_t pdg_id = gen_particles.pruned_pdg_id[i];
    uint16_t a_pdg_id = std::abs(pdg_id);

    // We only care of particles with PDG id <= 16 (16 is neutrino tau) and weak vector bosons
    if ( (a_pdg_id > 16) && ( ( a_pdg_id < 23 ) || ( a_pdg_id > 24 ) ) )
      continue;

    GenStatusFlags flags(gen_particles.pruned_status_flags[i]);

    if (! flags.isLastCopy() && ! flags.isFirstCopy())
      continue;

    if (! flags.fromHardProcess())
      continue;

#if TT_GEN_DEBUG
    std::cout << "---" << std::endl;
    std::cout << "Gen particle #" << i << ": PDG id: " << gen_particles.pruned_pdg_id[i];
    print_mother_chain(i);
    // flags.dump(); // TODO put back, but this is a bit too verbose
#endif

    if (pdg_id == 6) {
      FILL_GEN_COLL(t);
      continue;
    } else if (pdg_id == -6) {
      FILL_GEN_COLL(tbar);
      continue;
    } else if ( (gen_t_beforeFSR == -1) && (gen_tbar_beforeFSR == -1) && (gen_W_beforeFSR == -1) && flags.fromHardProcess() && flags.isHardProcess() && flags.isFirstCopy() ) {
      // rely on the fact that the W (or its products) will be inserted before the top quarks (and additional hard partons after)
#if TT_GEN_DEBUG
      std::cout << "W chain decay" << std::endl;
#endif
      if (a_pdg_id ==24) {
        FILL_GEN_COLL(W);
      } else if (a_pdg_id >= 1 && a_pdg_id <= 5) {
        FILL_GEN_COLL2(jet1_W, jet2_W, "Error: more than two quarks coming from W decay");
      } else if (a_pdg_id == 11 || a_pdg_id == 13 || a_pdg_id == 15) {
        FILL_GEN_COLL(lepton_W);
      } else if (a_pdg_id == 12 || a_pdg_id == 14 || a_pdg_id == 16) {
        FILL_GEN_COLL(neutrino_W);
      } else {
        std::cout << "Error: unknown particle coming from W decay - #" << i << " ; PDG Id: " << pdg_id << std::endl;
      }
      continue;
    }

    if ( ( gen_W_beforeFSR != -1 )
      || ( ( gen_lepton_W_beforeFSR   != -1 ) && ( gen_neutrino_W_beforeFSR != -1 ) )
      || ( ( gen_jet1_W_beforeFSR     != -1 ) && ( gen_jet2_W_beforeFSR     != -1 ) ) )
    { // try to find the before-FSR (if W itself found) and after-FSR W decay products
      bool from_W_decay{false};
      if ( gen_W_beforeFSR != -1 ) {
        from_W_decay = pruned_decays_from(i, genParticles[gen_W_beforeFSR].pruned_idx);
      } else if ( ( gen_lepton_W_beforeFSR   != -1 ) && ( gen_neutrino_W_beforeFSR != -1 ) ) {
        from_W_decay = ( pruned_decays_from(i, genParticles[gen_lepton_W_beforeFSR  ].pruned_idx)
                      || pruned_decays_from(i, genParticles[gen_neutrino_W_beforeFSR].pruned_idx) );
      } else {
        from_W_decay = ( pruned_decays_from(i, genParticles[gen_jet1_W_beforeFSR    ].pruned_idx)
                      || pruned_decays_from(i, genParticles[gen_jet2_W_beforeFSR    ].pruned_idx) );
      }
      if ( from_W_decay ) { // either first or last -> these are the two we want
#if TT_GEN_DEBUG
        std::cout << "Coming from the W chain decay" << std::endl;
#endif
        if (a_pdg_id ==24) {
          FILL_GEN_COLL(W);
        } else if (a_pdg_id >= 2 && a_pdg_id <= 5) {
          FILL_GEN_COLL2(jet1_W, jet2_W, "Error: more than two quarks coming from W decay");
        } else if (a_pdg_id == 11 || a_pdg_id == 13 || a_pdg_id == 15) {
          FILL_GEN_COLL(lepton_W);
        } else if (a_pdg_id == 12 || a_pdg_id == 14 || a_pdg_id == 16) {
          FILL_GEN_COLL(neutrino_W);
        } else {
          std::cout << "Error: unknown particle coming from W decay - #" << i << " ; PDG Id: " << pdg_id << std::endl;
        }
        continue;
      }
    }

    if ( ( gen_t != -1 ) && ( gen_tbar != -1 ) ) { // Don't bother if we don't have found the tops
      bool from_t_decay = pruned_decays_from(i, genParticles[gen_t].pruned_idx);
      bool from_tbar_decay = pruned_decays_from(i, genParticles[gen_tbar].pruned_idx);
      // Only keep particles coming from the top decays
      if ( from_t_decay || from_tbar_decay ) {
        if (pdg_id == 5) {
          // Maybe it's a b coming from the W decay
          if (!flags.isFirstCopy() && flags.isLastCopy() && gen_b == -1) {

            // This can be a B decaying from a W
            // However, we can't rely on the presence of the W in the decay chain, as it may be generator specific
            // Since it's the last copy (ie, after FSR), we can check if this B comes from the B assigned to the W decay (ie, gen_jet1_t_beforeFSR, gen_jet2_t_beforeFSR)
            // If yes, then it's not the B coming directly from the top decay
            if ((gen_jet1_t_beforeFSR != -1 && std::abs(genParticles[gen_jet1_t_beforeFSR].pdg_id) == 5) ||
                (gen_jet2_t_beforeFSR != -1 && std::abs(genParticles[gen_jet2_t_beforeFSR].pdg_id) == 5) ||
                (gen_jet1_tbar_beforeFSR != -1 && std::abs(genParticles[gen_jet1_tbar_beforeFSR].pdg_id) == 5) ||
                (gen_jet2_tbar_beforeFSR != -1 && std::abs(genParticles[gen_jet2_tbar_beforeFSR].pdg_id) == 5)) {

#if TT_GEN_DEBUG
              std::cout << "A quark coming from t-W decay is a b" << std::endl;
#endif

              if (! (gen_jet1_tbar_beforeFSR != -1 && pruned_decays_from(i, genParticles[gen_jet1_tbar_beforeFSR].pruned_idx)) &&
                  ! (gen_jet2_tbar_beforeFSR != -1 && pruned_decays_from(i, genParticles[gen_jet2_tbar_beforeFSR].pruned_idx)) &&
                  ! (gen_jet1_t_beforeFSR != -1 && pruned_decays_from(i, genParticles[gen_jet1_t_beforeFSR].pruned_idx)) &&
                  ! (gen_jet2_t_beforeFSR != -1 && pruned_decays_from(i, genParticles[gen_jet2_t_beforeFSR].pruned_idx))) {
#if TT_GEN_DEBUG
                std::cout << "This after-FSR b quark is not coming from a t-W decay" << std::endl;
#endif
                FILL_GEN_COLL(b);
                continue;
              }
#if TT_GEN_DEBUG
              else {
                std::cout << "This after-FSR b quark comes from a t-W decay" << std::endl;
              }
#endif
            } else {
#if TT_GEN_DEBUG
              std::cout << "Assigning gen_b" << std::endl;
#endif
              FILL_GEN_COLL(b);
              continue;
            }
          } else if (flags.isFirstCopy() && gen_b_beforeFSR == -1) {
            FILL_GEN_COLL(b);
            continue;
          } else {
#if TT_GEN_DEBUG
            std::cout << "This should not happen!" << std::endl;
#endif
          }
        } else if (pdg_id == -5) {
          if (!flags.isFirstCopy() && flags.isLastCopy() && gen_bbar == -1) {

            // This can be a B decaying from a W
            // However, we can't rely on the presence of the W in the decay chain, as it may be generator specific
            // Since it's the last copy (ie, after FSR), we can check if this B comes from the B assigned to the W decay (ie, gen_jet1_t_beforeFSR, gen_jet2_t_beforeFSR)
            // If yes, then it's not the B coming directly from the top decay
            if ((gen_jet1_t_beforeFSR != -1 && std::abs(genParticles[gen_jet1_t_beforeFSR].pdg_id) == 5) ||
                (gen_jet2_t_beforeFSR != -1 && std::abs(genParticles[gen_jet2_t_beforeFSR].pdg_id) == 5) ||
                (gen_jet1_tbar_beforeFSR != -1 && std::abs(genParticles[gen_jet1_tbar_beforeFSR].pdg_id) == 5) ||
                (gen_jet2_tbar_beforeFSR != -1 && std::abs(genParticles[gen_jet2_tbar_beforeFSR].pdg_id) == 5)) {

#if TT_GEN_DEBUG
              std::cout << "A quark coming from W decay is a bbar" << std::endl;
#endif

              if (! (gen_jet1_tbar_beforeFSR != -1 && pruned_decays_from(i, genParticles[gen_jet1_tbar_beforeFSR].pruned_idx)) &&
                  ! (gen_jet2_tbar_beforeFSR != -1 && pruned_decays_from(i, genParticles[gen_jet2_tbar_beforeFSR].pruned_idx)) &&
                  ! (gen_jet1_t_beforeFSR != -1 && pruned_decays_from(i, genParticles[gen_jet1_t_beforeFSR].pruned_idx)) &&
                  ! (gen_jet2_t_beforeFSR != -1 && pruned_decays_from(i, genParticles[gen_jet2_t_beforeFSR].pruned_idx))) {
#if TT_GEN_DEBUG
                std::cout << "This after-fsr b anti-quark is not coming from a W decay" << std::endl;
#endif
                FILL_GEN_COLL(bbar);
                continue;
              }
#if TT_GEN_DEBUG
              else {
                std::cout << "This after-fsr b anti-quark comes from a W decay" << std::endl;
              }
#endif
            } else {
#if TT_GEN_DEBUG
              std::cout << "Assigning gen_bbar" << std::endl;
#endif
              FILL_GEN_COLL(bbar);
              continue;
            }
          } else if (flags.isFirstCopy() && gen_bbar_beforeFSR == -1) {
              FILL_GEN_COLL(bbar);
              continue;
          }
        }

        if ((gen_tbar == -1) || (gen_t == -1))
          continue;

        if (gen_t != -1 && from_t_decay) {
#if TT_GEN_DEBUG
          std::cout << "Coming from the top chain decay" << std::endl;
#endif
          if (a_pdg_id >= 1 && a_pdg_id <= 5) {
            FILL_GEN_COLL2(jet1_t, jet2_t, "Error: more than two quarks coming from top decay");
          } else if (a_pdg_id == 11 || a_pdg_id == 13 || a_pdg_id == 15) {
            FILL_GEN_COLL(lepton_t);
          } else if (a_pdg_id == 12 || a_pdg_id == 14 || a_pdg_id == 16) {
            FILL_GEN_COLL(neutrino_t);
          } else if (a_pdg_id == 24) {
            FILL_GEN_COLL(W_t);
          } else {
            std::cout << "Error: unknown particle coming from top decay - #" << i << " ; PDG Id: " << pdg_id << std::endl;
          }
        } else if (gen_tbar != -1 && from_tbar_decay) {
#if TT_GEN_DEBUG
          std::cout << "Coming from the anti-top chain decay" << std::endl;
#endif
          if (a_pdg_id >= 1 && a_pdg_id <= 5) {
            FILL_GEN_COLL2(jet1_tbar, jet2_tbar, "Error: more than two quarks coming from anti-top decay");
          } else if (a_pdg_id == 11 || a_pdg_id == 13 || a_pdg_id == 15) {
            FILL_GEN_COLL(lepton_tbar);
          } else if (a_pdg_id == 12 || a_pdg_id == 14 || a_pdg_id == 16) {
            FILL_GEN_COLL(neutrino_tbar);
          } else if (a_pdg_id == 24) {
            FILL_GEN_COLL(W_tbar);
          } else {
            std::cout << "Error: unknown particle coming from anti-top decay - #" << i << " ; PDG Id: " << pdg_id << std::endl;
          }
        }
      }
    }
  }

  //
  // CLASSIFY TTBAR DECAY
  //
  if ( ( gen_t == -1 ) || ( gen_tbar == -1 ) ) {
#if TT_GEN_DEBUG
    std::cout << "This is not a ttbar event" << std::endl;
#endif
    gen_ttbar_decay_type = NotTT;
  } else {
    if ((gen_jet1_t != -1) && (gen_jet2_t != -1) && (gen_jet1_tbar != -1) && (gen_jet2_tbar != -1)) {
#if TT_GEN_DEBUG
      std::cout << "Hadronic ttbar decay" << std::endl;
#endif
      gen_ttbar_decay_type = Hadronic;
    } else if (
            ((gen_lepton_t != -1) && (gen_lepton_tbar == -1)) ||
            ((gen_lepton_t == -1) && (gen_lepton_tbar != -1))
            ) {

#if TT_GEN_DEBUG
      std::cout << "Semileptonic ttbar decay" << std::endl;
#endif

      uint16_t lepton_pdg_id;
      if (gen_lepton_t != -1)
        lepton_pdg_id = std::abs(genParticles[gen_lepton_t].pdg_id);
      else
        lepton_pdg_id = std::abs(genParticles[gen_lepton_tbar].pdg_id);

      if (lepton_pdg_id == 11)
        gen_ttbar_decay_type = Semileptonic_e;
      else if (lepton_pdg_id == 13)
        gen_ttbar_decay_type = Semileptonic_mu;
      else
        gen_ttbar_decay_type = Semileptonic_tau;
    } else if (gen_lepton_t != -1 && gen_lepton_tbar != -1) {
      uint16_t lepton_t_pdg_id = std::abs(genParticles[gen_lepton_t].pdg_id);
      uint16_t lepton_tbar_pdg_id = std::abs(genParticles[gen_lepton_tbar].pdg_id);

#if TT_GEN_DEBUG
      std::cout << "Dileptonic ttbar decay" << std::endl;
#endif

      if (lepton_t_pdg_id == 11 && lepton_tbar_pdg_id == 11)
        gen_ttbar_decay_type = Dileptonic_ee;
      else if (lepton_t_pdg_id == 13 && lepton_tbar_pdg_id == 13)
        gen_ttbar_decay_type = Dileptonic_mumu;
      else if (lepton_t_pdg_id == 15 && lepton_tbar_pdg_id == 15)
        gen_ttbar_decay_type = Dileptonic_tautau;
      else if (
              (lepton_t_pdg_id == 11 && lepton_tbar_pdg_id == 13) ||
              (lepton_t_pdg_id == 13 && lepton_tbar_pdg_id == 11)
              ) {
        gen_ttbar_decay_type = Dileptonic_mue;
      }
      else if (
              (lepton_t_pdg_id == 11 && lepton_tbar_pdg_id == 15) ||
              (lepton_t_pdg_id == 15 && lepton_tbar_pdg_id == 11)
              ) {
        gen_ttbar_decay_type = Dileptonic_etau;
      }
      else if (
              (lepton_t_pdg_id == 13 && lepton_tbar_pdg_id == 15) ||
              (lepton_t_pdg_id == 15 && lepton_tbar_pdg_id == 13)
              ) {
        gen_ttbar_decay_type = Dileptonic_mutau;
      } else {
        std::cout << "Error: unknown dileptonic ttbar decay." << std::endl;
        gen_ttbar_decay_type = NotTT;
        return;
      }
    } else {
      std::cout << "Error: unknown ttbar decay." << std::endl;
      gen_ttbar_decay_type = UnknownTT;
    }

    //
    // ttbar system
    //
    gen_ttbar_p4 = genParticles[gen_t].p4 + genParticles[gen_tbar].p4;
    if (gen_t_beforeFSR != -1 && gen_tbar_beforeFSR != -1)
        gen_ttbar_beforeFSR_p4 = genParticles[gen_t_beforeFSR].p4 + genParticles[gen_tbar_beforeFSR].p4;

    gen_t_tbar_deltaR = VectorUtil::DeltaR(genParticles[gen_t].p4, genParticles[gen_tbar].p4);
    gen_t_tbar_deltaEta = DeltaEta(genParticles[gen_t].p4, genParticles[gen_tbar].p4);
    gen_t_tbar_deltaPhi = VectorUtil::DeltaPhi(genParticles[gen_t].p4, genParticles[gen_tbar].p4);

    gen_b_bbar_deltaR = VectorUtil::DeltaR(genParticles[gen_b].p4, genParticles[gen_bbar].p4);

    //
    // Lepton matching
    //
    if (gen_ttbar_decay_type > Hadronic) {

      float min_dr_lepton_t = std::numeric_limits<float>::max();
      float min_dr_lepton_tbar = std::numeric_limits<float>::max();

      for ( size_t iL{0}; leptons.size() != iL; ++iL ) {
        const auto& lepton = leptons[iL];
        if (gen_lepton_t != -1) {
          float dr = VectorUtil::DeltaR(genParticles[gen_lepton_t].p4, lepton.p4());
          if (dr < min_dr_lepton_t &&
                  (std::abs(lepton.recoCand()->pdgId()) == std::abs(genParticles[gen_lepton_t].pdg_id))) {
            min_dr_lepton_t = dr;
            gen_matched_lepton_t = iL;
          }
          gen_lepton_t_deltaR.push_back(dr);
        }

        if (gen_lepton_tbar != -1) {
          float dr = VectorUtil::DeltaR(genParticles[gen_lepton_tbar].p4, lepton.p4());
          if (dr < min_dr_lepton_tbar &&
                  (std::abs(lepton.recoCand()->pdgId()) == std::abs(genParticles[gen_lepton_tbar].pdg_id))) {
            min_dr_lepton_tbar = dr;
            gen_matched_lepton_tbar = iL;
          }
          gen_lepton_tbar_deltaR.push_back(dr);
        }
      }
    }
    //
    // Match b quarks to jets
    //
    const float MIN_DR_JETS = 0.8;
    for ( const auto& jWP : m_jetWPs ) {
      float min_dr_b = MIN_DR_JETS;
      float min_dr_bbar = MIN_DR_JETS;
      float min_dr_b_beforeFSR = MIN_DR_JETS;
      float min_dr_bbar_beforeFSR = MIN_DR_JETS;
      size_t jet_index = 0;

      int16_t local_gen_matched_b = -1;
      int16_t local_gen_matched_bbar = -1;
      int16_t local_gen_matched_b_beforeFSR = -1;
      int16_t local_gen_matched_bbar_beforeFSR = -1;
      for ( auto iJ: ttW_ana.selectedJets(jWP)) {
        auto jMom = jets[iJ]->p4();
        float dr = VectorUtil::DeltaR(genParticles[gen_b].p4, jMom);
        if (dr < min_dr_b) {
            min_dr_b = dr;
          local_gen_matched_b = jet_index;
        }
        gen_b_deltaR[jWP]->push_back(dr);

        dr = VectorUtil::DeltaR(genParticles[gen_b_beforeFSR].p4, jMom);
        if (dr < min_dr_b_beforeFSR) {
            min_dr_b_beforeFSR = dr;
          local_gen_matched_b_beforeFSR = jet_index;
        }
        gen_b_beforeFSR_deltaR[jWP]->push_back(dr);

        dr = VectorUtil::DeltaR(genParticles[gen_bbar].p4, jMom);
        if (dr < min_dr_bbar) {
            min_dr_bbar = dr;
          local_gen_matched_bbar = jet_index;
        }
        gen_bbar_deltaR[jWP]->push_back(dr);

        dr = VectorUtil::DeltaR(genParticles[gen_bbar_beforeFSR].p4, jMom);
        if (dr < min_dr_bbar_beforeFSR) {
            min_dr_bbar_beforeFSR = dr;
          local_gen_matched_bbar_beforeFSR = jet_index;
        }
        gen_bbar_beforeFSR_deltaR[jWP]->push_back(dr);

        jet_index++;
      }

      *(gen_matched_b[jWP]) = local_gen_matched_b;
      *(gen_matched_bbar[jWP]) = local_gen_matched_bbar;
      *(gen_matched_b_beforeFSR[jWP]) = local_gen_matched_b_beforeFSR;
      *(gen_matched_bbar_beforeFSR[jWP]) = local_gen_matched_bbar_beforeFSR;
    }

    if (gen_b > -1 && gen_lepton_t > -1) {
      gen_b_lepton_t_deltaR = VectorUtil::DeltaR(genParticles[gen_b].p4, genParticles[gen_lepton_t].p4);
    }

    if (gen_bbar > -1 && gen_lepton_tbar > -1) {
      gen_bbar_lepton_tbar_deltaR = VectorUtil::DeltaR(genParticles[gen_bbar].p4, genParticles[gen_lepton_tbar].p4);
    }
  }

  //
  // CLASSIFY W DECAY
  //
  if ( ( gen_W == -1 ) && ( ( gen_lepton_W == -1 ) || ( gen_neutrino_W == -1 ) ) && ( ( gen_jet1_W == -1 ) || ( gen_jet2_W == -1 ) ) ) {
#if TT_GEN_DEBUG
    std::cout << "This is not a W event" << std::endl;
#endif
    gen_W_decay_type = NotW;
  } else {
    if ( ( gen_jet1_W != -1 ) && ( gen_jet2_W != -1 ) ) {
#if TT_GEN_DEBUG
      std::cout << "Hadronic W decay" << std::endl;
#endif
      gen_W_decay_type = HadronicW;
    } else if ( ( gen_lepton_W != -1 ) && ( gen_neutrino_W != -1 ) ) {
#if TT_GEN_DEBUG
      std::cout << "Leptonic W decay" << std::endl;
#endif
      uint16_t lepton_pdg_id = std::abs(genParticles[gen_lepton_W].pdg_id);
      if ( lepton_pdg_id == 11 ) {
        gen_W_decay_type = LeptonicW_e;
      } else if ( lepton_pdg_id == 13 ) {
        gen_W_decay_type = LeptonicW_mu;
      } else {
        gen_W_decay_type = LeptonicW_tau;
      }
    } else {
      std::cout << "Error: unknown W decay." << std::endl;
      gen_W_decay_type = UnknownW;
    }

    //
    // W system
    //
    if ( gen_W != -1 ) {
      gen_W_p4 = genParticles[gen_W].p4;
    } else if ( ( gen_lepton_W_beforeFSR != -1 ) && ( gen_neutrino_W_beforeFSR != -1 ) ) {
      gen_W_p4 = genParticles[gen_lepton_W_beforeFSR].p4+genParticles[gen_neutrino_W_beforeFSR].p4;
    } else if ( ( gen_jet1_W_beforeFSR != -1 ) && ( gen_jet2_W_beforeFSR != -1 ) ) {
      gen_W_p4 = genParticles[gen_jet1_W_beforeFSR  ].p4+genParticles[gen_jet2_W_beforeFSR    ].p4;
    } else {
      gen_W_p4 = reco::Candidate::LorentzVector();
    }

    //
    // Lepton matching
    //
    if ( gen_W_decay_type > HadronicW ) {
      float min_dr_lepton_W = std::numeric_limits<float>::max();
      for ( size_t iL{0}; leptons.size() != iL; ++iL ) {
        const auto& lepton = leptons[iL];
        if (gen_lepton_W != -1) {
          float dr = VectorUtil::DeltaR(genParticles[gen_lepton_W].p4, lepton.p4());
          if (dr < min_dr_lepton_W &&
                  (std::abs(lepton.recoCand()->pdgId()) == std::abs(genParticles[gen_lepton_W].pdg_id))) {
            min_dr_lepton_W = dr;
            gen_matched_lepton_W = iL;
          }
          gen_lepton_W_deltaR.push_back(dr);
        }
      }
    }
  }

  auto get2lSignCombination = [this] ( int lW, int lt, int ltbar ) -> TTWDecayType
  {
    if ( ( lW != -1 ) && ( lt != -1 ) ) {
      if ( ltbar != -1 ) {
        std::cout << "ERROR: Called 2-lepton sign classification for a 3-lepton decay" << std::endl;
        return UnknownTTW;
      }
      return ( genParticles[lW].pdg_id*genParticles[lt].pdg_id > 0 ) ? DiLeptonWSST : DiLeptonWOST;
    } else if ( ( lW != -1 ) && ( ltbar != -1 ) ) {
      if ( lt != -1 ) {
        std::cout << "ERROR: Called 2-lepton sign classification for a 3-lepton decay" << std::endl;
        return UnknownTTW;
      }
      return ( genParticles[lW].pdg_id*genParticles[ltbar].pdg_id > 0 ) ? DiLeptonWSST : DiLeptonWOST;
    } else {
      std::cout << "ERROR: Invalid lepton combination for 2-lepton sign classification" << std::endl;
      return UnknownTTW;
    }
  };

  if ( ( gen_ttbar_decay_type > NotTT ) && ( gen_W_decay_type > NotW ) ) {
    switch (gen_ttbar_decay_type) {
      case Hadronic: switch(gen_W_decay_type) {
          case HadronicW    : gen_ttbarW_decay_type = HadronicTTW   ; gen_ttbarW_flavour_type = ZeroLepton   ; break;
          case LeptonicW_e  : gen_ttbarW_decay_type = SingleLeptonW ; gen_ttbarW_flavour_type = Single_e     ; break;
          case LeptonicW_mu : gen_ttbarW_decay_type = SingleLeptonW ; gen_ttbarW_flavour_type = Single_mu    ; break;
          case LeptonicW_tau: gen_ttbarW_decay_type = SingleLeptonW ; gen_ttbarW_flavour_type = Single_tau   ; break;
        } break;
      case Semileptonic_e: switch(gen_W_decay_type) {
          case HadronicW    : gen_ttbarW_decay_type = SingleLeptonTT; gen_ttbarW_flavour_type = Single_e     ; break;
          case LeptonicW_e  : gen_ttbarW_decay_type = get2lSignCombination(gen_lepton_W, gen_lepton_t, gen_lepton_tbar); gen_ttbarW_flavour_type = Di_ee    ; break;
          case LeptonicW_mu : gen_ttbarW_decay_type = get2lSignCombination(gen_lepton_W, gen_lepton_t, gen_lepton_tbar); gen_ttbarW_flavour_type = Di_emu   ; break;
          case LeptonicW_tau: gen_ttbarW_decay_type = get2lSignCombination(gen_lepton_W, gen_lepton_t, gen_lepton_tbar); gen_ttbarW_flavour_type = Di_etau  ; break;
        } break;
      case Semileptonic_mu: switch(gen_W_decay_type) {
          case HadronicW    : gen_ttbarW_decay_type = SingleLeptonTT; gen_ttbarW_flavour_type = Single_mu    ; break;
          case LeptonicW_e  : gen_ttbarW_decay_type = get2lSignCombination(gen_lepton_W, gen_lepton_t, gen_lepton_tbar); gen_ttbarW_flavour_type = Di_emu    ; break;
          case LeptonicW_mu : gen_ttbarW_decay_type = get2lSignCombination(gen_lepton_W, gen_lepton_t, gen_lepton_tbar); gen_ttbarW_flavour_type = Di_mumu   ; break;
          case LeptonicW_tau: gen_ttbarW_decay_type = get2lSignCombination(gen_lepton_W, gen_lepton_t, gen_lepton_tbar); gen_ttbarW_flavour_type = Di_mutau  ; break;
        } break;
      case Dileptonic_mumu: switch(gen_W_decay_type) {
          case HadronicW    : gen_ttbarW_decay_type = DiLeptonTT    ; gen_ttbarW_flavour_type = Di_mumu      ; break;
          case LeptonicW_e  : gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_emumu    ; break;
          case LeptonicW_mu : gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_mumumu   ; break;
          case LeptonicW_tau: gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_mumutau  ; break;
        } break;
      case Dileptonic_ee: switch(gen_W_decay_type) {
          case HadronicW    : gen_ttbarW_decay_type = DiLeptonTT    ; gen_ttbarW_flavour_type = Di_ee        ; break;
          case LeptonicW_e  : gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_eee      ; break;
          case LeptonicW_mu : gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_eemu     ; break;
          case LeptonicW_tau: gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_eetau    ; break;
        } break;
      case Dileptonic_mue: switch(gen_W_decay_type) {
          case HadronicW    : gen_ttbarW_decay_type = DiLeptonTT    ; gen_ttbarW_flavour_type = Di_emu       ; break;
          case LeptonicW_e  : gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_eemu     ; break;
          case LeptonicW_mu : gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_emumu    ; break;
          case LeptonicW_tau: gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_emutau   ; break;
        } break;
      case Semileptonic_tau: switch(gen_W_decay_type) {
          case HadronicW    : gen_ttbarW_decay_type = SingleLeptonTT; gen_ttbarW_flavour_type = Single_mu    ; break;
          case LeptonicW_e  : gen_ttbarW_decay_type = get2lSignCombination(gen_lepton_W, gen_lepton_t, gen_lepton_tbar); gen_ttbarW_flavour_type = Di_etau    ; break;
          case LeptonicW_mu : gen_ttbarW_decay_type = get2lSignCombination(gen_lepton_W, gen_lepton_t, gen_lepton_tbar); gen_ttbarW_flavour_type = Di_mutau   ; break;
          case LeptonicW_tau: gen_ttbarW_decay_type = get2lSignCombination(gen_lepton_W, gen_lepton_t, gen_lepton_tbar); gen_ttbarW_flavour_type = Di_tautau  ; break;
        } break;
      case Dileptonic_tautau: switch(gen_W_decay_type) {
          case HadronicW    : gen_ttbarW_decay_type = DiLeptonTT    ; gen_ttbarW_flavour_type = Di_tautau    ; break;
          case LeptonicW_e  : gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_etautau  ; break;
          case LeptonicW_mu : gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_mutautau ; break;
          case LeptonicW_tau: gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_tautautau; break;
        } break;
      case Dileptonic_mutau: switch(gen_W_decay_type) {
          case HadronicW    : gen_ttbarW_decay_type = DiLeptonTT    ; gen_ttbarW_flavour_type = Di_mutau     ; break;
          case LeptonicW_e  : gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_emutau   ; break;
          case LeptonicW_mu : gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_mumutau  ; break;
          case LeptonicW_tau: gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_mutautau ; break;
        } break;
      case Dileptonic_etau: switch(gen_W_decay_type) {
          case Hadronic     : gen_ttbarW_decay_type = DiLeptonTT    ; gen_ttbarW_flavour_type = Di_etau      ; break;
          case LeptonicW_e  : gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_eetau    ; break;
          case LeptonicW_mu : gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_emutau   ; break;
          case LeptonicW_tau: gen_ttbarW_decay_type = TriLepton     ; gen_ttbarW_flavour_type = Tri_etautau  ; break;
        } break;
    }
    //
    // Forward-backward categories
    //
    gen_forwardbackward_top    = ( std::abs(genParticles[gen_t].p4.eta()) > std::abs(std::abs(genParticles[gen_tbar].p4.eta())) ) ? 1 : -1;
    gen_forwardbackward_beauty = ( std::abs(genParticles[gen_b].p4.eta()) > std::abs(std::abs(genParticles[gen_bbar].p4.eta())) ) ? 1 : -1;
    gen_forwardbackward_lepton = ( ( gen_lepton_t != -1 ) && ( gen_lepton_tbar != -1 ) ) ? ( ( std::abs(genParticles[gen_lepton_t].p4.eta()) > std::abs(std::abs(genParticles[gen_lepton_tbar].p4.eta())) ) ? 1 : -1 ) : 0;
    // TODO ttW system (?)
  } else {
    gen_ttbarW_decay_type = NotTTW;
  }
}

void TTWTruthAnalyzer::registerCategories(CategoryManager& manager, const edm::ParameterSet& config) {}

#include "FWCore/PluginManager/interface/PluginFactory.h"
DEFINE_EDM_PLUGIN(ExTreeMakerAnalyzerFactory, TTWTruthAnalyzer, "ttwtruth_analyzer");
