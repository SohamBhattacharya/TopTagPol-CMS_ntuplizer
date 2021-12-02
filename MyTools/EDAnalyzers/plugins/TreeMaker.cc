// -*- C++ -*-
//
// Package:    EDAnalyzers/TreeMaker
// Class:      TreeMaker
//
/**\class TreeMaker TreeMaker.cc EDAnalyzers/TreeMaker/plugins/TreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Sat, 11 May 2019 13:14:55 GMT
//
//


// system include files
# include <algorithm>
# include <memory>

// user include files
# include "CommonTools/UtilAlgos/interface/TFileService.h"
# include "DataFormats/CaloTowers/interface/CaloTowerDefs.h"
# include "DataFormats/Common/interface/MapOfVectors.h"
# include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
# include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
# include "DataFormats/EgammaCandidates/interface/Photon.h"
# include "DataFormats/FWLite/interface/ESHandle.h"
# include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
# include "DataFormats/HepMCCandidate/interface/GenParticle.h"
# include "DataFormats/JetReco/interface/PFJet.h"
# include "DataFormats/Math/interface/LorentzVector.h"
# include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
# include "DataFormats/PatCandidates/interface/Jet.h"
# include "DataFormats/PatCandidates/interface/MET.h"
# include "DataFormats/PatCandidates/interface/Muon.h"
# include "DataFormats/PatCandidates/interface/Electron.h"
# include "DataFormats/PatCandidates/interface/Photon.h"
# include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
# include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
# include "DataFormats/TrackReco/interface/Track.h"
# include "DataFormats/TrackReco/interface/TrackFwd.h"
# include "DataFormats/VertexReco/interface/Vertex.h"
# include "FWCore/Framework/interface/ConsumesCollector.h"
# include "FWCore/Framework/interface/Event.h"
# include "FWCore/Framework/interface/ESHandle.h"
# include "FWCore/Framework/interface/Frameworkfwd.h"
# include "FWCore/Framework/interface/MakerMacros.h"
# include "FWCore/Framework/interface/one/EDAnalyzer.h"
# include "FWCore/ParameterSet/interface/ParameterSet.h"
# include "FWCore/ServiceRegistry/interface/Service.h"
# include "FWCore/Utilities/interface/InputTag.h"
# include "Geometry/Records/interface/CaloGeometryRecord.h"
# include "Geometry/Records/interface/IdealGeometryRecord.h"
# include "RecoEgamma/EgammaTools/interface/MVAVariableHelper.h"
# include "RecoEgamma/EgammaTools/interface/MVAVariableManager.h"
# include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
# include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
# include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
# include "SimDataFormats/CaloHit/interface/PCaloHit.h"
# include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
# include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

# include <Compression.h>
# include <Math/VectorUtil.h>
# include <TH1F.h>
# include <TH2F.h>
# include <TLorentzVector.h>
# include <TMatrixD.h>
# include <TTree.h> 
# include <TVector2.h> 
# include <TVectorD.h> 


# include <CLHEP/Matrix/Matrix.h>
# include <CLHEP/Vector/LorentzVector.h>
# include <CLHEP/Vector/ThreeVector.h>
# include <CLHEP/Vector/ThreeVector.h>

# include "fastjet/contrib/AxesDefinition.hh"
# include "fastjet/contrib/EnergyCorrelator.hh"
# include "fastjet/contrib/MeasureDefinition.hh"
# include "fastjet/contrib/Nsubjettiness.hh"
# include "fastjet/contrib/SoftDrop.hh"
# include "fastjet/PseudoJet.hh"
# include "fastjet/ClusterSequence.hh"

# include "MyTools/EDAnalyzers/interface/Common.h"
# include "MyTools/EDAnalyzers/interface/Constants.h"
# include "MyTools/EDAnalyzers/interface/TreeOutputInfo.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



class TreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
    public:
    
    explicit TreeMaker(const edm::ParameterSet&);
    ~TreeMaker();
    
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    
    private:
    
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    
    
    TreeOutputInfo::TreeOutput *treeOutput;
    
    
    // My stuff //
    bool debug;
    bool isGunSample;
    
    fastjet::Strategy fj_strategy;
    fastjet::RecombinationScheme fj_recombScheme;
    
    fastjet::JetDefinition *fj_fatJetExcSubJetDef;
    fastjet::JetDefinition *fj_akJetReclusterDef;
    
    
    // GenEventInfoProduct //
    edm::EDGetTokenT <GenEventInfoProduct> tok_generator;
    
    
    // Gen particles //
    edm::EDGetTokenT <std::vector <reco::GenParticle> > tok_genParticle;
    
    
    // PV //
    edm::EDGetTokenT <std::vector <reco::Vertex> > tok_primaryVertex;
    
    
    // SV //
    edm::EDGetTokenT <std::vector <reco::VertexCompositePtrCandidate> > tok_secondaryVertex;
    
    
    // Pileup //
    edm::EDGetTokenT <std::vector <PileupSummaryInfo> > tok_pileup;
    
    
    // Rho //
    edm::EDGetTokenT <double> tok_rho;
    
    
    // Electrons //
    edm::EDGetTokenT <std::vector <pat::Electron> > tok_electron_reco;
    
    std::string eleMvaVariablesFile;
    //MVAVariableManager <pat::Electron> eleMvaVarManager;
    MVAVariableManager <reco::GsfElectron> eleMvaVarManager;
    
    //MVAVariableHelper <pat::Electron> eleMvaVarHelper;
    MVAVariableHelper <reco::GsfElectron> eleMvaVarHelper;
    
    
    // Muons //
    edm::EDGetTokenT <std::vector <pat::Muon> > tok_muon_reco;
    
    
    // Jets //
    std::vector <edm::InputTag> v_jetCollectionTag;
    std::vector <std::string> v_jetCollectionName;
    std::vector <edm::EDGetTokenT <std::vector <pat::Jet> > > v_tok_jet_reco;
    
    std::vector <edm::ParameterSet> v_recoJetPSet;
    
    
    // Other stuff //
    
    
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TreeMaker::TreeMaker(const edm::ParameterSet& iConfig) :
    eleMvaVarManager(iConfig.getParameter <std::string>("eleMvaVariablesFile")),
    eleMvaVarHelper(consumesCollector())
{
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    
    // Compression
    //fs->file().SetCompressionAlgorithm(ROOT::kLZMA);
    //fs->file().SetCompressionLevel(8);
    
    
    //now do what ever initialization is needed
    
    treeOutput = new TreeOutputInfo::TreeOutput("tree", fs);
    
    
    // My stuff //
    debug = iConfig.getParameter <bool>("debug");
    isGunSample = iConfig.getParameter <bool>("isGunSample");
    
    fj_strategy = fastjet::Best;
    fj_recombScheme = fastjet::E_scheme;
    
    //fj_jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, _jetR, fj_recombScheme, fj_strategy);
    //fj_fatJetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm, _fatJetR, fj_recombScheme, fj_strategy);
    fj_fatJetExcSubJetDef = new fastjet::JetDefinition(fastjet::kt_algorithm, 1.0, fj_recombScheme, fj_strategy);
    fj_akJetReclusterDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, 1000.0, fj_recombScheme, fj_strategy);
    
    
    // GenEventInfoProduct //
    tok_generator = consumes <GenEventInfoProduct>(iConfig.getParameter <edm::InputTag>("label_generator"));

    
    // Gen particles //
    tok_genParticle = consumes <std::vector <reco::GenParticle> >(iConfig.getParameter <edm::InputTag>("label_genParticle"));
    
    
    // PV //
    tok_primaryVertex = consumes <std::vector <reco::Vertex> >(iConfig.getParameter <edm::InputTag>("label_primaryVertex"));
    
    
    // SV //
    tok_secondaryVertex = consumes <std::vector <reco::VertexCompositePtrCandidate> >(iConfig.getParameter <edm::InputTag>("label_secondaryVertex"));
    
    
    // Pileup //
    tok_pileup = consumes <std::vector <PileupSummaryInfo> >(iConfig.getParameter <edm::InputTag>("label_pileup"));
    
    
    // Rho //
    tok_rho = consumes <double>(iConfig.getParameter <edm::InputTag>("label_rho"));
    
    
    // Electrons //
    tok_electron_reco = consumes <std::vector <pat::Electron> >(iConfig.getParameter <edm::InputTag>("label_slimmedElectron"));
    
    //eleMvaVariablesFile = iConfig.getParameter <std::string>("eleMvaVariablesFile");
    //eleMvaVarManager(eleMvaVariablesFile, MVAVariableHelper::indexMap());
    
    
    // Muons //
    tok_muon_reco = consumes <std::vector <pat::Muon> >(iConfig.getParameter <edm::InputTag>("label_slimmedMuon"));
    
    
    // Jets //
    v_recoJetPSet = iConfig.getParameter <std::vector <edm::ParameterSet> >("v_recoJetPSet");
    
    for(const auto &jetPSet : v_recoJetPSet)
    {
        edm::ConsumesCollector ccollector = consumesCollector();
        
        std::string jetName = treeOutput->addJetInfo(jetPSet, ccollector, treeOutput->tree, eleMvaVarManager);
        
        v_jetCollectionName.push_back(jetName);
    }
    
    
    // Other stuff
    //mvaVarHelper(consumesCollector());
}


TreeMaker::~TreeMaker()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    
    delete treeOutput;
}


//
// member functions
//


// ------------ method called for each event  ------------
void TreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    
    long long eventNumber = iEvent.id().event();
    //printf("Event %llu \n", eventNumber);
    
    
    treeOutput->clear();
    
    //recHitTools.getEventSetup(iSetup);
    
    //////////////////// Run info ////////////////////
    treeOutput->runNumber = iEvent.id().run();
    treeOutput->eventNumber = iEvent.id().event();
    treeOutput->luminosityNumber = iEvent.id().luminosityBlock();
    treeOutput->bunchCrossingNumber = iEvent.bunchCrossing();
    
    
    
    //////////////////// GenEventInfoProduct ////////////////////
    edm::Handle <GenEventInfoProduct> generatorHandle;
    iEvent.getByToken(tok_generator, generatorHandle);
    GenEventInfoProduct generator = *generatorHandle;
    
    if(debug)
    {
        printf("[%llu] Gen. evt. wt. %0.4g \n", eventNumber, generator.weight());
    }
    
    treeOutput->genEventWeight = generator.weight();
    
    
    //////////////////// Gen particle ////////////////////
    edm::Handle <std::vector <reco::GenParticle> > v_genParticle;
    iEvent.getByToken(tok_genParticle, v_genParticle);
    
    
    std::vector <CLHEP::HepLorentzVector> v_genTopVis_4mom;
    std::vector <CLHEP::HepLorentzVector> v_genWvis_4mom;
    std::vector <CLHEP::HepLorentzVector> v_genZ_4mom;
    
    for(int iPart = 0; iPart < (int) v_genParticle->size(); iPart++)
    {
        reco::GenParticle part = v_genParticle->at(iPart);
        
        int pdgId = part.pdgId();
        int status = part.status();
        
        if(std::abs(pdgId) == 6 && part.isLastCopy())
        {
            treeOutput->v_genTop_id.push_back(pdgId);
            treeOutput->v_genTop_E.push_back(part.energy());
            treeOutput->v_genTop_px.push_back(part.px());
            treeOutput->v_genTop_py.push_back(part.py());
            treeOutput->v_genTop_pz.push_back(part.pz());
            treeOutput->v_genTop_pT.push_back(part.pt());
            treeOutput->v_genTop_eta.push_back(part.eta());
            treeOutput->v_genTop_y.push_back(part.y());
            treeOutput->v_genTop_phi.push_back(part.phi());
            treeOutput->v_genTop_m.push_back(part.mass());
            
            int isLeptonicTop = Common::isLeptonicTopWZ(&part);
            
            treeOutput->v_genTop_isLeptonic.push_back(isLeptonicTop);
            
            CLHEP::HepLorentzVector genTopVis_4mom = Common::lorentzVector2clhep(Common::getVisibleCompnent(&part));
            v_genTopVis_4mom.push_back(genTopVis_4mom);
            
            if(debug)
            {
                printf(
                    "[%llu] "
                    "Gen t (isLeptonic %d) found (orig/vis): E %0.2f/%0.2f, pT %0.2f/%0.2f, eta %+0.2f/%+0.2f, phi %+0.2f/%+0.2f, "
                    "\n",
                    eventNumber,
                    (int) isLeptonicTop,
                    part.energy(), genTopVis_4mom.e(),
                    part.pt(), genTopVis_4mom.perp(),
                    part.eta(), genTopVis_4mom.eta(),
                    part.phi(), genTopVis_4mom.phi()
                );
            }
            
            treeOutput->v_genTopVis_E.push_back(genTopVis_4mom.e());
            treeOutput->v_genTopVis_px.push_back(genTopVis_4mom.px());
            treeOutput->v_genTopVis_py.push_back(genTopVis_4mom.py());
            treeOutput->v_genTopVis_pz.push_back(genTopVis_4mom.pz());
            treeOutput->v_genTopVis_pT.push_back(genTopVis_4mom.perp());
            treeOutput->v_genTopVis_eta.push_back(genTopVis_4mom.eta());
            treeOutput->v_genTopVis_y.push_back(genTopVis_4mom.rapidity());
            treeOutput->v_genTopVis_phi.push_back(genTopVis_4mom.phi());
            treeOutput->v_genTopVis_m.push_back(genTopVis_4mom.m());
        }
        
        
        if(std::abs(pdgId) == 24 && part.isLastCopy())
        {
            treeOutput->v_genW_id.push_back(pdgId);
            treeOutput->v_genW_E.push_back(part.energy());
            treeOutput->v_genW_px.push_back(part.px());
            treeOutput->v_genW_py.push_back(part.py());
            treeOutput->v_genW_pz.push_back(part.pz());
            treeOutput->v_genW_pT.push_back(part.pt());
            treeOutput->v_genW_eta.push_back(part.eta());
            treeOutput->v_genW_y.push_back(part.y());
            treeOutput->v_genW_phi.push_back(part.phi());
            treeOutput->v_genW_m.push_back(part.mass());
            
            int isLeptonic = Common::isLeptonicTopWZ(&part);
            
            treeOutput->v_genW_isLeptonic.push_back(isLeptonic);
            
            CLHEP::HepLorentzVector genWvis_4mom = Common::lorentzVector2clhep(Common::getVisibleCompnent(&part));
            v_genWvis_4mom.push_back(genWvis_4mom);
            
            
            if(debug)
            {
                printf(
                    "[%llu] "
                    "Gen W (isLeptonic %d) found (orig/vis): E %0.2f/%0.2f, pT %0.2f/%0.2f, eta %+0.2f/%+0.2f, phi %+0.2f/%+0.2f, "
                    "\n",
                    eventNumber,
                    (int) isLeptonic,
                    part.energy(), genWvis_4mom.e(),
                    part.pt(), genWvis_4mom.perp(),
                    part.eta(), genWvis_4mom.eta(),
                    part.phi(), genWvis_4mom.phi()
                );
            }
            
            
            treeOutput->v_genWvis_E.push_back(genWvis_4mom.e());
            treeOutput->v_genWvis_px.push_back(genWvis_4mom.px());
            treeOutput->v_genWvis_py.push_back(genWvis_4mom.py());
            treeOutput->v_genWvis_pz.push_back(genWvis_4mom.pz());
            treeOutput->v_genWvis_pT.push_back(genWvis_4mom.perp());
            treeOutput->v_genWvis_eta.push_back(genWvis_4mom.eta());
            treeOutput->v_genWvis_y.push_back(genWvis_4mom.rapidity());
            treeOutput->v_genWvis_phi.push_back(genWvis_4mom.phi());
            treeOutput->v_genWvis_m.push_back(genWvis_4mom.m());
        }
        
        
        if(std::abs(pdgId) == 23 && part.isLastCopy())
        {
            treeOutput->v_genZ_id.push_back(pdgId);
            treeOutput->v_genZ_E.push_back(part.energy());
            treeOutput->v_genZ_px.push_back(part.px());
            treeOutput->v_genZ_py.push_back(part.py());
            treeOutput->v_genZ_pz.push_back(part.pz());
            treeOutput->v_genZ_pT.push_back(part.pt());
            treeOutput->v_genZ_eta.push_back(part.eta());
            treeOutput->v_genZ_y.push_back(part.y());
            treeOutput->v_genZ_phi.push_back(part.phi());
            treeOutput->v_genZ_m.push_back(part.mass());
            
            int isLeptonic = Common::isLeptonicTopWZ(&part);
            
            treeOutput->v_genZ_isLeptonic.push_back(isLeptonic);
            
            CLHEP::HepLorentzVector genZ_4mom = Common::lorentzVector2clhep(part.p4());
            v_genZ_4mom.push_back(genZ_4mom);
            
            if(debug)
            {
                printf(
                    "[%llu] "
                    "Gen Z (isLeptonic %d) found: E %0.2f, pT %0.2f, eta %+0.2f, phi %+0.2f, "
                    "\n",
                    eventNumber,
                    (int) isLeptonic,
                    part.energy(),
                    part.pt(),
                    part.eta(),
                    part.phi()
                );
            }
        }
    }
    
    treeOutput->genTop_n = treeOutput->v_genTop_id.size();
    treeOutput->genW_n = treeOutput->v_genW_id.size();
    treeOutput->genZ_n = treeOutput->v_genZ_id.size();
    
    
    //// Sort the electrons
    //std::sort(
    //    treeOutput->v_genEl_HGCalEEP_EsortedIndex.begin(), treeOutput->v_genEl_HGCalEEP_EsortedIndex.end(),
    //    [&](int iEle1, int iEle2)
    //    {
    //        return (treeOutput->v_genEl_E[iEle1] > treeOutput->v_genEl_E[iEle2]);
    //    }
    //);
    
    
    // PV
    edm::Handle <std::vector <reco::Vertex> > handle_primaryVertex;
    iEvent.getByToken(tok_primaryVertex, handle_primaryVertex);
    const auto &v_primaryVertex = *handle_primaryVertex;
    
    int iPV = -1;
    int pmVtx_idx = -1;
    int nGoodVertex = 0;
    reco::Vertex pmVtx;
    
    //for(int iVtx = 0; iVtx < (int) v_primaryVertex->size(); iVtx++)
    for(const auto &vertex : v_primaryVertex)
    {
        iPV++;
        
        bool isGoodVertex = (
            !vertex.isFake() &&
            vertex.ndof() >= 4 &&
            fabs(vertex.z()) <= 24.0 &&
            fabs(vertex.position().rho()) <= 2.0
        );
        
        nGoodVertex += (int) isGoodVertex;
        
        // The first good vertex
        if(pmVtx_idx < 0 && nGoodVertex == 1)
        {
            pmVtx_idx = iPV;
            pmVtx = vertex;
        }
    }
    
    
    // SV
    edm::Handle <std::vector <reco::VertexCompositePtrCandidate> > handle_secondaryVertex;
    iEvent.getByToken(tok_secondaryVertex, handle_secondaryVertex);
    const auto &v_secondaryVertex = *handle_secondaryVertex;
    
    int iSV = -1;
    int nSecVtx = 0;
    
    std::vector <math::XYZVectorD> v_secVtxDir;
    
    for(const auto &sv : v_secondaryVertex)
    {
        iSV++;
        
        //VertexDistance3D vertTool;
        //double distance = vertTool.distance(sv, pmVtx).value();
        //double dist_err = vertTool.distance(sv, pmVtx).error();
        
        //if(nGoodVertex && sv.vertexNdof() >= 4)
        if(nGoodVertex)
        {
            math::XYZVectorD xyz_secVtxDir = sv.position() - pmVtx.position();
            
            v_secVtxDir.push_back(xyz_secVtxDir);
            
            if(debug)
            {
                printf(
                    "[%llu] "
                    "SV (%d) found: "
                    "(x, y, z) (%0.4f, %0.4f, %0.4f), "
                    "(vx, vy, vz) (%0.4f, %0.4f, %0.4f), "
                    "flight dist %f, ndof %f, chi2/ndof %f, "
                    "\n",
                    eventNumber,
                    iSV,
                    sv.position().x(), sv.position().y(), sv.position().z(),
                    sv.vx(), sv.vy(), sv.vz(),
                    std::sqrt(xyz_secVtxDir.mag2()), sv.vertexNdof(), sv.vertexNormalizedChi2()
                );
            }
        }
        
        nSecVtx++;
    }
    
    treeOutput->sv_n = nSecVtx;
    
    
    // Pileup
    edm::Handle <std::vector <PileupSummaryInfo> > pileUps_reco;
    iEvent.getByToken(tok_pileup, pileUps_reco);
    treeOutput->pileup_n = Common::getPileup(pileUps_reco);
    
    
    // Rho
    edm::Handle <double> handle_rho;
    iEvent.getByToken(tok_rho, handle_rho);
    double rho = *handle_rho;
    
    treeOutput->rho = rho;
    
    
    // Electrons
    edm::Handle <std::vector <pat::Electron> > handle_electron_reco;
    iEvent.getByToken(tok_electron_reco, handle_electron_reco);
    auto electrons_reco = *handle_electron_reco;
    
    int nEle = 0;
    
    std::vector <CLHEP::HepLorentzVector> v_ele_4mom;
    
    for(const auto &ele : electrons_reco)
    {
        edm::Ptr <pat::Electron> elePtr(handle_electron_reco, nEle);
        
        nEle++;
        
        CLHEP::HepLorentzVector ele_4mom = Common::lorentzVector2clhep(ele.p4());
        v_ele_4mom.push_back(ele_4mom);
        
        auto sc = ele.superCluster();
        
        if(debug)
        {
            printf(
                "[%llu] "
                "ele (%d) found: E %0.2f, pT %0.2f, eta %+0.2f, phi %+0.2f, "
                "SC E %0.2f, eta %+0.2f, phi %+0.2f"
                "\n",
                eventNumber,
                nEle,
                ele.energy(), ele.pt(), ele.eta(), ele.phi(),
                sc->energy(), sc->eta(), sc->phi()
            );
        }
        
        
        //std::vector <float> extraVariables = eleMvaVarHelper.getAuxVariables(elePtr, iEvent);
        //
        //for(int iVar = 0; iVar < eleMvaVarManager.getNVars(); iVar++)
        //{
        //    printf("%s %0.4f, ", eleMvaVarManager.getName(iVar).c_str(), eleMvaVarManager.getValue(iVar, ele, extraVariables));
        //}
        //printf("\n");
    }
    
    
    // Muons
    edm::Handle <std::vector <pat::Muon> > handle_muon_reco;
    iEvent.getByToken(tok_muon_reco, handle_muon_reco);
    auto muons_reco = *handle_muon_reco;
    
    int nMu = 0;
    
    std::vector <CLHEP::HepLorentzVector> v_mu_4mom;
    
    for(const auto &mu : muons_reco)
    {
        nMu++;
        
        if(debug)
        {
            printf(
                "[%llu] "
                "mu (%d) found: E %0.2f, pT %0.2f, eta %+0.2f, phi %+0.2f, "
                "\n",
                eventNumber,
                nMu,
                mu.energy(), mu.pt(), mu.eta(), mu.phi()
            );
        }
    }
    
    
    int iJetCollection = -1;
    
    for(const auto &jetName : v_jetCollectionName)
    {
        iJetCollection++;
        
        if(treeOutput->m_jetInfo.find(jetName) == treeOutput->m_jetInfo.end())
        {
            printf("Error in analyze(...): Key \"%s\" not found in JetInfo map. Add it to the map first. Exiting...", jetName.c_str());
            exit(EXIT_FAILURE);
        }
        
        const auto &jetInfo = treeOutput->m_jetInfo[jetName];
        
        edm::Handle <std::vector <pat::Jet> > handle_jet_reco;
        //iEvent.getByToken(v_tok_jet_reco.at(iJetCollection), handle_jet_reco);
        iEvent.getByToken(jetInfo->tok_jet, handle_jet_reco);
        auto jets_reco = *handle_jet_reco;
        
        int nJet = 0;
        
        std::vector <CLHEP::HepLorentzVector> v_jet_4mom;
        
        for(const auto &jet : jets_reco)
        {
            if(jet.pt() < jetInfo->minPt)
            {
                continue;
            }
            
            nJet++;
            
            if(debug)
            {
                printf(
                    "[%llu] "
                    "%s jet (%d) found: E %0.2f, pT %0.2f, eta %+0.2f, phi %+0.2f, "
                    "isPF %d, "
                    "nConsti %d, nDaughter %d, "
                    "\n",
                    eventNumber,
                    jetName.c_str(),
                    nJet,
                    jet.energy(), jet.pt(), jet.eta(), jet.phi(),
                    (int) jet.isPFJet(),
                    (int) jet.getJetConstituents().size(),
                    //(jet.isPFJet()? (int) jet.getPFConstituents().size() : -1),
                    //jet.isPFJet()
                    //(int) jet.pfCandidatesFwdPtr().size(),
                    (int) jet.numberOfDaughters()
                );
            }
            
            std::vector <fastjet::PseudoJet> fj_input_jet;
            
            auto const &v_jet_consti = jet.getJetConstituents();
            
            int iConsti = -1;
            int consti_n = 0;
            double constiE_sum = 0;
            
            for(auto const &consti : v_jet_consti)
            {
                iConsti++;
                
                fastjet::PseudoJet fj_pseudoJet(
                    consti->px(),
                    consti->py(),
                    consti->pz(),
                    consti->energy()
                );
                
                fj_pseudoJet.set_user_index(iConsti);
                
                fj_input_jet.push_back(fj_pseudoJet);
                
                consti_n++;
                constiE_sum += consti->energy();
            }
            
            fastjet::ClusterSequence fj_jet_clustSeq(fj_input_jet, *fj_akJetReclusterDef);
            std::vector <fastjet::PseudoJet> fj_jets = fj_jet_clustSeq.inclusive_jets();
            fastjet::PseudoJet fj_jet = fj_jets.at(0);
            
            
            // For soft-drop
            fastjet::PseudoJet fj_jet_softDrop = fj_jet;
            
            if(jetInfo->apply_sd)
            {
                double sd_zcut = jetInfo->sd_zcut;
                double sd_beta = jetInfo->sd_beta;
                double sd_R0   = jetInfo->sd_R0;
                fastjet::contrib::SoftDrop fj_softDrop(sd_beta, sd_zcut, sd_R0);
                
                fj_jet_softDrop = fj_softDrop(fj_jet_softDrop);
            }
            
            
            fastjet::PseudoJet fj_jet_raw = fj_jet;
            fj_jet = fj_jet_softDrop;
            
            CLHEP::HepLorentzVector jet_raw_4mom = Common::PseudoJetToHepLorentzVector(fj_jet_raw);
            
            jetInfo->v_jet_raw_E_reco.push_back(jet_raw_4mom.e());
            jetInfo->v_jet_raw_px_reco.push_back(jet_raw_4mom.px());
            jetInfo->v_jet_raw_py_reco.push_back(jet_raw_4mom.py());
            jetInfo->v_jet_raw_pz_reco.push_back(jet_raw_4mom.pz());
            jetInfo->v_jet_raw_pT_reco.push_back(jet_raw_4mom.perp());
            jetInfo->v_jet_raw_eta_reco.push_back(jet_raw_4mom.eta());
            jetInfo->v_jet_raw_y_reco.push_back(jet_raw_4mom.rapidity());
            jetInfo->v_jet_raw_phi_reco.push_back(jet_raw_4mom.phi());
            jetInfo->v_jet_raw_m_reco.push_back(jet_raw_4mom.m());
            
            CLHEP::HepLorentzVector jet_4mom = Common::PseudoJetToHepLorentzVector(fj_jet);
            v_jet_4mom.push_back(jet_4mom);
            
            jetInfo->v_jet_E_reco.push_back(jet_4mom.e());
            jetInfo->v_jet_px_reco.push_back(jet_4mom.px());
            jetInfo->v_jet_py_reco.push_back(jet_4mom.py());
            jetInfo->v_jet_pz_reco.push_back(jet_4mom.pz());
            jetInfo->v_jet_pT_reco.push_back(jet_4mom.perp());
            jetInfo->v_jet_eta_reco.push_back(jet_4mom.eta());
            jetInfo->v_jet_y_reco.push_back(jet_4mom.rapidity());
            jetInfo->v_jet_phi_reco.push_back(jet_4mom.phi());
            jetInfo->v_jet_m_reco.push_back(jet_4mom.m());
            
            
            // Secondary vertex stuff
            int iSV = -1;
            std::vector <const reco::VertexCompositePtrCandidate*> v_secVtx_inJet;
            
            for(const auto &sv : v_secondaryVertex)
            {
                iSV++;
                
                math::XYZTLorentzVectorD lv_temp(jet_4mom.px(), jet_4mom.py(), jet_4mom.pz(), jet_4mom.e());
                
                double secVtxDR = ROOT::Math::VectorUtil::DeltaR(v_secVtxDir.at(iSV), lv_temp);
                
                //printf("dR(jet %d, sv %d) %f \n", nJet, iSV+1, secVtxDR);
                
                // A somewhat relaxed threshold
                if(secVtxDR < jet.maxDistance() * 1.5)
                {
                    v_secVtx_inJet.push_back(&sv);
                }
            }
            
            jetInfo->v_jet_nSecVtxInJet_reco.push_back(v_secVtx_inJet.size());
            
            
            fastjet::PseudoJet fj_subStruc = fj_jet;
            
            // N-subjettiness
            double tauNm1 = 1;
            
            for(int iTauN = 0; iTauN <= jetInfo->maxTauN; iTauN++)
            {
                double tauN = 1;
                double tauNratio = 0;
                
                if(iTauN)
                {
                    fastjet::contrib::Nsubjettiness nSubjettiness(
                        iTauN,
                        fastjet::contrib::OnePass_KT_Axes(),
                        fastjet::contrib::UnnormalizedMeasure(1.0)
                    );
                    
                    tauN = nSubjettiness.result(fj_subStruc);
                    
                    if(tauNm1)
                    {
                        tauNratio = tauN / tauNm1;
                    }
                }
                
                jetInfo->vv_jet_tauN_reco.at(iTauN).push_back(tauN);
                jetInfo->vv_jet_tauNratio_reco.at(iTauN).push_back(tauNratio);
                
                //printf(
                //    "Fat jet %d: "
                //    "tau%d (ratio) %0.2f (%0.2f) "
                //    "\n",
                //    
                //    nJet,
                //    iTauN, tauN, tauNratio
                //);
                
                tauNm1 = tauN;
            }
            
            //printf("Starting image stuff... \n");
            
            
            // For jet image formation
            fastjet::PseudoJet fj_image = fj_jet;
            
            // Cluster into exactly jets
            // If there are fewer than 3, returns the constituents
            fastjet::ClusterSequence fj_jetExcSubJet_clustSeq(fj_image.constituents(), *fj_fatJetExcSubJetDef);
            std::vector <fastjet::PseudoJet> fj_jetExcSubJets = fj_jetExcSubJet_clustSeq.exclusive_jets_up_to(3);
            fj_jetExcSubJets = sorted_by_E(fj_jetExcSubJets);
            
            int nConsti = fj_image.constituents().size();
            int nFatJetExcSubJet = fj_jetExcSubJets.size();
            
            jetInfo->v_jet_nConsti_reco.push_back(nConsti);
            
            TVectorD direc1(2);
            TVectorD direc2(2);
            TVectorD finalTranslation(2);
            
            finalTranslation(0) = fj_image.eta() - fj_jetExcSubJets.at(0).eta();
            finalTranslation(1) = fj_jetExcSubJets.at(0).delta_phi_to(fj_image);
            
            std::vector <fastjet::PseudoJet> v_GSaxis;
            
            v_GSaxis.push_back(fj_image);
            
            if(nFatJetExcSubJet >= 2)
            {
                direc1(0) = fj_jetExcSubJets.at(1).eta() - fj_jetExcSubJets.at(0).eta();
                direc1(1) = fj_jetExcSubJets.at(0).delta_phi_to(fj_jetExcSubJets.at(1));
                
                v_GSaxis.push_back(fj_jetExcSubJets.at(0));
            }
            
            else
            {
                v_GSaxis.push_back(fastjet::PseudoJet(0, 0, 0, 0));
            }
            
            
            if(nFatJetExcSubJet >= 3)
            {
                direc2(0) = fj_jetExcSubJets.at(2).eta() - fj_jetExcSubJets.at(0).eta();
                direc2(1) = fj_jetExcSubJets.at(0).delta_phi_to(fj_jetExcSubJets.at(2));
                
                v_GSaxis.push_back(fj_jetExcSubJets.at(1));
            }
            
            else
            {
                v_GSaxis.push_back(fastjet::PseudoJet(0, 0, 0, 0));
            }
            
            
            // Rotation in pt-eta plane
            std::vector <std::vector <double> > v_consti_dEta_dPhi_transformed = Common::getTransformedDetaDphi(
                fj_image.constituents(),
                
                fj_jetExcSubJets.at(0),
                //fj_image,
                
                direc1,
                direc2,
                finalTranslation
            );
            
            
            // GS transform
            std::vector <CLHEP::HepLorentzVector> v_consti_boosted = Common::getGStranformed4mom(
                v_GSaxis,
                fj_image.constituents()
            );
            
            CLHEP::HepLorentzVector jet_4mom_boosted = Common::getHepLorentzVectorSum(v_consti_boosted);
            
            // Rescale
            // NOTE: Do not use jet_4mom_boosted.m() as it can occasionally be zero (due to precision) if jet_4mom.m() is very small
            // Also, very rarely, the mass can be a very small -ve value. Hence use fabs().
            double rescaleFactor = jetInfo->jetRescale_m0 / std::fabs(jet_4mom.m());
            jet_4mom_boosted *= rescaleFactor;
            
            double boostDir = -1;
            
            if(jet_4mom_boosted.e() < jetInfo->jetLorentzBoost_e0)
            {
                boostDir *= -1;
            }
            
            // Boost
            double boostGamma = 1.0/(jetInfo->jetRescale_m0*jetInfo->jetRescale_m0) * (
                jet_4mom_boosted.e() * jetInfo->jetLorentzBoost_e0 -
                jetInfo->jetLorentzBoost_p0 * jet_4mom_boosted.v().mag()
            );
            
            double boostBeta = boostDir * std::sqrt(1.0 - 1.0/(boostGamma*boostGamma));
            jet_4mom_boosted.boostX(boostBeta);
            
            
            // Rescale and boost the constituents
            //#pragma omp parallel for
            for(int iConsti = 0; iConsti < nConsti; iConsti++)
            {
                v_consti_boosted.at(iConsti) *= rescaleFactor;
                v_consti_boosted.at(iConsti).boostX(boostBeta);
            }
            
            CLHEP::HepLorentzVector jet_4mom_boosted_sumConsti = Common::getHepLorentzVectorSum(v_consti_boosted);
            
            std::vector <double> v_jet_consti_E_reco;
            std::vector <double> v_jet_consti_px_reco;
            std::vector <double> v_jet_consti_py_reco;
            std::vector <double> v_jet_consti_pz_reco;
            std::vector <double> v_jet_consti_pT_reco;
            std::vector <double> v_jet_consti_eta_reco;
            std::vector <double> v_jet_consti_phi_reco;
            std::vector <double> v_jet_consti_m_reco;
            
            std::vector <double> v_jet_consti_id_reco;
            
            std::vector <double> v_jet_consti_vx_reco;
            std::vector <double> v_jet_consti_vy_reco;
            std::vector <double> v_jet_consti_vz_reco;
            std::vector <double> v_jet_consti_v2d_reco;
            std::vector <double> v_jet_consti_v3d_reco;
            
            std::vector <double> v_jet_consti_pvdxy_reco;
            std::vector <double> v_jet_consti_pvdz_reco;
            
            std::vector <double> v_jet_consti_svdxy_reco;
            std::vector <double> v_jet_consti_svdz_reco;
            
            std::vector <double> v_jet_consti_PtEtaRot_dEta_reco;
            std::vector <double> v_jet_consti_PtEtaRot_dPhi_reco;
            
            std::vector <double> v_jet_consti_LBGS_x_reco;
            std::vector <double> v_jet_consti_LBGS_y_reco;
            
            std::vector <double> v_jet_consti_enFrac_reco;
            
            std::vector <double> v_jet_consti_pTwrtJet_reco;
            std::vector <double> v_jet_consti_dRwrtJet_reco;
            
            
            std::unordered_map <std::string, std::vector <double> > m_jet_consti_electronInfo_reco;
            
            for(int iVar = 0; iVar < eleMvaVarManager.getNVars(); iVar++)
            {
                std::string varName = eleMvaVarManager.getName(iVar);
                
                std::vector <double> v_temp(nConsti, -99);
                m_jet_consti_electronInfo_reco[varName] = v_temp;
            }
            
            
            for(int iConsti = 0; iConsti < nConsti; iConsti++)
            {
                fastjet::PseudoJet pseudoJet_consti = fj_image.constituents().at(iConsti);
                
                double x_LBGS = v_consti_boosted.at(iConsti).py() / v_consti_boosted.at(iConsti).e();
                double y_LBGS = v_consti_boosted.at(iConsti).pz() / v_consti_boosted.at(iConsti).e();
                double enFrac = v_consti_boosted.at(iConsti).e() / jetInfo->jetLorentzBoost_e0;
                
                int idx = pseudoJet_consti.user_index();
                
                const auto &consti = v_jet_consti.at(idx);
                
                //In the rare cases that a jet as < 3 constituents, the ratio can slightly exceed 1 due to precision, etc.
                enFrac = std::min(enFrac, 1.0);
                
                if(fabs(x_LBGS) > 1 || fabs(y_LBGS) > 1 || fabs(enFrac) > 1)
                {
                    std::cout << x_LBGS << ", " << y_LBGS << ", " << enFrac << "\n";
                }
                
                v_jet_consti_PtEtaRot_dEta_reco.push_back(v_consti_dEta_dPhi_transformed.at(iConsti).at(0));
                v_jet_consti_PtEtaRot_dPhi_reco.push_back(v_consti_dEta_dPhi_transformed.at(iConsti).at(1));
                
                v_jet_consti_LBGS_x_reco.push_back(x_LBGS);
                v_jet_consti_LBGS_y_reco.push_back(y_LBGS);
                
                v_jet_consti_enFrac_reco.push_back(enFrac);
                
                v_jet_consti_E_reco.push_back(consti->energy());
                v_jet_consti_px_reco.push_back(consti->px());
                v_jet_consti_py_reco.push_back(consti->py());
                v_jet_consti_pz_reco.push_back(consti->pz());
                v_jet_consti_pT_reco.push_back(consti->pt());
                v_jet_consti_eta_reco.push_back(consti->eta());
                v_jet_consti_phi_reco.push_back(consti->phi());
                v_jet_consti_m_reco.push_back(consti->mass());
                v_jet_consti_id_reco.push_back(consti->pdgId());
                
                v_jet_consti_vx_reco.push_back(consti->vx());
                v_jet_consti_vy_reco.push_back(consti->vy());
                v_jet_consti_vz_reco.push_back(consti->vz());
                v_jet_consti_v2d_reco.push_back(std::sqrt(consti->vertex().perp2()));
                v_jet_consti_v3d_reco.push_back(std::sqrt(consti->vertex().mag2()));
                
                // wrt PV
                double pvdxy = -1;
                double pvdz = -1;
                
                if(consti->bestTrack() && nGoodVertex)
                {
                    pvdxy = std::fabs(consti->bestTrack()->dxy(pmVtx.position()));
                    pvdz = std::fabs(consti->bestTrack()->dz(pmVtx.position()));
                }
                
                v_jet_consti_pvdxy_reco.push_back(pvdxy);
                v_jet_consti_pvdz_reco.push_back(pvdz);
                
                
                // wrt SV
                double svdxy_min = -1;
                double svdz_min = -1;
                
                if(consti->bestTrack())
                {
                    for(const auto svptr : v_secVtx_inJet)
                    {
                        const auto &sv = *svptr;
                        
                        double svdxy = std::fabs(consti->bestTrack()->dxy(sv.position()));
                        double svdz = std::fabs(consti->bestTrack()->dz(sv.position()));
                        
                        if(svdxy_min < 0 || svdxy < svdxy_min)
                        {
                            svdxy_min = svdxy;
                            svdz_min = svdz;
                        }
                    }
                }
                
                v_jet_consti_svdxy_reco.push_back(svdxy_min);
                v_jet_consti_svdz_reco.push_back(svdz_min);
                
                
                // 2D isolation
               TLorentzVector lv_jet;
               TLorentzVector lv_consti;
                
                lv_jet.SetXYZT(
                    jet.px(),
                    jet.py(),
                    jet.pz(),
                    jet.energy()
                );
                
                lv_consti.SetXYZT(
                    consti->px(),
                    consti->py(),
                    consti->pz(),
                    consti->energy()
                );
                
                double pTwrtJet = lv_consti.Vect().Perp(lv_jet.Vect());
                double dRwrtJet = lv_consti.DeltaR(lv_jet);
                
                v_jet_consti_pTwrtJet_reco.push_back(pTwrtJet);
                v_jet_consti_dRwrtJet_reco.push_back(dRwrtJet);
                
                
                int iEle = -1;
                int nearestEleIdx = -1;
                double minDR = 9999;
                
                for(const auto &ele : electrons_reco)
                {
                    iEle++;
                    double dR = ROOT::Math::VectorUtil::DeltaR(ele.p4(), consti->p4());
                    
                    if(dR < minDR)
                    {
                        nearestEleIdx = iEle;
                        minDR = dR;
                    }
                }
                
                if(nearestEleIdx >= 0 && minDR < 0.05)
                {
                    edm::Ptr <pat::Electron> elePtr(handle_electron_reco, nearestEleIdx);
                    
                    std::vector <float> extraVariables = eleMvaVarHelper.getAuxVariables(elePtr, iEvent);
                    
                    for(int iVar = 0; iVar < eleMvaVarManager.getNVars(); iVar++)
                    {
                        std::string varName = eleMvaVarManager.getName(iVar);
                        
                        m_jet_consti_electronInfo_reco[varName].at(iConsti) = eleMvaVarManager.getValue(iVar, *elePtr, extraVariables);
                    }
                }
            }
            
            jetInfo->vv_jet_consti_E_reco.push_back(v_jet_consti_E_reco);
            jetInfo->vv_jet_consti_px_reco.push_back(v_jet_consti_px_reco);
            jetInfo->vv_jet_consti_py_reco.push_back(v_jet_consti_py_reco);
            jetInfo->vv_jet_consti_pz_reco.push_back(v_jet_consti_pz_reco);
            jetInfo->vv_jet_consti_pT_reco.push_back(v_jet_consti_pT_reco);
            jetInfo->vv_jet_consti_eta_reco.push_back(v_jet_consti_eta_reco);
            jetInfo->vv_jet_consti_phi_reco.push_back(v_jet_consti_phi_reco);
            jetInfo->vv_jet_consti_m_reco.push_back(v_jet_consti_m_reco);
            
            jetInfo->vv_jet_consti_id_reco.push_back(v_jet_consti_id_reco);
            
            jetInfo->vv_jet_consti_vx_reco.push_back(v_jet_consti_vx_reco);
            jetInfo->vv_jet_consti_vy_reco.push_back(v_jet_consti_vy_reco);
            jetInfo->vv_jet_consti_vz_reco.push_back(v_jet_consti_vz_reco);
            jetInfo->vv_jet_consti_v2d_reco.push_back(v_jet_consti_v2d_reco);
            jetInfo->vv_jet_consti_v3d_reco.push_back(v_jet_consti_v3d_reco);
            
            jetInfo->vv_jet_consti_pvdxy_reco.push_back(v_jet_consti_pvdxy_reco);
            jetInfo->vv_jet_consti_pvdz_reco.push_back(v_jet_consti_pvdz_reco);
            
            jetInfo->vv_jet_consti_svdxy_reco.push_back(v_jet_consti_svdxy_reco);
            jetInfo->vv_jet_consti_svdz_reco.push_back(v_jet_consti_svdz_reco);
            
            jetInfo->vv_jet_consti_PtEtaRot_dEta_reco.push_back(v_jet_consti_PtEtaRot_dEta_reco);
            jetInfo->vv_jet_consti_PtEtaRot_dPhi_reco.push_back(v_jet_consti_PtEtaRot_dPhi_reco);
            
            jetInfo->vv_jet_consti_LBGS_x_reco.push_back(v_jet_consti_LBGS_x_reco);
            jetInfo->vv_jet_consti_LBGS_y_reco.push_back(v_jet_consti_LBGS_y_reco);
            
            jetInfo->vv_jet_consti_enFrac_reco.push_back(v_jet_consti_enFrac_reco);
            
            jetInfo->vv_jet_consti_pTwrtJet_reco.push_back(v_jet_consti_pTwrtJet_reco);
            jetInfo->vv_jet_consti_dRwrtJet_reco.push_back(v_jet_consti_dRwrtJet_reco);
            
            
            for(int iVar = 0; iVar < eleMvaVarManager.getNVars(); iVar++)
            {
                std::string varName = eleMvaVarManager.getName(iVar);
                
                jetInfo->m_jet_consti_electronInfo_reco[varName].push_back(m_jet_consti_electronInfo_reco[varName]);
            }
        }
        
        
        jetInfo->jet_n_reco = jetInfo->v_jet_E_reco.size();
        
        
        // Matching to gen top
        TMatrixD mat_genTopMatch;
        std::vector <int> v_matchedGenTop_idx;
        
        std::vector <double> v_genTop_deltaR = Common::getMinDeltaR(
            v_jet_4mom,
            v_genTopVis_4mom,
            mat_genTopMatch,
            v_matchedGenTop_idx,
            true
        );
        
        for(int iJet = 0; iJet < jetInfo->jet_n_reco; iJet++)
        {
            int genTop_idx = v_matchedGenTop_idx.at(iJet);
            double genTop_minDR = v_genTop_deltaR.at(iJet);
            int genTop_isLeptonic = (genTop_idx >= 0)? (int) treeOutput->v_genTop_isLeptonic.at(genTop_idx) : 0;
            
            jetInfo->v_jet_nearestGenTopIdx_reco.push_back(genTop_idx);
            jetInfo->v_jet_nearestGenTopDR_reco.push_back(genTop_minDR);
            jetInfo->v_jet_nearestGenTopIsLeptonic_reco.push_back(genTop_isLeptonic);
        }
        
        
        // Matching to gen W
        TMatrixD mat_genWMatch;
        std::vector <int> v_matchedGenW_idx;
        
        std::vector <double> v_genW_deltaR = Common::getMinDeltaR(
            v_jet_4mom,
            v_genWvis_4mom,
            mat_genWMatch,
            v_matchedGenW_idx,
            true
        );
        
        for(int iJet = 0; iJet < jetInfo->jet_n_reco; iJet++)
        {
            int genW_idx = v_matchedGenW_idx.at(iJet);
            double genW_minDR = v_genW_deltaR.at(iJet);
            int genW_isLeptonic = (genW_idx >= 0)? (int) treeOutput->v_genW_isLeptonic.at(genW_idx) : 0;
            
            jetInfo->v_jet_nearestGenWIdx_reco.push_back(genW_idx);
            jetInfo->v_jet_nearestGenWDR_reco.push_back(genW_minDR);
            jetInfo->v_jet_nearestGenWIsLeptonic_reco.push_back(genW_isLeptonic);
        }
        
        
        // Matching to gen Z
        TMatrixD mat_genZMatch;
        std::vector <int> v_matchedGenZ_idx;
        
        std::vector <double> v_genZ_deltaR = Common::getMinDeltaR(
            v_jet_4mom,
            v_genZ_4mom,
            mat_genZMatch,
            v_matchedGenZ_idx,
            true
        );
        
        for(int iJet = 0; iJet < jetInfo->jet_n_reco; iJet++)
        {
            int genZ_idx = v_matchedGenZ_idx.at(iJet);
            double genZ_minDR = v_genZ_deltaR.at(iJet);
            int genZ_isLeptonic = (genZ_idx >= 0)? (int) treeOutput->v_genZ_isLeptonic.at(genZ_idx) : 0;
            
            jetInfo->v_jet_nearestGenZIdx_reco.push_back(genZ_idx);
            jetInfo->v_jet_nearestGenZDR_reco.push_back(genZ_minDR);
            jetInfo->v_jet_nearestGenZIsLeptonic_reco.push_back(genZ_isLeptonic);
        }
    }
    
    
    // Fill tree
    treeOutput->fill();
    
    //#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    //ESHandle<SetupData> pSetup;
    //iSetup.get<SetupRecord>().get(pSetup);
    //#endif
    
    if(debug)
    {
        printf("\n\n");
    }
    
    fflush(stdout);
    fflush(stderr);
}


// ------------ method called once each job just before starting event loop  ------------
void
TreeMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TreeMaker::endJob()
{
    fflush(stdout);
    fflush(stderr);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeMaker);
