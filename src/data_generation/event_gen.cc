/*
 * In order for this to run, you need to set
 * ROOT_INCLUDE_PATH=~/Programs/Delphes/Delphes-3.5.0/external/ Solution from
 * https://root-forum.cern.ch/t/error-in-cling-insertintoautoloadingstate/29347.
 */

// normal includes
#include <cassert>
#include <cmath>
#include <iostream>
// pythia includes
#include "Pythia8/Pythia.h"
// ROOT includes
#include "TApplication.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TParticlePDG.h"
#include "TROOT.h"
// Delphes includes
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "modules/Delphes.h"

// FastJet include
#include "fastjet/ClusterSequence.hh"

using namespace Pythia8;
using namespace fastjet;

/**
 * @brief Write out parton jet 4-momenta as a text file with spaces as the
 * delimiters, and each jet on its own line
 *
 * @param partonJet list of floats in the order pt, eta, phi, m^2
 */
void write_momenta(std::ofstream &stream, std::vector<float> &partonJet)
{
    for (int i = 0; i < partonJet.size(); i++)
    {
        int index = i % 4;

        stream << partonJet[i] << " ";
        if (index == 3)
        {
            stream << std::endl;
        }
        else
        {
            stream << " ";
        }
    }
}

/**
 * @brief Convert a vector of jets (of size n) defined by pt, eta, phi, and e,
 * to a length 4*n vector of pt, eta, phi, and m^2. Only include jets with pt >
 * 20 GeV and a positive m^2.
 */
std::vector<float> convertKinematicVariables(std::vector<PseudoJet> &jets)
{
    std::vector<float> partonJetMomenta;
    for (int i = 0; i < jets.size(); i++)
    {
        float partonJetPt = jets[i].pt();
        float partonJetEta = jets[i].rap();
        float partonJetPhi = jets[i].phi_std();
        float partonJetPSquared =
            pow(jets[i].px(), 2) + pow(jets[i].py(), 2) + pow(jets[i].pz(), 2);
        float partonJetMSquared = pow(jets[i].e(), 2) - partonJetPSquared;
        if (partonJetPt > 20 && partonJetMSquared > 0)
        {
            partonJetMomenta.push_back(partonJetPt);
            partonJetMomenta.push_back(partonJetEta);
            partonJetMomenta.push_back(partonJetPhi);
            partonJetMomenta.push_back(partonJetMSquared);
        }
    }
    return partonJetMomenta;
}

/**
 * @brief Fill attributes of HepMCEvent with info from pythia event
 */
void assignElementInfo(HepMCEvent *element, Pythia *pythia)
{
    element->ProcessID = pythia->info.code();
    element->MPI = 1;
    element->Weight = pythia->info.weight();
    element->Scale = pythia->info.QRen();
    element->AlphaQED = pythia->info.alphaEM();
    element->AlphaQCD = pythia->info.alphaS();
    element->ID1 = pythia->info.id1();
    element->ID2 = pythia->info.id2();
    element->X1 = pythia->info.x1();
    element->X2 = pythia->info.x2();
    element->ScalePDF = pythia->info.QFac();
    element->PDF1 = pythia->info.pdf1();
    element->PDF2 = pythia->info.pdf2();
    return;
}

/**
 * @brief Fill delphes candidate attributes with info from pythia particle
 */
void assignCandidateInfo(Candidate *candidate, Particle &particle,
                         TDatabasePDG *pdgData, TParticlePDG *pdgParticle)
{
    candidate->PID = particle.id();
    candidate->Status = particle.statusHepMC();
    candidate->M1 = particle.mother1() - 1;
    candidate->M2 = particle.mother2() - 1;
    candidate->D1 = particle.daughter1() - 1;
    candidate->D2 = particle.daughter2() - 1;
    pdgParticle = pdgData->GetParticle(particle.id());
    candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge() / 3.0) : -999;
    candidate->Mass = particle.m();
    candidate->Momentum.SetPxPyPzE(particle.px(), particle.py(), particle.pz(),
                                   particle.e());
    candidate->Position.SetXYZT(particle.xProd(), particle.yProd(),
                                particle.zProd(), particle.tProd());
    return;
}

/**
 * @brief Convert every particle in the event into a Delphes candidate
 */
void convertToDelphesCandidate(int eventNumber, Pythia *pythia,
                               DelphesFactory *factory,
                               TObjArray *allParticleOutputArray,
                               TObjArray *stableParticleOutputArray,
                               TObjArray *partonOutputArray,
                               ExRootTreeBranch *branch)
{
    HepMCEvent *element = static_cast<HepMCEvent *>(branch->NewEntry());
    Candidate *candidate;
    TDatabasePDG *pdgData = TDatabasePDG::Instance();
    TParticlePDG *pdgParticle = NULL;
    element->Number = eventNumber;
    std::cout << "Event number = " << element->Number << std::endl;

    assignElementInfo(element, pythia);

    std::cout << "There are " << pythia->event.size() << " particles in event "
              << eventNumber << std::endl;
    for (int i = 0; i < pythia->event.size(); i++)
    {
        Particle &particle = pythia->event[i];
        candidate = factory->NewCandidate();
        assignCandidateInfo(candidate, particle, pdgData, pdgParticle);
        int pdgCode = TMath::Abs(candidate->PID);

        allParticleOutputArray->Add(candidate);
        if (!pdgParticle)
            continue;

        if (candidate->Status == 1)
        {
            stableParticleOutputArray->Add(candidate);
        }
        else if (pdgCode <= 5 || pdgCode == 21 || pdgCode == 15)
        {
            partonOutputArray->Add(candidate);
        }
    }
    return;
}

/**
 * @brief Create and initialize a pythia object
 */
Pythia *create_pythia_object()
{
    Pythia *pythia = new Pythia;
    // process settings
    pythia->readString("Top:gg2ttbar = on");
    pythia->readString("Top:qqbar2ttbar = on");
    pythia->readString("Beams:eCM = 13000.");
    pythia->readString("Tune:pp 14");
    pythia->readString("Tune:ee 7");
    pythia->readString("Random:seed = 0");
    // pythia->init();

    return pythia;
}

/**
 * @brief Determine if a parton should be used for jet clustering
 */
bool isValidParton(int status, int id)
{
    if (status > 0)
    {
        if (
            // quarks
            abs(id) == 1 || abs(id) == 2 || abs(id) == 3 || abs(id) == 4 ||
            abs(id) == 5 || abs(id) == 6 || abs(id) == 7 || abs(id) == 8 ||
            // gluons
            id == 21 ||
            // diquarks
            id == 1103 || id == 2101 || id == 2103 || id == 2203 ||
            id == 3101 || id == 3103 || id == 3201 || id == 3203 ||
            id == 3303 || id == 4101 || id == 4103 || id == 4201 ||
            id == 4203 || id == 4301 || id == 4303 || id == 4403 ||
            id == 5101 || id == 5103 || id == 5201 || id == 5203 ||
            id == 5301 || id == 5303 || id == 5401 || id == 5403 || id == 5503)
        {
            return true;
        }
    }
    return false;
}

/**
 * @brief Process the hard scatter + parton shower part of an event. This
 * involves two things: 1) Appending all of the particles from the hard process
 * pythia object to the hadronization pythia object. 2) Clustering jets out of
 * the partons
 *
 * @return std::vector<PseudoJet>
 */
std::vector<PseudoJet> processHardScatter(Pythia *hard_process,
                                          Pythia *hadronization)
{
    hadronization->event.reset();
    std::vector<PseudoJet> particles;
    for (int i = 0; i < hard_process->event.size(); i++)
    {
        Particle &particle = hard_process->event[i];
        hadronization->event.append(particle);
        int status = particle.status();
        int id = particle.id();

        if (isValidParton(status, id))
        {
            particles.push_back(PseudoJet(particle.px(), particle.py(),
                                          particle.pz(), particle.e()));
        }
        else
        {
            if (status > 0)
            {
                std::cerr << "Not including particle " << id
                          << " with pt = " << particle.pT() << " and status "
                          << status << std::endl;
            }
        }
    }
    double R = 0.5; // value used in delphes cms card
    JetDefinition jet_def(antikt_algorithm, R);
    ClusterSequence cs(particles, jet_def);

    return cs.inclusive_jets();
}

int main(int argc, char const *argv[])
{
    if (argc != 3)
    {
        std::cerr << "Arguments are num events and path to delphes card"
                  << std::endl;
        exit(1);
    }

    // create output files
    std::ofstream partonJetOutFile("test.txt");
    TFile *delphesJetOutFile = TFile::Open("test.root", "CREATE");
    assert(partonJetOutFile);
    partonJetOutFile << "pt eta phi m**2" << std::endl;

    // Create some Delphes Objects
    Delphes *modularDelphes = new Delphes("Delphes");
    DelphesFactory *factory = NULL;
    ExRootConfReader *confReader = new ExRootConfReader;
    ExRootTreeWriter *treeWriter =
        new ExRootTreeWriter(delphesJetOutFile, "Delphes");
    ExRootTreeBranch *branchEvent =
        treeWriter->NewBranch("Event", HepMCEvent::Class());
    std::cout << "Trying to open configuration file " << argv[2] << std::endl;
    confReader->ReadFile(argv[2]);
    modularDelphes->SetConfReader(confReader);
    modularDelphes->SetTreeWriter(treeWriter);
    factory = modularDelphes->GetFactory();

    // Create some ROOT objects
    TObjArray *allParticleOutputArray =
        modularDelphes->ExportArray("allParticles");
    TObjArray *stableParticleOutputArray =
        modularDelphes->ExportArray("stableParticles");
    TObjArray *partonOutputArray = modularDelphes->ExportArray("partons");

    modularDelphes->InitTask();

    // create pythia object
    // we use two pythia objects: 1 to simulate the hard scatter and parton
    // shower, and one for hadronization. This way, we can easily cluster parton
    // level jets and save the 4-momenta for training.
    // See https://pythia.org/manuals/pythia8307/HadronLevelStandalone.html for
    // more info
    Pythia *pythia_hard_process = create_pythia_object();
    pythia_hard_process->readString("HadronLevel:all = off");
    pythia_hard_process->init();
    Pythia *pythia_hadronization = create_pythia_object();
    pythia_hadronization->readString("ProcessLevel:all = off");
    pythia_hadronization->init();

    // gROOT->SetBatch();

    int numEvents = atoi(argv[1]);
    for (int i = 0; i < numEvents; i++)
    {
        std::cout << "Event " << i << std::endl;
        pythia_hard_process->next();
        std::vector<PseudoJet> jets =
            processHardScatter(pythia_hard_process, pythia_hadronization);

        std::vector<float> convertedJets = convertKinematicVariables(jets);
        write_momenta(partonJetOutFile, convertedJets);
        pythia_hadronization->next();
        convertToDelphesCandidate(
            i, pythia_hadronization, factory, allParticleOutputArray,
            stableParticleOutputArray, partonOutputArray, branchEvent);
        modularDelphes->ProcessTask();

        treeWriter->Fill();
        treeWriter->Clear();
        modularDelphes->Clear();
    }

    pythia_hard_process->stat();

    modularDelphes->FinishTask();
    treeWriter->Write();

    partonJetOutFile.close();

    // delete pointers
    delete pythia_hard_process;
    delete pythia_hadronization;
    delete modularDelphes;
    delete confReader;
    delete treeWriter;
    exit(0);
}
