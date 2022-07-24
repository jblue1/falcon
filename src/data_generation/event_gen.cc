/*
 * In order for this to run, you need to set
 * ROOT_INCLUDE_PATH=~/Programs/Delphes/Delphes-3.5.0/external/ Solution from
 * https://root-forum.cern.ch/t/error-in-cling-insertintoautoloadingstate/29347.
 */

// normal includes
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
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "modules/Delphes.h"

#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

using namespace Pythia8;

/**
 * @brief Fill attributes of HepMCEvent with info from pythia event
 */
void assignElementInfo(HepMCEvent *element, Pythia *pythia) {
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
                         TDatabasePDG *pdgData, TParticlePDG *pdgParticle) {
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
void ConvertInput(int eventNumber, Pythia *pythia, DelphesFactory *factory,
                  TObjArray *allParticleOutputArray,
                  TObjArray *stableParticleOutputArray,
                  TObjArray *partonOutputArray, ExRootTreeBranch *branch) {

  HepMCEvent *element = static_cast<HepMCEvent *>(branch->NewEntry());
  Candidate *candidate;
  TDatabasePDG *pdgData = TDatabasePDG::Instance();
  TParticlePDG *pdgParticle = NULL;
  element->Number = eventNumber;
  std::cout << "Event number = " << element->Number << std::endl;

  assignElementInfo(element, pythia);

  std::cout << "There are " << pythia->event.size() << " particles in event "
            << eventNumber << std::endl;
  for (int i = 0; i < pythia->event.size(); i++) {
    Particle &particle = pythia->event[i];
    candidate = factory->NewCandidate();
    assignCandidateInfo(candidate, particle, pdgData, pdgParticle);
    int pdgCode = TMath::Abs(candidate->PID);

    allParticleOutputArray->Add(candidate);
    if (!pdgParticle) continue;

    if (candidate->Status == 1) {
        stableParticleOutputArray->Add(candidate);
    }
    else if (pdgCode <= 5 || pdgCode == 21 || pdgCode == 15) {
        partonOutputArray->Add(candidate);
    }
  }
  return;
}

/**
 * @brief Create and initialize a pythia object
 */
Pythia *create_pythia_object() {
  Pythia *pythia = new Pythia;
  // process settings
  pythia->readString("Top:gg2ttbar = on");
  pythia->readString("Top:qqbar2ttbar = on");
  pythia->readString("Beams:eCM = 13000.");
  pythia->readString("Tune:pp 14");
  pythia->readString("Tune:ee 7");
  // pythia->init();

  return pythia;
}

int main(int argc, char const *argv[]) {
  if (argc != 3) {
    std::cerr << "Arguments are num events and path to delphes card"
              << std::endl;
    exit(1);
  }

  TFile *outputFile = TFile::Open("test.root", "CREATE");
  // Create some Delphes Objects
  Delphes *modularDelphes = new Delphes("Delphes");
  DelphesFactory *factory = NULL;
  ExRootConfReader *confReader = new ExRootConfReader;
  ExRootTreeWriter *treeWriter = new ExRootTreeWriter(outputFile, "Delphes");
  ExRootTreeBranch *branchEvent = treeWriter->NewBranch("Event", HepMCEvent::Class());
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
  Pythia *pythia = create_pythia_object();
  pythia->init();

  // gROOT->SetBatch();

  int numEvents = atoi(argv[1]);
  for (int i = 0; i < numEvents; i++) {
    std::cout << "Event " << i << std::endl;
    pythia->next();
    ConvertInput(i, pythia, factory, allParticleOutputArray, stableParticleOutputArray, partonOutputArray, branchEvent);
    modularDelphes->ProcessTask();

    treeWriter->Fill();
    treeWriter->Clear();
    modularDelphes->Clear();
  }

  pythia->stat();

  modularDelphes->FinishTask();
  treeWriter->Write();

  // delete pointers
  std::cout << "About to delete the pointers" << std::endl;
  delete pythia;
  delete modularDelphes;
  delete confReader;
  delete treeWriter;
  exit(0);
}
