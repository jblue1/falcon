/**
 * Script used to test the amount of time taken to generate ttbar events with
 * hadronization turned on.
 */
#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"

using namespace Pythia8;
using namespace fastjet;

void write_momenta(std::ofstream &stream, std::vector<float> &parton_jet) {

  for (int i = 0; i < parton_jet.size(); i++) {
    int index = i % 4;

    stream << parton_jet[i] << " ";

    if (index == 3) {
      stream << std::endl;
    } else {
      stream << " ";
    }
  }
}

int main(int argc, char const *argv[]) {
  if (argc != 2) {
    std::cout << "Incorrect number of arguments" << std::endl;
    return 1;
  }
  std::ofstream write_out("test-noML.txt");
  // std::ofstream energy_out("total_energy-cppflow.txt");
  assert(write_out.is_open());
  // assert(energy_out.is_open());

  // LHEF file initialization

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;

  // LHEF file initialization

  LHEF3FromPythia8 myLHEF3(&pythia.event, &pythia.settings, &pythia.info,
                           &pythia.particleData);
  myLHEF3.openLHEF("test.lhe");

  myLHEF3.setInit();

  int numEvents = atoi(argv[1]);
  // process settings
  pythia.readString("Top:gg2ttbar = on");
  pythia.readString("Top:qqbar2ttbar = on");
  pythia.readString("6:m0 = 175");
  pythia.readString("Beams:eCM = 13000.");

  // Tune Settings
  pythia.readString("Tune:pp 14");
  pythia.readString("Tune:ee 7");
  pythia.readString("MultipartonInteractions:expPow=1.6");
  pythia.readString("MultipartonInteractions:ecmPow=0.25208'");
  pythia.readString("MultipartonInteractions:pT0Ref = 2.4024");

  // Common Settings
  pythia.readString("Tune:preferLHAPDF = 2");
  pythia.readString("Next:numberShowEvent = 2");
  pythia.readString("Check:epTolErr = 0.010");
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10");
  pythia.readString("ParticleDecays:allowPhotonRadiation = on");
  pythia.readString("SLHA:minMassSM = 1000.");
  // pythia.readString("SLHA:keepSM = on");
  pythia.readString("Beams:setProductionScalesFromLHEF = off");

  // other settings

  // Turn off hadronization
  // pythia.readString("HadronLevel:all = off");
  pythia.init();

  int num_jets = 0;

  std::vector<float> parton_jet_momenta;

  // Loop over events.
  for (int i = 0; i < numEvents; ++i) {

    // Generate an event.
    pythia.next();

    myLHEF3.setEvent();

    int numParticles = pythia.event.size();

    std::vector<PseudoJet> particles;
    for (int j = 0; j < numParticles; j++) {
      Particle &p = pythia.event[j];
      int status = p.status();
      if (status > 0) {
        particles.push_back(PseudoJet(p.px(), p.py(), p.pz(), p.e()));
      }
    }
    double R = 0.4;
    JetDefinition jet_def(antikt_algorithm, R);
    ClusterSequence cs(particles, jet_def);

    std::vector<PseudoJet> jets = cs.inclusive_jets();

    for (int j = 0; j < jets.size(); j++) {
      float partonJetPt = jets[j].pt();
      float partonJetEta = jets[j].rap();
      float partonJetPhi = jets[j].phi_std();
      float partonJetE = jets[j].E();
      if (partonJetPt > 20) {
        parton_jet_momenta.push_back(partonJetPt);
        parton_jet_momenta.push_back(partonJetEta);
        parton_jet_momenta.push_back(partonJetPhi);
        parton_jet_momenta.push_back(partonJetE);
        num_jets++;
      }
    }
  }

  write_momenta(write_out, parton_jet_momenta);

  // Statistics: full printout.
  pythia.stat();

  myLHEF3.closeLHEF(true);

  write_out.close();
  return 0;
}
