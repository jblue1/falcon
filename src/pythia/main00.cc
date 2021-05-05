#include "Pythia8/Pythia.h"
#include "cppflow/cppflow.h"
#include "fastjet/ClusterSequence.hh"
#include <chrono>
#include <cmath>
#include <random>
#include <cassert>

using namespace Pythia8;
using namespace fastjet;

// constants used for normalizing and un-normalizing data
// To see where these come from look at notebook evaluatecWGAN.ipynb

std::vector<float> PARTON_MEANS = {1.7325336127869422, -4.3696245126315283e-03,
                                   -8.1989632208202834e-04, 2.0738273199892845};
std::vector<float> PARTON_STD_DEVS = {0.2664415599646997, 1.6049471051941415,
                                      1.813214019844585, 0.4007043722622637};
std::vector<float> RECO_MEANS = {1.6574159582302348, -4.4009885497273684e-03,
                                 -5.8530799456628130e-04, 1.9991621188618929};
std::vector<float> RECO_STD_DEVS = {0.3107208060660628, 1.609712071592146,
                                    1.8132994467877084, 0.4030490622744889};

/**
 * Normalize a set of four-momenta so each component has
 * a mean of zero a variance of one
 */
std::vector<float> normalize(std::vector<float> four_vec,
                             std::vector<float> mean,
                             std::vector<float> std_dev) {
  for (int i = 0; i < four_vec.size(); i++) {
    int index = i % 4;

    if (index == 0 || index == 3) {
      // take base 10 log of pT and E components
      four_vec[i] = log10(four_vec[i]);
    }
    // subtract mean and divide by std dev for each component
    four_vec[i] = (four_vec[i] - mean[index]) / std_dev[index];
  }

  return four_vec;
}

/**
 * Take data set of four-momenta with mean of zero and variance 
 * of one and scale them back to original distribution
 */
std::vector<float> unnormalize(std::vector<float> four_vec,
                               std::vector<float> mean,
                               std::vector<float> std_dev) {
  for (int i = 0; i < four_vec.size(); i++) {
    int index = i % 4;
    // multiply by std dev and add mean for each component
    four_vec[i] = four_vec[i] * std_dev[index] + mean[index];

    if (index == 0 || index == 3) {
      // take base 10 log of pT and E components
      four_vec[i] = pow(10, four_vec[i]);
    }
  }

  return four_vec;
}

/**
 * Write matched parton and reco jet four momenta to a text file
 */
void write_momenta(std::ofstream &stream, std::vector<int> &events, std::vector<float> &parton_jet,
                   std::vector<float> &reco_jet) {

  int jet_number = 0;
  for (int i = 0; i < reco_jet.size(); i++) {
    int index = i % 4;
    

    if (index == 0) {
      stream << events[jet_number] << " ";
    }

    
    stream << parton_jet[i] << " " << reco_jet[i];

    if (index == 3) {
      stream << std::endl;
      jet_number++;
    } else {
      stream << " ";
      
    }
  }
}

// Set up random number generator
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);

/**
 * Sample a vector of a given length from a uniform distribution
 */
std::vector<float> sample_rand_vec(int size) {
  std::vector<float> output;

  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  for (int i = 0; i < size; i++) {
    double number = distribution(generator);
    output.push_back(number);
  }
  return output;
}

int main(int argc, char const *argv[]) {
  if (argc != 2) {
    std::cout << "Incorrect number of arguments" << std::endl;
    return 1;
  }
  std::ofstream write_out("test-cppflow.txt");
  // std::ofstream energy_out("total_energy-cppflow.txt");
  assert(write_out.is_open());
  // assert(energy_out.is_open());

  // load model
  cppflow::model model(
      "/home/DAVIDSON/joblue/phy2/falcon/models/cWGAN/Run_2021-03-09_0/model");
  auto ops = model.get_operations();
  for (int kk = 0; kk < ops.size(); kk++) {
    std::cout << ops[kk] << std::endl;
  }

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;

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

  // Turn off hadronization
  pythia.readString("HadronLevel:all = off");
  pythia.init();

  int num_jets = 0;

  std::vector<float> parton_jet_momenta;

  std::vector<int> events;

  // Loop over events.
  for (int i = 0; i < numEvents; ++i) {

    // Generate an event.
    pythia.next();

    int numParticles = pythia.event.size();

    float totalE = 0.0;
    std::vector<PseudoJet> particles;
    for (int j = 0; j < numParticles; j++) {
      Particle &p = pythia.event[j];
      int status = p.status();
      if (status > 0) {
        int id = p.id();
        if (
            // quarks
            id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id == 6 ||
            id == 7 || id == 8 || id == -1 || id == -2 || id == -3 ||
            id == -4 || id == -5 || id == -6 || id == -7 || id == -8 ||
            // gluons
            id == 21 ||
            // diquarks
            id == 1103 || id == 2101 || id == 2103 || id == 2203 ||
            id == 3101 || id == 3103 || id == 3201 || id == 3203 ||
            id == 3303 || id == 4101 || id == 4103 || id == 4201 ||
            id == 4203 || id == 4301 || id == 4303 || id == 4403 ||
            id == 5101 || id == 5103 || id == 5201 || id == 5203 ||
            id == 5301 || id == 5303 || id == 5401 || id == 5403 ||
            id == 5503) {

          totalE += p.e();
          particles.push_back(PseudoJet(p.px(), p.py(), p.pz(), p.e()));
        } else if (id > 22 || id < -22) {
          std::cout << "Did not include particle: " << id << " with Pt "
                    << p.pT() << std::endl;
        }
      }
    }
    // energy_out << totalE << std::endl;
    // std::cout << "Total Energy: " << totalE << std::endl;
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
        events.push_back(i);
        num_jets++;
      }
    }
  }
  const std::vector<float> rand = sample_rand_vec(10 * num_jets);
  const std::vector<int64_t> rand_shape = {num_jets, 10};
  const std::vector<int64_t> jets_shape = {num_jets, 4};

  std::vector<float> normalized_parton_jet_momenta =
      normalize(parton_jet_momenta, PARTON_MEANS, PARTON_STD_DEVS);

  cppflow::tensor noise = cppflow::tensor(rand, rand_shape);
  cppflow::tensor jets_tensor =
      cppflow::tensor(normalized_parton_jet_momenta, jets_shape);
  auto output = model({{"serving_default_noiseIn", noise},
                       {"serving_default_pjetIn", jets_tensor}},
                      {"StatefulPartitionedCall:0"});
  cppflow::tensor *data = output.data();
  std::vector<float> data_vec = data->get_data<float>();
  std::vector<float> reco_jet_momenta =
      unnormalize(data_vec, RECO_MEANS, RECO_STD_DEVS);

  write_momenta(write_out, events, parton_jet_momenta, reco_jet_momenta);

  // Statistics: full printout.
  pythia.stat();

  write_out.close();
  // energy_out.close();
  return 0;
}
