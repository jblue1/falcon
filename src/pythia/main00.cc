/**
 * This code is used to create an executable generate reco level jet 4-vectors using pythia, 
 * cppflow, and a saved version of a trained cWGAN. timeTestModel.py is an example of a 
 * python scrip that uses the executable. The executable takes in 
 *
 * 
 * NOTE: To compile this script, you need to not have the falcon conda env activated, and 
 * you need some additional packages. First, you need to install 
 * the tensorflow 2 C api: https://www.tensorflow.org/install/lang_c. Next, you need the cppflow
 * library: https://github.com/serizba/cppflow. Finally, you need to have a version of pythia
 * installed: https://pythia.org/. 
 * ADDITIONALLY: You need to set the following environmental variables: 
 * export LIBRARY_PATH=$LIBRARY_PATH:/path/to/where/you/installed/libtensorflow2/lib
 * export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/where/you/installed/libtensorflow2/lib
 * export PYTHIAROOT=/path/to/where/you/installed/pythia
 * FINALLY: The makefile in this directory has TENSORFLOWFLAGS and CPPFLOWFLAGS. You need 
 * to set these to where you have installed libtensorflow2 and cppflow. 
 */
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

std::vector<float> normalize(std::vector<float> four_vec,
                             std::vector<float> mean,
                             std::vector<float> std_dev)
{
  for (int i = 0; i < four_vec.size(); i++)
  {
    int index = i % 4;

    if (index == 0 || index == 3)
    {
      // take base 10 log of pT and E components
      four_vec[i] = log10(four_vec[i]);
    }
    // subtract mean and divide by std dev for each component
    four_vec[i] = (four_vec[i] - mean[index]) / std_dev[index];
  }

  return four_vec;
}

std::vector<float> unnormalize(std::vector<float> four_vec,
                               std::vector<float> mean,
                               std::vector<float> std_dev)
{
  for (int i = 0; i < four_vec.size(); i++)
  {
    int index = i % 4;
    // multiply by std dev and add mean for each component
    four_vec[i] = four_vec[i] * std_dev[index] + mean[index];

    if (index == 0 || index == 3)
    {
      // take base 10 log of pT and E components
      four_vec[i] = pow(10, four_vec[i]);
    }
  }

  return four_vec;
}

void write_momenta(std::ofstream &stream, std::vector<int> &events, std::vector<float> &parton_jet,
                   std::vector<float> &reco_jet)
{


  int jet_number = 0;
  for (int i = 0; i < reco_jet.size(); i++)
  {
    int index = i % 4;

    if (index == 0)
    {
      stream << events[jet_number] << " ";
    }

    stream << parton_jet[i] << " " << reco_jet[i];

    if (index == 3)
    {
      stream << std::endl;
      jet_number++;
    }
    else
    {
      stream << " ";
    }
  }
}

// Set up random number generator
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);

// Create a 1-D vector of random floats sampled from UNI[0,1]
std::vector<float> sample_rand_vec(int size)
{
  std::vector<float> output;

  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  for (int i = 0; i < size; i++)
  {
    double number = distribution(generator);
    output.push_back(number);
  }
  return output;
}

int main(int argc, char const *argv[])
{
  if (argc != 18)
  {
    std::cout << "Incorrect number of arguments" << std::endl;
    return 1;
  }
  auto start = std::chrono::high_resolution_clock::now();
  // load the means and std_devs for normalization.
  std::vector<float> parton_means;
  std::vector<float> parton_std_devs;
  std::vector<float> reco_means;
  std::vector<float> reco_std_devs;
  for (int i=0; i<4; i++) {
    parton_means.push_back(atof(argv[i+2])); //argv[1] is the number of events
  }
  for (int i=0; i<4; i++) {
    parton_std_devs.push_back(atof(argv[i+6]));
  }
  for (int i=0; i<4; i++) {
    reco_means.push_back(atof(argv[i+10]));
  }
  for (int i=0; i<4; i++) {
    reco_std_devs.push_back(atof(argv[i+14]));
  }
  auto end_load_constants = std::chrono::high_resolution_clock::now();

  std::ofstream write_out("test-cppflow.txt");
  // std::ofstream energy_out("total_energy-cppflow.txt");
  assert(write_out.is_open());
  // assert(energy_out.is_open());

  auto start_load_model = std::chrono::high_resolution_clock::now();
  // load model
  cppflow::model model(
      "/home/DAVIDSON/joblue/phy2/falcon/models/cWGAN/Run_2021-07-23_2/model");
  auto end_load_model = std::chrono::high_resolution_clock::now();
  //auto ops = model.get_operations();
  //for (int kk = 0; kk < ops.size(); kk++)
  //{
  //  std::cout << ops[kk] << std::endl;
  //}

  // Generator. Process selection. LHC initialization. Histogram.
  auto start_generate_events = std::chrono::high_resolution_clock::now();
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
  for (int i = 0; i < numEvents; ++i)
  {

    // Generate an event.
    pythia.next();

    int numParticles = pythia.event.size();

    float totalE = 0.0;
    std::vector<PseudoJet> particles;
    for (int j = 0; j < numParticles; j++)
    {
      Particle &p = pythia.event[j];
      int status = p.status();
      if (status > 0)
      {
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
            id == 5503)
        {

          totalE += p.e();
          particles.push_back(PseudoJet(p.px(), p.py(), p.pz(), p.e()));
        }
        else if (id > 22 || id < -22)
        {
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

    for (int j = 0; j < jets.size(); j++)
    {
      float partonJetPt = jets[j].pt();
      float partonJetEta = jets[j].rap();
      float partonJetPhi = jets[j].phi_std();
      float partonJetE = jets[j].E();
      float partonJetPSquared = pow(jets[j].px(), 2) + pow(jets[j].py(), 2) + pow(jets[j].pz(), 2);
      float partonJetMSquared = pow(jets[j].e(), 2) - partonJetPSquared;
      if (partonJetPt > 20 && partonJetMSquared > 0)
      {
        parton_jet_momenta.push_back(partonJetPt);
        parton_jet_momenta.push_back(partonJetEta);
        parton_jet_momenta.push_back(partonJetPhi);
        parton_jet_momenta.push_back(sqrt(partonJetMSquared));
        events.push_back(i);
        num_jets++;
      }
    }
  }
  auto end_generate_events = std::chrono::high_resolution_clock::now();

  auto start_use_model = std::chrono::high_resolution_clock::now();
  const std::vector<float> rand = sample_rand_vec(10 * num_jets);
  const std::vector<int64_t> rand_shape = {num_jets, 10};
  const std::vector<int64_t> jets_shape = {num_jets, 4};

  std::vector<float> normalized_parton_jet_momenta =
      normalize(parton_jet_momenta, parton_means, parton_std_devs);

  cppflow::tensor noise = cppflow::tensor(rand, rand_shape);
  cppflow::tensor jets_tensor =
      cppflow::tensor(normalized_parton_jet_momenta, jets_shape);
  auto output = model({{"serving_default_noiseIn", noise},
                       {"serving_default_pjetIn", jets_tensor}},
                      {"StatefulPartitionedCall:0"});
  cppflow::tensor *data = output.data();
  auto end_use_model = std::chrono::high_resolution_clock::now();
  auto start_unnormalize = std::chrono::high_resolution_clock::now();
  std::vector<float> data_vec = data->get_data<float>();
  std::vector<float> reco_jet_momenta =
      unnormalize(data_vec, reco_means, reco_std_devs);
  auto end_unnormalize = std::chrono::high_resolution_clock::now();
  auto start_write_out = std::chrono::high_resolution_clock::now();
  write_momenta(write_out, events, parton_jet_momenta, reco_jet_momenta);
  auto end_write_out = std::chrono::high_resolution_clock::now();

  // Statistics: full printout.
  //pythia.stat();

  write_out.close();
  auto end = std::chrono::high_resolution_clock::now();
  // energy_out.close();
  std::chrono::duration<double> load_constants = end_load_constants - start;
  std::chrono::duration<double> load_model = end_load_model - start_load_model;
  std::chrono::duration<double> generate_events = end_generate_events - start_generate_events;
  std::chrono::duration<double> use_model = end_use_model - start_use_model;
  std::chrono::duration<double> write_out_data = end_write_out - start_write_out;
  std::chrono::duration<double> unnormlize_data = end_unnormalize - start_unnormalize;
  std::chrono::duration<double> total = end - start;
  std::cout << "TIMING RESULTS: " << std::endl;
  std::cout << "load constants: " << load_constants.count() << std::endl;
  std::cout << "load model: " << load_model.count() << std::endl;
  std::cout << "generate events: " << generate_events.count() << std::endl;
  std::cout << "use model: " << use_model.count() << std::endl;
  std::cout << "write out: " << write_out_data.count() << std::endl;
  std::cout << "unnormalize: " << unnormlize_data.count() << std::endl;
  std::cout << "total: " << total.count() << std::endl;
  return 0;
}
