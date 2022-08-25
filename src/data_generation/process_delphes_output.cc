#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "classes/DelphesClasses.h"

#include "../helpers.h"

void parseLine(std::string line, std::vector<float> &jet,
               int &currentEventNumber)
{
    jet.clear();
    std::string entry;
    std::stringstream s(line);
    for (int i = 0; i < 5; i++)
    {
        getline(s, entry, ' ');
        if (i < 4)
        {
            jet.push_back(std::stof(entry));
        }
        else
        {
            currentEventNumber = std::stoi(entry);
        }
    }
}

void getJetsFromFile(std::ifstream &file, int eventNumber,
                     std::vector<std::vector<float>> &currentEventJets,
                     std::vector<float> &firstJetNextEvent)
{
    std::string line;
    int currentEventNumber = eventNumber;
    while (currentEventNumber == eventNumber)
    {
        if (getline(file, line))
        {
            std::vector<float> jet;
            parseLine(line, jet, currentEventNumber);
            if (currentEventNumber == eventNumber)
            {
                currentEventJets.push_back(jet);
            }
            else
            {
                firstJetNextEvent = jet;
            }
        }
        else
        {
            currentEventNumber++;
        }
    }
}

int main(int argc, char const *argv[])
{
    std::string data_dir(argv[1]);
    std::string delphesJetFileString =
        data_dir + std::string("/delphesOut.root");
    std::string partonJetFileString = data_dir + std::string("/partonJets.txt");
    std::string matchedJetsFileString =
        data_dir + std::string("/matchedJets.txt");
    TChain chain("Delphes");
    chain.Add(delphesJetFileString.c_str());
    std::ifstream readFile;
    readFile.open(partonJetFileString);
    std::string headerLine;
    getline(readFile, headerLine); // first line says what kind of jets
    getline(readFile, headerLine); // second line is just column headers

    std::ofstream writeFile(matchedJetsFileString);

    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numDelphesEvents = treeReader->GetEntries();
    TClonesArray *branchRecoJet = treeReader->UseBranch("Jet");

    std::vector<std::vector<float>> partonJetsCurrentEvent;
    int numTotalPartonJets = 0;
    int numMatchedPartonJets = 0;
    for (int i = 0; i < numDelphesEvents; i++)
    {
        std::vector<float> firstPartonJetNextEvent;
        getJetsFromFile(readFile, i, partonJetsCurrentEvent,
                        firstPartonJetNextEvent);
        int numPartonJetsInEvent = partonJetsCurrentEvent.size();
        std::cout << "Number of parton jets in event " << i << ": "
                  << numPartonJetsInEvent << std::endl;
        treeReader->ReadEntry(i);
        int numRecoJets = branchRecoJet->GetEntries();
        std::cout << "Number of reco jets in event " << i << ": " << numRecoJets
                  << std::endl;

        // loop through parton jets looking for reco jet matches
        if (partonJetsCurrentEvent.size() > 0)
        {
            for (int j = 0; j < partonJetsCurrentEvent.size(); j++)
            {
                if (partonJetsCurrentEvent[j][0] > 20)
                {
                    numTotalPartonJets++;
                    float minDR = 10.0;
                    int recoJetIndex = 0;

                    for (int k = 0; k < numRecoJets; k++)
                    {
                        Jet *jet = (Jet *)branchRecoJet->At(k);
                        float dR = deltaR(jet->Eta, jet->Phi,
                                          partonJetsCurrentEvent[j][1],
                                          partonJetsCurrentEvent[j][2]);
                        if (dR < minDR)
                        {
                            minDR = dR;
                            recoJetIndex = k;
                        }
                    }
                    if (minDR < 0.35 && partonJetsCurrentEvent[j][3] > 0)
                    {
                        numMatchedPartonJets++;
                        Jet *matchedRecoJet =
                            (Jet *)branchRecoJet->At(recoJetIndex);
                        if (matchedRecoJet->Mass > 0)
                        {
                            writeFile << partonJetsCurrentEvent[j][0] << " "
                                      << partonJetsCurrentEvent[j][1] << " "
                                      << partonJetsCurrentEvent[j][2] << " "
                                      << partonJetsCurrentEvent[j][3] << " "
                                      << matchedRecoJet->PT << " "
                                      << matchedRecoJet->Eta << " "
                                      << matchedRecoJet->Phi << " "
                                      << pow(matchedRecoJet->Mass, 2)
                                      << std::endl;
                        }
                    }
                }
            }
        }

        partonJetsCurrentEvent.clear();
        if (i < numDelphesEvents - 1 && firstPartonJetNextEvent.size() > 0)
        {
            partonJetsCurrentEvent.push_back(firstPartonJetNextEvent);
        }
        firstPartonJetNextEvent.clear();
    }

    float fractionMatched = (float)numMatchedPartonJets / numTotalPartonJets;
    std::cout << "Fraction of matched jets = " << fractionMatched << std::endl;
    readFile.close();
    writeFile.close();
    delete treeReader;
    return 0;
}
