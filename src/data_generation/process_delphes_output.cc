#include <iostream>

#include "TClonesArray.h"
#include "TChain.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
int main(int argc, char const *argv[]) {
    TChain chain("Delphes");
    chain.Add("test.root");

    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numDelphesEvents = treeReader->GetEntries();
    TClonesArray *branchRecoJet = treeReader->UseBranch("Jet");

    for (int i = 0; i < numDelphesEvents; i++) {
        treeReader->ReadEntry(i);

        int numJets = branchRecoJet->GetEntries();
        std::cout << "There were " << numJets << " jets in the event" << std::endl;
        for (int j = 0; j < numJets; j++ ) {
            Jet *jet = (Jet*) branchRecoJet->At(i);
            std::cout << "Jet pt: " << jet->PT << std::endl;
        }

    }

    delete treeReader;
    return 0;
}
