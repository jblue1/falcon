#include <iostream>

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
int main(int argc, char const *argv[]) {
    TChain chain("Delphes");
    chain.Add("test.root");

    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numDelphesEvents = treeReader->GetEntries();
    TClonesArray *branchRecoJet = treeReader->UseBranch("Jet");

    delete treeReader;
    return 0;
}