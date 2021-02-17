CXXFLAGS = -Wall
ROOTFLAGS = `$(ROOTSYS)/bin/root-config --cflags --libs`
FASTJETFLAGS = `$(FASTJETSYS)/bin/fastjet-config --cxxflags --libs --plugins`
DELPHESFLAGS = -I ~/CMSSW_10_2_24/src/Delphes-3.4.2/ -I ~/CMSSW_10_2_24/src/Delphes-3.4.2/external/ -L ~/CMSSW_10_2_24/src/Delphes-3.4.2/ -lDelphes

makeHistos.out : src/analyzeJets.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/makeHistos.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)

writeJetMomenta.out : src/writeJetMomenta.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/writeJetMomenta.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)

analyzeDelphesJets.out : src/analyzeDelphesJets.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/analyzeDelphesJets.out $^ $(ROOTFLAGS) $(FASTJETFLAGS) $(DELPHESFLAGS)

histos : 
	./bin/makeHistos.out newPartonEvents.root histos.root momentaInfo

data : 
	./bin/writeJetMomenta.out newPartonEvents.root newPartonMatchedJetsNoRecoPtCut.txt
	
delphes :
	./bin/analyzeDelphesJets.out ~/CMSSW_10_2_24/src/analyzer/MyAnalyzer/Delphes/test_CMS_PhaseII_0PU_card.root ./data/raw/newPartonEvents.root

.PHONY: data histos delphes
