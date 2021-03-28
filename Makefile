CXXFLAGS = -Wall
ROOTFLAGS = `$(ROOTSYS)/bin/root-config --cflags --libs`
FASTJETFLAGS = `$(FASTJETSYS)/bin/fastjet-config --cxxflags --libs --plugins`
DELPHESFLAGS = -I ~/CMSSW_10_2_24/src/Delphes-3.4.2/ -I ~/CMSSW_10_2_24/src/Delphes-3.4.2/external/ -L ~/CMSSW_10_2_24/src/Delphes-3.4.2/ -lDelphes

makeHistos.out : src/analyzeJets.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/makeHistos.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)

writeJetMomentaAngular.out : src/writeJetMomentaAngular.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/writeJetMomentaAngular.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)

writeJetMomentaCartesian.out : src/writeJetMomentaCartesian.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/writeJetMomentaCartesian.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)

analyzeDelphesJets.out : src/analyzeDelphesJets.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/analyzeDelphesJets.out $^ $(ROOTFLAGS) $(FASTJETFLAGS) $(DELPHESFLAGS)

histos : 
	./bin/makeHistos.out newPartonEvents.root histos.root momentaInfo

angularData : 
	./bin/writeJetMomentaAngular.out newPartonEvents.root newPartonMatchedJetsNoRecoPtCut.txt

cartesianData : 
	./bin/writeJetMomentaCartesian.out newPartonEvents.root newPartonMatchedJetsNoRecoPtCutCartesian.txt
	
delphes :
	./bin/analyzeDelphesJets.out ~/test/Delphes/default_CMS_card.root ./data/raw/newPartonEventswPFCands.root

.PHONY: angularData cartesianData histos delphes
