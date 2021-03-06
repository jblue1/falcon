CXXFLAGS = -Wall
ROOTFLAGS = `$(ROOTSYS)/bin/root-config --cflags --libs`
FASTJETFLAGS = `$(FASTJETSYS)/bin/fastjet-config --cxxflags --libs --plugins`
DELPHESFLAGS = -I ~/CMSSW_10_2_24/src/Delphes-3.4.2/ -I ~/CMSSW_10_2_24/src/Delphes-3.4.2/external/ -L ~/CMSSW_10_2_24/src/Delphes-3.4.2/ -lDelphes

makeHistos.out : src/analyzeJets.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/makeHistos.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)

makeCheckMass.out : src/checkMass.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/makeCheckMass.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)

writeJetMomentaAngular.out : src/writeJetMomentaAngular.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/writeJetMomentaAngular.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)

writeJetMomentaCartesian.out : src/writeJetMomentaCartesian.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/writeJetMomentaCartesian.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)

analyzeDelphesJets.out : src/analyzeDelphesJets.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/analyzeDelphesJets.out $^ $(ROOTFLAGS) $(FASTJETFLAGS) $(DELPHESFLAGS)

makeKDTreeBins.out : src/makeKDTreeBins.cc 
	$(CXX) $(CXXFLAGS) -o bin/makeKDTreeBins.out $^ $(ROOTFLAGS) 

histos : 
	./bin/makeHistos.out newPartonEvents.root histos.root momentaInfo

checkMass : 
	./bin/makeCheckMass.out newPartonEventsFixedZCoord.root 

angularData : 
	./bin/writeJetMomentaAngular.out newPartonEventsFixedZCoord.root newPartonMatchedJetsNoRecoPtCutFixRapMass.txt

cartesianData : 
	./bin/writeJetMomentaCartesian.out newPartonEventsFixedZCoord.root testingforPt25Jets.txt
	
delphes :
	./bin/analyzeDelphesJets.out /home/DAVIDSON/joblue/CMSSW_10_2_24/src/analyzer/MyAnalyzer/Delphes/adjusted_delphes_card_CMS.root ./data/raw/newPartonEventswPFCands.root

KDTreeBins : 
	./bin/makeKDTreeBins.out
.PHONY: angularData cartesianData histos delphes
