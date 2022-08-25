CXXFLAGS = -Wall -O3
ROOTFLAGS = `$(ROOTSYS)/bin/root-config --cflags --libs`
FASTJETFLAGS = `$(FASTJETSYS)/bin/fastjet-config --cxxflags --libs --plugins`
DELPHESFLAGS = -I $(DELPHES_PATH) -I $(DELPHES_PATH)/external -L $(DELPHES_PATH) -lDelphes


makeHistos.out : src/analyzeJets.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/makeHistos.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)

makeCheckMass.out : src/checkMass.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/makeCheckMass.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)

writeJetMomentaAngular.out : src/writeJetMomentaAngular.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/writeJetMomentaAngular.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)

writeJetMomentaJetTagged.out : src/writeJetMomentaJetTagged.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/writeJetMomentaJetTagged.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)

writeJetMomentaCartesian.out : src/writeJetMomentaCartesian.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/writeJetMomentaCartesian.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)

analyzeDelphesJets.out : src/analyzeDelphesJets.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/analyzeDelphesJets.out $^ $(ROOTFLAGS) $(FASTJETFLAGS) $(DELPHESFLAGS)

event_gen.out : src/data_generation/event_gen.cc
	$(CXX) $(CXXFLAGS) -o bin/event_generation.out $^ $(DELPHESFLAGS) $(ROOTFLAGS) '-L$(CONDA_PREFIX)/lib/' -lpythia8 -lEG

process_delphes_output.out : src/data_generation/process_delphes_output.cc
	$(CXX) $(CXXFLAGS) -o bin/process_delphes_output.out $^ $(DELPHESFLAGS) $(ROOTFLAGS)

makeKDTreeBins.out : src/analyze/makeKDTreeBins.cc
	$(CXX) $(CXXFLAGS) -o bin/makeKDTreeBins.out $^ $(ROOTFLAGS)

makeRecoKDTreeBins.out : src/analyze/makeRecoKDTreeBins.cc
	$(CXX) $(CXXFLAGS) -o bin/makeRecoKDTreeBins.out $^ $(ROOTFLAGS)

makeGetNumEvents.out : src/getNumberofEvents.cc
	$(CXX) $(CXXFLAGS) -o bin/getNumberofEvents.out $^ $(ROOTFLAGS)

histos :
	./bin/makeHistos.out newPartonEvents.root histos.root momentaInfo

checkMass :
	./bin/makeCheckMass.out newPartonEventsFixedZCoord.root

angularData :
	./bin/writeJetMomentaAngular.out newPartonEventsJetTagged.root newPartonMatchedJetsNoRecoPtCutFixRapMassEventNumberJetTagged.txt

jetTaggedData :
	./bin/writeJetMomentaJetTagged.out newPartonEventsJetTagged.root newPartonMatchedJetsNoRecoPtCutJetTagged.txt

cartesianData :
	./bin/writeJetMomentaCartesian.out newPartonEventsJetTagged.root testing.txt

delphes :
	./bin/analyzeDelphesJets.out /home/DAVIDSON/joblue/CMSSW_10_2_24/src/analyzer/MyAnalyzer/Delphes/adjusted_delphes_card_CMS.root ./data/raw/newPartonEventswPFCands.root

KDTreeBins :
	./bin/makeKDTreeBins.out
.PHONY: angularData cartesianData histos delphes
