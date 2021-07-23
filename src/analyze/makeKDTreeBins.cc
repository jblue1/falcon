#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TKDTreeBinning.h"

/**
 * Append a line of numbers seperated by spaces into the corresponding 
 * number of vectors. Will not work if the number of numbers in line 
 * and the number of vectors in vector are different.
 */
void append_numbers(std::vector<std::vector<double> > &vectors, std::string &line) {
    int num_spaces = 0;
    std::string number;
    for (int i=0; i < line.size(); i++) {
        for (int j = 0; j < vectors.size(); j++) {
            if (num_spaces == j) {
                number.push_back(line[i]);
            }
        }
        if (line[i] == ' ' || i == line.size() - 1) {
            vectors[num_spaces].push_back(std::stod(number));
            num_spaces++;
            number.clear();
        }
    }
}

int main(int arc, char const *argv[]) {
    std::ifstream myfile ("./data/processed/newPartonMatchedJetsNoRecoPtCutFixRapMass.txt"); 

    int num_data_points = 10000;
    std::vector<double> partonPt;
    std::vector<double> partonEta;
    std::vector<double> partonPhi;
    std::vector<double> partonE;
    std::vector<double> recoPt;
    std::vector<double> recoEta;
    std::vector<double> recoPhi;
    std::vector<double> recoE;
    std::vector<std::vector<double> > jets;
    jets.push_back(partonPt);
    jets.push_back(partonEta);
    jets.push_back(partonPhi);
    jets.push_back(partonE);
    jets.push_back(recoPt);
    jets.push_back(recoEta);
    jets.push_back(recoPhi);
    jets.push_back(recoE);

    int loop = 0;
    std::string line;
    if (myfile.is_open()) {
        // skip past first two lines
        std::getline(myfile,line); 
        std::getline(myfile,line); 
            while (loop < num_data_points) {
                std::getline(myfile,line); 
                append_numbers(jets, line);
                loop++;
            }
            myfile.close(); 
        }

    std::vector<double> data_vec;
    for (int i = 0; i < 8; i++) {
        data_vec.insert(data_vec.end(), jets[i].begin(), jets[i].end());
    }

    int num_bins = 10;
    TKDTreeBinning test = TKDTreeBinning(num_data_points, 8, data_vec, num_bins);

    
    std::ofstream write_out;
    write_out.open("./data/processed/testbins.txt");
    for (int i = 0; i < num_bins; i++) {
        std::pair<const Double_t *, const Double_t * > edges =  test.GetBinEdges(i);
        for (int j = 0; j < 8; j++){
            write_out << *(edges.first + j) << " ";
        }
        for (int j = 0; j < 8; j++){
            write_out << *(edges.second + j) << " ";
        }
        write_out << std::endl;
    }
    write_out.close();
    exit(0);
}
