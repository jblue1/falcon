#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TKDTreeBinning.h"

void append_numbers(std::vector<double> &first, std::vector<double> &second, 
                    std::vector<double> &third, std::vector<double> &fourth,
                   std::string &line) {
    int num_spaces = 0;
    std::string first_number;
    std::string second_number;
    std::string third_number;
    std::string fourth_number;
    for (int i=0; i < line.size(); i++) {
        if (num_spaces == 4) {
            first_number.push_back(line[i]);
        }
        if (num_spaces == 5) {
            second_number.push_back(line[i]);
        }
        if (num_spaces == 6) {
            third_number.push_back(line[i]);
        }
        if (num_spaces == 7) {
            fourth_number.push_back(line[i]);
        }

        if (line[i] == ' ') {
            num_spaces++;
        }
    }
    
    first.push_back(std::stod(first_number));
    second.push_back(std::stod(second_number));
    third.push_back(std::stod(third_number));
    fourth.push_back(std::stod(fourth_number));
}

int main(int argc, char* argv[]) {
    int NUM_FEATURES = 8;
    std::string dataFile(argv[1]);
    std::string binsFile(argv[2]);
    int num_data_points = atoi(argv[3]);
    int num_bins = atoi(argv[4]);
    std::ifstream myfile (dataFile); 

    std::vector<double> Pt;
    std::vector<double> eta;
    std::vector<double> phi;
    std::vector<double> E;
    int loop = 0;
    std::string line;
    if (myfile.is_open()) {
        // skip past first two lines
        getline (myfile,line); 
        getline (myfile,line); 
            while (loop < num_data_points) {
                getline (myfile,line); 
                append_numbers(Pt, eta, phi, E, line);
                loop++;
            }
            myfile.close(); 
    }
    std::vector<double> data_vec;
    data_vec.insert(data_vec.end(), Pt.begin(), Pt.end());
    data_vec.insert(data_vec.end(), eta.begin(), eta.end());
    data_vec.insert(data_vec.end(), phi.begin(), phi.end());
    data_vec.insert(data_vec.end(), E.begin(), E.end());

    TKDTreeBinning test = TKDTreeBinning(num_data_points, 4, data_vec, num_bins);

    std::ofstream write_out;
    write_out.open(binsFile);
    for (int i = 0; i < num_bins; i++) {
        std::pair<const Double_t *, const Double_t * > edges =  test.GetBinEdges(i);
        for (int j = 0; j < 4; j++){
            write_out << *(edges.first + j) << " ";
        }
        for (int j = 0; j < 4; j++){
            write_out << *(edges.second + j) << " ";
        }
        write_out << std::endl;
    }
    write_out.close();
    return 0; 
}
