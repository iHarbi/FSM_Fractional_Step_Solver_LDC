#include <iostream>
#include <fstream>
#include <vector>
#include "mainMesh.h"

using namespace std;

class CSVWriter {
public:
    CSVWriter(const string& filename) : filename(filename) {}

    bool writeMatrix(const vector<vector<double>>& Velocity, const Position& Position_Struct) {
        ofstream outputFile(filename);
        if (!outputFile.is_open()) {
            cerr << "Error: Unable to open the file." << endl;
            return false;
        }
        // Write header cells
        outputFile << "X,Y,Z,U,V,W" << endl;

        for (size_t j = 0; j < Velocity[0].size(); ++j) {
            for (size_t i = 0; i < Velocity.size(); ++i) {
                outputFile << Position_Struct.CellCentroidX[i] << "," << Position_Struct.CellCentroidY[j] << ","<<0 << "," << Velocity[i][j] << ",";
                if (i < Velocity.size() - 1) {
                    outputFile << endl;
                }
            }
            outputFile << endl; // Add a blank line after each column
        }

        outputFile.close();
        return true;
    }

private:
    string filename;
};
