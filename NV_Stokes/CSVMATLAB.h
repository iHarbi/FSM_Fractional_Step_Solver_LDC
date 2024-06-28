#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class CSVMATLAB {
public:
    CSVMATLAB(const string& filename) : filename(filename) {}

    bool writeMatrix(const vector<vector<double>>& matrix) {
        ofstream outputFile(filename);
        if (!outputFile.is_open()) {
            cerr << "Error: Unable to open the file." << endl;
            return false;
        }

        for (size_t i = 0; i < matrix.size(); ++i) {
            for (size_t j = 0; j < matrix[i].size(); ++j) {
                outputFile << matrix[i][j];
                if (j < matrix[i].size() - 1) {
                    outputFile << ",";
                }
            }
            outputFile << endl;
        }

        outputFile.close();
        return true;
    }

private:
    string filename;
};