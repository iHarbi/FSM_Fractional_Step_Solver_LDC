#pragma once

#include <iostream>
#include <vector>

struct LinePlots
{
	std::vector<std::vector<double>> u_y;
	std::vector<std::vector<double>> v_x;
};

class Export
{
private:
	LinePlots Lineplot;

public:
    // Function to add u_y data to LinePlots, extracting specific rows and columns
    void addUyData(const std::vector<std::vector<double>>& uCurrent, int row)
    {
        std::vector<double> extractedRow;
        for (size_t i = 0; i < uCurrent[row].size(); ++i) {
            extractedRow.push_back(uCurrent[row][i]);
        }
        Lineplot.u_y.push_back(extractedRow);
    }

    // Function to add v_x data to LinePlots, extracting specific rows and columns
    void addVxData(const std::vector<std::vector<double>>& vCurrent, int col)
    {
        std::vector<double> extractedCol;
        for (size_t i = 0; i < vCurrent.size(); ++i) {
            extractedCol.push_back(vCurrent[i][col]);
        }
        Lineplot.v_x.push_back(extractedCol);
    }

    LinePlots& getLinePlots() { return Lineplot; }
};

