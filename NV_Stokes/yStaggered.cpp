#include "yStaggered.h"
#include <iostream>
#include <vector>
#include <string>
#include "mainMesh.h"
#include <iomanip>

StaggeredMesh_Y yStaggered::structGetter_Y(Cell Cell_, Position Position_, Area Area_)
{
	/*Object Initiation*/
	StaggeredMesh_Y yStaggered;
	/*Fields[Cell, Area, Position] Augmentation*/
	yStaggered.Cell = Cell_;
	yStaggered.Area = Area_;
	yStaggered.Position = Position_;

	return yStaggered;
}

CellY yStaggered::yStaggeringOPs(StaggeredMesh_Y mainMesh)
{

	/*Recall StagY Size is [MeshSize_X + 2, MeshSize_Y + 1]*/
	
	mesh2Modify_Y.CellCentroid_X_Staggered.resize(MeshSize_X + 2); /*StagY CellCentroids in X-Direction*/
	mesh2Modify_Y.CellCentroid_Y_Staggered.resize(MeshSize_Y + 1); /*StagY CellCentroids in Y-Direction*/
	mesh2Modify_Y.dPE.resize(MeshSize_X + 2); /*StagY dPE[Cell- East] Distance*/
	mesh2Modify_Y.dPW.resize(MeshSize_X + 2); /*StagY dPW[Cell- West] Distance*/
	mesh2Modify_Y.dPS.resize(MeshSize_Y + 1); /*StagY dPS[Cell-South] Distance*/
	mesh2Modify_Y.dPN.resize(MeshSize_Y + 1); /*StagY dPN[Cell-North] Distance*/
	mesh2Modify_Y.AreaEast.resize(MeshSize_X , vector<double>(MeshSize_Y + 1 ));
	mesh2Modify_Y.AreaWest.resize(MeshSize_X , vector<double>(MeshSize_Y + 1 ));
	mesh2Modify_Y.AreaNorth.resize(MeshSize_X, vector<double>(MeshSize_Y + 1 ));
	mesh2Modify_Y.AreaSouth.resize(MeshSize_X, vector<double>(MeshSize_Y + 1 ));
	mesh2Modify_Y.Volume.resize(MeshSize_X, vector<double>(MeshSize_Y + 1));
	mesh2Modify_Y.D_East.resize(MeshSize_X  , vector<double>(MeshSize_Y + 1));
	mesh2Modify_Y.D_West.resize(MeshSize_X  , vector<double>(MeshSize_Y + 1));
	mesh2Modify_Y.D_North.resize(MeshSize_X , vector<double>(MeshSize_Y + 1));
	mesh2Modify_Y.D_South.resize(MeshSize_X , vector<double>(MeshSize_Y + 1));


	for (int i = 0; i <= MeshSize_X + 1; ++i) { mesh2Modify_Y.CellCentroid_X_Staggered[i] = mainMesh.Position.CellCentroidX[i]; }

	for (int i = 0; i <= MeshSize_Y; ++i) { mesh2Modify_Y.CellCentroid_Y_Staggered[i] = mainMesh.Position.FaceY[i]; }

	/*dPE*/
	for (int i = 0; i < MeshSize_X + 1; ++i)
		//mesh2Modify_Y.dPE[i] = mesh2Modify_Y.CellCentroid_X_Staggered[i + 1] - mesh2Modify_Y.CellCentroid_X_Staggered[i];
	    mesh2Modify_Y.dPE[i] = mainMesh.Cell.dPE[i];
	/*dPW*/
	for (int i = 1; i <= MeshSize_X + 1; ++i)
		//mesh2Modify_Y.dPW[i] = mesh2Modify_Y.CellCentroid_X_Staggered[i] - mesh2Modify_Y.CellCentroid_X_Staggered[i - 1];
		mesh2Modify_Y.dPW[i] = mainMesh.Cell.dPW[i];
	/*dPN*/
	for (int i = 0; i < MeshSize_Y; ++i)
		mesh2Modify_Y.dPN[i] = mesh2Modify_Y.CellCentroid_Y_Staggered[i + 1] - mesh2Modify_Y.CellCentroid_Y_Staggered[i];
	/*dPS*/
	for (int i = 1; i <= MeshSize_Y ; ++i)
		mesh2Modify_Y.dPS[i] = mesh2Modify_Y.CellCentroid_Y_Staggered[i] - mesh2Modify_Y.CellCentroid_Y_Staggered[i - 1];
	/*---------------------Areas---------------------*/
	for (int i = 0; i < MeshSize_X; ++i)
		for (int j = 0; j <= MeshSize_Y; ++j)
			mesh2Modify_Y.AreaEast[i][j] = mainMesh.Area.East[i];
	for (int i = 0; i < MeshSize_X; ++i)
		for (int j = 0; j <= MeshSize_Y; ++j)
			mesh2Modify_Y.AreaWest[i][j] = mainMesh.Area.West[i];
	for (int i = 0; i < MeshSize_X; ++i)
		for (int j = 0; j <= MeshSize_Y; ++j)
			mesh2Modify_Y.AreaNorth[i][j] = mainMesh.Area.North[i];
	for (int i = 0; i < MeshSize_X; ++i)
		for (int j = 0; j <= MeshSize_Y; ++j)
			mesh2Modify_Y.AreaSouth[i][j] = mainMesh.Area.South[i];
	/*---------------------Volume---------------------*/
	for (int i = 0; i < MeshSize_X; ++i)
		for (int j = 0; j < MeshSize_Y + 1; ++j)
			mesh2Modify_Y.Volume[i][j] = mesh2Modify_Y.AreaEast[i][j] * mesh2Modify_Y.AreaNorth[i][j];
	/*D_East | D_West | D_North | D_South | [MeshSize_X,  MeshSize_Y + 1]*/
	for (int i = 0; i < MeshSize_X; ++i)
		for (int j = 0; j < MeshSize_Y + 1 ; ++j)
		{
			mesh2Modify_Y.D_East[i][j] =  mesh2Modify_Y.AreaEast[i][j]  /  mesh2Modify_Y.dPE[i + 1];
			mesh2Modify_Y.D_West[i][j] =  mesh2Modify_Y.AreaWest[i][j]  /  mesh2Modify_Y.dPW[i + 1];
			mesh2Modify_Y.D_North[i][j] = mesh2Modify_Y.AreaNorth[i][j] * MeshSize_Y;/*+1 is added here */
			mesh2Modify_Y.D_South[i][j] = mesh2Modify_Y.AreaSouth[i][j] * MeshSize_Y;/*+1 is added here */
		
			/* mesh2Modify_Y.D_East[i][j] = mesh2Modify_Y.AreaEast[i][j]  * MeshSize_X; 
			mesh2Modify_Y.D_West[i][j] = mesh2Modify_Y.AreaWest[i][j]  * MeshSize_X; 
			mesh2Modify_Y.D_North[i][j] = mesh2Modify_Y.AreaNorth[i][j] * MeshSize_Y;
			mesh2Modify_Y.D_South[i][j] = mesh2Modify_Y.AreaSouth[i][j] * MeshSize_Y; */
		}



	return mesh2Modify_Y;
}

void yStaggered::checkMesh(CellY cell)
{
	cout << "\033[1;32m-+-+-+-+-+-+-+-+-+-+-+-+checkMesh | Staggered-Y-Mesh | [N + 2 * M + 1] -+-+-+-+-+-+-+-+-+-+-+-+\033[0m" << endl;
	cout << "\033[1;31mdPE\t" << "\033[1;31mdPW\t" << "\033[1;31mdPN\t" << "\033[1;31mdPS\033[0m" << endl;
	size_t maxSize = max({ cell.dPE.size(), cell.dPW.size(), cell.dPN.size(), cell.dPS.size() });
	for (size_t i = 0; i < maxSize; ++i) {
		// Print dPE if available
		if (i < cell.dPE.size()) {
			cout << cell.dPE[i] << "\t";
		}
		else {
			cout << "N/A\t"; // Print N/A if (dPE) doesn't have an element at index i
		}

		// Print dPW if available
		if (i < cell.dPW.size()) {
			cout << cell.dPW[i] << "\t";
		}
		else {
			cout << "N/A\t"; // Print N/A if dPW doesn't have an element at index i
		}

		// Print dPN if available
		if (i < cell.dPN.size()) {
			cout << cell.dPN[i] << "\t";
		}
		else {
			cout << "N/A\t"; // Print N/A if dPN doesn't have an element at index i
		}

		// Print dPS if available
		if (i < cell.dPS.size()) {
			cout << cell.dPS[i] << endl;
		}
		else {
			cout << "N/A" << endl; // Print N/A if dPS doesn't have an element at index i
		}
	}

}

void yStaggered::massFlowRate_StaggeredY_Initializer()
{
	
	massFlowRate_StaggeredY_.massFlowRate_East.resize(MeshSize_X  , vector<double>(MeshSize_Y + 1));
	massFlowRate_StaggeredY_.massFlowRate_West.resize(MeshSize_X  , vector<double>(MeshSize_Y + 1));
	massFlowRate_StaggeredY_.massFlowRate_North.resize(MeshSize_X , vector<double>(MeshSize_Y + 1));
	massFlowRate_StaggeredY_.massFlowRate_South.resize(MeshSize_X , vector<double>(MeshSize_Y + 1));
}

void yStaggered::massFlowRate_StaggeredY_Setter(v_Staggered_Y v_Staggered_Y_, CellY CellY_, u_Staggered_X u_Staggered_X_, double roh)
{
	/*massFlowRate Required Variables| East and West: uA, uB and A_An, A_Bn */
	for (int i = 0; i < MeshSize_X ; ++i)
		for (int j = 0; j < MeshSize_Y + 1; ++j)
		{
			massFlowRate_StaggeredY_.massFlowRate_East[i][j] = (u_Staggered_X_.u_P[i + 1][j + 1] * CellY_.AreaEast[i][j] + u_Staggered_X_.u_P[i + 1][j] * CellY_.AreaEast[i][j]) * 0.5 * roh;
			massFlowRate_StaggeredY_.massFlowRate_West[i][j] = (u_Staggered_X_.u_P[i][j + 1] * CellY_.AreaWest[i][j] + u_Staggered_X_.u_P[i][j] * CellY_.AreaWest[i][j]) * 0.5 * roh;
		}
		
	/*Ghost Cell Concept Would have been an Appropriate Implementation instead of Bifurcation*/
	for (int i = 0; i < MeshSize_X; ++i) {
		for (int j = 0; j < MeshSize_Y + 1; ++j) {
			if (j < MeshSize_Y) {
				massFlowRate_StaggeredY_.massFlowRate_North[i][j] = 0.5 * (v_Staggered_Y_.v_P[i + 1][j] + v_Staggered_Y_.v_P[i + 1][j + 1]) * CellY_.AreaNorth[i][j] * roh;
			}
			else {
				// Handling B.C when j = MeshSize_Y
				massFlowRate_StaggeredY_.massFlowRate_North[i][j] = (v_Staggered_Y_.v_P[i + 1][j]) * CellY_.AreaNorth[i][j] * roh;
			}
			if (j > 0) {
				massFlowRate_StaggeredY_.massFlowRate_South[i][j] = 0.5 * (v_Staggered_Y_.v_P[i + 1][j] + v_Staggered_Y_.v_P[i + 1][j - 1]) * CellY_.AreaSouth[i][j] * roh;
			}
			else {
				// Handling B.C when j = 0
				massFlowRate_StaggeredY_.massFlowRate_South[i][j] = (v_Staggered_Y_.v_P[i + 1][j]) * CellY_.AreaSouth[i][j] * roh;
			}
		}
	}
}

void yStaggered::v_Staggered_Y_Initializer()
{
	
	v_Staggered_Y_.v_P.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 1));
	v_Staggered_Y_.v_E.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 1));
	v_Staggered_Y_.v_W.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 1));
	v_Staggered_Y_.v_N.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 1));
	v_Staggered_Y_.v_S.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 1));
}




void yStaggered::v_Staggered_Y_OPs(Cell Cell_Struct, vector<vector<double>> V_Pred, double Delta_t, double roh, vector<vector<double>> PressureField_New)
{

	for (int i = 1; i < MeshSize_X + 1; ++i) {
		for (int j = 1; j < MeshSize_Y; ++j) {
			v_Staggered_Y_.v_P[i][j] = V_Pred[i][j];
			v_Staggered_Y_.v_P[i][j] -= (Delta_t / roh) * (PressureField_New[i][j+1] - PressureField_New[i][j]) / Cell_Struct.dPN[j];
		}
	}
}

void yStaggered::printMassFlowRates() {
	// Print massFlowRate_East
	std::cout << "\033[1;32myStaggeredMass Flow Rates:\n";
	std::cout << std::setw(15) << "East" << std::setw(15) << "West" << std::setw(15) << "North" << std::setw(15) << "South" << std::endl;
	for (int i = 0; i < MeshSize_X; ++i) {
		for (int j = 0; j < MeshSize_Y + 1; ++j) {
			std::cout << std::setw(15) << massFlowRate_StaggeredY_.massFlowRate_East[i][j]
				<< std::setw(15) << massFlowRate_StaggeredY_.massFlowRate_West[i][j]
				<< std::setw(15) << massFlowRate_StaggeredY_.massFlowRate_North[i][j]
				<< std::setw(15) << massFlowRate_StaggeredY_.massFlowRate_South[i][j]
				<< std::endl;
		}
		std::cout << std::endl;
	}
}














yStaggered::yStaggered(float L, float H, float N, float M) :L(L), H(H), MeshSize_X(N), MeshSize_Y(M) {}
yStaggered::~yStaggered() {}
