#include "xStaggered.h"

StaggeredMesh_X xStaggered::structGetter_X(Cell Cell, Position Position, Area Area)
{
	/*Object Initiation*/
	StaggeredMesh_X xStaggered;
	/*Fields[Cell, Area, Position] Augmentation*/
	xStaggered.Cell = Cell;
	xStaggered.Area = Area;
	xStaggered.Position = Position;

		return xStaggered;
}

CellX xStaggered::xStaggeringOPs(StaggeredMesh_X mainMesh)
{
	/*mesh2Modify --> Staggered Mesh Information | mainMesh Struct is the Encapsulated mainMesh Information*/
	/*Recall StagX Size is [MeshSize_X + 1, MeshSize_Y + 2]*/
	mesh2Modify.CellCentroid_X_Staggered.resize(MeshSize_X + 1); /*StagX CellCentroids in X-Direction*/
	mesh2Modify.CellCentroid_Y_Staggered.resize(MeshSize_Y + 2); /*StagX CellCentroids in Y-Direction*/
	mesh2Modify.dPE.resize(MeshSize_X + 1); /*StagX dPE[Cell-East] Distance*/
	mesh2Modify.dPW.resize(MeshSize_X + 1); /*StagX dPW[Cell-West] Distance*/
	mesh2Modify.dPN.resize(MeshSize_Y + 2); /*StagX dPN[Cell-North] Distance*/
	mesh2Modify.dPS.resize(MeshSize_Y + 2); /*StagX dPS[Cell-South] Distance*/
	mesh2Modify.AreaEast.resize(MeshSize_X + 1, vector<double> (MeshSize_Y) );
	mesh2Modify.AreaWest.resize(MeshSize_X + 1, vector<double>(MeshSize_Y));
	mesh2Modify.AreaNorth.resize(MeshSize_X + 1, vector<double>(MeshSize_Y));
	mesh2Modify.AreaSouth.resize(MeshSize_X + 1, vector<double>(MeshSize_Y));
	mesh2Modify.Volume.resize(MeshSize_X + 1, vector<double>(MeshSize_Y));
	mesh2Modify.D_East.resize(MeshSize_X + 1, vector<double>(MeshSize_Y));
	mesh2Modify.D_West.resize(MeshSize_X + 1, vector<double>(MeshSize_Y));
	mesh2Modify.D_North.resize(MeshSize_X + 1, vector<double>(MeshSize_Y));
	mesh2Modify.D_South.resize(MeshSize_X + 1, vector<double>(MeshSize_Y));

	for (int i = 0; i <= MeshSize_X; ++i) { mesh2Modify.CellCentroid_X_Staggered[i] = mainMesh.Position.FaceX[i]; }
		
	for (int i = 0; i <= MeshSize_Y + 1; ++i) { mesh2Modify.CellCentroid_Y_Staggered[i] = mainMesh.Position.CellCentroidY[i];}
	
	/*dPE*/
	for (int i = 0; i < MeshSize_X; ++i){
	mesh2Modify.dPE[i] = mesh2Modify.CellCentroid_X_Staggered[i + 1] - mesh2Modify.CellCentroid_X_Staggered[i];
	}
	/*dPW*/
	for (int i = 1; i <= MeshSize_X; ++i)
		mesh2Modify.dPW[i] = mesh2Modify.CellCentroid_X_Staggered[i] - mesh2Modify.CellCentroid_X_Staggered[i-1];
	/*dPN*/
	for (int i = 0; i <= MeshSize_Y + 1; ++i)
		mesh2Modify.dPN[i] = mainMesh.Cell.dPN[i];
	/*dPS*/
	for (int i = 0; i <= MeshSize_Y + 1; ++i)
		mesh2Modify.dPS[i] = mainMesh.Cell.dPS[i];
	/*---------------------Areas---------------------*/
	for (int i = 0; i <= MeshSize_X; ++i)
		for (int j = 0; j < MeshSize_Y; ++j)
			mesh2Modify.AreaEast[i][j] = mainMesh.Area.East[i];
	for (int i = 0; i <= MeshSize_X; ++i)
		for (int j = 0; j < MeshSize_Y; ++j)
			mesh2Modify.AreaWest[i][j] = mainMesh.Area.West[i];
	for (int i = 0; i <= MeshSize_X; ++i)
		for (int j = 0; j < MeshSize_Y; ++j)
			mesh2Modify.AreaNorth[i][j] = mainMesh.Area.North[j];
	for (int i = 0; i <= MeshSize_X; ++i)
		for (int j = 0; j < MeshSize_Y; ++j)
			mesh2Modify.AreaSouth[i][j] = mainMesh.Area.South[j];
	/*---------------------Volume---------------------*/
	for (int i = 0; i <= MeshSize_X; ++i)
		for (int j = 0; j < MeshSize_Y; ++j)
			mesh2Modify.Volume[i][j] = mesh2Modify.AreaEast[i][j] * mesh2Modify.AreaNorth[i][j];
	/*D_East | D_West | D_North | D_South | [MeshSize_X + 1,  MeshSize_Y]*/
	for (int i = 1; i < MeshSize_X ; ++i)
		for (int j = 0; j < MeshSize_Y; ++j)
		{
			mesh2Modify.D_East[i][j] = mesh2Modify.AreaEast[i][j]   / mesh2Modify.dPE[i];
			mesh2Modify.D_West[i][j] = mesh2Modify.AreaWest[i][j]   / mesh2Modify.dPW[i];
			mesh2Modify.D_North[i][j] = mesh2Modify.AreaNorth[i][j] / mesh2Modify.dPN[j+1];
			mesh2Modify.D_South[i][j] = mesh2Modify.AreaSouth[i][j] / mesh2Modify.dPS[j+1];
			/*This is not RIGHT; GET BACK TO THE PRESVIOUS ONE*/
			/*mesh2Modify.D_East[i][j] = mesh2Modify.AreaEast[i][j] * MeshSize_X;
			mesh2Modify.D_West[i][j] = mesh2Modify.AreaWest[i][j] * MeshSize_X;
			mesh2Modify.D_North[i][j] = mesh2Modify.AreaNorth[i][j] * MeshSize_Y;
			mesh2Modify.D_South[i][j] = mesh2Modify.AreaSouth[i][j] * MeshSize_Y;
			*/

		}
	cout << " TRIAL";
	return mesh2Modify;
}

void xStaggered::massFlowRate_StaggeredX_Initializer()
{
	massFlowRate_StaggeredX_.massFlowRate_East.resize(MeshSize_X + 1, vector<double>(MeshSize_Y));
	massFlowRate_StaggeredX_.massFlowRate_West.resize(MeshSize_X + 1, vector<double>(MeshSize_Y));
	massFlowRate_StaggeredX_.massFlowRate_North.resize(MeshSize_X + 1, vector<double>(MeshSize_Y));
	massFlowRate_StaggeredX_.massFlowRate_South.resize(MeshSize_X + 1, vector<double>(MeshSize_Y));
}

void xStaggered::u_Staggered_X_Initializer()
{
	u_Staggered_X_.u_P.resize(MeshSize_X +1, vector<double>( MeshSize_Y + 2 ));
	u_Staggered_X_.u_E.resize(MeshSize_X + 1, vector<double>(MeshSize_Y + 2));
	u_Staggered_X_.u_W.resize(MeshSize_X + 1, vector<double>(MeshSize_Y + 2));
	u_Staggered_X_.u_N.resize(MeshSize_X + 1, vector<double>(MeshSize_Y + 2));
	u_Staggered_X_.u_S.resize(MeshSize_X + 1, vector<double>(MeshSize_Y + 2));
}

void xStaggered::u_Staggered_X_BC(double Re, double mu, double roh)
{
	for (int i = 0; i < MeshSize_X + 1; ++i) {
			u_Staggered_X_.u_P[i][MeshSize_Y + 1] = Re * mu / roh;
	}
}

void xStaggered::massFlowRate_StaggeredX_Setter(u_Staggered_X u_Staggered_X_, CellX CellX_, v_Staggered_Y v_Staggered_Y_, double roh)
{
	/*Ghost Cell Concept Would have been an Appropriate Implementation instead of Bifurcation*/
	for (int i = 0; i < MeshSize_X + 1; ++i) {
		for (int j = 0; j < MeshSize_Y; ++j) {
			if (i > 0) {
				massFlowRate_StaggeredX_.massFlowRate_West[i][j] = 0.5 * (u_Staggered_X_.u_P[i][j + 1] + u_Staggered_X_.u_P[i - 1][j + 1]) * CellX_.AreaWest[i][j] * roh;
			}
			else {
				// Handling B.C When i = 0
				//massFlowRate_StaggeredX_.massFlowRate_West[i][j] = 0.5 * (u_Staggered_X_.u_P[i][j + 1]) * CellX_.AreaWest[i][j] * roh;
				massFlowRate_StaggeredX_.massFlowRate_West[i][j] = (u_Staggered_X_.u_P[i][j + 1]) * CellX_.AreaWest[i][j] * roh; /* Newly Modified */
			}
			if (i < MeshSize_X) {
				massFlowRate_StaggeredX_.massFlowRate_East[i][j] = 0.5 * (u_Staggered_X_.u_P[i][j + 1] + u_Staggered_X_.u_P[i + 1][j + 1]) * CellX_.AreaEast[i][j] * roh;
			}
			else {
				// Handling B.C When i = MeshSize_X
				//massFlowRate_StaggeredX_.massFlowRate_East[i][j] = 0.5 * (u_Staggered_X_.u_P[i][j + 1]) * CellX_.AreaEast[i][j] * roh;
				massFlowRate_StaggeredX_.massFlowRate_East[i][j] = (u_Staggered_X_.u_P[i][j + 1]) * CellX_.AreaEast[i][j] * roh;
			}
		}
	}

	for (int i = 0; i < MeshSize_X + 1; ++i)
		for (int j = 0; j < MeshSize_Y; ++j)
		{
			if (i == 0) {
				// Handling B.C When i = 0
				massFlowRate_StaggeredX_.massFlowRate_North[i][j] = v_Staggered_Y_.v_P[i + 1][j + 1] * 0.5 * CellX_.AreaNorth[i][j] * roh;
				massFlowRate_StaggeredX_.massFlowRate_South[i][j] = v_Staggered_Y_.v_P[i + 1][j] * 0.5 * CellX_.AreaSouth[i][j] * roh;
			}
			else if (i == MeshSize_X) {
				// Handling B.C When i = MeshSize_X
				massFlowRate_StaggeredX_.massFlowRate_North[i][j] = v_Staggered_Y_.v_P[i][j + 1] * 0.5 * CellX_.AreaNorth[i][j] * roh;
				massFlowRate_StaggeredX_.massFlowRate_South[i][j] = v_Staggered_Y_.v_P[i][j] * 0.5 * CellX_.AreaSouth[i][j] * roh;
			}
			else if (i > 0 && i < MeshSize_X) {
				massFlowRate_StaggeredX_.massFlowRate_North[i][j] = (v_Staggered_Y_.v_P[i + 1][j + 1] * CellX_.AreaNorth[i][j] + v_Staggered_Y_.v_P[i][j + 1]* CellX_.AreaNorth[i][j]) * 0.5 * roh;
				massFlowRate_StaggeredX_.massFlowRate_South[i][j] = (v_Staggered_Y_.v_P[i + 1][j] * CellX_.AreaSouth[i][j] + v_Staggered_Y_.v_P[i][j] * CellX_.AreaSouth[i][j]) * 0.5 * roh;
			}
		}
}

// Function to print the Cell Distances: [dPE, dPW, dPN, dPS]
void xStaggered::checkMesh(CellX cell)
{
	cout << "\033[1;32m-+-+-+-+-+-+-+-+-+-+-+-+ checkMesh | Staggered-X-Mesh | [N + 1 * M + 2] -+-+-+-+-+-+-+-+-+-+-+-+\033[0m" << endl;
	cout << "\033[1;31mdPE\t" << "\033[1;31mdPW\t" << "\033[1;31mdPN\t" << "\033[1;31mdPS\033[0m" << endl;
	size_t maxSize = max({ cell.dPE.size(), cell.dPW.size(), cell.dPN.size(), cell.dPS.size() });
	for (size_t i = 0; i < maxSize; ++i) {
		// Print dPE if available
		if (i < cell.dPE.size()) {
			cout << cell.dPE[i] << "\t";
		}
		else {
			cout << "N/A\t"; // Print N/A if dPE doesn't have an element at index i
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




void xStaggered::u_Staggered_X_OPs(Cell Cell_Struct, vector<vector<double>> U_Pred, double Delta_t, double roh, vector<vector<double>> PressureField_New)
{

	for (int i = 1; i < MeshSize_X; ++i) {
		for (int j = 1; j < MeshSize_Y + 1; ++j) {	
			
			u_Staggered_X_.u_P[i][j] = U_Pred[i][j];
			u_Staggered_X_.u_P[i][j] -= (Delta_t / roh) * (PressureField_New[i + 1][j] - PressureField_New[i][j]) / Cell_Struct.dPE[i];
		}
	}


}

void xStaggered::printMassFlowRates() {
	// Print massFlowRate_East
	std::cout << "\033[1;32mxStaggeredMass Flow Rates:\n";
	std::cout << std::setw(15) << "East" << std::setw(15) << "West" << std::setw(15) << "North" << std::setw(15) << "South" << std::endl;
	for (int i = 0; i < MeshSize_X + 1; ++i) {
		for (int j = 0; j < MeshSize_Y ; ++j) {
			std::cout << std::setw(15) << massFlowRate_StaggeredX_.massFlowRate_East[i][j]
				<< std::setw(15) << massFlowRate_StaggeredX_.massFlowRate_West[i][j]
				<< std::setw(15) << massFlowRate_StaggeredX_.massFlowRate_North[i][j]
				<< std::setw(15) << massFlowRate_StaggeredX_.massFlowRate_South[i][j]
				<< std::endl;
		}
		std::cout << std::endl;
	}
}




xStaggered::xStaggered(float L, float H, float N, float M):L(L), H(H), MeshSize_X(N), MeshSize_Y(M) {}
xStaggered::~xStaggered() {}