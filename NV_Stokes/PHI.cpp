#include "PHI.h"

void PHI::PHI_Initializer()
{
	
	u.u_e.resize(MeshSize_X + 1, vector<double>(MeshSize_Y));
	u.u_w.resize(MeshSize_X + 1, vector<double>(MeshSize_Y));
	u.u_n.resize(MeshSize_X + 1, vector<double>(MeshSize_Y));
	u.u_s.resize(MeshSize_X + 1, vector<double>(MeshSize_Y));
	
	v.v_e.resize(MeshSize_X, vector<double>(MeshSize_Y + 1));
	v.v_w.resize(MeshSize_X, vector<double>(MeshSize_Y + 1));
	v.v_n.resize(MeshSize_X, vector<double>(MeshSize_Y + 1));
	v.v_s.resize(MeshSize_X, vector<double>(MeshSize_Y + 1));
}

void PHI::CDS_Scheme(u_Staggered_X u_Staggered_X, v_Staggered_Y v_Staggered_Y)
{
	/* CDS Implmentation Only till Now | For UDS: Create another Method that takes both massFlowRates as Inputs too!*/
	/*-+-+-+-+-+-+-+-+-+-+-+-+-+ <u> Internal -+-+-+-+-+-+-+-+-+-+-+-+-+*/
	for (int i = 1; i < MeshSize_X; ++i)
		for (int j = 1; j < MeshSize_Y - 1; ++j)
		{
			u.u_e[i][j] = 0.5 *  (u_Staggered_X.u_P[i][j + 1] + u_Staggered_X.u_P[i + 1][j + 1] );
			u.u_n[i][j] = 0.5 *  (u_Staggered_X.u_P[i][j + 1] + u_Staggered_X.u_P[i][j + 2] );
			u.u_w[i][j] = 0.5 *  (u_Staggered_X.u_P[i][j + 1] + u_Staggered_X.u_P[i - 1][j + 1] );
			u.u_s[i][j] = 0.5 *  (u_Staggered_X.u_P[i][j + 1] + u_Staggered_X.u_P[i][j] );
		}
	/*-+-+-+-+-+-+-+-+-+-+-+-+-+ <v> Internal -+-+-+-+-+-+-+-+-+-+-+-+-+*/
	for (int i = 1; i < MeshSize_X - 1; ++i)
		for (int j = 1; j < MeshSize_Y ; ++j)
		{
			v.v_e[i][j] = 0.5 * (v_Staggered_Y.v_P[i + 1][j] + v_Staggered_Y.v_P[i + 2][j]);
			v.v_n[i][j] = 0.5 * (v_Staggered_Y.v_P[i + 1][j] + v_Staggered_Y.v_P[i + 1][j + 1]);
			v.v_w[i][j] = 0.5 * (v_Staggered_Y.v_P[i + 1][j] + v_Staggered_Y.v_P[i][j]);
			v.v_s[i][j] = 0.5 * (v_Staggered_Y.v_P[i + 1][j] + v_Staggered_Y.v_P[i + 1][j - 1]);
		}
}


void PHI::PHI_Setter_BC(float U, u_Staggered_X u_Staggered_X, v_Staggered_Y  v_Staggered_Y)
{
	/*-+-+-+-+-+-+-+-+-+-+-+-+-+ <u> B.Cs -+-+-+-+-+-+-+-+-+-+-+-+-+*/
	/* No Slip B.C -- > lowerPatch */
	for (int i = 0; i <= MeshSize_X; ++i)
	{
		u.u_n[i][0] = 0.5 * (u_Staggered_X.u_P[i][1] + u_Staggered_X.u_P[i][2]);
		u.u_s[i][0] = 0.5 * (u_Staggered_X.u_P[i][1] + u_Staggered_X.u_P[i][0]); /* Was Set to ZERO*/
		
		if (i == 0){
			//u.u_w[i][0] = 0.5 * (u_Staggered_X.u_P[i][1]);
			u.u_w[i][0] = 0; /*Newly Set | Conrner CELL*/
		}
		else {
			u.u_w[i][0] = 0.5 * (u_Staggered_X.u_P[i][1] + u_Staggered_X.u_P[i - 1][1]);
		}
		
		if (i == MeshSize_X) {
			//u.u_e[i][0] = 0.5 * (u_Staggered_X.u_P[i][0]);
			u.u_e[i][0] = 0; /*Newly Set | Conrner CELL*/
		}
		else {
			u.u_e[i][0] = 0.5 * (u_Staggered_X.u_P[i][1] + u_Staggered_X.u_P[i + 1][1]);
		}
	}
		
	/* Lid Velocity B.C -- > upperPatch */
	for (int i = 0; i <= MeshSize_X; ++i)
	{
		u.u_n[i][MeshSize_Y - 1] = U;
		u.u_s[i][MeshSize_Y - 1] = 0.5 * (u_Staggered_X.u_P[i][MeshSize_Y] + u_Staggered_X.u_P[i][MeshSize_Y - 1]);
		
		if (i == 0) {
			//u.u_w[i][MeshSize_Y - 1] = 0.5 * (u_Staggered_X.u_P[i][MeshSize_Y + 1]);
			u.u_w[i][MeshSize_Y - 1] = 0; /*Newly Set | Corner CELL*/
		}
		else {
			u.u_w[i][MeshSize_Y - 1] = 0.5 * (u_Staggered_X.u_P[i][MeshSize_Y ] + u_Staggered_X.u_P[i - 1][MeshSize_Y ]); /* Newly Modified [J] Index*/
		}

		if (i == MeshSize_X) {
			//u.u_e[i][MeshSize_Y - 1] = 0.5 * (u_Staggered_X.u_P[i][MeshSize_Y + 1]);
			u.u_e[i][MeshSize_Y - 1] = 0; /*Newly Set | Corner CELL*/
		}
		else {
			u.u_e[i][MeshSize_Y - 1] = 0.5 * (u_Staggered_X.u_P[i][MeshSize_Y] + u_Staggered_X.u_P[i + 1][MeshSize_Y]); /* Newly Modified [J] Index*/
		}
	}
	/* No Slip B.C -- > leftPatch */
	for (int j = 0; j < MeshSize_Y; ++j)
	{
		u.u_w[0][j] = 0.0;
		u.u_e[0][j] = 0.5 * (u_Staggered_X.u_P[0][j + 1] + u_Staggered_X.u_P[1][j + 1]);

		if (j == 0) {
			//u.u_s[0][j] = 0.5 * (u_Staggered_X.u_P[0][1]);
			u.u_s[0][j] = 0; /*Newly Set | CORNER CELL*/
		}
		else {
			u.u_s[0][j] = 0.5 * (u_Staggered_X.u_P[0][j + 1] + u_Staggered_X.u_P[0][j]);
		}

		if (j == MeshSize_Y - 1) {
			//u.u_n[0][j] = 0.5 * (u_Staggered_X.u_P[0][MeshSize_Y + 1]);
			u.u_n[0][j] = U; /*Newly Set | CORNER CELL*/
		}
		else {
			u.u_n[0][j] = 0.5 * (u_Staggered_X.u_P[0][j + 1] + u_Staggered_X.u_P[0][j + 2]);
		}
	
	}
	/* No Slip B.C -- > rightPatch */
	for (int j = 0; j < MeshSize_Y; ++j)
	{
		u.u_e[MeshSize_X][j] = 0.0;
		u.u_w[MeshSize_X][j] = 0.5 * (u_Staggered_X.u_P[MeshSize_X][j + 1] + u_Staggered_X.u_P[MeshSize_X - 1][j + 1]);

		if (j == 0) {
			//u.u_s[MeshSize_X][j] = 0.5 * (u_Staggered_X.u_P[MeshSize_X][j + 1]);
			u.u_s[MeshSize_X][j] = 0; /* Newly Set  | Corner CELL*/
		}
		else {
			u.u_s[MeshSize_X][j] = 0.5 * (u_Staggered_X.u_P[MeshSize_X][j + 1] + u_Staggered_X.u_P[MeshSize_X][j]);
		}

		if (j == MeshSize_Y - 1) {
			//u.u_n[MeshSize_X][j] = 0.5 * (u_Staggered_X.u_P[MeshSize_X][MeshSize_Y + 1]);
			u.u_n[MeshSize_X][j] = U; /* Newly Set | Corner CELL*/
		}
		else {
			u.u_n[MeshSize_X][j] = 0.5 * (u_Staggered_X.u_P[MeshSize_X][j + 1] + u_Staggered_X.u_P[MeshSize_X][j + 2]);
		}
	}
	



	
	/*-+-+-+-+-+-+-+-+-+-+-+-+-+ <v> B.Cs-+-+-+-+-+-+-+-+-+-+-+-+-+*/
	/* No Slip B.C -- > lowerPatch */
	for (int i = 0; i < MeshSize_X; ++i)
	{
		v.v_s[i][0] = 0.0;
		v.v_n[i][0] = 0.5 * (v_Staggered_Y.v_P[i + 1][0] + v_Staggered_Y.v_P[i + 1][1]);

		if (i == 0) {
			//v.v_w[i][0] = 0.5 * (v_Staggered_Y.v_P[i + 1][0]);
			v.v_w[i][0] = 0; /* Newly Set | Corner CELL*/
		}
		else {
			v.v_w[i][0] = 0.5 * (v_Staggered_Y.v_P[i + 1][0] + v_Staggered_Y.v_P[i][0]);
		}

		if (i == MeshSize_X - 1) {
			//v.v_e[i][0] = 0.5 * (v_Staggered_Y.v_P[i + 1][0]);
			v.v_e[i][0] = 0;/* Newly Set | Corner CELL*/
		}
		else {
			v.v_e[i][0] = 0.5 * (v_Staggered_Y.v_P[i + 1][0] + v_Staggered_Y.v_P[i + 2][0]);
		}



	}
	/* Lid Velocity B.C -- > UpperPatch */
	for (int i = 0; i < MeshSize_X; ++i)
	{
		v.v_n[i][MeshSize_Y] = 0.0;
		v.v_s[i][MeshSize_Y] = 0.5 * (v_Staggered_Y.v_P[i + 1][MeshSize_Y] + v_Staggered_Y.v_P[i + 1][MeshSize_Y -  1]);

		if (i == 0) {
			v.v_w[i][MeshSize_Y] = 0.5 * (v_Staggered_Y.v_P[i + 1][MeshSize_Y]);
			v.v_w[i][MeshSize_Y] = 0; /* Newly Set | Corner CELL */
		}
		else {
			v.v_w[i][MeshSize_Y] = 0.5 * (v_Staggered_Y.v_P[i + 1][MeshSize_Y] + v_Staggered_Y.v_P[i][MeshSize_Y]);
		}

		if (i == MeshSize_X - 1){
			//v.v_e[i][MeshSize_Y] = 0.5 * (v_Staggered_Y.v_P[i + 1][MeshSize_Y]);
			v.v_e[i][MeshSize_Y] = 0; /* Newly Set | Corner CELL */
		}
		else {
			v.v_e[i][MeshSize_Y] = 0.5 * (v_Staggered_Y.v_P[i + 1][MeshSize_Y] + v_Staggered_Y.v_P[i + 2][MeshSize_Y]);
		}
	}

	/* No Slip B.C --> leftPatch */
	for (int j = 0; j <= MeshSize_Y; ++j)
	{ 
		v.v_w[0][j] = 0.0;
		v.v_e[0][j] = 0.5 * (v_Staggered_Y.v_P[1][j] + v_Staggered_Y.v_P[2][j]);

		if (j == 0)
		{
			v.v_s[0][j] = 0.5 * (v_Staggered_Y.v_P[1][j]);
			v.v_s[0][j] = 0;  /* Newly Set | Corner CELL */
		}
		else {
			v.v_s[0][j] = 0.5 * (v_Staggered_Y.v_P[1][j] + v_Staggered_Y.v_P[1][j - 1]);
		}

		if (j == MeshSize_Y)
		{
			v.v_n[0][j] = 0.5 * (v_Staggered_Y.v_P[1][j]);
		}
		else {
			v.v_n[0][j] = 0.5 * (v_Staggered_Y.v_P[1][j] + v_Staggered_Y.v_P[1][j + 1]);
		}
	}
	/* No Slip B.C --> rightPatch */
	for (int j = 0; j <= MeshSize_Y; ++j)
	{
		v.v_e[MeshSize_X - 1][j] = 0.0;
		v.v_w[MeshSize_X - 1][j] = 0.5 * (v_Staggered_Y.v_P[MeshSize_X - 1][j] + v_Staggered_Y.v_P[MeshSize_X][j]);

		if (j == 0)
		{
			//v.v_s[MeshSize_X - 1][j] = 0.5 * (v_Staggered_Y.v_P[MeshSize_X - 1][j]);
			v.v_s[MeshSize_X - 1][j] = 0; /* Newly Set | Corner CELL */
		}
		else {
			v.v_s[MeshSize_X - 1][j] = 0.5 * (v_Staggered_Y.v_P[MeshSize_X][j] + v_Staggered_Y.v_P[MeshSize_X ][j - 1]);
		}

		if (j == MeshSize_Y)
		{
			//v.v_n[MeshSize_X - 1][j] = 0.5 * (v_Staggered_Y.v_P[MeshSize_X - 1][j]);
			v.v_n[MeshSize_X - 1][j] = 0; /* Newly Set | Corner CELL */
		}
		else {
			v.v_n[MeshSize_X - 1][j] = 0.5 * (v_Staggered_Y.v_P[MeshSize_X][j] + v_Staggered_Y.v_P[MeshSize_X][j + 1]);
		}
	}
}



void PHI::printAndCheck() {
	// Print and check the result in a tabular format
	std::cout << std::setw(15) << "u_e" << std::setw(15) << "u_w" << std::setw(15) << "u_n" << std::setw(15) << "u_s" << std::endl;
	for (int j = 0; j < MeshSize_Y; ++j) {
		for (int i = 0; i <= MeshSize_X; ++i) {
			std::cout << std::setw(15) << u.u_e[i][j]
				<< std::setw(15) << u.u_w[i][j]
				<< std::setw(15) << u.u_n[i][j]
				<< std::setw(15) << u.u_s[i][j]
				<< std::endl;
		}
		std::cout << std::endl;
	}
}

void PHI::HRS(u_Staggered_X u_Staggered_X, v_Staggered_Y v_Staggered_Y, massFlowRate_StaggeredX massFlowRate_StaggeredX_, massFlowRate_StaggeredY massFlowRate_StaggeredY_, CellX xStaggeredMesh, CellY yStaggeredMesh)
{
	HRS_Schemes HRS;

	/*-+-+-+-+-+-+-+-+-+-+-+-+-+ <u> Internal -+-+-+-+-+-+-+-+-+-+-+-+-+*/
	for (int i = 1; i < MeshSize_X; ++i)
		for (int j = 1; j < MeshSize_Y - 1; ++j)
			/*Not Completed: Initial Condition might Cuase Numerical Isssues, causing instabilities*/
		{	
			/* I/P(m_f, x_f, x_P, x_E, x_W, x_EE, phi_P, phi_E,phi_W, phi_EE, "Scheme", "Face" )*/

			if (massFlowRate_StaggeredX_.massFlowRate_East[i][j] >= 0) {
				double u_e = HRS.computeConvectiveValue(massFlowRate_StaggeredX_.massFlowRate_East[i][j],  0.5 * (xStaggeredMesh.CellCentroid_X_Staggered[i] + xStaggeredMesh.CellCentroid_X_Staggered[i + 1]), xStaggeredMesh.CellCentroid_X_Staggered[i], xStaggeredMesh.CellCentroid_X_Staggered[i + 1], xStaggeredMesh.CellCentroid_X_Staggered[i - 1], xStaggeredMesh.CellCentroid_X_Staggered[i + 2], u_Staggered_X.u_P[i][j + 1], u_Staggered_X.u_P[i + 1][j + 1], u_Staggered_X.u_P[i - 1][j + 1], u_Staggered_X.u_P[i + 2][j + 1], "UDS", "East");
			}																							  
			else {																						  
				double u_e = HRS.computeConvectiveValue(massFlowRate_StaggeredX_.massFlowRate_East[i][j],  0.5 * (xStaggeredMesh.CellCentroid_X_Staggered[i] + xStaggeredMesh.CellCentroid_X_Staggered[i + 1]), xStaggeredMesh.CellCentroid_X_Staggered[i], xStaggeredMesh.CellCentroid_X_Staggered[i + 1], xStaggeredMesh.CellCentroid_X_Staggered[i - 1], xStaggeredMesh.CellCentroid_X_Staggered[i + 2], u_Staggered_X.u_P[i][j + 1], u_Staggered_X.u_P[i + 1][j + 1], u_Staggered_X.u_P[i - 1][j + 1], u_Staggered_X.u_P[i + 2][j + 1], "UDS", "East");
			}
			if (massFlowRate_StaggeredX_.massFlowRate_North[i][j] >= 0) {
				double u_n = HRS.computeConvectiveValue(massFlowRate_StaggeredX_.massFlowRate_North[i][j], 0.5 * (yStaggeredMesh.CellCentroid_Y_Staggered[j] + yStaggeredMesh.CellCentroid_Y_Staggered[j + 1]), yStaggeredMesh.CellCentroid_Y_Staggered[j], yStaggeredMesh.CellCentroid_Y_Staggered[j + 1], yStaggeredMesh.CellCentroid_Y_Staggered[j - 1], 0, u_Staggered_X.u_P[i][j + 1], u_Staggered_X.u_P[i][j + 2], u_Staggered_X.u_P[i][j], 0, "UDS", "North");
			}
			else {
				double u_n = HRS.computeConvectiveValue(massFlowRate_StaggeredX_.massFlowRate_North[i][j], 0.5 * (yStaggeredMesh.CellCentroid_Y_Staggered[j] + yStaggeredMesh.CellCentroid_Y_Staggered[j + 1]), yStaggeredMesh.CellCentroid_Y_Staggered[j], yStaggeredMesh.CellCentroid_Y_Staggered[j + 1], 0, yStaggeredMesh.CellCentroid_Y_Staggered[j + 2], u_Staggered_X.u_P[i][j + 1], u_Staggered_X.u_P[i][j + 2], 0, u_Staggered_X.u_P[i][j + 3], "UDS", "North");
			}
			if (massFlowRate_StaggeredX_.massFlowRate_West[i][j] >= 0) {
				double u_w = HRS.computeConvectiveValue(massFlowRate_StaggeredX_.massFlowRate_West[i][j],  0.5 * (xStaggeredMesh.CellCentroid_X_Staggered[i] + xStaggeredMesh.CellCentroid_X_Staggered[i - 1]), xStaggeredMesh.CellCentroid_X_Staggered[i], 0, xStaggeredMesh.CellCentroid_X_Staggered[i - 1], xStaggeredMesh.CellCentroid_X_Staggered[i - 2], u_Staggered_X.u_P[i][j + 1],0 , u_Staggered_X.u_P[i - 1][j + 1], u_Staggered_X.u_P[i - 2][j + 1], "UDS", "West");
			}																							   
			else {																						   
				double u_w = HRS.computeConvectiveValue(massFlowRate_StaggeredX_.massFlowRate_West[i][j],  0.5 * (xStaggeredMesh.CellCentroid_X_Staggered[i] + xStaggeredMesh.CellCentroid_X_Staggered[i - 1]), xStaggeredMesh.CellCentroid_X_Staggered[i], xStaggeredMesh.CellCentroid_X_Staggered[i + 1], xStaggeredMesh.CellCentroid_X_Staggered[i - 1], 0, u_Staggered_X.u_P[i][j + 1], 0, u_Staggered_X.u_P[i - 1][j + 1], 0, "UDS", "West");
			}
			if (massFlowRate_StaggeredX_.massFlowRate_South[i][j] >= 0) {
				double u_s = HRS.computeConvectiveValue(massFlowRate_StaggeredX_.massFlowRate_South[i][j], 0.5 * (yStaggeredMesh.CellCentroid_Y_Staggered[j] + yStaggeredMesh.CellCentroid_Y_Staggered[j - 1]), yStaggeredMesh.CellCentroid_Y_Staggered[j], 0, yStaggeredMesh.CellCentroid_Y_Staggered[j - 1], yStaggeredMesh.CellCentroid_Y_Staggered[j - 2], u_Staggered_X.u_P[i][j + 1], 0, u_Staggered_X.u_P[i][j], u_Staggered_X.u_P[i][j - 1], "UDS", "South");
			}
			else {
				double u_s = HRS.computeConvectiveValue(massFlowRate_StaggeredX_.massFlowRate_South[i][j], 0.5 * (yStaggeredMesh.CellCentroid_Y_Staggered[j] + yStaggeredMesh.CellCentroid_Y_Staggered[j - 1]), yStaggeredMesh.CellCentroid_Y_Staggered[j], yStaggeredMesh.CellCentroid_Y_Staggered[j + 1], yStaggeredMesh.CellCentroid_Y_Staggered[j - 1], 0, u_Staggered_X.u_P[i][j + 1], u_Staggered_X.u_P[i][j + 2], u_Staggered_X.u_P[i][j], 0, "UDS", "South");
			}

		}
	/*-+-+-+-+-+-+-+-+-+-+-+-+-+ <v> Internal -+-+-+-+-+-+-+-+-+-+-+-+-+*/
	for (int i = 1; i < MeshSize_X - 1; ++i)
		for (int j = 1; j < MeshSize_Y; ++j)
		{

		}



}

void PHI::Upwind_Scheme(u_Staggered_X u_Staggered_X, v_Staggered_Y v_Staggered_Y, massFlowRate_StaggeredX massFlowRate_StaggeredX_, massFlowRate_StaggeredY massFlowRate_StaggeredY_)
{

	/*-+-+-+-+-+-+-+-+-+-+-+-+-+ <u> Internal -+-+-+-+-+-+-+-+-+-+-+-+-+*/
	for (int i = 1; i < MeshSize_X; ++i)
		for (int j = 1; j < MeshSize_Y - 1; ++j)
		{

			if (massFlowRate_StaggeredX_.massFlowRate_East[i][j] >= 0) {u.u_e[i][j] = u_Staggered_X.u_P[i][j + 1];}
			else { u.u_e[i][j] = u_Staggered_X.u_P[i + 1][j + 1];}
			if (massFlowRate_StaggeredX_.massFlowRate_West[i][j] <= 0) { u.u_w[i][j] = u_Staggered_X.u_P[i][j + 1]; }
			else { u.u_w[i][j] = u_Staggered_X.u_P[i - 1][j + 1]; }
			if (massFlowRate_StaggeredX_.massFlowRate_North[i][j] >= 0) { u.u_n[i][j] = u_Staggered_X.u_P[i][j + 1]; }
			else { u.u_n[i][j] = u_Staggered_X.u_P[i][j + 2]; }
			if (massFlowRate_StaggeredX_.massFlowRate_South[i][j] <= 0) { u.u_s[i][j] = u_Staggered_X.u_P[i][j + 1]; }
			else { u.u_s[i][j] = u_Staggered_X.u_P[i][j]; }

		}
	/*-+-+-+-+-+-+-+-+-+-+-+-+-+ <v> Internal -+-+-+-+-+-+-+-+-+-+-+-+-+*/
	for (int i = 1; i < MeshSize_X - 1; ++i)
		for (int j = 1; j < MeshSize_Y; ++j)
		{

			if (massFlowRate_StaggeredY_.massFlowRate_East[i][j] >= 0)  { v.v_e[i][j] = v_Staggered_Y.v_P[i + 1][j]; }
			else { v.v_e[i][j] = v_Staggered_Y.v_P[i + 2][j]; }		    
			if (massFlowRate_StaggeredY_.massFlowRate_West[i][j] <= 0) { v.v_w[i][j] = v_Staggered_Y.v_P[i + 1][j]; }
			else { v.v_w[i][j] = v_Staggered_Y.v_P[i][j]; }
			if (massFlowRate_StaggeredY_.massFlowRate_North[i][j] >= 0) { v.v_n[i][j] = v_Staggered_Y.v_P[i + 1][j]; }
			else { v.v_n[i][j] = v_Staggered_Y.v_P[i + 1][j + 1 ]; }
			if (massFlowRate_StaggeredY_.massFlowRate_South[i][j] <= 0) { v.v_s[i][j] = v_Staggered_Y.v_P[i + 1][j]; }
			else { v.v_s[i][j] = v_Staggered_Y.v_P[i + 1][j - 1]; }

		}

}






PHI::PHI(float L, float H, float N, float M): L(L), H(H), MeshSize_X(N), MeshSize_Y(M) {}
PHI::~PHI() {}

void PHI::QUICK(u_Staggered_X u_Staggered_X, v_Staggered_Y v_Staggered_Y, massFlowRate_StaggeredX massFlowRate_StaggeredX_, massFlowRate_StaggeredY massFlowRate_StaggeredY_)
{
	
	/*-+-+-+-+-+-+-+-+-+-+-+-+-+ <u> Internal -+-+-+-+-+-+-+-+-+-+-+-+-+*/
	for (int i = 2; i < MeshSize_X - 1; ++i)
		for (int j = 1; j < MeshSize_Y - 1; ++j)
		{

			if (massFlowRate_StaggeredX_.massFlowRate_East[i][j] >= 0) { u.u_e[i][j] = 0.75 * u_Staggered_X.u_P[i][j + 1] - 0.125 * u_Staggered_X.u_P[i - 1][j + 1] + 0.375 * u_Staggered_X.u_P[i + 1][j + 1]; }
			else { u.u_e[i][j] = 0.75 * u_Staggered_X.u_P[i + 1][j + 1] - 0.125 * u_Staggered_X.u_P[i + 2][j + 1] + 0.375 * u_Staggered_X.u_P[i][j + 1]; }
			if (massFlowRate_StaggeredX_.massFlowRate_West[i][j] <= 0) { u.u_w[i][j] = 0.75 * u_Staggered_X.u_P[i][j + 1] - 0.125 * u_Staggered_X.u_P[i + 1][j + 1] + 0.375 * u_Staggered_X.u_P[i - 1][j + 1]; }
			else { u.u_w[i][j] = 0.75 * u_Staggered_X.u_P[i - 1][j + 1] - 0.125 * u_Staggered_X.u_P[i - 2][j + 1] + 0.375 * u_Staggered_X.u_P[i][j + 1]; }
			if (massFlowRate_StaggeredX_.massFlowRate_North[i][j] >= 0) { u.u_n[i][j] = 0.75 * u_Staggered_X.u_P[i][j + 1] - 0.125 * u_Staggered_X.u_P[i][j] + 0.375 * u_Staggered_X.u_P[i][j + 2]; }
			else { u.u_n[i][j] = 0.75 * u_Staggered_X.u_P[i][j + 2] - 0.125 * u_Staggered_X.u_P[i][j + 3] + 0.375 * u_Staggered_X.u_P[i][j + 1]; }
			if (massFlowRate_StaggeredX_.massFlowRate_South[i][j] <= 0) { u.u_s[i][j] = 0.75 * u_Staggered_X.u_P[i][j + 1] - 0.125 * u_Staggered_X.u_P[i][j + 2] + 0.375 * u_Staggered_X.u_P[i][j]; }
			else { u.u_s[i][j] = 0.75 * u_Staggered_X.u_P[i][j] - 0.125 * u_Staggered_X.u_P[i][j - 1] + 0.375 * u_Staggered_X.u_P[i][j + 1]; }

		}
	/*-+-+-+-+-+-+-+-+-+-+-+-+-+ <v> Internal -+-+-+-+-+-+-+-+-+-+-+-+-+*/
	for (int i = 1; i < MeshSize_X - 1; ++i)
		for (int j = 2; j < MeshSize_Y - 1; ++j)
		{

			if (massFlowRate_StaggeredY_.massFlowRate_East[i][j] >= 0) { v.v_e[i][j] = 0.75 * v_Staggered_Y.v_P[i + 1][j] - 0.125 * v_Staggered_Y.v_P[i][j] + 0.375 * v_Staggered_Y.v_P[i + 2][j]; }
			else { v.v_e[i][j] = 0.75 * v_Staggered_Y.v_P[i + 2][j] - 0.125 * v_Staggered_Y.v_P[i + 3][j] + 0.375 * v_Staggered_Y.v_P[i + 1][j]; }
			if (massFlowRate_StaggeredY_.massFlowRate_West[i][j] <= 0) { v.v_w[i][j] = 0.75 * v_Staggered_Y.v_P[i + 1][j] - 0.125 * v_Staggered_Y.v_P[i + 2][j] + 0.375 * v_Staggered_Y.v_P[i][j]; }
			else { v.v_w[i][j] = 0.75 * v_Staggered_Y.v_P[i][j] - 0.125 * v_Staggered_Y.v_P[i - 1][j] + 0.375 * v_Staggered_Y.v_P[i + 1][j]; }
			if (massFlowRate_StaggeredY_.massFlowRate_North[i][j] >= 0) { v.v_n[i][j] = 0.75 * v_Staggered_Y.v_P[i + 1][j] - 0.125 * v_Staggered_Y.v_P[i + 1][j - 1] + 0.375 * v_Staggered_Y.v_P[i + 1][j + 1]; }
			else { v.v_n[i][j] = 0.75 * v_Staggered_Y.v_P[i + 1][j + 1] - 0.125 * v_Staggered_Y.v_P[i + 1][j + 2] + 0.375 * v_Staggered_Y.v_P[i + 1][j]; }
			if (massFlowRate_StaggeredY_.massFlowRate_South[i][j] <= 0) { v.v_s[i][j] = 0.75 * v_Staggered_Y.v_P[i + 1][j] - 0.125 * v_Staggered_Y.v_P[i + 1][j + 1] + 0.375 * v_Staggered_Y.v_P[i + 1][j - 1]; }
			else { v.v_s[i][j] = 0.75 * v_Staggered_Y.v_P[i + 1][j - 1] - 0.125 * v_Staggered_Y.v_P[i + 1][j - 2] + 0.375 * v_Staggered_Y.v_P[i + 1][j]; }

		}
}

void PHI::SOU(u_Staggered_X u_Staggered_X, v_Staggered_Y v_Staggered_Y, massFlowRate_StaggeredX massFlowRate_StaggeredX_, massFlowRate_StaggeredY massFlowRate_StaggeredY_)
{
	/*-+-+-+-+-+-+-+-+-+-+-+-+-+ <u> Internal -+-+-+-+-+-+-+-+-+-+-+-+-+*/
	for (int i = 2; i < MeshSize_X - 1; ++i)
		for (int j = 1; j < MeshSize_Y - 1; ++j)
		{

			if (massFlowRate_StaggeredX_.massFlowRate_East[i][j] >= 0) { u.u_e[i][j] = 1.5 * u_Staggered_X.u_P[i][j + 1] - 0.5 * u_Staggered_X.u_P[i - 1][j + 1] ; }
			else { u.u_e[i][j] = 1.5 * u_Staggered_X.u_P[i + 1][j + 1] - 0.5 * u_Staggered_X.u_P[i + 2][j + 1] ; }
			if (massFlowRate_StaggeredX_.massFlowRate_West[i][j] <= 0) { u.u_w[i][j] = 1.5 * u_Staggered_X.u_P[i][j + 1] - 0.5 * u_Staggered_X.u_P[i + 1][j + 1] ; }
			else { u.u_w[i][j] = 1.5 * u_Staggered_X.u_P[i - 1][j + 1] - 0.5 * u_Staggered_X.u_P[i - 2][j + 1] ; }
			if (massFlowRate_StaggeredX_.massFlowRate_North[i][j] >= 0) { u.u_n[i][j] = 1.5 * u_Staggered_X.u_P[i][j + 1] - 0.5 * u_Staggered_X.u_P[i][j] ; }
			else { u.u_n[i][j] = 1.5 * u_Staggered_X.u_P[i][j + 2] - 0.5 * u_Staggered_X.u_P[i][j + 3] ; }
			if (massFlowRate_StaggeredX_.massFlowRate_South[i][j] <= 0) { u.u_s[i][j] = 1.5 * u_Staggered_X.u_P[i][j + 1] - 0.5 * u_Staggered_X.u_P[i][j + 2] ; }
			else { u.u_s[i][j] = 1.5 * u_Staggered_X.u_P[i][j] - 0.5 * u_Staggered_X.u_P[i][j - 1] ; }
		}
	/*-+-+-+-+-+-+-+-+-+-+-+-+-+ <v> Internal -+-+-+-+-+-+-+-+-+-+-+-+-+*/
	for (int i = 1; i < MeshSize_X - 1; ++i)
		for (int j = 2; j < MeshSize_Y - 1; ++j)
		{

			if (massFlowRate_StaggeredY_.massFlowRate_East[i][j] >= 0) { v.v_e[i][j] = 1.5 * v_Staggered_Y.v_P[i + 1][j] - 0.5 * v_Staggered_Y.v_P[i][j] ; }
			else { v.v_e[i][j] = 1.5 * v_Staggered_Y.v_P[i + 2][j] - 0.5 * v_Staggered_Y.v_P[i + 3][j] ; }
			if (massFlowRate_StaggeredY_.massFlowRate_West[i][j] <= 0) { v.v_w[i][j] = 1.5 * v_Staggered_Y.v_P[i + 1][j] - 0.5 * v_Staggered_Y.v_P[i + 2][j] ; }
			else { v.v_w[i][j] = 1.5 * v_Staggered_Y.v_P[i][j] - 0.5 * v_Staggered_Y.v_P[i - 1][j] ; }
			if (massFlowRate_StaggeredY_.massFlowRate_North[i][j] >= 0) { v.v_n[i][j] = 1.5 * v_Staggered_Y.v_P[i + 1][j] - 0.5 * v_Staggered_Y.v_P[i + 1][j - 1] ; }
			else { v.v_n[i][j] = 1.5 * v_Staggered_Y.v_P[i + 1][j + 1] - 0.5 * v_Staggered_Y.v_P[i + 1][j + 2] ; }
			if (massFlowRate_StaggeredY_.massFlowRate_South[i][j] <= 0) { v.v_s[i][j] = 1.5 * v_Staggered_Y.v_P[i + 1][j] - 0.5 * v_Staggered_Y.v_P[i + 1][j + 1] ; }
			else { v.v_s[i][j] = 1.5 * v_Staggered_Y.v_P[i + 1][j - 1] - 0.5 * v_Staggered_Y.v_P[i + 1][j - 2]; }

		}
}
