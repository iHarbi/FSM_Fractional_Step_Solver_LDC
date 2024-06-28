#include "Pressure.h"

void Pressure::PressureField_Initialization()
{
	//Pressure_Field_.P_P.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 2));
	//Pressure_Field_.P_E.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 2));
	//Pressure_Field_.P_W.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 2));
	//Pressure_Field_.P_N.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 2));
	//Pressure_Field_.P_S.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 2));
	Pressure_Field_.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 2));
	Pressure_Field_New.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 2));
}

void Pressure::Coeff_Vectors_Initialization()
{
	Coeff_Vectors_.aE.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 2));
	Coeff_Vectors_.aW.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 2));
	Coeff_Vectors_.aN.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 2));
	Coeff_Vectors_.aS.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 2));
	Coeff_Vectors_.aP.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 2));
	Coeff_Vectors_.bP.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 2));

}

/* Check Indicies --> 11/03/2024 */
void Pressure::Coeff_Vectors_InternalNodes(Cell Cell_, Area Area, vector<vector<double>>& U_Pred, vector<vector<double>>& V_Pred, double Delta_t, double roh)
{
	for (int i = 1; i < MeshSize_X + 1; ++i) {
		for (int j = 1; j < MeshSize_Y + 1; ++j) {
			/*Check the Const. Coeff*/
				Coeff_Vectors_.aE[i][j] = Area.East[i] / Cell_.dPE[i];
				Coeff_Vectors_.aW[i][j] = Area.West[i] / Cell_.dPW[i];
				Coeff_Vectors_.aN[i][j] = Area.North[j] / Cell_.dPN[j];
				Coeff_Vectors_.aS[i][j] = Area.South[j] / Cell_.dPS[j];
				Coeff_Vectors_.aP[i][j] = Coeff_Vectors_.aE[i][j] + Coeff_Vectors_.aW[i][j] + Coeff_Vectors_.aN[i][j] + Coeff_Vectors_.aS[i][j];
				Coeff_Vectors_.bP[i][j] = -roh * (1.0 / Delta_t);
				Coeff_Vectors_.bP[i][j] *= (U_Pred[i][j] * Area.East[i] - U_Pred[i - 1][j] * Area.West[i - 1] + V_Pred[i][j] * Area.North[j] - V_Pred[i][j - 1] * Area.South[j-1]);
		}
	}
}
void Pressure::Coeff_Vectors_BCs()
{
	/*Discretization Matrices | Boundary Nodes | LeftBoundary*/
	for (int j = 0; j < MeshSize_Y + 2; ++j) {
		Coeff_Vectors_.aE[0][j] = 1; // aE
		Coeff_Vectors_.aW[0][j] = 0; // aW
		Coeff_Vectors_.aN[0][j] = 0; // aN
		Coeff_Vectors_.aS[0][j] = 0; // aS
		Coeff_Vectors_.aP[0][j] = 1; // aP
		Coeff_Vectors_.bP[0][j] = 0; // bP
	}
	/*Discretization Matrices | Boundary Nodes | RightBoundary*/
	for (int j = 0; j < MeshSize_Y + 2; ++j) {
		Coeff_Vectors_.aE[MeshSize_X + 1][j] = 0; // aE
		Coeff_Vectors_.aW[MeshSize_X + 1][j] = 1; // aW
		Coeff_Vectors_.aN[MeshSize_X + 1][j] = 0; // aN
		Coeff_Vectors_.aS[MeshSize_X + 1][j] = 0; // aS
		Coeff_Vectors_.aP[MeshSize_X + 1][j] = 1; // aP
		Coeff_Vectors_.bP[MeshSize_X + 1][j] = 0; // bP
	}

	/*Discretization Matrices | Boundary Nodes | LowerWall*/
	for (int i = 0; i < MeshSize_X + 2; ++i) {
		Coeff_Vectors_.aE[i][0] = 0; // aE
		Coeff_Vectors_.aW[i][0] = 0; // aW
		Coeff_Vectors_.aN[i][0] = 1; // aN
		Coeff_Vectors_.aS[i][0] = 0; // aS
		Coeff_Vectors_.aP[i][0] = 1; // aP
		Coeff_Vectors_.bP[i][0] = 0; // bP
	}

	/* Discretization Matrices | Boundary Nodes | UpperWall*/
	for (int i = 0; i < MeshSize_X + 2; ++i) {
		Coeff_Vectors_.aE[i][MeshSize_Y + 1] = 0; // aE
		Coeff_Vectors_.aW[i][MeshSize_Y + 1] = 0; // aW
		Coeff_Vectors_.aN[i][MeshSize_Y + 1] = 0; // aN
		Coeff_Vectors_.aS[i][MeshSize_Y + 1] = 1; // aS
		Coeff_Vectors_.aP[i][MeshSize_Y + 1] = 1; // aP
		Coeff_Vectors_.bP[i][MeshSize_Y + 1] = 0; // bP 
	}
}

void Pressure::LIN_EQ_Solver(Coeff_Vectors& Coeff_Vectors_, vector<vector<double>>& Pressure_Field_, double Conv_Crit, int Itrmax)
{
	vector<vector<double>> phi_new = Pressure_Field_;
	vector<vector<double>> phi_estimated = Pressure_Field_;

	int Length_X = phi_new.size() - 2;
	int Length_Y = phi_new[0].size() - 2;
	
	double Delta = 1E3;
	int iter = 0;
	int print_interval = 5000;
	
	while (Delta >= Conv_Crit && iter < Itrmax) {

		/*Inner Cells*/
		for (int i = 1; i < Length_X + 1; ++i) {
			for (int j = 1; j < Length_Y + 1; ++j) {
				phi_new[i][j] = (Coeff_Vectors_.aE[i][j] * phi_new[i + 1][j] + Coeff_Vectors_.aW[i][j] * phi_new[i - 1][j] + Coeff_Vectors_.aN[i][j] * phi_new[i][j + 1] + Coeff_Vectors_.aS[i][j] * phi_new[i][j - 1] + Coeff_Vectors_.bP[i][j]) / Coeff_Vectors_.aP[i][j];
			}
		}
		/*UpperWall B.C.*/
		for (int i = 0; i <=  Length_X + 1; ++i) {
			//phi_new[i][Length_Y + 1] = (Coeff_Vectors_.aE[i][Length_Y + 1] * phi_new[i + 1][Length_Y + 1] + Coeff_Vectors_.aW[i][Length_Y + 1] * phi_new[i - 1][Length_Y + 1] + Coeff_Vectors_.aS[i][Length_Y + 1] * phi_new[i][Length_Y ] + Coeff_Vectors_.bP[i][Length_Y - 1]) / Coeff_Vectors_.aP[i][Length_Y + 1];
			phi_new[i][Length_Y + 1] = (Coeff_Vectors_.aS[i][Length_Y + 1] * phi_new[i][Length_Y]) / Coeff_Vectors_.aP[i][Length_Y + 1];

		}
		/*Right B.C.*/
		
		for (int j = 0; j <= Length_Y + 1; ++j) {
			//phi_new[Length_X + 1][j] = (Coeff_Vectors_.aW[Length_X + 1][j] * phi_new[Length_X][j] + Coeff_Vectors_.aN[Length_X + 1][j] * phi_new[Length_X + 1][j + 1] + Coeff_Vectors_.aS[Length_X + 1][j - 1] * phi_new[0][j - 1] + Coeff_Vectors_.bP[Length_X + 1][j]) / Coeff_Vectors_.aP[Length_X + 1][j];
			phi_new[Length_X + 1][j] = (Coeff_Vectors_.aW[Length_X + 1][j] * phi_new[Length_X][j]) / Coeff_Vectors_.aP[Length_X + 1][j];

		}
		/*Left B.C.*/
		for (int j = 0; j <= Length_Y + 1; ++j) {
			//phi_new[0][j] = (Coeff_Vectors_.aE[0][j] * phi_new[1][j] + Coeff_Vectors_.aN[0][j] * phi_new[0][j + 1] + Coeff_Vectors_.aS[0][j] * phi_new[0][j - 1] + Coeff_Vectors_.bP[0][j]) / Coeff_Vectors_.aP[0][j];
			phi_new[0][j] = (Coeff_Vectors_.aE[0][j] * phi_new[1][j]) / Coeff_Vectors_.aP[0][j];

		}

		/*Lower Patch | First Row*/
		for (int i = 0; i <= Length_X + 1; ++i) {
			//phi_new[i][0] = (Coeff_Vectors_.aE[i][0] * phi_new[i + 1][0] + Coeff_Vectors_.aW[i][0] * phi_new[i - 1][0] + Coeff_Vectors_.aN[i][0] * phi_new[i][1] + Coeff_Vectors_.bP[i][0]) / Coeff_Vectors_.aP[i][0];
			phi_new[i][0] = (Coeff_Vectors_.aN[i][0] * phi_new[i][1]) / Coeff_Vectors_.aP[i][0];

		}

		// Initialize max_diff with the smallest possible value for double
		double max_diff = numeric_limits<double>::lowest();

		// Iterate through each element and find the maximum absolute difference
		for (size_t i = 0; i < phi_new.size(); ++i) {
			for (size_t j = 0; j < phi_new[i].size(); ++j) {
				double diff = abs(phi_new[i][j] - phi_estimated[i][j]);
				if (diff > max_diff) {
					max_diff = diff;
				}
			}
		}

		Delta = max_diff;
		iter++;
		if (iter % print_interval == 0) {
			cout << "Iteration " << iter << ": delta = " << Delta << endl;
		}
		phi_estimated = phi_new;
	}
	cout << "Pressure Solver: GS Solver Convergence Criterion = " << Delta << endl;
	Pressure_Field_New = phi_new;

}

//void Pressure::InitialField_Setter(){}

Pressure::Pressure(float L, float H, float N, float M) : L(L), H(H), MeshSize_X(N), MeshSize_Y(M)
{
}

Pressure::~Pressure()
{

	//cout << "PressureField Instance was Created Successfully!";
}
