#include <vector>
#include "mainMesh.h"

/*struct Pressure_Field
{
	vector<vector<double>> P_P;
	vector<vector<double>> P_E;
	vector<vector<double>> P_W;
	vector<vector<double>> P_N;
	vector<vector<double>> P_S;
};*/


struct Coeff_Vectors
{
	vector<vector<double>> aE;
	vector<vector<double>> aW;
	vector<vector<double>> aN;
	vector<vector<double>> aS;
	vector<vector<double>> aP;
	vector<vector<double>> bP;
};

class Pressure
{
private:
	float L, H, MeshSize_X, MeshSize_Y;
	Coeff_Vectors Coeff_Vectors_;
	vector<vector<double>> Pressure_Field_, Pressure_Field_New;
	vector<vector<double>> InitialField;

public:
	void PressureField_Initialization();
	void Coeff_Vectors_Initialization();
	void Coeff_Vectors_InternalNodes(Cell Cell_,Area Area, vector<vector<double>>& U_Pred, vector<vector<double>>& V_Pred, double Delta_t, double roh);
	void Coeff_Vectors_BCs();
	void LIN_EQ_Solver(Coeff_Vectors& Coeff_Vectors_, vector<vector<double>>& Pressure_Field_,  double Conv_Crit, int Itrmax);
	Coeff_Vectors& Coffs_Vector_Getter() { return Coeff_Vectors_; }
	//void InitialField_Setter();
	vector<vector<double>>& PressureField_Getter() { return Pressure_Field_; }
	vector<vector<double>>& PressureFieldNew_Getter() { return Pressure_Field_New; }


	Pressure(float L, float H, float N, float M);
	~Pressure();

};

