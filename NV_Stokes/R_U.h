#include <iostream>
#include <vector>
#include <string>



// Forward declarations
struct CellX;
struct massFlowRate_StaggeredX;
struct u_Staggered_X;
struct PHI_u;


using namespace std;
class R_U
{

private:
	vector<vector<double>> _R_U_; 	
	float L, H, MeshSize_X, MeshSize_Y;
public:
	R_U();
	~R_U();
	void _R_U_Initializer();
	void _R_U_OPS(CellX xStaggeredMesh,massFlowRate_StaggeredX massFlowRate_StaggeredX_, u_Staggered_X u_Staggered_X_, PHI_u u, double mu);
	vector<vector<double>>& _R_U_Getter() { return _R_U_; }
};

