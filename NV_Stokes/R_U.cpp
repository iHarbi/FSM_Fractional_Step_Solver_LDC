#include "R_U.h"

R_U::R_U()
{}
R_U::~R_U()
{}

void R_U::_R_U_Initializer()
{
	_R_U_.resize(MeshSize_X + 1, vector<double> (MeshSize_Y) );

}

void R_U::_R_U_OPS(CellX xStaggeredMesh,massFlowRate_StaggeredX massFlowRate_StaggeredX_, u_Staggered_X u_Staggered_X_, PHI_u u, double mu)
{


	for (int i = 1; i < MeshSize_X; ++i)
		for (int j = 0; j < MeshSize_Y; ++j)
		{
			{
			_R_U_[i][j] = -(massFlowRate_StaggeredX_.massFlowRate_East[i][j] * u.u_e[i][j] - massFlowRate_StaggeredX_.massFlowRate_West[i][j] * u.u_w[i][j] + massFlowRate_StaggeredX_.massFlowRate_North[i][j] * u.u_n[i][j] - massFlowRate_StaggeredX_.massFlowRate_South[i][j] * u.u_s[i][j])
				+ mu * (u_Staggered_X_.u_E[i][j] - u_Staggered_X_.u_P[i][j]) * xStaggeredMesh.D_East[i][j]
				- mu * (u_Staggered_X_.u_P[i][j] - u_Staggered_X_.u_W[i][j]) * xStaggeredMesh.D_West[i][j]
				+ mu * (u_Staggered_X_.u_N[i][j] - u_Staggered_X_.u_P[i][j]) * xStaggeredMesh.D_North[i][j]
				- mu * (u_Staggered_X_.u_P[i][j] - u_Staggered_X_.u_S[i][j]) * xStaggeredMesh.D_South[i][j];
			_R_U_[i][j] = _R_U_[i][j] / xStaggeredMesh.Volume[i][j];
			} 
		}



}

