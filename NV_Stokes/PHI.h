// PHI.h
#ifndef PHI_H
#define PHI_H
#include "PHI.h"
#include <iostream>
#include <vector>
#include <string>
#include "mainMesh.h"
#include <iomanip>
#include "xStaggered.h"
#include "yStaggered.h"

#include "HRS_Schemes.h"

using namespace std;


struct PHI_u
{
    vector<vector<double>> u_e;
    vector<vector<double>> u_w;
    vector<vector<double>> u_n;
    vector<vector<double>> u_s;

};

struct PHI_v
{
    vector<vector<double>> v_e;
    vector<vector<double>> v_w;
    vector<vector<double>> v_n;
    vector<vector<double>> v_s;

};


class PHI
{

private:
    PHI_u u;
    PHI_v v;
    float L, H, MeshSize_X, MeshSize_Y;

public:
    void PHI_Initializer();

    void PHI_Setter_BC(float U, u_Staggered_X u_Staggered_X, v_Staggered_Y  v_Staggered_Y);
    void CDS_Scheme(u_Staggered_X u_Staggered_X, v_Staggered_Y  v_Staggered_Y);
    void Upwind_Scheme(u_Staggered_X u_Staggered_X, v_Staggered_Y  v_Staggered_Y, massFlowRate_StaggeredX massFlowRate_StaggeredX_, massFlowRate_StaggeredY massFlowRate_StaggeredY_);
    void HRS(u_Staggered_X u_Staggered_X, v_Staggered_Y  v_Staggered_Y, massFlowRate_StaggeredX massFlowRate_StaggeredX_, massFlowRate_StaggeredY massFlowRate_StaggeredY_, CellX xStaggeredMesh, CellY yStaggeredMesh);
    void QUICK(u_Staggered_X u_Staggered_X, v_Staggered_Y  v_Staggered_Y, massFlowRate_StaggeredX massFlowRate_StaggeredX_, massFlowRate_StaggeredY massFlowRate_StaggeredY_);
    void SOU(u_Staggered_X u_Staggered_X, v_Staggered_Y  v_Staggered_Y, massFlowRate_StaggeredX massFlowRate_StaggeredX_, massFlowRate_StaggeredY massFlowRate_StaggeredY_);
// Getter functions for PHI_u and PHI_v
    PHI_u& getPHI_u() { return u; }
    PHI_v& getPHI_v() { return v; }


    void printAndCheck();

	PHI(float L, float H, float N, float M);
	~PHI();
};

#endif // PHI_H