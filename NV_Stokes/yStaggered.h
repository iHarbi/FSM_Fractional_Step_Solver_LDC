// yStaggered.h
#ifndef YSTAGGERED_H
#define YSTAGGERED_H
#include <iostream>
#include <vector>
#include <string>
#include "mainMesh.h"
using namespace std;

struct StaggeredMesh_Y
{
    Cell Cell;
    Area Area;
    Position Position;
};

struct CellY {
    vector<vector<double>> Volume;
    vector<vector<double>> D_East;
    vector<vector<double>> D_West;
    vector<vector<double>> D_North;
    vector<vector<double>> D_South;
    vector<double> dPE;
    vector<double> dPW;
    vector<double> dPN;
    vector<double> dPS;
    vector<double> CellCentroid_X_Staggered;
    vector<double> CellCentroid_Y_Staggered;
    vector<vector<double>> AreaEast;
    vector<vector<double>> AreaWest;
    vector<vector<double>> AreaNorth;
    vector<vector<double>> AreaSouth;

};


struct massFlowRate_StaggeredY
{
    vector<vector<double>> massFlowRate_East;
    vector<vector<double>> massFlowRate_West;
    vector<vector<double>> massFlowRate_North;
    vector<vector<double>> massFlowRate_South;
};
/* Velocities at the Faces of the mainMesh, or the Centroids of the Staggered One */
/*struct v_Staggered_Y
{
    vector<vector<double>> v_P;
    vector<vector<double>> v_E;
    vector<vector<double>> v_W;
    vector<vector<double>> v_N;
    vector<vector<double>> v_S;
};*/


class yStaggered
{
private:
    float L, H, MeshSize_X, MeshSize_Y;
    massFlowRate_StaggeredY massFlowRate_StaggeredY_;
    v_Staggered_Y v_Staggered_Y_;
    CellY mesh2Modify_Y;

public:
    /*Staggered_Y Mesh Struct: [Cell, Area, Position] : to Return a Struct with Encapsulated Structs*/
    StaggeredMesh_Y structGetter_Y(Cell Cell_, Position Position_, Area Area_);
    
    
    /*Staggering Processes [Method/Function] */
    CellY yStaggeringOPs(StaggeredMesh_Y mainMesh);
    void checkMesh(CellY cell);
    void massFlowRate_StaggeredY_Initializer();
    void massFlowRate_StaggeredY_Setter(v_Staggered_Y v_Staggered_Y_, CellY CellY_, u_Staggered_X u_Staggered_X_, double roh);
    massFlowRate_StaggeredY& getMassFlowRate_StaggeredY() {return massFlowRate_StaggeredY_;}
    void v_Staggered_Y_Initializer();
    v_Staggered_Y& get_v_Staggered_Y() {return v_Staggered_Y_;}
    void v_Staggered_Y_OPs(Cell Cell_Struct, vector<vector<double>> V_Pred, double Delta_t, double roh, vector<vector<double>> PressureField_New);
    void printMassFlowRates();



    yStaggered(float L, float H, float N, float M);
    ~yStaggered();
};
#endif // YSTAGGERED_H