// xStaggered.h
#ifndef XSTAGGERED_H
#define XSTAGGERED_H

#include <iostream>
#include <vector>
#include <string>
#include "mainMesh.h"
#include <iomanip>

using namespace std;

struct StaggeredMesh_X
{
    Cell Cell;
    Area Area; 
    Position Position;
};

struct CellX {
    vector<vector<double>> Volume;
    vector<vector<double>> D_East ;
    vector<vector<double>> D_West ;
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

struct massFlowRate_StaggeredX
{
    vector<vector<double>> massFlowRate_East;
    vector<vector<double>> massFlowRate_West;
    vector<vector<double>> massFlowRate_North;
    vector<vector<double>> massFlowRate_South;
};

/* Velocities at the Faces of the mainMesh, or the Centroids of the Staggered One */
/*struct u_Staggered_X
{
    vector<vector<double>> u_P;
    vector<vector<double>> u_E;
    vector<vector<double>> u_W;
    vector<vector<double>> u_N;
    vector<vector<double>> u_S;
};*/





class xStaggered
{
private:
    float L, H, MeshSize_X, MeshSize_Y;
    massFlowRate_StaggeredX massFlowRate_StaggeredX_;
    CellX mesh2Modify;
    u_Staggered_X u_Staggered_X_;

public:
    /*Staggered Mesh Struct: [Cell, Area, Position] : ro Return a Struct with Encapsulated Structs*/
    StaggeredMesh_X structGetter_X(Cell Cell, Position Position, Area Area);
    /*Staggering Processes [Method] */
    CellX xStaggeringOPs(StaggeredMesh_X mainMesh);
    /*Initialization/Resizing of massFlowRate_StaggeredX*/
    void massFlowRate_StaggeredX_Initializer();
    void u_Staggered_X_Initializer();
    void u_Staggered_X_BC(double Re, double mu, double roh);
    void massFlowRate_StaggeredX_Setter(u_Staggered_X u_Staggered_X_, CellX CellX_, v_Staggered_Y v_Staggered_Y_, double roh);
    /*massFlowRate_StaggeredX Operations*/
    massFlowRate_StaggeredX& getMassFlowRate_StaggeredX() {return massFlowRate_StaggeredX_;}
    u_Staggered_X& get_u_Staggered_X() {return u_Staggered_X_;}
    void u_Staggered_X_OPs(Cell Cell_Struct,vector<vector<double>> U_Pred, double Delta_t, double roh, vector<vector<double>> PressureField_New);
    void printMassFlowRates();
    











    /*Constructor*/
    xStaggered(float L, float H, float N, float M);
    void checkMesh(CellX cell);
    /*Destructor */
    ~xStaggered();

};

#endif // XSTAGGERED_H