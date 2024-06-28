#ifndef MAINMESH_H
#define MAINMESH_H
#include <iostream>
#include <vector>
#include <string>
using namespace std;

/*-+-+-+-+-+-+-+-+-+-+-+-+ Structs Accessible by the Class-+-+-+-+-+-+-+-+-+-+-+-+*/

/*Struct Declaration | XY Positions for Faces and Cell Centroids*/
struct Position {
    vector<double> FaceX;
    vector<double> FaceY;
    vector<double> East;
    vector<double> West;
    vector<double> North;
    vector<double> South;
    vector<double> CellCentroidX;
    vector<double> CellCentroidY;
};
/*Area Struct | Allocating the Positions of Faces: Ease, West, North, and South*/
struct Area {
    vector<double> East;
    vector<double> West;
    vector<double> North;
    vector<double> South;
};
/*Cell Struct | Essential, Passed as an Arg, for the Numerical Solver */
struct Cell {
    vector<vector<double>> Volume;
    vector<vector<double>> D;
    vector<double> dPE;
    vector<double> dPW;
    vector<double> dPN;
    vector<double> dPS;
    vector<double> dPe;
    vector<double> dPw;
    vector<double> dPn;
    vector<double> dPs;
    vector<double> dEe;
    // StaggeredX staggeredX; // Staggered mesh data for X direction
    // StaggeredY staggeredY; // Staggered mesh data for Y direction

};


struct v_Staggered_Y
{
    vector<vector<double>> v_P;
    vector<vector<double>> v_E;
    vector<vector<double>> v_W;
    vector<vector<double>> v_N;
    vector<vector<double>> v_S;
};
/* Velocities at the Faces of the mainMesh, or the Centroids of the Staggered One */
struct u_Staggered_X
{
    vector<vector<double>> u_P;
    vector<vector<double>> u_E;
    vector<vector<double>> u_W;
    vector<vector<double>> u_N;
    vector<vector<double>> u_S;
};

class mainMesh
{
private:
	float L, H, MeshSize_X, MeshSize_Y;

public:
	/*OverLoaded Constructors*/
	mainMesh(float L, float H, float N, float M);
	mainMesh();
    /*Destructor*/
	~mainMesh();

    float dx();
    float dy();
    
    Cell Cell_Setter();
    Position Position_Getter();
    Area Area_Getter(Position position);
    void checkMesh(Cell cell);
    void checkPosition(Position position);
    void checkArea(Area area);
};
#endif // MAINMESH_H
