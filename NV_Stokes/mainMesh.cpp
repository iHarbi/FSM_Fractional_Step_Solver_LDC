#include "mainMesh.h"
float mainMesh::dx()
{
	return (L / MeshSize_X);
}

float mainMesh::dy()
{
	return (H / MeshSize_Y);
}

Cell mainMesh::Cell_Setter()
{
    /* 'Position' Struct Instantiation , Fields Addition, and  Resizing | [N+2, M+2]*/
    Position Position;
    Position.FaceX.resize(MeshSize_X + 2);
    Position.FaceY.resize(MeshSize_Y + 2);

    // Node Centroids Positions | Just to be used for the Sizing Operations
    vector<double> x_P(MeshSize_X + 2); // Cell Centroids X-Direction
    vector<double> y_P(MeshSize_Y + 2); // Cell Centroids Y-Direction

    //
    Position.East.resize(x_P.size());
    Position.West.resize(x_P.size());
    Position.North.resize(y_P.size());
    Position.South.resize(y_P.size());

    // Centroid Position - Struct
    Position.CellCentroidX.resize(x_P.size());
    Position.CellCentroidY.resize(y_P.size());

    /* 'Area' Struct Instantiation , Fields Addition, and  Resizing*/
    Area Area;
    Area.East.resize(x_P.size());
    Area.West.resize(x_P.size());
    Area.North.resize(y_P.size());
    Area.South.resize(y_P.size());

    /* 'Cell' Struct Instantiation , Fields Addition, and  Resizing*/
    Cell Cell;
    Cell.Volume.resize(x_P.size(), vector<double>(y_P.size()));
    Cell.D.resize(x_P.size(), vector<double>(y_P.size()));
    Cell.dPE.resize(x_P.size());
    Cell.dPW.resize(x_P.size());
    Cell.dPN.resize(y_P.size());
    Cell.dPS.resize(y_P.size());
    Cell.dPe.resize(x_P.size());
    Cell.dPw.resize(x_P.size());
    Cell.dPn.resize(y_P.size());
    Cell.dPs.resize(y_P.size());
    Cell.dEe.resize(x_P.size());

    /*Mesh Generation | Face Positioning*/

    Position.FaceX[0] = 0;    Position.FaceX[Position.FaceX.size() - 1] = L; //First and Last Faces [East and West]
    Position.FaceY[0] = 0;    Position.FaceY[Position.FaceY.size() - 1] = H; //First and Last Faces [North and South] 
    for (int i = 1; i <= MeshSize_X; ++i)
        Position.FaceX[i] = (i)*dx();
    for (int i = 1; i <= MeshSize_Y; ++i)
        Position.FaceY[i] = (i)*dy();

    /*Mesh Generation | Cell Centroid Positioning*/
    Position.CellCentroidX[0] = 0;     Position.CellCentroidX[Position.CellCentroidX.size() - 1] = L;
    Position.CellCentroidY[0] = 0;     Position.CellCentroidY[Position.CellCentroidY.size() - 1] = H;

    for (int i = 1; i <= MeshSize_X; ++i)
        Position.CellCentroidX[i] = (Position.FaceX[i] + Position.FaceX[i - 1]) / 2;
    for (int i = 1; i <= MeshSize_Y; ++i)
        Position.CellCentroidY[i] = (Position.FaceY[i] + Position.FaceY[i - 1]) / 2;
    /* Mesh Generation | Eastern & Western Cells | Interior Nodes */
    Position.East[0] = Position.FaceX[0];     Position.West[Position.West.size() - 1] = Position.FaceX[Position.FaceX.size() - 1];
    Position.North[0] = Position.FaceY[0];    Position.South[Position.South.size() - 1] = Position.FaceY[Position.FaceY.size() - 1];
    // Interior Nodes | X/Y Directions
    for (int i = 1; i <= MeshSize_X; ++i) {
        Position.East[i] = Position.FaceX[i];
        Position.West[i] = Position.FaceX[i - 1];
    }
    Position.East[Position.East.size() - 1] = L;

    for (int i = 1; i <= MeshSize_Y + 1; ++i) {
        Position.North[i] = Position.FaceY[i];
        Position.South[i] = Position.FaceY[i - 1];
    }
    Position.South[Position.South.size() - 1] = H;

    // Mesh Generation | Volume
    for (int i = 1; i <= MeshSize_X; ++i)
        for (int j = 1; j <= MeshSize_Y; ++j)
            Cell.Volume[i][j] = (Position.East[i] - Position.West[i]) * (Position.North[j] - Position.South[j]);

    /*Mesh Generation | Cell Distancing */
    for (int i = 0; i <= MeshSize_X; ++i) {
        Cell.dPE[i] = Position.CellCentroidX[i + 1] - Position.CellCentroidX[i];
        Cell.dPW[i + 1] = Position.CellCentroidX[i + 1] - Position.CellCentroidX[i];
        Cell.dPw[i + 1] = Position.CellCentroidX[i + 1] - Position.FaceX[i];
        Cell.dEe[i] = Position.CellCentroidX[i + 1] - Position.FaceX[i];
    }

    for (int i = 0; i <= MeshSize_Y; ++i) {
        Cell.dPN[i] = Position.CellCentroidY[i + 1] - Position.CellCentroidY[i];
        Cell.dPS[i + 1] = Position.CellCentroidY[i + 1] - Position.CellCentroidY[i];
        // Cell.dPn(i+1) = Position.CellCentroidX(i+1) - Position.FaceX(i);
        // Cell.dEe(i) =   Position.CellCentroidX(i+1) - Position.FaceX(i);
    }

    for (int i = 0; i <= MeshSize_X + 1; ++i) {
        Cell.dPe[i] = Position.FaceX[i] - Position.CellCentroidX[i];
    }

    for (int i = 0; i <= MeshSize_Y + 1; ++i) {
        Cell.dPn[i] = Position.FaceY[i] - Position.CellCentroidY[i];
    }

    /* Mesh Generation | Eastern & Western Area */
    for (int i = 0; i < MeshSize_X + 2; ++i) {
        if (i == MeshSize_X) {
            Area.East[i] = Position.FaceY[i] - Position.FaceY[i - 1];
            Area.West[i] = Position.FaceY[i - 1] - Position.FaceY[i - 2];
        }
        else if (i == 0) {
            Area.East[i] = Position.FaceY[i + 1] - Position.FaceY[i];
            Area.West[i] = Position.FaceY[i + 1] - Position.FaceY[i];
        }
        else if (i == MeshSize_X + 1) {
            Area.East[i] = Position.FaceY[i - 1] - Position.FaceY[i - 2];
            Area.West[i] = Position.FaceY[i - 1] - Position.FaceY[i - 2];
        }
        else {
            Area.East[i] = Position.FaceY[i + 1] - Position.FaceY[i];
            Area.West[i] = Position.FaceY[i] - Position.FaceY[i - 1];
        }
    }

    /*Mesh Generation | Northern & Southern Area*/
    for (int i = 0; i < MeshSize_Y + 2; ++i) {
        if (i == MeshSize_Y) {
            Area.North[i] = Position.FaceX[i] - Position.FaceX[i - 1];
            Area.South[i] = Position.FaceX[i - 1] - Position.FaceX[i - 2];
        }
        else if (i == 0) {
            Area.North[i] = Position.FaceX[i + 1] - Position.FaceX[i];
            Area.South[i] = Position.FaceX[i + 1] - Position.FaceX[i];
        }
        else if (i == MeshSize_Y + 1) {
            Area.North[i] = Position.FaceX[i - 2] - Position.FaceX[i - 3];
            Area.South[i] = Position.FaceX[i - 2] - Position.FaceX[i - 3];
        }
        else {
            Area.North[i] = Position.FaceX[i + 1] - Position.FaceX[i];
            Area.South[i] = Position.FaceX[i] - Position.FaceX[i - 1];
        }
    }

    return Cell;
}

Position mainMesh::Position_Getter()
{
    /* 'Position' Struct Instantiation , Fields Addition, and  Resizing | [N+2, M+2]*/
    Position Position;
    Position.FaceX.resize(MeshSize_X + 2);
    Position.FaceY.resize(MeshSize_Y + 2);

    // Node Centroids Positions | Just to be used for the Sizing Operations
    vector<double> x_P(MeshSize_X + 2); // Cell Centroids X-Direction
    vector<double> y_P(MeshSize_Y + 2); // Cell Centroids Y-Direction

    //
    Position.East.resize(x_P.size());
    Position.West.resize(x_P.size());
    Position.North.resize(y_P.size());
    Position.South.resize(y_P.size());

    // Centroid Position - Struct
    Position.CellCentroidX.resize(x_P.size());
    Position.CellCentroidY.resize(y_P.size());
    /*Mesh Generation | Face Positioning*/

    Position.FaceX[0] = 0;    Position.FaceX[Position.FaceX.size() - 1] = L; //First and Last Faces [East and West]
    Position.FaceY[0] = 0;    Position.FaceY[Position.FaceY.size() - 1] = H; //First and Last Faces [North and South] 
    for (int i = 1; i <= MeshSize_X; ++i)
        Position.FaceX[i] = (i)*dx();
    for (int i = 1; i <= MeshSize_Y; ++i)
        Position.FaceY[i] = (i)*dy();

    /*Mesh Generation | Cell Centroid Positioning*/
    Position.CellCentroidX[0] = 0;     Position.CellCentroidX[Position.CellCentroidX.size() - 1] = L;
    Position.CellCentroidY[0] = 0;     Position.CellCentroidY[Position.CellCentroidY.size() - 1] = H;

    for (int i = 1; i <= MeshSize_X; ++i)
        Position.CellCentroidX[i] = (Position.FaceX[i] + Position.FaceX[i - 1]) / 2;
    for (int i = 1; i <= MeshSize_Y; ++i)
        Position.CellCentroidY[i] = (Position.FaceY[i] + Position.FaceY[i - 1]) / 2;
    /* Mesh Generation | Eastern & Western Cells | Interior Nodes */
    Position.East[0] = Position.FaceX[0];     Position.West[Position.West.size() - 1] = Position.FaceX[Position.FaceX.size() - 1];
    Position.North[0] = Position.FaceY[0];    Position.South[Position.South.size() - 1] = Position.FaceY[Position.FaceY.size() - 1];
    // Interior Nodes | X/Y Directions
    for (int i = 1; i <= MeshSize_X; ++i) {
        Position.East[i] = Position.FaceX[i];
        Position.West[i] = Position.FaceX[i - 1];
    }
    Position.East[Position.East.size() - 1] = L;

    for (int i = 1; i <= MeshSize_Y + 1; ++i) {
        Position.North[i] = Position.FaceY[i];
        Position.South[i] = Position.FaceY[i - 1];
    }
    Position.South[Position.South.size() - 1] = H;
    return Position;
}

Area mainMesh::Area_Getter(Position Position)
{
    // Node Centroids Positions | Just to be used for the Sizing Operations
    vector<double> x_P(MeshSize_X + 2); // Cell Centroids X-Direction
    vector<double> y_P(MeshSize_Y + 2); // Cell Centroids Y-Direction
    Area Area;
    Area.East.resize(x_P.size());
    Area.West.resize(x_P.size());
    Area.North.resize(y_P.size());
    Area.South.resize(y_P.size());
    /* Mesh Generation | Eastern & Western Area */
    for (int i = 0; i < MeshSize_X + 2; ++i) {
        if (i == MeshSize_X) {
            Area.East[i] = Position.FaceY[i] - Position.FaceY[i - 1];
            Area.West[i] = Position.FaceY[i - 1] - Position.FaceY[i - 2];
        }
        else if (i == 0) {
            Area.East[i] = Position.FaceY[i + 1] - Position.FaceY[i];
            Area.West[i] = Position.FaceY[i + 1] - Position.FaceY[i];
        }
        else if (i == MeshSize_X + 1) {
            Area.East[i] = Position.FaceY[i - 1] - Position.FaceY[i - 2];
            Area.West[i] = Position.FaceY[i - 1] - Position.FaceY[i - 2];
        }
        else {
            Area.East[i] = Position.FaceY[i + 1] - Position.FaceY[i];
            Area.West[i] = Position.FaceY[i] - Position.FaceY[i - 1];
        }
    }

    /*Mesh Generation | Northern & Southern Area*/
    for (int i = 0; i < MeshSize_Y + 2; ++i) {
        if (i == MeshSize_Y) {
            Area.North[i] = Position.FaceX[i] - Position.FaceX[i - 1];
            Area.South[i] = Position.FaceX[i - 1] - Position.FaceX[i - 2];
        }
        else if (i == 0) {
            Area.North[i] = Position.FaceX[i + 1] - Position.FaceX[i];
            Area.South[i] = Position.FaceX[i + 1] - Position.FaceX[i];
        }
        else if (i == MeshSize_Y + 1) {
            Area.North[i] = Position.FaceX[i - 2] - Position.FaceX[i - 3];
            Area.South[i] = Position.FaceX[i - 2] - Position.FaceX[i - 3];
        }
        else {
            Area.North[i] = Position.FaceX[i + 1] - Position.FaceX[i];
            Area.South[i] = Position.FaceX[i] - Position.FaceX[i - 1];
        }
    }
    return Area;
}

// Function to print the Cell Distances: [dPE, dPW, dPN, dPS]
void mainMesh::checkMesh(Cell cell) {
    cout << "\033[1;32m-+-+-+-+-+-+-+-+-+-+-+-+checkMesh-+-+-+-+-+-+-+-+-+-+-+-+\033[0m" << endl;
    cout << "\033[1;31mdPE\t" << "\033[1;31mdPW\t" << "\033[1;31mdPN\t" << "\033[1;31mdPS\033[0m" << endl;
    for (size_t i = 0; i < cell.dPE.size(); ++i) {
        cout << cell.dPE[i] << "\t" << cell.dPW[i] << "\t" << cell.dPN[i] << "\t" << cell.dPS[i] << endl;
    }
}

// Function to print the <<Position>> Struct: [FaceX, FaceY, East, West, North, South, CentroidX, CentroidY] 
void mainMesh::checkPosition(Position position)
{
        cout << "\033[1;32m-+-+-+-+-+-+-+-+-+-+-+-+checkFaces-+-+-+-+-+-+-+-+-+-+-+-+\033[0m" << endl;  
        cout << "\033[1;31mFaceX\t\t" << "\033[1;31mFaceY\t\t" << "\033[1;31mEast\t\t" << "\033[1;31mWest\t\t" << "\033[1;31mNorth\t\t" << "\033[1;31mSouth\t\t" << "\033[1;31mCentroidX\t" << "\033[1;31mCentroidY\033[0m" << endl;
        size_t size = max({ position.FaceX.size(), position.FaceY.size(), position.East.size(), position.West.size(), position.North.size(), position.South.size(), position.CellCentroidX.size(), position.CellCentroidY.size() });
        for (size_t i = 0; i < size; ++i) {
            cout << (i < position.FaceX.size() ? to_string(position.FaceX[i]) : " ") << "\t";
            cout << (i < position.FaceY.size() ? to_string(position.FaceY[i]) : " ") << "\t";
            cout << (i < position.East.size() ?  to_string(position.East[i]) : " ") << "\t";
            cout << (i < position.West.size() ?  to_string(position.West[i]) : " ") << "\t";
            cout << (i < position.North.size() ? to_string(position.North[i]) : " ") << "\t";
            cout << (i < position.South.size() ? to_string(position.South[i]) : " ") << "\t";
            cout << (i < position.CellCentroidX.size() ? to_string(position.CellCentroidX[i]) : " ") << "\t";
            cout << (i < position.CellCentroidY.size() ? to_string(position.CellCentroidY[i]) : " ") << endl;
        }
}
// Function to print the <<Area>> Struct: [East, West, North, South] Areas 

void mainMesh::checkArea(Area area)
{
    cout << "\033[1;32m-+-+-+-+-+-+-+-+-+-+-+-+checkAreas-+-+-+-+-+-+-+-+-+-+-+-+\033[0m" << endl;
    cout << "\033[1;31mEast\t\t" << "West\t\t" << "North\t\t" << "South\033[0m" << endl;
    size_t size = max({ area.East.size(), area.West.size(), area.North.size(), area.South.size() });
    for (size_t i = 0; i < size; ++i) {
        cout << (i < area.East.size() ? to_string(area.East[i]) : " ") << "\t";
        cout << (i < area.West.size() ? to_string(area.West[i]) : " ") << "\t";
        cout << (i < area.North.size() ? to_string(area.North[i]) : " ") << "\t";
        cout << (i < area.South.size() ? to_string(area.South[i]) : " ") << endl;
    }
}


mainMesh::mainMesh(float L, float H, float N, float M) :L(L), H(H), MeshSize_X(N), MeshSize_Y(M)
{
}

mainMesh::mainMesh() :L(0.0), H(0.0), MeshSize_X(0.0), MeshSize_Y(0.0)
{
}

mainMesh::~mainMesh()
{
    // Change text color to green
    cout << "\033[1;32mDestructor FLAG: Works fine tell Now\033[0m" << endl;
}
