// NV_Stokes.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include <iostream>
using namespace std;
#include "mainMesh.h"
#include "xStaggered.h"
#include "yStaggered.h"
#include "PHI.h"
#include "Pressure.h"
#include "CSVWriter.h"
#include "CSVMATLAB.h"
#include <chrono>
#include "Export.h"

/*Functions Declaration | _R_U_ */
void _R_U_Initializer(vector<vector<double>>& _R_U_,  vector<vector<double>>& U_Pred, float MeshSize_X, float MeshSize_Y);
vector<vector<double>> _R_U_OPS(CellX xStaggeredMesh, massFlowRate_StaggeredX massFlowRate_StaggeredX_, u_Staggered_X u_Staggered_X_, PHI_u u, double mu, float MeshSize_X, float MeshSize_Y, vector<vector<double>> _R_U_);
/*Functions Declaration | _R_V_ */
void _R_V_Initializer(vector<vector<double>>& _R_V_,  vector<vector<double>>& V_Pred, float MeshSize_X, float MeshSize_Y);
vector<vector<double>> _R_V_OPS(CellY yStaggeredMesh, massFlowRate_StaggeredY massFlowRate_StaggeredY_, v_Staggered_Y v_Staggered_Y_, PHI_v v, double mu, float MeshSize_X, float MeshSize_Y, vector<vector<double>> _R_V_);
/*Functions Declaration | Predictor uP*/
vector<vector<double>> U_Predictor(vector<vector<double>>& _R_U_, vector<vector<double>>& _R_U_Prev, vector<vector<double>>& U_Pred, u_Staggered_X u_Staggered_X_, double Delta_t, double roh, float MeshSize_X, float MeshSize_Y, const string& spatialOperator);
/*Functions Declaration | Predictor vP*/
vector<vector<double>> V_Predictor(vector<vector<double>>& _R_V_, vector<vector<double>>& _R_V_Prev, vector<vector<double>>& V_Pred, v_Staggered_Y v_Staggered_Y_, double Delta_t, double roh, float MeshSize_X, float MeshSize_Y, const string& spatialOperator);

/*CFL Condition | Δt */
double CFL(u_Staggered_X u_Staggered_X_, v_Staggered_Y v_Staggered_Y_, double MeshSize_X, double MeshSize_Y, double roh, double Viscosity);
bool checkSteadyState(const u_Staggered_X& uCurrent, const u_Staggered_X& uPrevious, const v_Staggered_Y& vCurrent, const v_Staggered_Y& vPrevious, double tolerance, double MeshSize_X, double MeshSize_Y);
/*Mass Balance | Divergence-free Field*/
double massBalance(u_Staggered_X u_Staggered_X_, v_Staggered_Y v_Staggered_Y_, double roh, Area Area_, float MeshSize_X, float MeshSize_Y);

int main()
{
    double L = 1.0, H = 1.0, Nx = 100, My = 100, mu = 1.8e-5;
    /*Define Re*/
    double Re = 100;
    string Scheme = "CDS";
    double roh = Re * mu, Conv_Crit = 1E-6, t = 0.0, timeElapsed = 0.0, U = 1.0;
    int itrMax = 30000, LoopNo = 1;
    vector<vector<double>> _R_U_, _R_V_, U_Pred, V_Pred, _R_U_Prev, _R_V_Prev, FinalPressureField;
    double tolerance = 1E-6; /*Temporal Loop*/
    bool isSteadyState = false;
    double Simulation_time = 15; /* a Rule of thumb: The fluid should pass through the domain 10 times to reach S.S. */
    auto start = std::chrono::high_resolution_clock::now();
    /* -+-+-+-+-+--+-+-+-+-+-Class: ONE-+-+-+-+-+--+-+-+-+-+- */
    /*<<mainMesh>> a Class that Generates the Mesh, with a Constructor taking [L, H, Nx, My] as I/P */
    mainMesh MeshTrial(L, H, Nx, My);
    /*-+-+-+-+-+-+-+-+-+-+-+-+ Position/Area Structs -+-+-+-+-+-+-+-+-+-+-+-+*/
    Cell Cell_Struct = MeshTrial.Cell_Setter();
    Position Position_Struct = MeshTrial.Position_Getter();
    Area Area_Struct = MeshTrial.Area_Getter(Position_Struct);
    /*-+-+-+-+-+-+-+-+-+-+-+-+checkMesh-+-+-+-+-+-+-+-+-+-+-+-+*/
    MeshTrial.checkPosition(Position_Struct);
    MeshTrial.checkArea(Area_Struct);
    MeshTrial.checkMesh(Cell_Struct);
    /*-+-+-+-+-+-+-+-+-+-+-+-+Staggering-+-+-+-+-+-+-+-+-+-+-+-+*/
    /* Class: TWO */
    /*Definition of Return Struct / Class Object*/

    StaggeredMesh_X mainMesh_X_arg; 
    StaggeredMesh_Y mainMesh_Y_arg;
    
    xStaggered xStaggered(L, H, Nx, My);
    yStaggered yStaggered(L, H, Nx, My);
    
    mainMesh_X_arg = xStaggered.structGetter_X(Cell_Struct, Position_Struct, Area_Struct);
    mainMesh_Y_arg = yStaggered.structGetter_Y(Cell_Struct, Position_Struct, Area_Struct);
    CellX xStaggeredMesh = xStaggered.xStaggeringOPs(mainMesh_X_arg);
    CellY yStaggeredMesh = yStaggered.yStaggeringOPs(mainMesh_Y_arg);

    /*-+-+-+-+-+-+-+-+-+-+-+-+staggeringCheck-+-+-+-+-+-+-+-+-+-+-+-+*/
    xStaggered.checkMesh(xStaggeredMesh);
    yStaggered.checkMesh(yStaggeredMesh);
    /*-+-+-+-+-+-+-+-+-+-+-+-+Advected FIELD: PHI | Velocity-+-+-+-+-+-+-+-+-+-+-+-+*/
    PHI PHI(L, H, Nx, My);
    PHI.PHI_Initializer(); /*-RESIZING*/
    /*/*-+-+-+-+-+-+-+-+-+-+-+-+Initialization, Resizing and & Reference Return/*-+-+-+-+-+-+-+-+-+-+-+-+*/
    xStaggered.massFlowRate_StaggeredX_Initializer();
    /*-+-+-+-+-+-+-+-+-+-+-+-+-+Initialization, Resizing and Reference Return-+-+-+-+-+-+-+-+-+-+-+-+*/
    yStaggered.massFlowRate_StaggeredY_Initializer();
    /*-+-+-+-+-+-+-+-+-+-+-+-+Initialization, Resizing and Reference Return-+-+-+-+-+-+-+-+-+-+-+-+*/
    xStaggered.u_Staggered_X_Initializer(); /*uP, uE, uW, uN, uS*/
    xStaggered.u_Staggered_X_BC(Re, mu, roh);
    u_Staggered_X u_Staggered_X_ = xStaggered.get_u_Staggered_X();
    /*-+-+-+-+-+-+-+-+-+-+-+-+Initialization, Resizing and Reference Return-+-+-+-+-+-+-+-+-+-+-+-+*/
    yStaggered.v_Staggered_Y_Initializer(); /*vP, vE, vW, vN, vS*/
    v_Staggered_Y v_Staggered_Y_ = yStaggered.get_v_Staggered_Y();

    /*PressureField Operations*/
    Pressure PressureField(L, H, Nx, My);
    //PressureField.PressureField_Initialization();
    PressureField.PressureField_Initialization();
    PressureField.Coeff_Vectors_Initialization();

    u_Staggered_X uPrevious = u_Staggered_X_, uCurrent = u_Staggered_X_;
    v_Staggered_Y vPrevious = v_Staggered_Y_, vCurrent = v_Staggered_Y_;

    double timeStep = CFL(u_Staggered_X_, v_Staggered_Y_, Nx, My, roh, mu);
    timeElapsed += timeStep;



    _R_U_Prev.resize(Nx + 1, vector<double>(My + 2));
    _R_V_Prev.resize(Nx + 2, vector<double>(My + 1));
    FinalPressureField.resize(Nx + 2, vector<double>(My + 2));

    Export Export_Obj;
    /* *-*-*-*-*-*-*-**-*-*- Temporal Loop *-*-*-*-*-*-*-**-*-*- */
    while (!isSteadyState && timeElapsed < Simulation_time) {

        LoopNo++;
        std::cout << "Iteration: " << LoopNo << endl;
       
        u_Staggered_X_ = uPrevious;
        v_Staggered_Y_ = vPrevious;


        /*massFlow Rate Calculations*/
        xStaggered.massFlowRate_StaggeredX_Setter(u_Staggered_X_, xStaggeredMesh, v_Staggered_Y_, roh);
        massFlowRate_StaggeredX massFlowRate_StaggeredX_ = xStaggered.getMassFlowRate_StaggeredX();
        yStaggered.massFlowRate_StaggeredY_Setter(v_Staggered_Y_, yStaggeredMesh, u_Staggered_X_, roh);
        massFlowRate_StaggeredY massFlowRate_StaggeredY_ = yStaggered.getMassFlowRate_StaggeredY();

        PHI.PHI_Setter_BC(U, u_Staggered_X_, v_Staggered_Y_); /*B. Convection*/

        if (Scheme == "UDS"){PHI.Upwind_Scheme(u_Staggered_X_, v_Staggered_Y_, massFlowRate_StaggeredX_, massFlowRate_StaggeredY_);}
        else if (Scheme == "CDS") { PHI.CDS_Scheme(u_Staggered_X_, v_Staggered_Y_); }
        else if (Scheme == "QUICK") { PHI.QUICK(u_Staggered_X_, v_Staggered_Y_, massFlowRate_StaggeredX_, massFlowRate_StaggeredY_); }
        else if (Scheme == "SOU") {PHI.SOU(u_Staggered_X_, v_Staggered_Y_, massFlowRate_StaggeredX_, massFlowRate_StaggeredY_);}

        PHI_u u = PHI.getPHI_u(); /*[u_e, u_w, u_n, u_s]*/
        PHI_v v = PHI.getPHI_v(); /*[v_e, v_w, v_n, v_s]*/

        /*R Operations | Predictor Velocity*/
        _R_U_Initializer(_R_U_, U_Pred, Nx, My);
        _R_U_ = _R_U_OPS(xStaggeredMesh, massFlowRate_StaggeredX_, u_Staggered_X_, u, mu, Nx, My, _R_U_);
        if (LoopNo == 2) { _R_U_Prev = _R_U_;  }
        U_Pred = U_Predictor(_R_U_, _R_U_Prev, U_Pred, u_Staggered_X_, timeStep, roh, Nx, My, "AdamBashforth");

        _R_V_Initializer(_R_V_, V_Pred, Nx, My);
        _R_V_ = _R_V_OPS(yStaggeredMesh, massFlowRate_StaggeredY_, v_Staggered_Y_, v, mu, Nx, My, _R_V_);
        if (LoopNo == 2) { _R_V_Prev = _R_V_; }
        V_Pred = V_Predictor(_R_V_, _R_V_Prev, V_Pred, v_Staggered_Y_, timeStep, roh, Nx, My, "AdamBashforth");


        PressureField.Coeff_Vectors_InternalNodes(Cell_Struct, Area_Struct, U_Pred, V_Pred, timeStep, roh);
        PressureField.Coeff_Vectors_BCs();
        Coeff_Vectors PressureLINSYS = PressureField.Coffs_Vector_Getter();
        vector<vector<double>> PressureField_ = PressureField.PressureField_Getter();

        PressureField.LIN_EQ_Solver(PressureLINSYS, PressureField_, Conv_Crit, itrMax);
        vector<vector<double>> PressureField_New = PressureField.PressureFieldNew_Getter();

        xStaggered.u_Staggered_X_OPs(Cell_Struct, U_Pred, timeStep, roh, PressureField_New);
        yStaggered.v_Staggered_Y_OPs(Cell_Struct, V_Pred, timeStep, roh, PressureField_New);
        /*Return u_P and v_P | Staggered X and Y*/
        uCurrent = xStaggered.get_u_Staggered_X();
        vCurrent = yStaggered.get_v_Staggered_Y();
        /* CFL Condition */
        timeStep = CFL(uCurrent, vCurrent, Nx, My, roh, mu);
        //timeStep = 0.0005;


        std::cout << "CFL Delta t:" << timeStep << endl;
        
        // Check for steady state condition
        // You need to implement a method to check convergence
        isSteadyState = checkSteadyState(uCurrent, uPrevious, vCurrent, vPrevious, tolerance, Nx, My);
        
        uPrevious = uCurrent;
        vPrevious = vCurrent;
        PressureField_ = PressureField_New;
        FinalPressureField = PressureField_New;
        _R_U_Prev = _R_U_;
        _R_V_Prev = _R_V_;
        // Increment time
        timeElapsed += timeStep;
        std::cout << "timeElapsed:" << timeElapsed << endl;
        /*CONVERGENCE FLAG*/
        if (isSteadyState) {
            std::cout << "Solution has converged!" << endl;
        }
        else {
            std::cout << "Solution has not converged yet." << endl;
        }

        /* --- Temporal Evolution of u(y), v(x)---*/
        Export_Obj.addUyData(uCurrent.u_P, Nx / 2 + 1 );
        Export_Obj.addVxData(vCurrent.v_P, My / 2 + 1);

        /*-+-+-+-+-+-+-+-+-+-+-+-+ Check Staggered massFlowRates -+-+-+-+-+-+-+-+-+-+-+-+*/
        //cout << "\033[1;32m massFlowRate Staggered X \t " << massFlowRate_StaggeredX_.massFlowRate_East.size() << "*" << massFlowRate_StaggeredX_.massFlowRate_North[0].size() << "\t" << "\033[1;31m PHI u: \t " << u.u_e.size() << "*" << u.u_e[0].size() << "\n";
        //cout << "\033[1;32m massFlowRate Staggered Y \t " << massFlowRate_StaggeredY_.massFlowRate_East.size() << "*" << massFlowRate_StaggeredY_.massFlowRate_North[0].size() << "\t" << "\033[1;31m PHI v: \t " << v.v_e.size() << "*" << v.v_e[0].size() << "\n";
        /*-+-+-+-+-+-+-+-+-+-+-+-+ */
        //PHI.printAndCheck();
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    //xStaggered.printMassFlowRates();
    //yStaggered.printMassFlowRates();

    /*Writer for ParaVIEW*/
       // CSVWriter uCurrentWriter("uCurrent.csv");
    //CSVWriter vCurrentWriter("vCurrent.csv");
    
/*
    if (uCurrentWriter.writeMatrix(uCurrent.u_P, Position_Struct)) {
        std::cout << "uCurrent successfully written to uCurrent.csv" << endl;
    }
    else {
        cerr << "Error writing uCurrent to uCurrent.csv" << endl;
    }

    // Write vCurrent to CSV file
    if (vCurrentWriter.writeMatrix(vCurrent.v_P, Position_Struct)) {
        std::cout << "vCurrent successfully written to vCurrent.csv" << endl;
    }
    else {
        cerr << "Error writing vCurrent to vCurrent.csv" << endl;
    }*/

    /*Writer for MATLAB*/
    CSVMATLAB uWriter("uCurrent_MATLAB_Re100_100_CDS.csv");
    CSVMATLAB vWriter("vCurrent_MATLAB_Re100_100_CDS.csv");
    CSVMATLAB PressureWriter("PressureField_Re100_100_CDS.csv");
    CSVMATLAB v_xWritter("v_x_Re100_100_CDS.csv");
    CSVMATLAB u_yWritter("u_y_Re100_100_CDS.csv");
    if (uWriter.writeMatrix(uCurrent.u_P)) {
        std::cout << "uCurrent successfully written to uCurrent_MATLAB.csv" << endl;
    }
    else {
        cerr << "Error writing uCurrent to vCurrent_MATLAB.csv" << endl;
    }

    // Write vCurrent to CSV file
    if (vWriter.writeMatrix(vCurrent.v_P)) {
        std::cout << "vCurrent successfully written to vCurrent.csv" << endl;
    }
    else {
        cerr << "Error writing vCurrent to vCurrent.csv" << endl;
    }
    if (PressureWriter.writeMatrix(FinalPressureField)) {
        std::cout << "Pressure Field successfully written to PressureField.csv" << endl;
    }
    else {
        cerr << "Error writing PressureField to PressureField.csv" << endl;
    }

    LinePlots Export_Struct = Export_Obj.getLinePlots();
    if (v_xWritter.writeMatrix(Export_Struct.v_x)) {
        std::cout << "Temporal Evolution of [v_x] successfully written to v_x.csv" << endl;
    }
    else {
        cerr << "Error writing uCurrent to v_x.csv" << endl;
    }
    if (u_yWritter.writeMatrix(Export_Struct.u_y)) {
        std::cout << "Temporal Evolution of [u_y] successfully written to u_y.csv" << endl;
    }
    else {
        cerr << "Error writing uCurrent to u_y.csv" << endl;
    }
    double massBalance_ = massBalance(uCurrent, vCurrent, roh, Area_Struct, Nx, My);
    std::cout << "massBalnce on the mainGrid is " << massBalance_ <<endl;
    std::cout << "Simulation took " << duration.count() << " seconds." << std::endl;
    //PHI.printAndCheck();
}


/*Functions Implementation: Use Utility Class Later*/
void _R_U_Initializer(vector<vector<double>>& _R_U_, vector<vector<double>>& U_Pred, float MeshSize_X, float MeshSize_Y)
{
    _R_U_.resize(MeshSize_X + 1, vector<double>(MeshSize_Y + 2 ));
    U_Pred.resize(MeshSize_X + 1, vector<double>(MeshSize_Y + 2));
    //_R_U_Prev.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 1));

}

vector<vector<double>> _R_U_OPS(CellX xStaggeredMesh, massFlowRate_StaggeredX massFlowRate_StaggeredX_, u_Staggered_X u_Staggered_X_, PHI_u u, double mu, float MeshSize_X, float MeshSize_Y, vector<vector<double>> _R_U_)
{
    for (int i = 1; i < MeshSize_X; ++i)
        for (int j = 1; j <= MeshSize_Y; ++j)
        {
            {
                double conv_term = -(massFlowRate_StaggeredX_.massFlowRate_East[i][j - 1] * u.u_e[i][j - 1] - massFlowRate_StaggeredX_.massFlowRate_West[i][j - 1] * u.u_w[i][j - 1] + massFlowRate_StaggeredX_.massFlowRate_North[i][j - 1] * u.u_n[i][j - 1] - massFlowRate_StaggeredX_.massFlowRate_South[i][j - 1] * u.u_s[i][j - 1]);
                double diff_term_x = mu * (u_Staggered_X_.u_P[i + 1][j] - u_Staggered_X_.u_P[i][j]) * xStaggeredMesh.D_East[i][j - 1] - mu * (u_Staggered_X_.u_P[i][j] - u_Staggered_X_.u_P[i - 1][j]) * xStaggeredMesh.D_West[i][j - 1];
                double diff_term_y = mu * (u_Staggered_X_.u_P[i][j + 1] - u_Staggered_X_.u_P[i][j]) * xStaggeredMesh.D_North[i][j - 1] - mu * (u_Staggered_X_.u_P[i][j] - u_Staggered_X_.u_P[i][j - 1]) * xStaggeredMesh.D_South[i][j - 1];

                _R_U_[i][j] = (conv_term + diff_term_x + diff_term_y) /  xStaggeredMesh.Volume[i][j - 1];
            }
        }
    return _R_U_;
}

void _R_V_Initializer(vector<vector<double>>& _R_V_, vector<vector<double>>& V_Pred, float MeshSize_X, float MeshSize_Y)
{
     _R_V_.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 1));
    V_Pred.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 1));
    //_R_V_Prev.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 1));
}

vector<vector<double>> _R_V_OPS(CellY yStaggeredMesh, massFlowRate_StaggeredY massFlowRate_StaggeredY_, v_Staggered_Y v_Staggered_Y_, PHI_v v, double mu, float MeshSize_X, float MeshSize_Y, vector<vector<double>> _R_V_)
{
    for (int i = 1; i <= MeshSize_X; ++i) {
        for (int j = 1; j < MeshSize_Y; ++j) {
            double conv_term = -( massFlowRate_StaggeredY_.massFlowRate_East[i - 1][j] * v.v_e[i - 1][j] -
                massFlowRate_StaggeredY_.massFlowRate_West[i - 1][j] * v.v_w[i - 1][j] +
                massFlowRate_StaggeredY_.massFlowRate_North[i - 1][j] * v.v_n[i - 1][j] -
                massFlowRate_StaggeredY_.massFlowRate_South[i - 1][j] * v.v_s[i - 1][j] );

            double diff_term_x = mu * (v_Staggered_Y_.v_P[i + 1][j] - v_Staggered_Y_.v_P[i][j]) * yStaggeredMesh.D_East[i - 1][j] -
                mu * (v_Staggered_Y_.v_P[i][j] - v_Staggered_Y_.v_P[i - 1][j]) * yStaggeredMesh.D_West[i - 1][j];

            double diff_term_y = mu * (v_Staggered_Y_.v_P[i][j + 1] - v_Staggered_Y_.v_P[i][j]) * yStaggeredMesh.D_North[i - 1][j] -
                mu * (v_Staggered_Y_.v_P[i][j] - v_Staggered_Y_.v_P[i][j - 1]) * yStaggeredMesh.D_South[i - 1][j];

            _R_V_[i][j] = (conv_term + diff_term_x + diff_term_y) / yStaggeredMesh.Volume[i - 1][j];
        }
    }
    return _R_V_;
}


vector<vector<double>> U_Predictor(vector<vector<double>>& _R_U_, vector<vector<double>>& _R_U_Prev, vector<vector<double>>& U_Pred, u_Staggered_X u_Staggered_X_, double Delta_t, double roh, float MeshSize_X, float MeshSize_Y, const string& spatialOperator)
{
    for (int i = 1; i < MeshSize_X; ++i)
        for (int j = 1; j <=  MeshSize_Y ; ++j)
        {
            {

                if (spatialOperator == "AdamBashforth") {
                    U_Pred[i][j] = u_Staggered_X_.u_P[i][j] + (Delta_t / roh) * (1.5 * _R_U_[i][j] - 0.5 * _R_U_Prev[i][j]);
                }
                else if (spatialOperator == "ForwardEuler") {
                    U_Pred[i][j] = u_Staggered_X_.u_P[i][j] + (Delta_t / roh) * (1.0 * _R_U_[i][j] - 0.0 * _R_U_Prev[i][j]);
                }
            }
        }
    return U_Pred;
}

vector<vector<double>> V_Predictor(vector<vector<double>>& _R_V_, vector<vector<double>>& _R_V_Prev, vector<vector<double>>& V_Pred, v_Staggered_Y v_Staggered_Y_, double Delta_t, double roh, float MeshSize_X, float MeshSize_Y, const string& spatialOperator)
{

    for (int i = 1; i <= MeshSize_X; ++i)
        for (int j = 1; j < MeshSize_Y; ++j)
        {
            {
                if (spatialOperator == "AdamBashforth") {
                    V_Pred[i][j] = v_Staggered_Y_.v_P[i][j] + (Delta_t / roh) * (1.5 * _R_V_[i][j] - 0.5 * _R_V_Prev[i][j]);
                }
                else if (spatialOperator == "ForwardEuler") {
                    V_Pred[i][j] = v_Staggered_Y_.v_P[i][j] + (Delta_t / roh) * (1.0 * _R_V_[i][j] - 0.0 * _R_V_Prev[i][j]);
                }

            }
        }
    return V_Pred;
}



double CFL(u_Staggered_X u_Staggered_X_, v_Staggered_Y v_Staggered_Y_, double MeshSize_X, double MeshSize_Y, double roh, double Viscosity)
{

    double dX = 1.0 / MeshSize_X;
    double dY = 1.0 / MeshSize_Y;
    vector<vector<double>> CFL_Cond, CFL_Conv;
    CFL_Cond.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 2));
    CFL_Conv.resize(MeshSize_X + 2, vector<double>(MeshSize_Y + 2));

    for (int i = 1; i <= MeshSize_X; ++i) {
        for (int j = 1; j <= MeshSize_Y; ++j) {
            double abs_U = max(u_Staggered_X_.u_P[i][j],  u_Staggered_X_.u_P[i - 1][j]);
            double abs_V = max(v_Staggered_Y_.v_P[i][j], v_Staggered_Y_.v_P[i][j - 1]);
            CFL_Conv[i][j] = 0.1 * dX / sqrt(abs_U * abs_U + abs_V * abs_V);
            CFL_Cond[i][j] = 0.1 * dX * dX / ( Viscosity / roh ); 
        }
    }

    // Finding Minimum Δt
    double min_delta_t = std::numeric_limits<double>::max();
    
    for (int i = 1; i <= MeshSize_X; ++i) {
        for (int j = 1; j <= MeshSize_Y; ++j) {
            double delta_t = std::min(CFL_Cond[i][j], CFL_Conv[i][j]);
            min_delta_t = std::min(min_delta_t, delta_t);
       }
    }
    return min_delta_t;
}

bool checkSteadyState(const u_Staggered_X& uCurrent, const u_Staggered_X& uPrevious,
    const v_Staggered_Y& vCurrent, const v_Staggered_Y& vPrevious,
    double tolerance, double MeshSize_X, double MeshSize_Y) {
    
    // Convergence Check for U
    for (int i = 0; i < MeshSize_X + 1; ++i) {
        for (int j = 0; j < MeshSize_Y; ++j) {
            if (abs(uCurrent.u_P[i][j] - uPrevious.u_P[i][j]) > tolerance) {
                return false; // Solution hasn't converged yet
            }
        }
    }

    // Convergence Check for V
    for (int i = 0; i < MeshSize_X + 2; ++i) {
        for (int j = 0; j < MeshSize_Y + 1; ++j) {
            if (abs(vCurrent.v_P[i][j] - vPrevious.v_P[i][j]) > tolerance) {
                return false; // Solution hasn't converged yet
            }
        }
    }

    return true; // Solution has converged
}

double massBalance(u_Staggered_X u_Staggered_X_, v_Staggered_Y v_Staggered_Y_, double roh, Area Area_, float MeshSize_X, float MeshSize_Y)
{
    double massBalance = 0.0;

    // Iterate over the interior cells
    for (int i = 1; i < MeshSize_X; ++i) {
        for (int j = 1; j < MeshSize_Y; ++j) {
            // Compute mass fluxes in the x-direction
            double Mass_X_u = roh * (Area_.East[i] * u_Staggered_X_.u_P[i][j] - Area_.West[i] * u_Staggered_X_.u_P[i - 1][j]);
            // Compute mass fluxes in the y-direction
            double Mass_Y_v = roh * (Area_.North[j] * v_Staggered_Y_.v_P[i][j] - Area_.South[j] * v_Staggered_Y_.v_P[i][j - 1]);

            // Accumulate the mass fluxes
            massBalance += Mass_X_u;
            massBalance += Mass_Y_v;
        }
    }
    return massBalance;
}