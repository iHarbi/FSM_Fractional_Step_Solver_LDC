#include <string>

class HRS_Schemes {
public:
    /*Public Accessible Method to Compute the PHI at the Faces*/
    double computeConvectiveValue(double faceMassFlowRate, double XY_e, double XY_P, double XY_EN, double XY_WS, double XY_EE_WW_SS_NN,
        double phi_P, double phi_EN, double phi_WS, double phi_EE_WW_SS_NN, const std::string& Scheme,
        const std::string& Face) {
        
        double Face_Phi;
        double x_D, phi_D, x_C, phi_C, x_U, phi_U, x_e;
        /*Double Check the Logic with Prof. CD*/
        if (faceMassFlowRate >= 0) {
            if (Face == "East") {
                x_D = XY_EN; phi_D = phi_EN; x_C = XY_P; phi_C = phi_P; x_U = XY_WS; phi_U = phi_WS; x_e = XY_e;
            }
            else if (Face == "West") {
                x_D = XY_P; phi_D = phi_P; x_C = XY_WS; phi_C = phi_WS; x_U = XY_EE_WW_SS_NN; phi_U = phi_EE_WW_SS_NN; x_e = XY_e;
            }
            else if (Face == "North") {
                x_D = XY_EN; phi_D = phi_EN; x_C = XY_P; phi_C = phi_P; x_U = XY_WS; phi_U = phi_WS; x_e = XY_e;
            }
            else if (Face == "South") {
                x_D = XY_P; phi_D = phi_P; x_C = XY_WS; phi_C = phi_WS; x_U = XY_EE_WW_SS_NN; phi_U = phi_EE_WW_SS_NN; x_e = XY_e;
            }
        }
        else if (faceMassFlowRate < 0) {
            if (Face == "East") {
                x_D = XY_P; phi_D = phi_P; x_C = XY_EN; phi_C = phi_EN; x_U = XY_EE_WW_SS_NN; phi_U = phi_EE_WW_SS_NN; x_e = XY_e;
            }
            else if (Face == "West") {
                x_D = XY_WS; phi_D = phi_WS; x_C = XY_P; phi_C = phi_P; x_U = XY_EN; phi_U = phi_EN; x_e = XY_e;
            }
            else if (Face == "North") {
                x_D = XY_P; phi_D = phi_P; x_C = XY_EN; phi_C = phi_EN; x_U = XY_EE_WW_SS_NN; phi_U = phi_EE_WW_SS_NN; x_e = XY_e;
            }
            else if (Face == "South") {
                x_D = XY_WS; phi_D = phi_WS; x_C = XY_P; phi_C = phi_P; x_U = XY_EN; phi_U = phi_EN; x_e = XY_e;
            }
        }
        /* Normalization Process */
        double phi_Normal_C = (phi_C - phi_U) / (phi_D - 0.99 * phi_U);
        double x_Normal_C = (x_C - x_U) / (x_D - x_U);
        double x_Normal_e = (x_e - x_U) / (x_D - x_U);
        
        double phi_Normal_e;
        /*Scheme of Choice */
        if (Scheme == "CDS") {
            phi_Normal_e = (x_Normal_e - x_Normal_C) / (1 - x_Normal_C) + (phi_Normal_C * (x_Normal_e - 1) / (x_Normal_C - 1));
        }
        else if (Scheme == "UDS") {
            phi_Normal_e = phi_Normal_C;
        }
        else if (Scheme == "SUDS") {
            phi_Normal_e = phi_Normal_C * (x_Normal_e / x_Normal_C);
        }
        else if (Scheme == "QUICK") {
            phi_Normal_e = x_Normal_e + (x_Normal_e / x_Normal_C) * ((x_Normal_e - 1) / (x_Normal_C - 1)) * (phi_Normal_C - x_Normal_C);
        }
        /*Don't recall why I addded this exclusion */
        if (faceMassFlowRate == 0) {
            phi_Normal_e = 0;
            Face_Phi = phi_U;
        }
        else {
            Face_Phi = phi_Normal_e * (phi_D - phi_U) + phi_U;
        }

        return Face_Phi;
    }
};
