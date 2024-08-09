#include "c3d4_ele.h"

double c3d4::getvolume(){
    double volu = 0.0;
    volu = (xyz0[2][0]*xyz0[1][1]*xyz0[0][2] - xyz0[3][0]*xyz0[1][1]*xyz0[0][2] - xyz0[1][0]*xyz0[2][1]*xyz0[0][2] + xyz0[3][0]*xyz0[2][1]*xyz0[0][2] + 
    xyz0[1][0]*xyz0[3][1]*xyz0[0][2] - xyz0[2][0]*xyz0[3][1]*xyz0[0][2] - xyz0[2][0]*xyz0[0][1]*xyz0[1][2] + xyz0[3][0]*xyz0[0][1]*xyz0[1][2] + 
    xyz0[0][0]*xyz0[2][1]*xyz0[1][2] - xyz0[3][0]*xyz0[2][1]*xyz0[1][2] - xyz0[0][0]*xyz0[3][1]*xyz0[1][2] + xyz0[2][0]*xyz0[3][1]*xyz0[1][2] + 
    xyz0[1][0]*xyz0[0][1]*xyz0[2][2] - xyz0[3][0]*xyz0[0][1]*xyz0[2][2] - xyz0[0][0]*xyz0[1][1]*xyz0[2][2] + xyz0[3][0]*xyz0[1][1]*xyz0[2][2] + 
    xyz0[0][0]*xyz0[3][1]*xyz0[2][2] - xyz0[1][0]*xyz0[3][1]*xyz0[2][2] - xyz0[1][0]*xyz0[0][1]*xyz0[3][2] + xyz0[2][0]*xyz0[0][1]*xyz0[3][2] + 
    xyz0[0][0]*xyz0[1][1]*xyz0[3][2] - xyz0[2][0]*xyz0[1][1]*xyz0[3][2] - xyz0[0][0]*xyz0[2][1]*xyz0[3][2] + xyz0[1][0]*xyz0[2][1]*xyz0[3][2])/6;
    volu = abs(volu);
    return volu;
};

void c3d4::getK_e(double K[][12], double KdE[][12]){
    double depsilon[3][3][12];
    double dsigdudE[3][3][12];
    double dsigdu[3][3][12];
    getdSigdu(depsilon,dsigdudE,dsigdu);

    for (int i = 0; i < 12; i++){
        for (int j = 0; j < 12; j++){
            double val = 0.0;
            double valdE = 0.0;
            for (int m = 0; m < 3; m++){
                for (int n = 0; n < 3; n++){
                    val += dsigdu[m][n][j] * depsilon[m][n][i];
                    valdE += dsigdudE[m][n][j] * depsilon[m][n][i];
                }
            }
            K[i][j] = val * vol;
            KdE[i][j] = valdE * vol;
        }
    }
};

void c3d4::getK_e(double K[][12]){
    double depsilon[3][3][12];
    double dsigdu[3][3][12];
    getdSigdu(depsilon, dsigdu);

    for (int i = 0; i < 12; i++){
        for (int j = 0; j < 12; j++){
            double val = 0.0;
            for (int m = 0; m < 3; m++){
                for (int n = 0; n < 3; n++){
                    val += dsigdu[m][n][j] * depsilon[m][n][i];
                }
            }
            K[i][j] = val * vol;
        }
    }
};

void c3d4::getdSigdu(double deps[][3][12], double sigdudE[][3][12], double sigdu[][3][12]){
    double lamda = mat->E * mat->mu / (1 + mat->mu) / (1 - 2*mat->mu);
    double Nu = mat->E / 2 / (1 + mat->mu);
    double lamdadE = mat->mu / (1 + mat->mu) / (1 - 2 * mat->mu);
    double NudE = 0.5 / (1 + mat->mu);

    double uldu[3][3][12] = {{{1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 
        0, -1, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0}}, {{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0}, 
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0}}, {{0, 0, 1, 0, 0, 0, 0, 
        0, 0, 0, 0, -1}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1}}};
    double dulddu[3][3][12] = {{{1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 
        0, -1, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0}}, {{0, 1, 0, 0,0, 0, 0, 0, 0, 0, -1, 0}, 
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0}}, {{0, 0, 1, 0, 0, 0, 0, 
        0, 0, 0, 0, -1}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1}}};

    double Xl[9] = {xyz0[0][0] - xyz0[3][0], xyz0[1][0] - xyz0[3][0], 
        xyz0[2][0] - xyz0[3][0], xyz0[0][1] - xyz0[3][1], xyz0[1][1] - xyz0[3][1], 
        xyz0[2][1] - xyz0[3][1], xyz0[0][2] - xyz0[3][2], xyz0[1][2] - xyz0[3][2], 
        xyz0[2][2] - xyz0[3][2]};
    int ipiv[3];
    // LU factorization
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 3, 3, Xl, 3, ipiv);
    // Compute the inverse
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, 3, Xl, 3, ipiv);

    // for (int i = 0; i < 3; i++){
    //     printf("%.4f ",Xl(0,i));
    // }
    double uxdu[3][3][12], duxddu[3][3][12];
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 12; k++){
                double val = 0.0;
                double dval = 0.0;
                for (int m = 0; m < 3; m++){
                    val += uldu[i][m][k] * Xl[3*m + j];
                    dval += dulddu[i][m][k] * Xl[3*m + j];
                }
                uxdu[i][j][k] = val;
                duxddu[i][j][k] = dval;
            }   
        }
    }
    double uddu[3][3][12];
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 12; k++){
                uddu[i][j][k] = (uxdu[i][j][k] + uxdu[j][i][k]) / 2.0;
                deps[i][j][k] = (duxddu[i][j][k] + duxddu[j][i][k]) / 2.0;
            }   
        }
    }

    for (int k = 0; k < 12; k++){
        double Truddu = 0.0;
        double TruddudE = 0.0;
        for (int i = 0; i < 3; i++){
            Truddu += uddu[i][i][k];
        }
        TruddudE = Truddu * lamdadE;
        Truddu = Truddu * lamda;
        
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                sigdu[i][j][k] = 2*Nu*uddu[i][j][k];
                sigdudE[i][j][k] = 2*NudE*uddu[i][j][k];
                if (i == j){
                    sigdu[i][j][k] += Truddu;
                    sigdudE[i][j][k] += TruddudE;
                }
            }   
        }
    }
};

void c3d4::getdSigdu(double deps[][3][12], double sigdu[][3][12]){
    double lamda = mat->E * mat->mu / (1 + mat->mu) / (1 - 2*mat->mu);
    double Nu = mat->E / 2 / (1 + mat->mu);

    double uldu[3][3][12] = {{{1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 
        0, -1, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0}}, {{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0}, 
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0}}, {{0, 0, 1, 0, 0, 0, 0, 
        0, 0, 0, 0, -1}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1}}};
    double dulddu[3][3][12] = {{{1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 
        0, -1, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0}}, {{0, 1, 0, 0,0, 0, 0, 0, 0, 0, -1, 0}, 
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0}}, {{0, 0, 1, 0, 0, 0, 0, 
        0, 0, 0, 0, -1}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1}}};
    
    double Xl[9] = {xyz0[0][0] - xyz0[3][0], xyz0[1][0] - xyz0[3][0], 
        xyz0[2][0] - xyz0[3][0], xyz0[0][1] - xyz0[3][1], xyz0[1][1] - xyz0[3][1], 
        xyz0[2][1] - xyz0[3][1], xyz0[0][2] - xyz0[3][2], xyz0[1][2] - xyz0[3][2], 
        xyz0[2][2] - xyz0[3][2]};
    int ipiv[3];
    // LU factorization
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 3, 3, Xl, 3, ipiv);
    // Compute the inverse
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, 3, Xl, 3, ipiv);
    // for (int i = 0; i < 3; i++){
    //     printf("%.4f ",Xl(0,i));
    // }
    double uxdu[3][3][12], duxddu[3][3][12];
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 12; k++){
                double val = 0.0;
                double dval = 0.0;
                for (int m = 0; m < 3; m++){
                    val += uldu[i][m][k] * Xl[3*m + j];
                    dval += dulddu[i][m][k] * Xl[3*m + j];
                }
                uxdu[i][j][k] = val;
                duxddu[i][j][k] = dval;
            }   
        }
    }
    double uddu[3][3][12];
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 12; k++){
                uddu[i][j][k] = (uxdu[i][j][k] + uxdu[j][i][k]) / 2.0;
                deps[i][j][k] = (duxddu[i][j][k] + duxddu[j][i][k]) / 2.0;
            }   
        }
    }

    for (int k = 0; k < 12; k++){
        double Truddu = 0.0;
        for (int i = 0; i < 3; i++){
            Truddu += uddu[i][i][k];
        }
        Truddu = Truddu * lamda;
        
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                sigdu[i][j][k] = 2*Nu*uddu[i][j][k];
                if (i == j){
                    sigdu[i][j][k] += Truddu;
                }
            }   
        }
    }
};