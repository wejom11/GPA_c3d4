#include "P9SF_ele.h"

bool P9SF::is_valid(){
    if(Nodetag[0] != 0 && XY.front() != nullptr && mat != nullptr){
        return true;
    }
    return false;
}

bool P9SF::is_f(){
    return is_fuild;
}

void P9SF::set_fs(bool is_f){
    this->is_fuild = is_f;
    fd_ele = is_f ? 31 : 9;
    if(is_fuild){
        T.resize(9, nullptr);
        UV.resize(9, nullptr);
        P.resize(4, nullptr);
    }
    else{
        T.resize(9, nullptr);
        UV.clear();
        UV.shrink_to_fit();
        P.clear();
        P.shrink_to_fit();
    }
}

void P9SF::set_nval(int* ntag, double* xyz, double* uvtp, int* tn_num){
    int i = 0;
    int fst_p = tn_num[0] * 4;
    int snd_p = tn_num[0] * 4 + tn_num[1] * 3;
    int snd = tn_num[0] + tn_num[1];

    for(i = 0; i < 9; i++){
        Nodetag[i] = ntag[i];
        XY[i] = &xyz[2 * ntag[i] - 2];
        
        if(is_fuild){
            if(i < 4){
                UV.at(i) = &uvtp[(ntag[i] - 1) * 4];
                T.at(i) = &uvtp[(ntag[i] - 1) * 4 + 2];
                P.at(i) = &uvtp[(ntag[i] - 1) * 4 + 3];
            }
            else{
                UV.at(i) = &uvtp[fst_p + (ntag[i] - tn_num[0] - 1) * 3];
                T.at(i) = &uvtp[fst_p + (ntag[i] - tn_num[0] - 1) * 3 + 2];
            }
        }
        else{
            if(ntag[i] <= tn_num[0]){
                T.at(i) = &uvtp[(ntag[i] - 1) * 4 + 3];
            }
            else if(ntag[i] > tn_num[0] && ntag[i] <= snd){
                T.at(i) = &uvtp[fst_p + (ntag[i] - tn_num[0] - 1) * 3 + 2];
            }
            else{
                T.at(i) = &uvtp[snd_p + ntag[i] - snd - 1];
            }
        }
    }
};

void P9SF::get_N_Np(double N[9 * 9], double Np[4 * 9], double weight[9]){
    double Gp[18] = {-sqrt(0.6),-sqrt(0.6),    -sqrt(0.6),0,   -sqrt(0.6),sqrt(0.6),
                     0         ,-sqrt(0.6),    0         ,0,   0         ,sqrt(0.6),
                     sqrt(0.6) ,-sqrt(0.6),    sqrt(0.6) ,0,   sqrt(0.6) ,sqrt(0.6)};
    int kk[9 * 2] = { 0, 0,  2, 0,  2,2, 
                      0, 2,  1, 0,  2,1,
                      1, 2,  0, 1,  1,1};
    int kk_2o[4 * 2] = { 0, 0,  1, 0,  1, 1,  0, 1};
    double wei1d[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};

    int ROW_I = 0, ROW_Ip = 0;
    for (int i = 0; i < 9; i++){
        double xi = Gp[2 * i], eta = Gp[2 * i + 1];
        double xi_pow2 = xi * xi;
        double eta_pow2 = eta * eta;
        ROW_I = 9 * i;
        ROW_Ip = 4 * i;

        double xi_1d_3o[3] = {(xi_pow2 - xi) / 2, 1 - xi_pow2, (xi_pow2 + xi) / 2};
        double eta_1d_3o[3] = {(eta_pow2 - eta) / 2, 1 - eta_pow2, (eta_pow2 + eta) / 2};

        double xi_1d_2o[2] = {(1 - xi)/2, (1 + xi)/2};
        double eta_1d_2o[2] = {(1 - eta)/2, (1 + eta)/2};

        for (int j = 0; j < 9; j++){
            N[ROW_I + j] = xi_1d_3o[kk[2*j]] * eta_1d_3o[kk[2*j + 1]];
        }
        for (int j = 0; j < 4; j++){
            Np[ROW_Ip + j] = xi_1d_2o[kk_2o[2*j]] * eta_1d_2o[kk_2o[2*j + 1]];
        }
    }

    for(int i = 0; i < 9; i++){
        weight[i] = wei1d[i/3] * wei1d[i%3];
    }
};

void P9SF::get_N(double N[9 * 9], double weight[9]){
    double Gp[18] = {-sqrt(0.6),-sqrt(0.6),    -sqrt(0.6),0,   -sqrt(0.6),sqrt(0.6),
                     0         ,-sqrt(0.6),    0         ,0,   0         ,sqrt(0.6),
                     sqrt(0.6) ,-sqrt(0.6),    sqrt(0.6) ,0,   sqrt(0.6) ,sqrt(0.6)};
    int kk[9 * 2] = { 0, 0,  2, 0,  2,2, 
                      0, 2,  1, 0,  2,1,
                      1, 2,  0, 1,  1,1};
    double wei1d[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};

    int ROW_I = 0;
    for (int i = 0; i < 9; i++){
        double xi = Gp[2 * i], eta = Gp[2 * i + 1];
        double xi_pow2 = xi * xi;
        double eta_pow2 = eta * eta;
        ROW_I = 9 * i;

        double xi_1d_3o[3] = {(xi_pow2 - xi) / 2, 1 - xi_pow2, (xi_pow2 + xi) / 2};
        double eta_1d_3o[3] = {(eta_pow2 - eta) / 2, 1 - eta_pow2, (eta_pow2 + eta) / 2};

        for (int j = 0; j < 9; j++){
            N[ROW_I + j] = xi_1d_3o[kk[2*j]] * eta_1d_3o[kk[2*j + 1]];
        }
    }

    for(int i = 0; i < 9; i++){
        weight[i] = wei1d[i/3] * wei1d[i%3];
    }
};

void P9SF::get_Ndl(double Ndl[9 * 9 * 2]){
    double Gp[18] = {-sqrt(0.6),-sqrt(0.6),    -sqrt(0.6),0,   -sqrt(0.6),sqrt(0.6),
                    0         ,-sqrt(0.6),    0         ,0,   0         ,sqrt(0.6),
                    sqrt(0.6) ,-sqrt(0.6),    sqrt(0.6) ,0,   sqrt(0.6) ,sqrt(0.6)};
    int kk[9 * 2] = { 0, 0,  2, 0,  2,2, 
                      0, 2,  1, 0,  2,1,
                      1, 2,  0, 1,  1,1};

    int ROW_I = 0, ROW_Ip = 0;
    for (int i = 0; i < 9; i++){
        double xi = Gp[2 * i], eta = Gp[2 * i + 1];
        double xi_pow2 = xi * xi;
        double eta_pow2 = eta * eta;
        ROW_I = 18 * i;
        ROW_Ip = 8 * i;

        double xi_1d_3o[3] = {(xi_pow2 - xi) / 2, 1 - xi_pow2, (xi_pow2 + xi) / 2};
        double xidl_3o[3] = {xi - 0.5, - 2*xi, xi + 0.5};
        double eta_1d_3o[3] = {(eta_pow2 - eta) / 2, 1 - eta_pow2, (eta_pow2 + eta) / 2};
        double etadl_3o[3] = {eta - 0.5, - 2*eta, eta + 0.5};

        for (int j = 0; j < 9; j++){
            Ndl[ROW_I + 2*j] = xidl_3o[kk[2*j]] * eta_1d_3o[kk[2*j + 1]];
            Ndl[ROW_I + 2*j + 1] = xi_1d_3o[kk[2*j]] * etadl_3o[kk[2*j + 1]];
        }
    }    
};

void P9SF::get_Ndx(double Ndx[9 * 9 * 2], double Det[9]){
    get_Ndl(Ndx);
    double J[4];
    double val = 0.0;

    int info = 0;
    int INTP = 0;
    int ipiv[2]{0};

    for(int i = 0; i < 9; i++){
        INTP = 18 * i;
        for(int j = 0; j < 2; j++){
            for(int k = 0; k < 2; k++){
                val = 0.0;
                for(int m = 0; m < 9; m++){
                    val += Ndx[INTP + 2*m + j] * XY.at(m)[k];
                }

                J[2*j + k] = val;
            }
        }

        Det[i] = abs(J[0] * J[3] - J[1] * J[2]);

        LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 2, 2, J, 2, ipiv);

        for(int j = 0; j < 9; j++){
            info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', 2, 1, J, 2, ipiv, &Ndx[INTP + 2*j], 1);

            if(info != 0){
                printf("ERROR: calculating Ndx failed!\n");
            }
        }
    }
};

void P9SF::get_UVTP_field(double UV_fld[9 * 2], double T_fld[9], double P_fld[9], const double N[9 * 9], const double Np[9 * 4]){
    int i = 0, j = 0, ROW_I = 0, ROW_NI = 0, ROW_NpI;
    double val_u = 0., val_v = 0., val_t = 0., val_p = 0.;

    for(i = 0; i < 9; i++){
        ROW_I = 2 * i;
        ROW_NI = 9 * i;
        ROW_NpI = 4 * i;

        val_u = 0.;
        val_v = 0.;
        val_t = 0.;
        for(j = 0; j < 9; j++){
            val_u += N[ROW_NI + j] * UV.at(j)[0];
            val_v += N[ROW_NI + j] * UV.at(j)[1];
            val_t += N[ROW_NI + j] * *(T.at(j));
        }
        UV_fld[ROW_I] = val_u;
        UV_fld[ROW_I + 1] = val_v;
        T_fld[i] = val_t;

        val_p = 0.;
        for(j = 0; j < 4; j++){
            val_p += Np[ROW_NpI + j] * *(P.at(j));
        }
        P_fld[i] = val_p;
    }
};

void P9SF::get_UVTP_field(double T_fld[9], const double N[9 * 9]){
    int i = 0, j = 0, ROW_I = 0;
    double val_t = 0.;

    for(i = 0; i < 9; i++){
        ROW_I = 9 * i;

        val_t = 0.;
        for(j = 0; j < 9; j++){
            val_t += N[ROW_I + j] * *(T.at(j));
        }
        T_fld[i] = val_t;
    }
}

void P9SF::get_UVTdx_field(double UVdx_fld[9 * 2 * 2], double Tdx_fld[9 * 2], const double Ndx[9 * 9 * 2]){
    int i = 0, j = 0, ROW_I = 0, ROW_NI = 0, ROW_NpI;
    double val_udx = 0., val_vdx = 0., val_tdx = 0., val_pdx = 0.,
           val_udy = 0., val_vdy = 0., val_tdy = 0., val_pdy = 0.;

    for(i = 0; i < 9; i++){
        ROW_I = 2 * 2 * i;
        ROW_NI = 2 * 9 * i;
        ROW_NpI = 2 * 4 * i;

        val_udx = 0.;
        val_udy = 0.;
        val_vdx = 0.;
        val_vdy = 0.;
        val_tdx = 0.;
        val_tdy = 0.;
        for(j = 0; j < 9; j++){
            val_udx += Ndx[ROW_NI + 2*j] * UV.at(j)[0];
            val_udy += Ndx[ROW_NI + 2*j + 1] * UV.at(j)[0];
            val_vdx += Ndx[ROW_NI + 2*j] * UV.at(j)[1];
            val_vdy += Ndx[ROW_NI + 2*j + 1] * UV.at(j)[1];
            val_tdx += Ndx[ROW_NI + 2*j] * *(T.at(j));
            val_tdy += Ndx[ROW_NI + 2*j + 1] * *(T.at(j));
        }
        UVdx_fld[ROW_I] = val_udx;
        UVdx_fld[ROW_I + 1] = val_udy;
        UVdx_fld[ROW_I + 2] = val_vdx;
        UVdx_fld[ROW_I + 3] = val_vdy;
        Tdx_fld[2*i] = val_tdx;
        Tdx_fld[2*i + 1] = val_tdy;
    }
};

void P9SF::get_UVTdx_field(double Tdx_fld[9 * 2], const double Ndx[9 * 9 * 2]){
    int i = 0, j = 0, ROW_I = 0, ROW_NI = 0;
    double val_tdx = 0., val_tdy = 0.;

    for(i = 0; i < 9; i++){
        ROW_I = 2 * i;
        ROW_NI = 2 * 9 * i;

        val_tdx = 0.;
        val_tdy = 0.;
        for(j = 0; j < 9; j++){
            val_tdx += Ndx[ROW_NI + 2*j] * *(T.at(j));
            val_tdx += Ndx[ROW_NI + 2*j + 1] * *(T.at(j));
        }
        Tdx_fld[ROW_I] = val_tdx;
        Tdx_fld[ROW_I + 1] = val_tdy;
    }
};

void P9SF::get_Fint_f(double F[31], const double N[9 * 9], const double UV_fld[9 * 2], const double UVdx_fld[9 * 2 * 2],
                      const double T_fld[9], const double Tdx_fld[9 * 2], const double P_fld[9], const double Np[9 * 4], 
                      const double Ndx[9 * 9 * 2], const double Det[9], const double weight[9]){
    int i = 0, j = 0, ROW_N = 0, /* ROW_Np = 0, */ I2, I4, ROW_Ndx;
    double rho = this->mat->rho;
    double alpha = this->mat->alpha;
    double c1 = this->mat->c1;
    double c2 = this->mat->c2;
    double beta = this->mat->beta;
    double T0 = this->mat->T0;

    double conv_Udu = 0., conv_Udv = 0., conv_Udt = 0.,
           div = 0., vis = 0., dU = 0.;

    for(i = 0; i < 31; i++){
        F[i] = 0.;
    }

    for(i = 0; i < 9; i++){ // Guass integration points
        ROW_N = 9 * i;
        ROW_Ndx = 18 * i;
        I2 = 2 * i;
        I4 = 4 * i;

        conv_Udu = UV_fld[I2]*UVdx_fld[I4] + UV_fld[I2 + 1]*UVdx_fld[I4 + 1];       // u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y}
        conv_Udv = UV_fld[I2]*UVdx_fld[I4 + 2] + UV_fld[I2 + 1]*UVdx_fld[I4 + 3];   // u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y}
        conv_Udt = UV_fld[I2]*Tdx_fld[I2] + UV_fld[I2 + 1]*Tdx_fld[I2 + 1];         // u\frac{\partial T}{\partial x} + v\frac{\partial T}{\partial y}
        div = UVdx_fld[I4 + 1] + UVdx_fld[I4 + 2];                                  // \frac{\partial u}{\partail y} + \frac{\partial v}{\partial x}
        dU = UVdx_fld[I4] + UVdx_fld[I4 + 3];                                       // \frac{\partial u}{\partail x} + \frac{\partial v}{\partial y}
        vis = c1 * T_fld[i] + c2;                                                   // viscosity
        for(j = 0; j < 4; j++){
            // U, V
            F[4*j] += (rho * N[ROW_N + j] * conv_Udu + vis * (2*UVdx_fld[I4]*Ndx[ROW_Ndx + 2*j] + Ndx[ROW_Ndx + 2*j + 1]*div) -
                       P_fld[i] * Ndx[ROW_Ndx + 2*j]) * weight[i] * Det[i];
            F[4*j + 1] += (rho * N[ROW_N + j] * conv_Udv + vis * (2*UVdx_fld[I4 + 3]*Ndx[ROW_Ndx + 2*j + 1] + Ndx[ROW_Ndx + 2*j]*div) -
                           P_fld[i] * Ndx[ROW_Ndx + 2*j + 1]) * weight[i] * Det[i];
            // T
            F[4*j + 2] += weight[i] * Det[i] * ((conv_Udt + beta*(T_fld[i] - T0)) * N[ROW_N + j] +
                          alpha*(Tdx_fld[I2]*Ndx[ROW_Ndx + 2*j] + Tdx_fld[I2 + 1]*Ndx[ROW_Ndx + 2*j + 1]));
            // P
            F[4*j + 3] += Np[I4 + j] * dU * weight[i] * Det[i];
        }

        for(j = 4; j < 9; j++){
            // U, V
            F[4 + 3*j] += (rho * N[ROW_N + j] * conv_Udu + vis * (2*UVdx_fld[I4]*Ndx[ROW_Ndx + 2*j] + Ndx[ROW_Ndx + 2*j + 1]*div) -
                       P_fld[i] * Ndx[ROW_Ndx + 2*j]) * weight[i] * Det[i];
            F[5 + 3*j] += (rho * N[ROW_N + j] * conv_Udv + vis * (2*UVdx_fld[I4 + 3]*Ndx[ROW_Ndx + 2*j + 1] + Ndx[ROW_Ndx + 2*j]*div) -
                       P_fld[i] * Ndx[ROW_Ndx + 2*j + 1]) * weight[i] * Det[i];
            // T
            F[6 + 3*j] += weight[i] * Det[i] * ((conv_Udt + beta*(T_fld[i] - T0)) * N[ROW_N + j] +
                           alpha*(Tdx_fld[I2]*Ndx[ROW_Ndx + 2*j] + Tdx_fld[I2 + 1]*Ndx[ROW_Ndx + 2*j + 1]));
        }
    }
};

void P9SF::get_K_f(double K[31][31], const double N[9 * 9], const double UV_fld[9 * 2], const double UVdx_fld[9 * 2 * 2],
                      const double T_fld[9], const double Tdx_fld[9 * 2], const double P_fld[9], const double Np[9 * 4], 
                      const double Ndx[9 * 9 * 2], const double Det[9], const double weight[9]){
    int i = 0, j = 0, k = 0, m = 0, n = 0,
        ROW_N = 0, /* ROW_Np = 0, */ I2, I4, ROW_Ndx,
        J4, K4, J3, K3;
    double rho = this->mat->rho;
    double alpha = this->mat->alpha;
    double c1 = this->mat->c1;  
    double c2 = this->mat->c2;
    double beta = this->mat->beta;

    double vis = 0.;
    for(i = 0; i < 31; i++){
        for(j = 0; j < 31; j++){
            K[i][j] = 0.;
        }
    }

    for(i = 0; i < 9; i++){ // Gauss integration points
        ROW_N = 9 * i;
        I4 = 4 * i;
        ROW_Ndx = 18 * i;
        I2 = 2 * i;

        vis = c1 * T_fld[i] + c2;                                                   // viscosity
        for(j = 0; j < 4; j++){
            J4 = 4 * j;

            for(k = 0; k < 4; k++){
                K4 = 4 * k;
                // K (j-th, k-th) block
                for(m = 0; m < 2; m++){
                    for(n = 0; n < 2; n++){
                        if(m == n){
                            K[J4+m][K4+n] += (rho * N[ROW_N + j] * (N[ROW_N + k]*UVdx_fld[I4 + 3*m] + UV_fld[I2]*Ndx[ROW_Ndx + 2*k] +
                                        UV_fld[I2 + 1] * Ndx[ROW_Ndx + 2*k + 1]) + vis*(2*Ndx[ROW_Ndx + 2*k + m]*Ndx[ROW_Ndx + 2*j + m] + 
                                        Ndx[ROW_Ndx + 2*k + 1 - m]*Ndx[ROW_Ndx + 2*j + 1 - m])) * weight[i] * Det[i];
                        }
                        else{
                            K[J4+m][K4+n] += (rho * N[ROW_N + j] * N[ROW_N + k] * UVdx_fld[I4 + m + 1] + 
                                              vis*Ndx[ROW_Ndx + 2*j + n]*Ndx[ROW_Ndx + 2*k + 1-n]) * weight[i] * Det[i];
                        }
                    }
                }

                K[J4+2][K4] += (N[ROW_N + j] * N[ROW_N + k] * Tdx_fld[I2]) * weight[i] * Det[i];
                K[J4+2][K4 + 1] += (N[ROW_N + j] * N[ROW_N + k] * Tdx_fld[I2 + 1]) * weight[i] * Det[i];
                K[J4+2][K4 + 2] += ((UV_fld[I2]*Ndx[ROW_Ndx + 2*k] + UV_fld[I2 + 1]*Ndx[ROW_Ndx + 2*k + 1] +
                                    beta*N[ROW_N + k]) * N[ROW_N + j] + alpha * (Ndx[ROW_Ndx + 2*j] * Ndx[ROW_Ndx + 2*k] +
                                    Ndx[ROW_Ndx + 2*j + 1] * Ndx[ROW_Ndx + 2*k + 1])) * weight[i] * Det[i];
                
                K[J4+3][K4] += (Np[I4 + j] * Ndx[ROW_Ndx + 2*k]) * weight[i] * Det[i];
                K[J4+3][K4 + 1] += (Np[I4 + j] * Ndx[ROW_Ndx + 2*k + 1]) * weight[i] * Det[i];
                K[J4][K4+3] += -(Np[I4 + k] * Ndx[ROW_Ndx + 2*j]) * weight[i] * Det[i];
                K[J4+1][K4+3] += -(Np[I4 + k] * Ndx[ROW_Ndx + 2*j + 1]) * weight[i] * Det[i];

            }

            for(k = 4; k < 9; k++){
                K3 = 3 * k;

                for(m = 0; m < 2; m++){
                    for(n = 0; n < 2; n++){
                        if(m == n){
                            K[J4+m][K3+n+4] += (rho * N[ROW_N + j] * (N[ROW_N + k]*UVdx_fld[I4 + 3*m] + UV_fld[I2]*Ndx[ROW_Ndx + 2*k] +
                                        UV_fld[I2 + 1] * Ndx[ROW_Ndx + 2*k + 1]) + vis*(2*Ndx[ROW_Ndx + 2*k + m]*Ndx[ROW_Ndx + 2*j + m] + 
                                        Ndx[ROW_Ndx + 2*k + 1 - m]*Ndx[ROW_Ndx + 2*j + 1 - m])) * weight[i] * Det[i];
                        }
                        else{
                            K[J4+m][K3+n+4] += (rho * N[ROW_N + j] * N[ROW_N + k] * UVdx_fld[I4 + m + 1] + 
                                              vis*Ndx[ROW_Ndx + 2*j + n]*Ndx[ROW_Ndx + 2*k + 1-n]) * weight[i] * Det[i];
                        }
                    }
                }

                K[J4+2][K3+4] += (N[ROW_N + j] * N[ROW_N + k] * Tdx_fld[I2]) * weight[i] * Det[i];
                K[J4+2][K3+5] += (N[ROW_N + j] * N[ROW_N + k] * Tdx_fld[I2 + 1]) * weight[i] * Det[i];
                K[J4+2][K3+6] += ((UV_fld[I2]*Ndx[ROW_Ndx + 2*k] + UV_fld[I2 + 1]*Ndx[ROW_Ndx + 2*k + 1] +
                                    beta*N[ROW_N + k]) * N[ROW_N + j] + alpha * (Ndx[ROW_Ndx + 2*j] * Ndx[ROW_Ndx + 2*k] +
                                    Ndx[ROW_Ndx + 2*j + 1] * Ndx[ROW_Ndx + 2*k + 1])) * weight[i] * Det[i];
                
                K[J4+3][K3+4] += (Np[I4 + j] * Ndx[ROW_Ndx + 2*k]) * weight[i] * Det[i];
                K[J4+3][K3+5] += (Np[I4 + j] * Ndx[ROW_Ndx + 2*k + 1]) * weight[i] * Det[i];

            }
        }

        for(j = 4; j < 9; j++){
            J3 = 3 * j;
            for(k = 0; k < 4; k++){
                K4 = 4 * k;
                for(m = 0; m < 2; m++){
                    for(n = 0; n < 2; n++){
                        if(m == n){
                            K[4+J3+m][K4+n] += (rho * N[ROW_N + j] * (N[ROW_N + k]*UVdx_fld[I4 + 3*m] + UV_fld[I2]*Ndx[ROW_Ndx + 2*k] +
                                        UV_fld[I2 + 1] * Ndx[ROW_Ndx + 2*k + 1]) + vis*(2*Ndx[ROW_Ndx + 2*k + m]*Ndx[ROW_Ndx + 2*j + m] + 
                                        Ndx[ROW_Ndx + 2*k + 1 - m]*Ndx[ROW_Ndx + 2*j + 1 - m])) * weight[i] * Det[i];
                        }
                        else{
                            K[4+J3+m][K4+n] += (rho * N[ROW_N + j] * N[ROW_N + k] * UVdx_fld[I4 + m + 1] + 
                                              vis*Ndx[ROW_Ndx + 2*j + n]*Ndx[ROW_Ndx + 2*k + 1-n]) * weight[i] * Det[i];
                        }
                    }
                }

                K[J3+6][K4] += (N[ROW_N + j] * N[ROW_N + k] * Tdx_fld[I2]) * weight[i] * Det[i];
                K[J3+6][K4+1] += (N[ROW_N + j] * N[ROW_N + k] * Tdx_fld[I2 + 1]) * weight[i] * Det[i];
                K[J3+6][K4+2] += ((UV_fld[I2]*Ndx[ROW_Ndx + 2*k] + UV_fld[I2 + 1]*Ndx[ROW_Ndx + 2*k + 1] +
                                    beta*N[ROW_N + k]) * N[ROW_N + j] + alpha * (Ndx[ROW_Ndx + 2*j] * Ndx[ROW_Ndx + 2*k] +
                                    Ndx[ROW_Ndx + 2*j + 1] * Ndx[ROW_Ndx + 2*k + 1])) * weight[i] * Det[i];                

                K[J3+4][K4 + 3] += -(Np[I4 + k] * Ndx[ROW_Ndx + 2*j]) * weight[i] * Det[i];
                K[J3+5][K4 + 3] += -(Np[I4 + k] * Ndx[ROW_Ndx + 2*j + 1]) * weight[i] * Det[i];
            }

            for(k = 4; k < 9; k++){
                K3 = 3 * k;
                // K (j-th, k-th) block
                for(m = 0; m < 2; m++){
                    for(n = 0; n < 2; n++){
                        if(m == n){
                            K[4+J3+m][4+K3+n] += (rho * N[ROW_N + j] * (N[ROW_N + k]*UVdx_fld[I4 + 3*m] + UV_fld[I2]*Ndx[ROW_Ndx + 2*k] +
                                                  UV_fld[I2 + 1] * Ndx[ROW_Ndx + 2*k + 1]) + vis*(2*Ndx[ROW_Ndx + 2*k + m]*Ndx[ROW_Ndx + 2*j + m] + 
                                                  Ndx[ROW_Ndx + 2*k + 1 - m]*Ndx[ROW_Ndx + 2*j + 1 - m])) * weight[i] * Det[i];
                        }
                        else{
                            K[4+J3+m][4+K3+n] += (rho * N[ROW_N + j] * N[ROW_N + k] * UVdx_fld[I4 + m + 1] + 
                                                  vis*Ndx[ROW_Ndx + 2*j + n]*Ndx[ROW_Ndx + 2*k + 1-n]) * weight[i] * Det[i];
                        }
                    }
                }

                K[J3+6][K3+4] += (N[ROW_N + j] * N[ROW_N + k] * Tdx_fld[I2]) * weight[i] * Det[i];
                K[J3+6][K3+5] += (N[ROW_N + j] * N[ROW_N + k] * Tdx_fld[I2 + 1]) * weight[i] * Det[i];
                K[J3+6][K3+6] += ((UV_fld[I2]*Ndx[ROW_Ndx + 2*k] + UV_fld[I2 + 1]*Ndx[ROW_Ndx + 2*k + 1] +
                                    beta*N[ROW_N + k]) * N[ROW_N + j] + alpha * (Ndx[ROW_Ndx + 2*j] * Ndx[ROW_Ndx + 2*k] +
                                    Ndx[ROW_Ndx + 2*j + 1] * Ndx[ROW_Ndx + 2*k + 1])) * weight[i] * Det[i];
            }
        }
    }

};

void P9SF::get_Fint_s(double F[9], const double N[9 * 9], const double T_fld[9], const double Tdx_fld[9 * 2], 
                   const double Ndx[9 * 9 * 2], const double Det[9], const double weight[9]){
    int i = 0, j = 0 ,ROW_N = 0, ROW_Ndx = 0, I2 = 0;

    double alpha = this->mat->alpha;
    double beta = this->mat->beta;
    double T0 = this->mat->T0;

    for(i = 0; i < 9; i++){
        F[i] = 0.;
    }

    for(i = 0; i < 9; i++){ // Gauss integration points
        ROW_N = 9 * i;
        ROW_Ndx = 18 * i;
        I2 = 2 * i;

        for(j = 0; j < 9; j++){
            F[j] += (beta * N[ROW_N + j] * (T_fld[i] - T0) + alpha * (Tdx_fld[I2]*Ndx[ROW_Ndx + 2*j] + 
                   Tdx_fld[I2 + 1]*Ndx[ROW_Ndx + 2*j + 1])) * weight[i] * Det[i];
        }
    }
}

void P9SF::get_K_s(double K[9][9], const double N[9 * 9], const double T_fld[9], const double Tdx_fld[9 * 2], 
                   const double Ndx[9 * 9 * 2], const double Det[9], const double weight[9]){
    int i = 0,     j = 0,       k = 0,
        ROW_N = 0, ROW_Ndx = 0, I2 = 0;

    double alpha = this->mat->alpha;
    double beta = this->mat->beta;
    double T0 = this->mat->T0;

    for(i = 0; i < 9; i++){
        for(j = 0; j < 9; j++){
            K[i][j] = 0.;
        }
    }

    for(i = 0; i < 9; i++){ // Gauss integration points
        ROW_N = 9 * i;
        ROW_Ndx = 18 * i;

        for(j = 0; j < 9; j++){
            for(k = 0; k < 9; k++){
                K[j][k] += weight[i] * Det[i] *(beta * N[ROW_N + j] * N[ROW_N + k] + 
                           alpha * (Ndx[ROW_Ndx + 2*j]*Ndx[ROW_Ndx + 2*k] + Ndx[ROW_Ndx + 2*j + 1]*Ndx[ROW_Ndx + 2*k + 1]));
            }
        }
    }

}

void P9SF::KF_f(double K[31][31], double F[31]){
    double N[9*9]{0},          Np[9*4]{0},      Ndx[9*9*2]{0},
           uv_fld[9*2]{0},     t_fld[9]{0},     p_fld[9]{0},
           uvdx_fld[9*2*2]{0}, tdx_fld[9*2]{0}, pdx_fld[9*2]{0},
           determinate[9]{0},  weightness[9]{0};

    get_N_Np(N,Np,weightness);
    get_Ndx(Ndx, determinate);
    get_UVTP_field(uv_fld, t_fld, p_fld, N, Np);
    get_UVTdx_field(uvdx_fld, tdx_fld, Ndx);

    get_Fint_f(F, N, uv_fld, uvdx_fld, t_fld, tdx_fld, p_fld, Np, Ndx,
               determinate, weightness);
    get_K_f(K, N, uv_fld, uvdx_fld, t_fld, tdx_fld, p_fld, Np, Ndx,
            determinate, weightness);

};

void P9SF::KF_s(double K[9][9], double F[9]){
    double N[9*9]{0},       Ndx[9*9*2]{0},      t_fld[9]{0},
           tdx_fld[9*2]{0}, determinate[9]{0},  weightness[9]{0};

    get_N(N, weightness);
    get_Ndx(Ndx, determinate);
    get_UVTP_field(t_fld, N);
    get_UVTdx_field(tdx_fld, Ndx);

    get_Fint_s(F ,N, t_fld, tdx_fld, Ndx, determinate, weightness);
    get_K_s(K, N, t_fld, tdx_fld, Ndx, determinate, weightness);

}

int P9SF::at_Nodetag(int index){
    return Nodetag[index];
};
