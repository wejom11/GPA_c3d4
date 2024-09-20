#pragma once
#include <Eigen/Eigen>

// 快速傅里叶变换及辅助函数
namespace EERAHelper {
constexpr double PI = 3.1415926535898;
using Eigen::RowVectorXcd;
using Eigen::RowVectorXd;
using Eigen::VectorXcd;
using Eigen::MatrixXcd;
inline void swap(double &a, double &b)
{
    double t;
    t = a;
    a = b;
    b = t;
}

void bitrp(double xreal[], double ximag[], int n)
{
    // 位反转置换 Bit-reversal Permutation
    int i, j, a, b, p;

    for (i = 1, p = 0; i < n; i *= 2)
    {
        p++;
    }
    for (i = 0; i < n; i++)
    {
        a = i;
        b = 0;
        for (j = 0; j < p; j++)
        {
            b = (b << 1) + (a & 1);    // b = b * 2 + a % 2;
            a >>= 1;        // a = a / 2;
        }
        if (b > i)
        {
            swap(xreal[i], xreal[b]);
            swap(ximag[i], ximag[b]);
        }
    }
}

void FFT(double xreal[], double ximag[], int n)
{
    // 快速傅立叶变换，将复数 x 变换后仍保存在 x 中，xreal, ximag 分别是 x 的实部和虚部
    std::vector<double> wreal, wimag;
    wreal.resize(n * 2);
    wimag.resize(n * 2);
    double treal, timag, ureal, uimag, arg;
    int m, k, j, t, index1, index2;

    bitrp(xreal, ximag, n);

    // 计算 1 的前 n / 2 个 n 次方根的共轭复数 W'j = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
    arg = -2 * PI / n;
    treal = cos(arg);
    timag = sin(arg);
    wreal[0] = 1.0;
    wimag[0] = 0.0;
    for (j = 1; j < n / 2; j++)
    {
        wreal[j] = wreal[j - 1] * treal - wimag[j - 1] * timag;
        wimag[j] = wreal[j - 1] * timag + wimag[j - 1] * treal;
    }

    for (m = 2; m <= n; m *= 2)
    {
        for (k = 0; k < n; k += m)
        {
            for (j = 0; j < m / 2; j++)
            {
                index1 = k + j;
                index2 = index1 + m / 2;
                t = n * j / m;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
                treal = wreal[t] * xreal[index2] - wimag[t] * ximag[index2];
                timag = wreal[t] * ximag[index2] + wimag[t] * xreal[index2];
                ureal = xreal[index1];
                uimag = ximag[index1];
                xreal[index1] = ureal + treal;
                ximag[index1] = uimag + timag;
                xreal[index2] = ureal - treal;
                ximag[index2] = uimag - timag;
            }
        }
    }
}

void  IFFT(double xreal[], double ximag[], int n)
{
    const int N = 1024;
    // 快速傅立叶逆变换
    std::vector<double> wreal, wimag;
    wreal.resize(n * 2);
    wimag.resize(n * 2);
    double treal, timag, ureal, uimag, arg;
    int m, k, j, t, index1, index2;

    bitrp(xreal, ximag, n);

    // 计算 1 的前 n / 2 个 n 次方根 Wj = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
    arg = 2 * PI / n;
    treal = cos(arg);
    timag = sin(arg);
    wreal[0] = 1.0;
    wimag[0] = 0.0;
    for (j = 1; j < n / 2; j++)
    {
        wreal[j] = wreal[j - 1] * treal - wimag[j - 1] * timag;
        wimag[j] = wreal[j - 1] * timag + wimag[j - 1] * treal;
    }

    for (m = 2; m <= n; m *= 2)
    {
        for (k = 0; k < n; k += m)
        {
            for (j = 0; j < m / 2; j++)
            {
                index1 = k + j;
                index2 = index1 + m / 2;
                t = n * j / m;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
                treal = wreal[t] * xreal[index2] - wimag[t] * ximag[index2];
                timag = wreal[t] * ximag[index2] + wimag[t] * xreal[index2];
                ureal = xreal[index1];
                uimag = ximag[index1];
                xreal[index1] = ureal + treal;
                ximag[index1] = uimag + timag;
                xreal[index2] = ureal - treal;
                ximag[index2] = uimag - timag;
            }
        }
    }

    for (j = 0; j < n; j++)
    {
        xreal[j] /= n;
        ximag[j] /= n;
    }
}

void  IFFT_Symmetric(double xreal[], double ximag[], int n)
{
    ximag[0] = 0;
    for (size_t i = 1; i < n / 2; i++)
    {
        xreal[n - i] = xreal[i];
        ximag[n - i] = -ximag[i];
    }
    if ((n - 1) % 2 == 1) //奇数需要处理
    {
        ximag[(n - 1) / 2 + 1] = 0;
    }
    const int N = 1024;
    // 快速傅立叶逆变换
    std::vector<double> wreal, wimag;
    wreal.resize(n * 2);
    wimag.resize(n * 2);
    double treal, timag, ureal, uimag, arg;
    int m, k, j, t, index1, index2;

    bitrp(xreal, ximag, n);

    // 计算 1 的前 n / 2 个 n 次方根 Wj = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
    arg = 2 * PI / n;
    treal = cos(arg);
    timag = sin(arg);
    wreal[0] = 1.0;
    wimag[0] = 0.0;
    for (j = 1; j < n / 2; j++)
    {
        wreal[j] = wreal[j - 1] * treal - wimag[j - 1] * timag;
        wimag[j] = wreal[j - 1] * timag + wimag[j - 1] * treal;
    }

    for (m = 2; m <= n; m *= 2)
    {
        for (k = 0; k < n; k += m)
        {
            for (j = 0; j < m / 2; j++)
            {
                index1 = k + j;
                index2 = index1 + m / 2;
                t = n * j / m;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
                treal = wreal[t] * xreal[index2] - wimag[t] * ximag[index2];
                timag = wreal[t] * ximag[index2] + wimag[t] * xreal[index2];
                ureal = xreal[index1];
                uimag = ximag[index1];
                xreal[index1] = ureal + treal;
                ximag[index1] = uimag + timag;
                xreal[index2] = ureal - treal;
                ximag[index2] = uimag - timag;
            }
        }
    }

    for (j = 0; j < n; j++)
    {
        xreal[j] /= n;
        ximag[j] /= n;
    }
}

// 共轭对称处理, 用于ifft
void ConjSym(RowVectorXcd& vin)
{
    auto n = vin.size();
    if(n == 0) return;
    vin(0) = {vin(0).real(), 0};
    for (int i = 1; i < n / 2; i++)
        vin(n-i) = {vin(i).real(), -vin(i).imag()};
    if ((n-1)%2 == 1) //奇数需要处理
    {
        vin((n-1)/2+1)=0;
    }
}

//MatrixXcd load(const char* filename)
//{
//    std::vector<double> cc;
//    FILE* fp = fopen(filename, "r");
//    while (!feof(fp))
//    {
//        float val;
//        fscanf(fp, "%f", &val);
//        cc.push_back(val);
//    }
//    fclose(fp);
//    int n = cc.size();
//    MatrixXcd tmp(1, n);
//    for (size_t i = 0; i < n; i++)
//    {
//        tmp(0, i) = std::complex<double>(cc[i], 0);
//    }
//    return tmp;
//}

MatrixXcd fft(MatrixXcd const& vin)
{
    std::vector<double> reals;
    std::vector<double> imgs;
    for (int i = 0; i < vin.size(); i++)
    {
        auto tmpdata = vin(0, i);
        reals.push_back(tmpdata.real());
        imgs.push_back(tmpdata.imag());
    }
    FFT(reals.data(), imgs.data(), vin.size());
    MatrixXcd vout(vin.rows(), vin.cols());
    for (int i = 0; i < vin.size(); i++)
    {
        vout(0, i) = std::complex<double>(reals[i], imgs[i]);
    }
    return vout;
}
MatrixXcd ifft(MatrixXcd const& vin)
{
    std::vector<double> reals;
    std::vector<double> imgs;
    for (int i = 0; i < vin.size(); i++)
    {
        auto tmpdata = vin(0, i);
        reals.push_back(tmpdata.real());
        imgs.push_back(tmpdata.imag());
    }
    IFFT_Symmetric(reals.data(), imgs.data(), vin.size());
    MatrixXcd vout(vin.rows(), vin.cols());
    for (int i = 0; i < vin.size(); i++)
    {
        vout(0, i) = std::complex<double>(reals[i], imgs[i]);
    }
    return vout;
}

RowVectorXcd Integration(const RowVectorXcd& vin, double dt)
{
    RowVectorXcd vout;
    vout.resize(vin.size());
    vout(0) = 0;
    for(int i = 1; i < vout.size(); i++)
        vout(i) = vout(i-1) + vin(i-1)*dt;
    return vout;
}

double Mod(const std::complex<double> v) { return sqrt(v.real()*v.real()+v.imag()*v.imag()); }

// 近似到零
void Chop(RowVectorXcd& v, double tol = 1e-16){
    for(int i = 0; i < v.size(); i++)
        if(Mod(v(i)) < tol)
            v(i) = 0;
}

void Chop(VectorXcd& v, double tol = 1e-16){
    for(int i = 0; i < v.size(); i++)
        if(Mod(v(i)) < tol)
            v(i) = 0;
}

void Chop(MatrixXcd& m, double tol = 1e-16) {
    for(int i = 0; i < m.size(); i++)
        if(Mod(m(i)) < tol)
            m(i) = 0;
}

}
