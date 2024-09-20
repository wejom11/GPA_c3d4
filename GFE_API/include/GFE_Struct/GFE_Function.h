#pragma once
#pragma warning(disable: 4251)
#include "../api_a.h"
#include <functional>
#include <map>

namespace GFE
{
struct GFE_API Function
{
    enum FuncType{
        TABULAR,
        TABULAR_XLOG,
        TABULAR_YLOG,
        TABULAR_XYLOG,
        SMOOTH_STEP,
        POWER_SERIES,
        SPECTRUM
    };
    static const std::map<int, std::string>& TypeStrMap();
    string TypeStr() const;

    int id;
    int type = TABULAR;
    string name;
    int internal_type = -1;		// for spectrum
    double gravity = 0;			// for spectrum
    vector<double> values;     //x0, y0, x1, y1, x2, y2 ...
};


template <typename Tp>
class Interp;
/**
 * @brief The Function2 class
 * Function因历史遗留问题(引用的地方很多), 不好改, 所以另外写了一个Function2
 */
class GFE_API Function2
{
public:
    using FuncType = Function::FuncType;
    Function2();
    //! @brief: 使用Function构造Function2(拷贝数据)
    /*explicit */Function2(const Function&);
    Function2(const Function2&);
    Function2(Function2&&);
    ~Function2();
    Function2& operator =(const Function2& rhs);
    Function2& operator =(Function2&& rhs);
    GFE::Function ToFunc1();
    void SetExtrap(bool f);     // 是否开启外插

    string Name() const { return m_name; }
    int Type() const { return m_type; }
    const vector<double>& DataArr() const { return m_data; }

    //! @warning 下面几个函数都是以type=TABULAR为前提的
    //! 如果之后增加了其他类型, 可能需要增加条件分支
    double Value(double x) const;
    std::tuple<std::size_t, double, double> Max() const;
    std::tuple<std::size_t, double, double> Min() const;
    GFE::Function2 Abs() const;
    std::tuple<std::size_t, double, double> AbsMax() const;
    std::tuple<std::size_t, double, double> AbsMin() const;


    void Update(const Function&);
    void Clear();
    //! 将绝对值远小于AbsMax的值置为0
    void Chop(double multiple = 1e7);

    //! IO
    //! @return last error line, 0: success, -1: file not exists
    long long Read(const std::string& path);
    //! @return last error line, 0: success, -1: file not exists
    long long Write(const std::string& path);

    //! 操作符重载
    Function2 operator -();
    GFE_API friend Function2 operator +(const Function2& a, const Function2& b);
    GFE_API friend Function2 operator -(const Function2& a, const Function2& b);
    GFE_API friend Function2 operator *(const Function2& a, const Function2& b);
    GFE_API friend Function2 operator /(const Function2& a, const Function2& b);

private:
    GFE_API friend Function2 OpHelper(const Function2& a, const Function2& b, char type);
    std::tuple<std::size_t, double, double> StatHelper(const std::function<bool(double,double)>& comp) const;

    int m_type = Function::TABULAR;
    string m_name;
    vector<double> m_data;
    Interp<double>* m_interp;
};
using FuncType = Function::FuncType;
}
