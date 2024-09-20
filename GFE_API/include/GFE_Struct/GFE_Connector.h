#ifndef GFE_CONNECTOR_H
#define GFE_CONNECTOR_H

#include <map>
#include <string>
#include <vector>
#include <array>

#include "api_a.h"

namespace GFE {


class GFE_API Orientation {
public:
	enum Type {
		Rectangular,
		Cylindrical,
		// 暂不支持
		//Spherical,
		//Z_Rectangular
	};
	enum Def {
		Coordinate,
	};

public:

	Orientation() : name(), type(0), definition(0), data() {}
	//! 反序列化的构造函数,
	Orientation(const char* binary);

	//! 序列化
	std::vector<char> ToBinary() const;

	std::array<double, 3> Origin() const;
	std::array<std::array<double, 3>, 3> BaseDirection() const;

public:

	std::string name;
	int type = 0;
	int definition = 0;
	std::vector<double> data;

};

class GFE_API ConnectorSub
{
public:
    enum TYPE {
        ELASTIC = 0,
        PLASTIC,
        HARDENING,
        DAMPING
    };
    explicit ConnectorSub(int _t) : type(_t) {}
    id_t id;
    int type;
    virtual int Type() const = 0;
};

class GFE_API ConnectorElastic : public ConnectorSub
{
public:
    ConnectorElastic() : ConnectorSub(ELASTIC) {}
    int Type() const { return ELASTIC; }

    int component = 0;
    std::string extrapolation = "CONSTANT";
    bool nonlinear = false;
    bool rigid = false;
    double stiffness = 0.0;

    // GFE新增
    std::string definition;
    double tensile_stiffness = 0.0;

    std::vector<std::array<double, 3>> nonlinear_table;

    static constexpr char strGfeEla[] = "GFE ELA";
};

class GFE_API ConnectorPlastic : public ConnectorSub
{
public:
    ConnectorPlastic() : ConnectorSub(PLASTIC) {}
    int Type() const { return PLASTIC; }

    int component = 0;
};

class GFE_API ConnectorHardening : public ConnectorSub
{
public:
    ConnectorHardening() : ConnectorSub(HARDENING) {}
    int Type() const { return HARDENING; }

    std::string type = strKinematic;
    std::string definition = strHalfCyle;
    std::vector<double> values;

    // type
    static constexpr char strIsotropic[] = "ISOTROPIC";
    static constexpr char strKinematic[] = "KINEMATIC";

    // definition
    static constexpr char strTabular[] = "TABULAR";
    static constexpr char strHalfCyle[] = "HALF CYCLE";
    static constexpr char strParam[] = "PARAMETERS";
    static constexpr char strStab[] = "STABILIZED";
    static constexpr char strExpo[] = "EXPONENTIAL LAW";
    static constexpr char strGfe2[] = "GFE HDN2";
    static constexpr char strGfeBW[] = "GFE BW";
    static constexpr char strGfePend[] = "GFE PEND";
};

class GFE_API ConnectorDamping : public ConnectorSub
{
public:
    ConnectorDamping() : ConnectorSub(DAMPING) {}
    int Type() const { return DAMPING; }

    std::string type = strViscous;
    int component = 0;
    std::string extrapolation = strConstant;
    bool nolinear = false;
    std::vector<double> values;

    // type
    static constexpr char strViscous[] = "VISCOUS";
    static constexpr char strStructual[] = "STRUCTURAL";
    static constexpr char strGfe2[] = "GFE DAMP2";

    // extrapolation
    static constexpr char strConstant[] = "CONSTANT";
    static constexpr char strLinear[] = "LINEAR";
};

class GFE_API ConnectorBehavior
{
public:
    static constexpr char key[] = "*CONNECTOR BEHAVIOR";
    std::string name;
    std::vector<std::shared_ptr<ConnectorSub>> behaviors;
};

struct InternalConnectorBehavior
{
    std::string name;
    std::vector<std::pair<int, int>> subtable_type_row;
};
struct PropertyConnector;

struct SpecialInteraction;
using SpecInt = SpecialInteraction;
struct SpecialInteraction
{
    string name;
    int type;
    vector<double> parameters;
    vector<string> namelist;
};

GFE_API void AddConnectorBehavior(shared_ptr<DB>, const ConnectorBehavior&);
GFE_API std::vector<std::string> GetAllConnectorBehaviorName(shared_ptr<DB>);
GFE_API std::vector<std::string> GetAllConnectorName(shared_ptr<DB>);
GFE_API shared_ptr<PropertyConnector> GetConnector(shared_ptr<DB>, const string& name);
GFE_API shared_ptr<ConnectorBehavior> GetConnectorBehavior(shared_ptr<DB>, const string& name);

}
#endif // GFE_CONNECTOR_H
