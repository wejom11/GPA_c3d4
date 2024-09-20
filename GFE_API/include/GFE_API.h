#pragma once
#include "GFE_API_2.h"
#include "api_st.h"
#include <array>
#include <functional>
#include <map>
#include <set>
#include <deprecated.h>

namespace GFE {

GFE_API void setLogPath(const string& path);
GFE_API bool activeLog(bool);
GFE_API void setBusyDelay(int ms);

/*
 * Common
 * @open: 打开数据库. 需要新建并打开, 或者数据库文件与api版本不匹配时请将sync设为true, 否则设为false(sync会降低效率)
 * @vacuum: 释放无用空间
 * @backup: 持久化内存数据库
 * @filepath: 获取数据库文件路径
 * @transaction: 执行一次事务
 * @clearModel: 清空数据库模型信息(不释放已占用的空间, 以提高下一次插入的效率, 使用vacuum释放空间)
 * @setBusyTimeout: default value = 5000ms
 */
GFE_API shared_ptr<DB> open(const string& path, bool sync = false);
GFE_API void vacuum(shared_ptr<DB>);
GFE_API void backup(shared_ptr<DB>, const string& path);
GFE_API string filepath(shared_ptr<DB>);
GFE_API bool transaction(shared_ptr<DB>, std::function<bool()> f);
GFE_API void clearModel(shared_ptr<DB>);
GFE_API void clearOutput(shared_ptr<DB>);
GFE_API void setBusyTimeout(shared_ptr<DB>, int ms);

/*
 * Model Information
 */
GFE_API vector<id_t> getTypeOffset(shared_ptr<DB>);
GFE_API id_t getNodeNum(shared_ptr<DB>);
GFE_API id_t getElementNum(shared_ptr<DB>);
GFE_API std::array<id_t, CT_NUM> getElementTypeNum(shared_ptr<DB>);

/*
 * Node
 */
GFE_API bool setNode(shared_ptr<DB>, const vector<data_t>& coordinates);
GFE_API bool addNodeAttribute(shared_ptr<DB>, const string& name, const vector<id_t>&);
GFE_API bool addNodeAttribute(shared_ptr<DB>, const string& name, const vector<data_t>&);
GFE_API vector<Node> getNode(shared_ptr<DB>, vector<id_t>&& subset = {-1});          // deprecated
GFE_API vector<data_t> getNodeCoordinate(shared_ptr<DB>, vector<id_t>&& subset = {-1});
GFE_API vector<id_t> getNodeAttribute(shared_ptr<DB>, const string& name, vector<id_t>&& subset = {-1});
GFE_API vector<data_t> getNodeAttributeF(shared_ptr<DB>, const string& name, vector<id_t>&& subset = {-1});

/*
 * Element
 */
GFE_API bool setElement(shared_ptr<DB>, const vector<id_t>& el_nodes, const std::array<id_t, CT_NUM>& each_type_num);
GFE_API bool addElementAttribute(shared_ptr<DB>, const string& name, const vector<id_t>&);
GFE_API bool addElementAttribute(shared_ptr<DB>, const string& name, const vector<data_t>&);
GFE_API vector<Element> getElement(shared_ptr<DB>, vector<id_t>&& subset = {-1});            // deprecated
GFE_API vector<std::tuple<int, vector<id_t>>> getElement2(shared_ptr<DB>, vector<id_t>&& subset = {-1});        // {{type, nodes}}
GFE_API vector<id_t> getElementAttribute(shared_ptr<DB>, const string& name, vector<id_t>&& subset = {-1});
GFE_API vector<data_t> getElementAttributeF(shared_ptr<DB>, const string& name, vector<id_t>&& subset = {-1});

/*
 * Node&Element
 */
inline vector<id_t> getAttribute(shared_ptr<DB> db, const string& name, vector<id_t>&& subset = {-1}) {
    if(name.substr(0, 4) == "Node") return getNodeAttribute(db, name.substr(4, string::npos), std::move(subset));
    else if(name.substr(0, 7) == "Element") return getElementAttribute(db, name.substr(7, string::npos), std::move(subset));
    else return vector<id_t>();
}
inline vector<data_t> getAttributeF(shared_ptr<DB> db, const string& name, vector<id_t>&& subset = {-1}) {
    if(name.substr(0, 4) == "Node") return getNodeAttributeF(db, name.substr(4, string::npos), std::move(subset));
    else if(name.substr(0, 7) == "Element") return getElementAttributeF(db, name.substr(7, string::npos), std::move(subset));
    else return vector<data_t>();
}
GFE_API bool hasAttribute(shared_ptr<DB> db, const string& name);

/*
 * NodeSet
 */
GFE_API void addNodeSet(shared_ptr<DB> db, const NodeSet& nset);
GFE_API std::unique_ptr<NodeSet> getNodeSet(shared_ptr<DB> db, const string& name, bool nocase = false);
GFE_API vector<NodeSet> getAllNodeSet(shared_ptr<DB> db);
GFE_API vector<string> getAllNodeSetName(shared_ptr<DB> db);

/*
 * ElementSet
 */
GFE_API void addElementSet(shared_ptr<DB>, const ElementSet& elset);
GFE_API std::unique_ptr<ElementSet> getElementSet(shared_ptr<DB>, const string& name, bool nocase = false);
GFE_API vector<ElementSet> getAllElementSet(shared_ptr<DB>);
GFE_API vector<string> getAllElementSetName(shared_ptr<DB>);
GFE_API vector<id_t> getElementSetNodes(shared_ptr<DB>, const string& name, bool nocase = false);

/*
 * Surface
 */
GFE_API void addSurface(shared_ptr<DB>, const Surface&);
GFE_API vector<Surface> getAllSurface(shared_ptr<DB>);
GFE_API std::unique_ptr<Surface> getSurface(shared_ptr<DB>, const string& name, bool nocase = false);
GFE_API vector<id_t> getSurfaceNodes(shared_ptr<DB>, const string& name, bool nocase = false);
GFE_API vector<string> getAllSurfaceName(shared_ptr<DB>);
GFE_API void addSurfacePair(shared_ptr<DB>, const SurfacePair&);
GFE_API vector<SurfacePair> getAllSurfacePair(shared_ptr<DB>);
GFE_API vector<string> getAllSurfacePairName(shared_ptr<DB>);
GFE_API std::unique_ptr<SurfacePair> getSurfacePair(shared_ptr<DB>, const string& name, bool nocase = false);

/*
 * Material
 */
GFE_API void addMaterial(shared_ptr<DB> db, const Material&);
GFE_API vector<string> getAllMaterialName(shared_ptr<DB>);
GFE_API vector<std::tuple<id_t, string>> getAllMatIdAndName(shared_ptr<DB>);
GFE_API std::unique_ptr<Material> getMaterial(shared_ptr<DB>, const string& name, bool nocase = false);
GFE_API string getMaterialName(shared_ptr<DB>, id_t mat_id);
GFE_API id_t getMaterialId(shared_ptr<DB>, const string& name);

/*
 * Property
 */
GFE_API void addProperty(shared_ptr<DB>, Property*);
GFE_API vector<string> getAllPropertyName(shared_ptr<DB>);
GFE_API shared_ptr<Property> getProperty(shared_ptr<DB>, const string& name, bool nocase = false);

/*
 * RigidBody
 */
GFE_API void addRigidBody(shared_ptr<DB>, const RigidBody&);
GFE_API vector<string> getAllRigidBodyName(shared_ptr<DB>);
GFE_API vector<RigidBody> getAllRigidBody(shared_ptr<DB>);
GFE_API std::unique_ptr<RigidBody> getRigidBody(shared_ptr<DB>, const string& name, bool nocase = false);

/*
 * MPC
 */
GFE_API void addMPC(shared_ptr<DB>, const MPC&);
GFE_API std::unique_ptr<MPC> getMPC(shared_ptr<DB>);

/*
 * Embed
 */
GFE_API void addEmbed(shared_ptr<DB>, const Embed&);
GFE_API vector<string> getAllEmbedName(shared_ptr<DB>);
GFE_API vector<Embed> getAllEmbed(shared_ptr<DB>);
GFE_API std::unique_ptr<Embed> getEmbed(shared_ptr<DB>, const string& name, bool nocase = false);

/*
 * Incident Wave Property
 */
GFE_API void addIncWaveProp(shared_ptr<DB>, const IncidentWaveProperty&);
GFE_API vector<string> getAllIncWavePropName(shared_ptr<DB>);
GFE_API vector<IncidentWaveProperty> getAllIncWaveProp(shared_ptr<DB>);
GFE_API std::unique_ptr<IncidentWaveProperty> getIncWaveProp(shared_ptr<DB>, const string& name, bool nocase = false);

/*
 * Incident Wave
 */
GFE_API void addIncWave(shared_ptr<DB>, const IncidentWave&);
GFE_API vector<string> getAllIncWaveName(shared_ptr<DB>);
GFE_API vector<IncidentWave> getAllIncWave(shared_ptr<DB>);
GFE_API std::unique_ptr<IncidentWave> getIncWave(shared_ptr<DB>, const string& name, bool nocase = false);

/*
 * Function
 */
GFE_API void addFunction(shared_ptr<DB>, const Function&);
GFE_API vector<string> getAllFunctionName(shared_ptr<DB>);
GFE_API shared_ptr<Function> getFunction(shared_ptr<DB>, const string&, bool nocase = false);

/*
 * FieldOutput
 */
GFE_API vector<FieldOutput> getFieldOutput(shared_ptr<DB>, const string& var, vector<int>&& frames = {-1});
inline vector<data_t> getFieldOutput(shared_ptr<DB> db, const string& var, int frame) {
    auto ret = getFieldOutput(db, var, vector<int>{frame});
    return ret.empty() ? vector<data_t>() : ret[0].data;
}
GFE_API int getFieldOutputFrameNum(shared_ptr<DB>);
GFE_API int getFieldOutputFrameNumHasOutput(shared_ptr<DB>);

/*
 * HistoryOutput
 */
GFE_API void addHistoryOutput(shared_ptr<DB>, const HistoryOutput& output, id_t request_id);
GFE_API void addHistoryOutputRequest(shared_ptr<DB>, const HistoryOutputRequest& request);
GFE_API std::unique_ptr<HistoryOutput> getHistoryOutput(shared_ptr<DB>, int id);
GFE_API vector<HistoryOutput> getHistoryOutput(shared_ptr<DB>, const string& reg = "%", const string& var = "%", const string& pid = "%", const string& pos = "%");
GFE_API vector<HistoryOutputRequest> getAllHistoryOutputRequest(shared_ptr<DB>);
GFE_API vector<string> getHistoryVariable(shared_ptr<DB>, int id);
GFE_API vector<std::tuple<id_t, id_t, int, string, string>> getAllHistoryOutputInfomation(shared_ptr<DB>);
GFE_API vector<double> getHistoryFrequency(shared_ptr<DB>, int id);
GFE_API vector<double> getHistoryFrequencyOfQuas(shared_ptr<DB>, const std::set<string>&, bool exclude);      // qua: quantity缩写
inline vector<double> getHistoryFrequencyByQuas(shared_ptr<DB> db, const std::set<string>& qs) { return getHistoryFrequencyOfQuas(db, qs, false); }
inline vector<double> getHistoryFrequencyExQuas(shared_ptr<DB> db, const std::set<string>& qs) { return getHistoryFrequencyOfQuas(db, qs, true); }

/*
 * Output(Common)
 */
GFE_API vector<double> getOutputFrequency(shared_ptr<DB>);

/*
 * Parameter
 */
GFE_API void addParameter(shared_ptr<DB>, const Parameter&);
GFE_API string getParameter(shared_ptr<DB>, const string& name);

/*
 * ArtBoundary
 */
GFE_API void addArtBoundary(shared_ptr<DB>, const ArtBoundary&);
GFE_API std::unique_ptr<ArtBoundary> getArtBoundary(shared_ptr<DB>);

/*
 * Soil
 */
GFE_API void addSoil(shared_ptr<DB>, const SoilLayer&);
GFE_API vector<string> getAllSoilName(shared_ptr<DB>);
GFE_API vector<SoilLayer> getAllSoil(shared_ptr<DB>);
GFE_API std::unique_ptr<SoilLayer> getSoil(shared_ptr<DB>, const string& name, bool nocase = false);

/*
 * Vibration Load
 */
GFE_API void addVibLoad(shared_ptr<DB>, const VibLoad&);
GFE_API vector<string> getAllVibLoadName(shared_ptr<DB>);
GFE_API vector<VibLoad> getAllVibLoad(shared_ptr<DB>);
GFE_API std::unique_ptr<VibLoad> getVibLoad(shared_ptr<DB>, const string& name, bool nocase = false);

/*
 * Spring & Dashpot
 */
GFE_API void addSpringDashpot(shared_ptr<DB>, const SpringDashpot&);
GFE_API vector<string> getAllSpringDashpotName(shared_ptr<DB>);
GFE_API vector<SpringDashpot> getAllSpringDashpot(shared_ptr<DB>);
GFE_API std::unique_ptr<SpringDashpot> getSpringDashpot(shared_ptr<DB>, const string& name, bool nocase = false);

}
