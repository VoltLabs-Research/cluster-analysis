#include <volt/cluster_analysis_service.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/utilities/json_utils.h>
#include <spdlog/spdlog.h>

namespace Volt{

using namespace Volt::Particles;

ClusterAnalysisService::ClusterAnalysisService()
    : _cutoff(3.2),
      _sortBySize(true),
      _unwrapParticleCoordinates(false),
      _computeCentersOfMass(false),
      _computeRadiusOfGyration(false){}

void ClusterAnalysisService::setCutoff(double cutoff){
    _cutoff = cutoff;
}

void ClusterAnalysisService::setOptions(
    bool sortBySize,
    bool unwrapParticleCoordinates,
    bool computeCentersOfMass,
    bool computeRadiusOfGyration
){
    _sortBySize = sortBySize;
    _unwrapParticleCoordinates = unwrapParticleCoordinates;
    _computeCentersOfMass = computeCentersOfMass;
    _computeRadiusOfGyration = computeRadiusOfGyration;
}

json ClusterAnalysisService::compute(const LammpsParser::Frame& frame, const std::string& outputFilename){
    auto startTime = std::chrono::high_resolution_clock::now();

    if(frame.natoms <= 0)
        return AnalysisResult::failure("Invalid number of atoms");

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if(!positions)
        return AnalysisResult::failure("Failed to create position property");

    spdlog::info("Starting cluster analysis (cutoff = {}, sort = {}, unwrap = {}, com = {}, rg = {})...",
        _cutoff, _sortBySize, _unwrapParticleCoordinates, _computeCentersOfMass, _computeRadiusOfGyration);

    ClusterAnalysisEngine engine(
        positions.get(),
        frame.simulationCell,
        ClusterAnalysis::CutoffRange,
        _cutoff,
        _sortBySize,
        _unwrapParticleCoordinates,
        _computeCentersOfMass,
        _computeRadiusOfGyration
    );

    engine.perform();

    auto clusters = engine.particleClusters();
    auto unwrapped = engine.unwrappedPositions();
    auto clusterSizes = engine.clusterSizes();
    auto clusterIds = engine.clusterIDs();
    auto centers = engine.centersOfMass();
    auto rg = engine.radiiOfGyration();
    auto gt = engine.gyrationTensors();

    const size_t k = engine.numClusters();

    json result;
    result["main_listing"] = {
        { "total_atoms", frame.natoms },
        { "clusters", static_cast<int>(k) },
        { "largest_cluster_size", engine.largestClusterSize() },
        { "has_zero_weight_cluster", engine.hasZeroWeightCluster() }
    };

    // sub_listings: cluster_list
    json clusterList = json::array();
    for(size_t ci = 0; ci < k; ci++){
        json c;
        c["cluster_id"] = static_cast<int64_t>(ci + 1);
        c["size"] = clusterSizes ? clusterSizes->getInt64(ci) : 0;

        if(centers){
            const Point3 p = centers->getPoint3(ci);
            c["center"] = {p.x(), p.y(), p.z()};
        }

        if(rg)  c["radius_of_gyration"] = rg->getDouble(ci);

        if(gt){
            json tensor = json::array();
            for(int comp = 0; comp < 6; comp++) tensor.push_back(gt->getDoubleComponent(ci, comp));
            c["gyration_tensor"] = tensor;
        }

        clusterList.push_back(c);
    }

    result["sub_listings"] = { { "clusters", clusterList } };

    // per-atom-properties
    json perAtom = json::array();
    for(int i = 0; i < frame.natoms; i++){
        json a;
        a["id"] = frame.ids[i];
        a["cluster"] = clusters ? clusters->getInt(i) : 0;

        if(unwrapped){
            const Point3 p = unwrapped->getPoint3(i);
            a["pos_unwrapped"] = {p.x(), p.y(), p.z()};
        }

        perAtom.push_back(a);
    }
    result["per-atom-properties"] = perAtom;

    if(!outputFilename.empty()){
        const std::string outputPath = outputFilename + "_cluster_analysis.msgpack";
        if(JsonUtils::writeJsonMsgpackToFile(result, outputPath, false)){
            spdlog::info("Cluster analysis msgpack written to {}", outputPath);
        }else{
            spdlog::warn("Could not write cluster analysis msgpack: {}", outputPath);
        }
    }

    spdlog::info("Cluster analysis completed. Clusters: {}, largest: {}", k, engine.largestClusterSize());
    return result;
}

}
