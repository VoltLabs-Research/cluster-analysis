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

    json result = AnalysisResult::success();
    AnalysisResult::addTiming(result, startTime);
    result["cutoff"] = _cutoff;
    result["cluster_count"] = engine.numClusters();
    result["largest_cluster_size"] = engine.largestClusterSize();
    result["has_zero_weight_cluster"] = engine.hasZeroWeightCluster();

    result["clusters"] = json::array();

    if(!outputFilename.empty()){
        const std::string outputPath = outputFilename + "_cluster_analysis.msgpack";
        if(JsonUtils::writeJsonMsgpackToFile(result, outputPath, false)){
            spdlog::info("Cluster analysis msgpack written to {}", outputPath);
        }else{
            spdlog::warn("Could not write cluster analysis msgpack: {}", outputPath);
        }
    }

    return result;
}

}
