#include <volt/cluster_analysis_service.h>
#include <volt/cluster_analysis_log.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/core/progress.h>
#include <volt/utilities/json_utils.h>
#include <volt/utilities/parquet_atom_writer.h>

#include <string>
#include <vector>

namespace Volt{

using namespace Volt::Particles;

const OutputFile<ClusterAnalysisService, ClusterAnalysisService::Outputs> ClusterAnalysisService::outputFiles[] = {
    {"Writing the cluster table", 1.0, nullptr, [](const Outputs& out){
        out.log.dump(JsonUtils::writeSummaryTable(out.clusterTable, out.base, "cluster_analysis", "Cluster analysis"));
    }},
    {"Writing the atoms table", 3.0, nullptr, [](const Outputs& out){
        const auto clusters = out.engine.particleClusters();
        const auto unwrapped = out.engine.unwrappedPositions();
        out.log.dump(writeAtomsTable(out.base, out.frame,
            [&](std::size_t i){
                const int cid = clusters ? clusters->getInt(i) : 0;
                return cid > 0 ? "Cluster_" + std::to_string(cid) : std::string("Unclustered");
            },
            [&](ColumnarRowWriter& w, std::size_t i){
                if(unwrapped){
                    w.field("pos_unwrapped", unwrapped->getPoint3(i));
                }
                w.field("cluster_id", clusters ? clusters->getInt(i) : 0);
            }
        ));
    }},
};

ClusterAnalysisService::ClusterAnalysisService()
    : _cutoff(3.2),
      _sortBySize(true),
      _unwrapParticleCoordinates(false),
      _computeCentersOfMass(false),
      _computeRadiusOfGyration(false){}

void ClusterAnalysisService::setCutoff(double cutoff){
    _cutoff = cutoff;
}

void ClusterAnalysisService::setSortBySize(bool sortBySize){
    _sortBySize = sortBySize;
}

void ClusterAnalysisService::setUnwrapParticleCoordinates(bool unwrapParticleCoordinates){
    _unwrapParticleCoordinates = unwrapParticleCoordinates;
}

void ClusterAnalysisService::setComputeCentersOfMass(bool computeCentersOfMass){
    _computeCentersOfMass = computeCentersOfMass;
}

void ClusterAnalysisService::setComputeRadiusOfGyration(bool computeRadiusOfGyration){
    _computeRadiusOfGyration = computeRadiusOfGyration;
}

json ClusterAnalysisService::compute(const LammpsParser::Frame& frame, const std::string& outputFilename){
    if(frame.natoms <= 0)
        return AnalysisResult::failure("Invalid number of atoms");

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if(!positions)
        return AnalysisResult::failure("Failed to create position property");

    ProgressSteps steps({
        {"Clustering atoms", 2.0},
        {"Writing outputs", totalWeight(outputSteps(outputFiles, *this))},
    });
    ClusterAnalysisLog log(static_cast<std::size_t>(frame.natoms), _cutoff);

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

    engine.perform(log);
    log.dump(engine);
    steps.next();

    auto clusterSizes = engine.clusterSizes();
    auto centers = engine.centersOfMass();
    auto rg = engine.radiiOfGyration();
    auto gt = engine.gyrationTensors();

    const size_t k = engine.numClusters();

    json clusterList = json::array();
    for(size_t ci = 0; ci < k; ci++){
        json c;
        c["cluster_id"] = ci + 1;
        c["size"] = clusterSizes ? clusterSizes->getInt64(ci) : 0;
        if(centers) c["center"] = centers->getPoint3(ci);
        if(rg) c["radius_of_gyration"] = rg->getDouble(ci);
        if(gt){
            json tensor = json::array();
            for(int comp = 0; comp < 6; comp++) tensor.push_back(gt->getDoubleComponent(ci, comp));
            c["gyration_tensor"] = tensor;
        }
        clusterList.push_back(c);
    }

    json result;
    result["main_listing"] = {
        { "total_atoms", frame.natoms },
        { "clusters", k },
        { "largest_cluster_size", engine.largestClusterSize() },
        { "has_zero_weight_cluster", engine.hasZeroWeightCluster() }
    };
    result["sub_listings"] = { { "clusters", clusterList } };

    if(!outputFilename.empty()){
        writeOutputFiles(outputFiles, *this, Outputs{outputFilename, frame, engine, result, log});
    }

    return result;
}

}
