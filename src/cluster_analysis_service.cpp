#include <volt/cluster_analysis_service.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/core/progress.h>
#include <volt/utilities/json_utils.h>
#include <volt/utilities/parquet_atom_writer.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <string>

namespace Volt{

using namespace Volt::Particles;

namespace{

constexpr std::size_t kSizeBuckets = 6;
constexpr const char* kSizeBucketLabels[kSizeBuckets] = {
    "1 atom", "2-9 atoms", "10-99 atoms", "100-999 atoms", "1000-9999 atoms", "10000+ atoms"
};

std::size_t sizeBucket(std::int64_t size){
    if(size <= 1) return 0;
    std::size_t bucket = 1;
    for(std::int64_t limit = 10; size >= limit && bucket + 1 < kSizeBuckets; limit *= 10) bucket++;
    return bucket;
}

// Size summary from the per-cluster sizes the engine already holds (one pass over the clusters).
void logClusterSizes(const ClusterAnalysisEngine& engine, std::size_t atomCount, double cutoff){
    const std::size_t k = engine.numClusters();
    const auto sizes = engine.clusterSizes();
    if(k == 0 || !sizes){
        spdlog::warn("No clusters found: all {} atoms are left unclustered", atomCount);
        return;
    }

    std::size_t largestIndex = 0;
    std::int64_t largest = 0;
    std::int64_t smallest = std::numeric_limits<std::int64_t>::max();
    std::int64_t clustered = 0;
    std::array<std::size_t, kSizeBuckets> histogram{};
    for(std::size_t ci = 0; ci < k; ci++){
        const std::int64_t size = sizes->getInt64(ci);
        if(size > largest){
            largest = size;
            largestIndex = ci;
        }
        smallest = std::min(smallest, size);
        clustered += size;
        histogram[sizeBucket(size)]++;
    }

    spdlog::info("Cluster sizes: largest {} atoms ({:.1f}%), smallest {}, mean {:.1f}",
        largest, percentOf(static_cast<std::uint64_t>(largest), atomCount), smallest,
        static_cast<double>(clustered) / static_cast<double>(k));
    if(k > 1){
        std::string histogramText;
        for(std::size_t b = 0; b < kSizeBuckets; b++){
            if(histogram[b] == 0) continue;
            if(!histogramText.empty()) histogramText += ", ";
            histogramText += kSizeBucketLabels[b];
            histogramText += ": ";
            histogramText += std::to_string(histogram[b]);
        }
        spdlog::info("Clusters with {}", histogramText);
    }
    if(k == atomCount && atomCount > 1){
        spdlog::warn("No two atoms are within the cutoff ({}): every atom is its own cluster", cutoff);
    }
    if(static_cast<std::size_t>(clustered) < atomCount){
        spdlog::warn("{} atoms are left unclustered", atomCount - static_cast<std::size_t>(clustered));
    }

    const auto centers = engine.centersOfMass();
    const auto rg = engine.radiiOfGyration();
    const bool hasCenter = centers && largestIndex < centers->size();
    const bool hasRg = rg && largestIndex < rg->size();
    if(hasCenter){
        const Point3 p = centers->getPoint3(largestIndex);
        if(hasRg){
            spdlog::info("Largest cluster (#{}): center ({:.3f}, {:.3f}, {:.3f}), radius of gyration {:.3f}",
                largestIndex + 1, p.x(), p.y(), p.z(), rg->getDouble(largestIndex));
        }else{
            spdlog::info("Largest cluster (#{}): center ({:.3f}, {:.3f}, {:.3f})",
                largestIndex + 1, p.x(), p.y(), p.z());
        }
    }else if(hasRg){
        spdlog::info("Largest cluster (#{}): radius of gyration {:.3f}", largestIndex + 1, rg->getDouble(largestIndex));
    }
    if(engine.hasZeroWeightCluster()){
        spdlog::warn("Some clusters have zero total weight: their center and gyration are reported at the origin");
    }
}

}

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
    if(frame.natoms <= 0)
        return AnalysisResult::failure("Invalid number of atoms");

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if(!positions)
        return AnalysisResult::failure("Failed to create position property");

    spdlog::info("Clustering {} atoms within cutoff {}...", frame.natoms, _cutoff);
    const Stopwatch timer;

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
    auto centers = engine.centersOfMass();
    auto rg = engine.radiiOfGyration();
    auto gt = engine.gyrationTensors();

    const size_t k = engine.numClusters();
    spdlog::info("Found {} clusters in {}", k, timer.elapsed());
    logClusterSizes(engine, static_cast<std::size_t>(frame.natoms), _cutoff);

    json clusterList = json::array();
    for(size_t ci = 0; ci < k; ci++){
        json c;
        c["cluster_id"] = static_cast<int64_t>(ci + 1);
        c["size"] = clusterSizes ? clusterSizes->getInt64(ci) : 0;
        if(centers){
            const Point3 p = centers->getPoint3(ci);
            c["center"] = {p.x(), p.y(), p.z()};
        }
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
        { "clusters", static_cast<int>(k) },
        { "largest_cluster_size", engine.largestClusterSize() },
        { "has_zero_weight_cluster", engine.hasZeroWeightCluster() }
    };
    result["sub_listings"] = { { "clusters", clusterList } };

    if(!outputFilename.empty()){
        const std::string outputPath = outputFilename + "_cluster_analysis.parquet";
        JsonUtils::writeJsonToParquet(result, outputPath, false);
        spdlog::info("Cluster analysis parquet written to {}", outputPath);

        const std::string atomsPath = outputFilename + "_atoms.parquet";
        streamAtomsToParquet(atomsPath, frame,
            [&](std::size_t i){
                const int cid = clusters ? clusters->getInt(i) : 0;
                return cid > 0 ? "Cluster_" + std::to_string(cid) : std::string("Unclustered");
            },
            [&](ColumnarAtomWriter& w, std::size_t i){
                if(unwrapped){
                    const Point3 pu = unwrapped->getPoint3(i);
                    w.field("pos_unwrapped", std::vector<double>{pu.x(), pu.y(), pu.z()});
                }
                w.field("cluster_id", clusters ? clusters->getInt(i) : 0);
            },
            {},
            false
        );
        spdlog::info("Exported atoms data to: {}", atomsPath);
    }

    return result;
}

}
