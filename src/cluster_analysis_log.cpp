#include <volt/cluster_analysis_log.h>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <limits>
#include <string>

namespace Volt{

using namespace Volt::Particles;

ClusterAnalysisLog::ClusterAnalysisLog(std::size_t atomCount, double cutoff)
    : _atomCount(atomCount),
      _cutoff(cutoff){}

std::size_t ClusterAnalysisLog::sizeBucket(std::int64_t size){
    if(size <= 1) return 0;
    std::size_t bucket = 1;
    for(std::int64_t limit = 10; size >= limit && bucket + 1 < kSizeBuckets; limit *= 10) bucket++;
    return bucket;
}

void ClusterAnalysisLog::dump(const ClusterAnalysisEngine::ConcurrentDisjointSet&){
    spdlog::info("  Union-find over neighbor pairs done in {}", phaseElapsed());
}

void ClusterAnalysisLog::dump(const ParticleProperty&){
    spdlog::info("  Flood fill with unwrapped coordinates done in {}", phaseElapsed());
}

void ClusterAnalysisLog::dump(const ClusterAnalysisEngine& engine){
    const std::size_t k = engine.numClusters();
    spdlog::info("Found {} clusters in {}", k, runElapsed());

    const auto sizes = engine.clusterSizes();
    if(k == 0 || !sizes){
        spdlog::warn("No clusters found: all {} atoms are left unclustered", _atomCount);
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
        largest, percentOf(static_cast<std::uint64_t>(largest), _atomCount), smallest,
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
    if(k == _atomCount && _atomCount > 1){
        spdlog::warn("No two atoms are within the cutoff ({}): every atom is its own cluster", _cutoff);
    }
    if(static_cast<std::size_t>(clustered) < _atomCount){
        spdlog::warn("{} atoms are left unclustered", _atomCount - static_cast<std::size_t>(clustered));
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
