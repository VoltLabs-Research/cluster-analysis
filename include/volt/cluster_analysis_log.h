#pragma once

#include <volt/cluster_analysis_engine.h>
#include <volt/cluster_analysis_service.h>
#include <volt/core/run_log.h>

#include <cstddef>
#include <cstdint>

namespace Volt{

// The execution log of one cluster analysis. Each phase hands it what it built, and the log writes what that holds
// and how long the phase took. Scripts parse these lines, so their wording stays as it is.
class ClusterAnalysisLog : public RunLog{
public:
    ClusterAnalysisLog(std::size_t atomCount, double cutoff);

    using RunLog::dump;

    void dump(const ClusterAnalysisEngine::ConcurrentDisjointSet& joinedNeighborPairs);
    // The flood fill's unwrapped coordinates.
    void dump(const Particles::ParticleProperty& unwrappedPositions);
    // The clusters once they are complete; the line times the whole clustering.
    void dump(const ClusterAnalysisEngine& engine);

private:
    static constexpr std::size_t kSizeBuckets = 6;
    static constexpr const char* kSizeBucketLabels[kSizeBuckets] = {
        "1 atom", "2-9 atoms", "10-99 atoms", "100-999 atoms", "1000-9999 atoms", "10000+ atoms"
    };

    static std::size_t sizeBucket(std::int64_t size);

    const std::size_t _atomCount;
    const double _cutoff;
};

}
