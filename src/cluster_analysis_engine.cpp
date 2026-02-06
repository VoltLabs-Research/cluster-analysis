#include <volt/cluster_analysis_engine.h>
#include <volt/analysis/cutoff_neighbor_finder.h>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <numeric>
#include <random>
#include <cmath>
#include <limits>

namespace Volt{

using namespace Volt::Particles;

ClusterAnalysisEngine::ClusterAnalysisEngine(
    ParticleProperty* positions,
    const SimulationCell& cell,
    ClusterAnalysis::NeighborMode neighborMode,
    double cutoff,
    bool sortBySize,
    bool unwrapParticleCoordinates,
    bool computeCentersOfMass,
    bool computeRadiusOfGyration
)
    : _positions(positions),
      _simCell(cell),
      _neighborMode(neighborMode),
      _cutoff(cutoff),
      _sortBySize(sortBySize),
      _unwrapParticleCoordinates(unwrapParticleCoordinates),
      _computeCentersOfMass(computeCentersOfMass),
      _computeRadiusOfGyration(computeRadiusOfGyration),
      _numClusters(0),
      _largestClusterSize(0),
      _hasZeroWeightCluster(false)
{
    if(!_positions){
        throw std::runtime_error("ClusterAnalysisEngine: positions is null");
    }
    if(_cutoff <= 0.0 && _neighborMode == ClusterAnalysis::CutoffRange){
        throw std::runtime_error("ClusterAnalysisEngine: cutoff must be > 0 for cutoff clustering");
    }
}


void ClusterAnalysisEngine::perform(){
    const size_t n = _positions->size();

    if(n == 0){
        _particleClusters.reset();
        _unwrappedPositions.reset();
        _clusterSizes.reset();
        _clusterIDs.reset();
        _centersOfMass.reset();
        _radiiOfGyration.reset();
        _gyrationTensors.reset();
        _numClusters = 0;
        _largestClusterSize = 0;
        return;
    }

    // per-particle cluster id (int).
    _particleClusters = std::make_shared<ParticleProperty>(n, ParticleProperty::ClusterProperty, 1, true);

    // -1 for "unassigned", 0 for "excluded"
    for(size_t i = 0; i < n; i++){
        _particleClusters->setInt(i, -1);
    }

    // Unwrapped positions output
    if(_unwrapParticleCoordinates || _computeCentersOfMass || _computeRadiusOfGyration){
        _unwrappedPositions = std::make_shared<ParticleProperty>(n, ParticleProperty::PositionProperty, 3, true);

        // Initialize with wrapped/current positions
        const Point3* pos = _positions->constDataPoint3();
        Point3* uw = _unwrappedPositions->dataPoint3();
        for(size_t i = 0; i < n; i++){
            uw[i] = pos[i];
        }
    }

    // Per-cluster outputs
    if(_computeCentersOfMass){
        _centersOfMass = std::make_shared<ParticleProperty>(0, ParticleProperty::PositionProperty, 3, false);
    }

    if(_computeRadiusOfGyration){
        _radiiOfGyration = std::make_shared<ParticleProperty>(0, DataType::Double, 1, 0, false);
        _gyrationTensors = std::make_shared<ParticleProperty>(0, DataType::Double, 6, 0, false);
    }

    if(_neighborMode == ClusterAnalysis::Bonding){
        throw std::runtime_error("ClusterAnalysisEngine: Bonding mode is not implemented in Volt CLI analyzer yet (missing bonds topology in frame).");
    }

    std::vector<Point3> centers;
    centers.reserve(128);

    doClusteringCutoff(centers);

    if(_computeRadiusOfGyration){
        computeGyration(centers);
    }

    buildClusterSizes();
    buildClusterIds();

    if(_sortBySize){
        sortClustersBySize();
    }
}

void ClusterAnalysisEngine::doClusteringCutoff(std::vector<Point3>& centerOfMass){
    const std::size_t n = _positions->size();

    CutoffNeighborFinder neighFinder;
    if(!neighFinder.prepare(_cutoff, _positions, _simCell)) return;

    const Point3* pos = _positions->constDataPoint3();
    Point3* uw = _unwrappedPositions ? _unwrappedPositions->dataPoint3() : nullptr;

    std::deque<size_t> queue;
    for(size_t seed = 0; seed < n; seed++){
        if(_particleClusters->getInt(seed) != -1){
            continue;
        }

        _numClusters++;
        const int clusterId = static_cast<int>(_numClusters);
        _particleClusters->setInt(seed, clusterId);

        Vector3 cm = Vector3::Zero();
        double totalW = 0.0;

        queue.clear();
        queue.push_back(seed);

        while(!queue.empty()){
            size_t curr = queue.front();
            queue.pop_front();

            // Iterate neighbors
            for(CutoffNeighborFinder::Query q(neighFinder, curr); !q.atEnd(); q.next()){
                const size_t nb = q.current();
                if(nb >= n) continue;

                if(_particleClusters->getInt(nb) != -1) continue;

                _particleClusters->setInt(nb, clusterId);
                queue.push_back(nb);

                if(uw){
                    Vector3 delta = q.delta();
                    uw[nb] = uw[curr] + delta;
                    if(_computeCentersOfMass || _computeRadiusOfGyration){
                        const double w = 1.0;
                        cm += w * (uw[nb] - Point3::Origin());
                        totalW += w;
                    }
                }
            }
        }

        if(_computeCentersOfMass || _computeRadiusOfGyration){
            if(!uw){
                centerOfMass.push_back(Point3::Origin());
                _hasZeroWeightCluster = true;
            }else{
                const double w0 = 1.0;
                cm += w0 * (uw[seed] - Point3::Origin());
                totalW += w0;

                if(totalW > 0.0){
                    centerOfMass.push_back(Point3::Origin() + (cm / totalW));
                }else{
                    centerOfMass.push_back(Point3::Origin());
                    _hasZeroWeightCluster = true;
                }
            }
        }
    }

    _largestClusterSize = 0;
    if(_centersOfMass){
        _centersOfMass->resize(centerOfMass.size(), false);
        for(size_t i = 0; i < centerOfMass.size(); i++){
            _centersOfMass->setPoint3(i, centerOfMass[i]);
        }
    }
}

void ClusterAnalysisEngine::computeGyration(const std::vector<Point3>& centersOfMass){
    if(!_radiiOfGyration || !_gyrationTensors) return;
    if(!_unwrappedPositions) return;

    const size_t n = _positions->size();
    const size_t k = centersOfMass.size();

    _radiiOfGyration->resize(k, true);
    _gyrationTensors->resize(k, true);

    double* rg = _radiiOfGyration->dataDouble();
    double* gt = _gyrationTensors->dataDouble();

    std::vector<double> clusterMass(k, 0.0);

    const Point3* uw = _unwrappedPositions->constDataPoint3();
    const double *massData = nullptr;

    for(size_t i = 0; i < n; i++){
        const int cid = _particleClusters->getInt(i);
        if(cid <= 0) continue;

        const size_t c = static_cast<size_t>(cid - 1);
        if(c >= k) continue;

        const double w = 1.0;
        clusterMass[c] += w;

        Vector3 d = uw[i] - centersOfMass[c];
        rg[c] += w * d.squaredLength();
        double* t = gt + c * 6;
        // Order: xx, yy, zz, xy, xz, yz (6 comps)
        t[0] += w * d.x() * d.x();
        t[1] += w * d.y() * d.y();
        t[2] += w * d.z() * d.z();
        t[3] += w * d.x() * d.y();
        t[4] += w * d.x() * d.z();
        t[5] += w * d.y() * d.z();
    }

    for(size_t c = 0; c < k; c++){
        double M = clusterMass[c];
        if(M <= 0.0){
            M = 1.0;
            _hasZeroWeightCluster = true;
        }

        rg[c] = std::sqrt(rg[c] / M);

        double* t = gt + c * 6;
        for(int j = 0; j < 6; j++){
            t[j] /= M;
        }
    }
}

void ClusterAnalysisEngine::buildClusterSizes(){
    _clusterSizes = std::make_shared<ParticleProperty>(_numClusters, DataType::Int64, 1, 0, true);

    for(size_t i = 0; i < _numClusters; i++){
        _clusterSizes->setInt64(i, 0);
    }

    const size_t n = _positions->size();
    for(size_t i = 0; i < n; i++){
        const int cid = _particleClusters->getInt(i);
        if(cid > 0){
            const size_t idx = static_cast<size_t>(cid - 1);
            _clusterSizes->setInt64(idx, _clusterSizes->getInt64(idx) + 1);
        }
    }

    _largestClusterSize = 0;
    for(size_t i = 0; i < _numClusters; i++){
        _largestClusterSize = std::max<size_t>(_largestClusterSize, static_cast<size_t>(_clusterSizes->getInt64(i)));
    }
}

void ClusterAnalysisEngine::buildClusterIds(){
    _clusterIDs = std::make_shared<ParticleProperty>(_numClusters, DataType::Int64, 1, 0, true);
    for(size_t i = 0; i < _numClusters; i++){
        _clusterIDs->setInt64(i, static_cast<int64_t>(i + 1));
    }
}

void ClusterAnalysisEngine::sortClustersBySize(){
    if(!_clusterSizes || _numClusters == 0) return;

    std::vector<size_t> mapping(_numClusters);
    std::iota(mapping.begin(), mapping.end(), size_t(0));

    // Sort mapping by clusterSizes desc
    std::sort(mapping.begin(), mapping.end(), [&](size_t a, size_t b){
        return _clusterSizes->getInt64(a) > _clusterSizes->getInt64(b);
    });

    // Build inverse mapping: oldClusterIndex -> newClusterIndex (0 ... k - 1)
    std::vector<size_t> inv(_numClusters + 1);
    inv[0] = 0;
    for(size_t newIdx = 0; newIdx < _numClusters; newIdx++){
        inv[mapping[newIdx] + 1] = newIdx + 1;
    }

    // Remap particle cluster IDs
    const size_t n = _positions->size();
    for(size_t i = 0; i < n; i++){
        int cid = _particleClusters->getInt(i);
        if(cid >= 0 && static_cast<size_t>(cid) < inv.size()){
            _particleClusters->setInt(i, static_cast<int>(inv[cid]));
        }
    }

    // Reorder clusterSizes/ids
    auto reorder_int64 = [&](std::shared_ptr<ParticleProperty>& prop){
        if(!prop) return;
        std::vector<std::int64_t> tmp(_numClusters);
        for(std::size_t i = 0; i < _numClusters; i++){
            tmp[i] = prop->getInt64(mapping[i]);
        }

        for(size_t i = 0; i < _numClusters; i++){
            prop->setInt64(i, tmp[i]);
        }
    };

    reorder_int64(_clusterSizes);
    reorder_int64(_clusterIDs);

    if(_centersOfMass){
        std::vector<Point3> tmp(_numClusters);
        for(std::size_t i = 0; i < _numClusters; i++){
            tmp[i] = _centersOfMass->getPoint3(mapping[i]);
        }

        for(size_t i = 0; i < _numClusters; i++){
            _centersOfMass->setPoint3(i, tmp[i]);
        }
    }

    if(_radiiOfGyration){
        std::vector<double> tmp(_numClusters);
        for(size_t i = 0; i < _numClusters; i++){
            tmp[i] = _radiiOfGyration->getDouble(mapping[i]);
        }
        for(size_t i = 0; i < _numClusters; i++){
            _radiiOfGyration->setDouble(i, tmp[i]);
        }
    }

    if(_gyrationTensors){
        std::vector<std::array<double, 6>> tmp(_numClusters);
        for(size_t i = 0; i < _numClusters; i++){
            std::array<double, 6> t{};
            for(int c = 0; c < 6; c++){
                t[c] = _gyrationTensors->getDoubleComponent(mapping[i], c);
            }
            tmp[i] = t;
        }

        for(size_t i = 0; i < _numClusters; i++){
            for(int c = 0; c < 6; c++){
                _gyrationTensors->setDoubleComponent(i, c, tmp[i][c]);
            }
        }
    }

    _largestClusterSize = (_clusterSizes && _numClusters > 0) ? static_cast<size_t>(_clusterSizes->getInt64(0)) : 0;
}

}
