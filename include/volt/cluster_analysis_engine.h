#pragma once

#include <volt/core/volt.h>
#include <volt/core/simulation_cell.h>
#include <volt/core/particle_property.h>

#include <memory>
#include <vector>
#include <deque>
#include <atomic>

namespace Volt{

using namespace Particles;

class CutoffNeighborFinder;

class ClusterAnalysis{
public:
    enum NeighborMode{
        CutoffRange = 0,
        Bonding = 1
    };
};

class ClusterAnalysisEngine{
public:
    ClusterAnalysisEngine(
        ParticleProperty* positions,
        const SimulationCell& cell,
        ClusterAnalysis::NeighborMode neighborMode,
        double cutoff,
        bool sortBySize,
        bool unwrapParticleCoordinates,
        bool computeCentersOfMass,
        bool computeRadiusOfGyration
    );
    
    void perform();

    std::shared_ptr<ParticleProperty> particleClusters() const{
        return _particleClusters;
    }

    std::shared_ptr<ParticleProperty> unwrappedPositions() const{
        return _unwrappedPositions;
    }

    std::shared_ptr<ParticleProperty> clusterSizes() const{
        return _clusterSizes;
    }

    std::shared_ptr<ParticleProperty> clusterIDs() const{
        return _clusterIDs;
    }

    std::shared_ptr<ParticleProperty> centersOfMass() const{
        return _centersOfMass;
    }

    std::shared_ptr<ParticleProperty> radiiOfGyration() const{
        return _radiiOfGyration;
    }

    std::shared_ptr<ParticleProperty> gyrationTensors() const{
        return _gyrationTensors;
    }

    size_t numClusters() const{
        return _numClusters;
    }

    size_t largestClusterSize() const{
        return _largestClusterSize;
    }

    bool hasZeroWeightCluster() const{
        return _hasZeroWeightCluster;
    }

private:
    void doClusteringCutoff(std::vector<Point3>& centersOfMass);
    void computeGyration(const std::vector<Point3>& centersOfMass);
    void buildClusterSizes();
    void buildClusterIds();
    void sortClustersBySize();
    void buildColors();

    ParticleProperty* _positions;
    SimulationCell _simCell;

    ClusterAnalysis::NeighborMode _neighborMode;
    double _cutoff;
    bool _onlySelectedParticles;
    bool _sortBySize;
    bool _unwrapParticleCoordinates;
    bool _computeCentersOfMass;
    bool _computeRadiusOfGyration;

    std::shared_ptr<Particles::ParticleProperty> _particleClusters;
    std::shared_ptr<Particles::ParticleProperty> _unwrappedPositions;

    std::shared_ptr<Particles::ParticleProperty> _clusterSizes;
    std::shared_ptr<Particles::ParticleProperty> _clusterIDs;

    std::shared_ptr<Particles::ParticleProperty> _centersOfMass;
    std::shared_ptr<Particles::ParticleProperty> _radiiOfGyration;
    std::shared_ptr<Particles::ParticleProperty> _gyrationTensors;

    std::size_t _numClusters;
    std::size_t _largestClusterSize;
    bool _hasZeroWeightCluster;
};

}
