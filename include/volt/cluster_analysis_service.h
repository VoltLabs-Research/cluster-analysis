#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <volt/core/particle_property.h>
#include <volt/cluster_analysis_engine.h>
#include <nlohmann/json.hpp>
#include <memory>
#include <string>

namespace Volt{
using json = nlohmann::json;

class ClusterAnalysisService{
public:
    ClusterAnalysisService();

    void setCutoff(double cutoff);
    void setOptions(
        bool sortBySize,
        bool unwrapParticleCoordinates,
        bool computeCenterOfMass,
        bool computeRadiusOfGyration
    );

    json compute(const LammpsParser::Frame& frame, const std::string &outputFilename);

private:
    double _cutoff;
    bool _sortBySize;
    bool _unwrapParticleCoordinates;
    bool _computeCentersOfMass;
    bool _computeRadiusOfGyration;
};

}
