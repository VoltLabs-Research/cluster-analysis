#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <volt/core/particle_property.h>
#include <volt/cluster_analysis_engine.h>
#include <volt/plugin/output_files.h>
#include <nlohmann/json.hpp>
#include <memory>
#include <string>

namespace Volt{
using json = nlohmann::json;

class ClusterAnalysisLog;

class ClusterAnalysisService{
public:
    // What the files are written from: where they go, the cluster table, the clusters the engine built, and the log
    // that names each file once it is written.
    struct Outputs{
        const std::string& base;
        const LammpsParser::Frame& frame;
        const ClusterAnalysisEngine& engine;
        const json& clusterTable;
        ClusterAnalysisLog& log;
    };

    ClusterAnalysisService();

    void setCutoff(double cutoff);
    void setSortBySize(bool sortBySize);
    void setUnwrapParticleCoordinates(bool unwrapParticleCoordinates);
    void setComputeCentersOfMass(bool computeCentersOfMass);
    void setComputeRadiusOfGyration(bool computeRadiusOfGyration);

    json compute(const LammpsParser::Frame& frame, const std::string &outputFilename);

private:
    // The files in the order they are written. None is optional, so no flag of the service asks for them.
    static const OutputFile<ClusterAnalysisService, Outputs> outputFiles[];

    double _cutoff;
    bool _sortBySize;
    bool _unwrapParticleCoordinates;
    bool _computeCentersOfMass;
    bool _computeRadiusOfGyration;
};

}
