#include <volt/plugin/plugin_entry.h>
#include <volt/cluster_analysis_service.h>

using namespace Volt;
using namespace Volt::Plugin;
using S = ClusterAnalysisService;

static const std::vector<OptionBinding<S>> bindings = {
    opt("--cutoff",             "Cutoff radius for neighbor search.",                     3.2,   &S::setCutoff),
    opt("--sort_by_size",       "Sort clusters by size (desc).",                          true,  &S::setSortBySize),
    opt("--unwrap",             "Unwrap particle coordinates inside clusters.",           false, &S::setUnwrapParticleCoordinates),
    opt("--centers_of_mass",    "Compute cluster centers (uniform weights).",             false, &S::setComputeCentersOfMass),
    opt("--radius_of_gyration", "Compute radii + tensors of gyration (uniform weights).", false, &S::setComputeRadiusOfGyration),
};

VOLT_SERVICE_PLUGIN("volt-cluster-analysis", "Cluster Analysis", S, bindings)
