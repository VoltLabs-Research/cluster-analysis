#include <volt/cli/common.h>
#include <volt/cluster_analysis_service.h>

using namespace Volt;
using namespace Volt::CLI;

void showUsage(const std::string& name) {
    printUsageHeader(name, "Volt - Cluster Analysis");
    std::cerr
        << "  --cutoff <float>              Cutoff radius for neighbor search. [default: 3.2]\n"
        << "  --sortBySize                  Sort clusters by size (desc). [default: true]\n"
        << "  --unwrap                      Unwrap particle coordinates inside clusters. [default: false]\n"
        << "  --centersOfMass               Compute cluster centers (uniform weights). [default: false]\n"
        << "  --radiusOfGyration            Compute radii + tensors of gyration (uniform weights). [default: false]\n"
        << "  --threads <int>               Max worker threads (TBB/OMP). [default: auto]\n";
    printHelpOption();
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        showUsage(argv[0]);
        return 1;
    }

    std::string filename, outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);

    if (hasOption(opts, "--help") || filename.empty()) {
        showUsage(argv[0]);
        return filename.empty() ? 1 : 0;
    }

    auto parallel = initParallelism(opts, false);
    initLogging("volt-cluster-analysis", parallel.threads);

    LammpsParser::Frame frame;
    if (!parseFrame(filename, frame)) return 1;

    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);

    ClusterAnalysisService analyzer;
    analyzer.setCutoff(getDouble(opts, "--cutoff", 3.2));

    analyzer.setOptions(
        getBool(opts, "--sortBySize", true),
        getBool(opts, "--unwrap", false),
        getBool(opts, "--centersOfMass", false),
        getBool(opts, "--radiusOfGyration", false)
    );

    spdlog::info("Starting cluster analysis...");
    json result = analyzer.compute(frame, outputBase);

    if (result.value("is_failed", false)) {
        spdlog::error("Analysis failed: {}", result.value("error", "Unknown error"));
        return 1;
    }

    spdlog::info("Cluster analysis completed.");
    spdlog::info("Clusters: {}, largest size: {}",
        result.value("cluster_count", 0),
        result.value("largest_cluster_size", 0));

    return 0;
}
