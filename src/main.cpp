#include <volt/cli/common.h>
#include <volt/cluster_analysis_service.h>
#include <volt/plugin/option_reader.h>
#include <oneapi/tbb/global_control.h>

using namespace Volt;
using namespace Volt::CLI;
using namespace Volt::Plugin;

static PluginDescriptor buildDescriptor() {
    return {
        "volt-cluster-analysis",
        "Cluster Analysis",
        {
            {"--cutoff", "float", "Cutoff radius for neighbor search.", "3.2", {}, ""},
            {"--sort_by_size", "bool", "Sort clusters by size (desc).", "true", {}, ""},
            {"--unwrap", "bool", "Unwrap particle coordinates inside clusters.", "false", {}, ""},
            {"--centers_of_mass", "bool", "Compute cluster centers (uniform weights).", "false", {}, ""},
            {"--radius_of_gyration", "bool",
             "Compute radii + tensors of gyration (uniform weights).", "false", {}, ""},
        }
    };
}

int main(int argc, char* argv[]) {
    const PluginDescriptor descriptor = buildDescriptor();

    if (argc < 2) {
        showPluginUsage(argv[0], descriptor);
        return 1;
    }

    std::string filename, outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);

    if (auto exitCode = handleIntrospection(argv[0], descriptor, opts, filename)) {
        return *exitCode;
    }

    const int requestedThreads = std::max(1, getInt(opts, "--threads", std::thread::hardware_concurrency() > 0
        ? static_cast<int>(std::thread::hardware_concurrency())
        : 1));
    oneapi::tbb::global_control parallelControl(
        oneapi::tbb::global_control::max_allowed_parallelism,
        static_cast<std::size_t>(requestedThreads)
    );
    initLogging("volt-cluster-analysis");
    spdlog::info("Using {} threads (OneTBB)", requestedThreads);

    LammpsParser::Frame frame;
    if (!parseFrame(filename, frame)) return 1;

    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);

    const OptionReader options(descriptor, opts);

    ClusterAnalysisService analyzer;
    analyzer.setCutoff(options.number("--cutoff"));

    analyzer.setOptions(
        options.boolean("--sort_by_size"),
        options.boolean("--unwrap"),
        options.boolean("--centers_of_mass"),
        options.boolean("--radius_of_gyration")
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
