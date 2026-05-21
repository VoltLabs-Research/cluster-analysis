# ClusterAnalysis

`ClusterAnalysis` groups atoms into geometric clusters using a cutoff neighbor graph.

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- ClusterAnalysis
```

## Build from source

Requires [Conan 2.x](https://docs.conan.io/2/installation.html), CMake 3.20+, and a C++23 compiler (GCC 14+ or Clang 17+).

### Prerequisites

The following Conan packages must be available in your local cache:

- `coretoolkit/1.0.0` (from the `CoreToolkit` repository)

For each dependency, clone its repository and create the package:

```bash
conan create <path-to-dependency-repo> --build=missing -o "hwloc/*:shared=True"
```

### Build

From the root of this repository:

```bash
conan install . -of build --build=missing -o "hwloc/*:shared=True"
cmake --preset conan-release
cmake --build build/build/Release -j
```

### Run

```bash
./build/build/Release/cluster-analysis --help
```

### Package as Conan recipe

To make this plugin available as a Conan package for other projects:

```bash
conan create . --build=missing -o "hwloc/*:shared=True"
```

## CLI

Usage:

```bash
cluster-analysis <lammps_file> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file. | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--cutoff <float>` | No | Cutoff radius for neighbor search. | `3.2` |
| `--sort_by_size` | No | Sort clusters by descending size. | `true` |
| `--unwrap` | No | Unwrap coordinates inside each cluster. | `false` |
| `--centers_of_mass` | No | Compute cluster centers of mass. | `false` |
| `--radius_of_gyration` | No | Compute radii and tensors of gyration. | `false` |
| `--threads <int>` | No | Maximum worker threads. | auto |
| `--help` | No | Print CLI help. | |
