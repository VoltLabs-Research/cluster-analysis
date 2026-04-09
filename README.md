# ClusterAnalysis

`ClusterAnalysis` groups atoms into geometric clusters using a cutoff neighbor graph.

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
| `--sortBySize` | No | Sort clusters by descending size. | `true` |
| `--unwrap` | No | Unwrap coordinates inside each cluster. | `false` |
| `--centersOfMass` | No | Compute cluster centers of mass. | `false` |
| `--radiusOfGyration` | No | Compute radii and tensors of gyration. | `false` |
| `--threads <int>` | No | Maximum worker threads. | auto |
| `--help` | No | Print CLI help. | |

## Build With CoreToolkit

```bash
cd /path/to/voltlabs-ecosystem/tools/CoreToolkit
conan create . -nr

cd /path/to/voltlabs-ecosystem/plugins/ClusterAnalysis
conan create . -nr
```
