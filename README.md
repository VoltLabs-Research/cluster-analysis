# Cluster Analysis

Groups atoms into clusters by distance (cutoff) and optionally computes centers of mass and radius of gyration.

## Install

```bash
vpm install @voltlabs/cluster-analysis
```

## CLI

```bash
cluster-analysis <input_dump> [output_base] [options]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `<input_dump>` | yes | — | Input LAMMPS dump. |
| `[output_base]` | no | derived from input | Base path for output files. |
| `--cutoff <float>` | no | `3.2` | Cutoff radius for neighbor search. |
| `--sort_by_size` | no | `true` | Sort clusters by descending size. |
| `--unwrap` | no | `false` | Unwrap coordinates inside each cluster. |
| `--centers_of_mass` | no | `false` | Compute cluster centers of mass. |
| `--radius_of_gyration` | no | `false` | Compute radii and tensors of gyration. |
| `--threads <int>` | no | auto | Maximum worker threads. |

## Exports

| Output file | Exposure | Exporter → artifact |
|---|---|---|
| `{output_base}_cluster_analysis.parquet` | Cluster Analysis | — (listing-only) |
| `{output_base}_atoms.parquet` | Cluster Model | AtomisticExporter → glb |

---

Full input contract and examples: https://docs.voltcloud.dev/docs/plugins/cluster-analysis
