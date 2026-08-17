# Run 8 TBS adjacency matrices

## What this is

Filtered connectivity graphs extracted from the Super-Selective ECR
pipeline output for the "Run 8 TBS" organoid MEA experiment
(`Run_8_TBS_Experiment_ecr_results/`, recorded 2023-06-01/02, processed
2024-05-21).

**6 wells x 36 timepoints = 216 directed graphs.**

- **6 wells**: `well000` ... `well005` (0-indexed; corresponds to wells 1-6
  as labeled on the physical MEA plate).
- **36 timepoints**: one whole-plate recording each, in chronological
  order, spanning two recording days:
  - timepoints 0-17: 2023-06-01, folder `230601 RUN 8 Wells 1-3`
  - timepoints 18-35: 2023-06-02, folder `230602 RUN 8 Wells 4-6`

  Every recording captures **all 6 wells simultaneously** -- the folder
  names ("Wells 1-3" / "Wells 4-6") and per-file labels
  (e.g. "well #4 post stim 3") only describe which well was being
  electrically stimulated and what elapsed time that implies for that
  well; they do not restrict which wells' data are present in the file.
  `metadata/source_filenames` in the HDF5 file gives the exact source
  file for each timepoint index if that stimulation context is needed.

## File format

Single HDF5 file: `Run8_TBS_adjacency_matrices.h5`

```
/                                          (root attrs: n_wells, n_timepoints,
                                             fs_hz, created, source_dir, description)
/metadata/
    source_filenames                       (36,) string - relative path of the
                                             source .pkl for each timepoint index
    time_since_baseline_min                (36,) float - minutes elapsed since
                                             the first (baseline) recording

/well000 .. /well005/
    num_vertices                           (36,) int - electrode/node count
                                             for that well at each timepoint
    adjacency/t00 .. t35                   (n_i, n_i) bool - directed adjacency
                                             matrix at timepoint i (see below)
    channel_numbers/t00 .. t35             (n_i,) int - electrode/channel ID
                                             for each row/column of that
                                             timepoint's adjacency matrix
```

**Why ragged, not one 3D array**: the number of active electrodes (and
which physical channels they are) differs per well and per timepoint, so
matrices are NOT padded to a common size and are NOT guaranteed to share
node identity across timepoints. Use `channel_numbers/tXX` to align nodes
across timepoints for the same well (e.g. before graph matching / MDS,
as done in the analysis notebooks) -- do not assume row i means the same
electrode at every timepoint.

## What "adjacency matrix" means here

- Square boolean matrix, `adj[i, j] = True` means a predicted directed
  functional link from electrode i to electrode j.
- Built from `adj_matrix_predicted` in the raw ECR output, then
  link-filtered: a link is removed if every correlation-peak delay
  between that channel pair was smaller than one sampling period
  (`1 / fs`), which flags links more likely explained by a shared
  global/simultaneous event than a direct connection. `fs` (sampling
  rate, Hz) is stored as a root-level HDF5 attribute.

## Minimal read example

```python
import h5py

with h5py.File("Run8_TBS_adjacency_matrices.h5", "r") as h5:
    adj = h5["well003/adjacency/t05"][:]      # (n, n) bool matrix, well003, timepoint 5
    ch = h5["well003/channel_numbers/t05"][:] # electrode IDs for that matrix's rows/cols
    t_min = h5["metadata/time_since_baseline_min"][5]  # minutes since baseline
```
