# raider-nisar-lookdir

Real-data evidence harness for [RAiDER issue #676](https://github.com/dbekaert/RAiDER/issues/676)
— *"Get look direction parameter ready for NISAR"*.

This is a **sibling project** to a [RAiDER](https://github.com/dbekaert/RAiDER)
fork: it lives inside the fork tree but is a fully independent git repository,
excluded from the fork history via `.git/info/exclude`. Nothing here is part of
the upstream contribution — it is the *evidence* that backs it.

## What it establishes

Issue #676 asked whether RAiDER's `look_dir` run-config parameter is wired up.
It isn't: `look_dir` is validated in `cli/raider.py`, stored on `RunConfig`, and
never read again. `Raytracing` is constructed in two places in
`cli/validators.py::get_los`, neither of which forwards it, so it always falls
back to its default `right`.

In 2024 the maintainer's reply was that this was latent — *"Right now everything
is right-looking, but very soon we will have a left-looking sensor (NISAR)"*.
That is now testable, and the answer is stronger than "soon": against products
live in the archive, this harness shows

1. Every NISAR GUNW carries `science/LSAR/identification/lookDirection`
   (`'Left'` / `'Right'`), so RAiDER never has to guess.
2. **NISAR is a left-looking mission.** All 221 (track, pass) groups in the
   `NISAR_L2_GUNW_BETA_V1` collection — 10,217 granules — are left-looking, as
   are sampled GSLC/GCOV/RSLC products. Not one right-looking acquisition
   exists. So the hardcoded `right` is wrong for *every* NISAR product.
3. The metadata is not taken on faith: the look side is re-derived from the
   state vectors and footprint embedded in the same product, and agrees on
   10/10 granules inspected in depth.
4. With the hardcoded `right`, RAiDER buffers the weather model on the wrong
   side of the scene for all of them.

Two reports are generated:

- [`reports/lookdir_verification.md`](reports/lookdir_verification.md) — depth:
  10 granules, metadata vs. re-derived geometry, and the buffer consequence.
- [`reports/archive_survey.md`](reports/archive_survey.md) — breadth: the
  look-direction census of the whole beta archive.

The one right-looking product in the verification set is the JPL sample-suite
granule, which is ALOS-1 PALSAR surrogate data rather than a NISAR acquisition.
It is kept deliberately, as the control that shows the reader and the geometric
check discriminate both values instead of always reporting `Left`.

Granules under test are declared in [`granules.json`](granules.json) — add
entries there rather than editing the script.

## How it works

Downloading the stack would mean ~15 GB. Instead the harness resolves the
Earthdata presigned URL and lets `h5py` read the remote file through an
`fsspec` HTTP handle, so only the HDF5 metadata blocks are actually fetched —
**~298 MB of 14.96 GB (2.0%), in about 90 seconds.**

The one piece of RAiDER logic under test, `llreader.AOI.calc_buffer_ray`, is
vendored verbatim (with its upstream provenance recorded in the docstring)
rather than imported, so this harness installs without isce3.

## Usage

```bash
make setup    # uv venv (.venv) + deps from pyproject.toml
make verify   # stream metadata, re-derive geometry, regenerate the report (~2 min)
make survey   # census the whole beta archive's look direction (~6 min)
```

`make` with no target lists the available targets.

Earthdata Login credentials are required for the beta-stack granules — put them
in `~/.netrc`:

```
machine urs.earthdata.nasa.gov login <user> password <pass>
```

The JPL sample-suite granule is public and needs no credentials.

If you already have the products on disk, skip the network entirely:

```bash
make verify-local DIR=/path/to/granules
```

## Layout

- `pyproject.toml` — deps for the isolated `.venv` (resolved with `uv`).
- `Makefile` — `setup` / `verify` / `verify-local` / `survey` / `clean`.
- `granules.json` — the granules under test, with per-entry `expect` and notes.
- `verify_lookdir.py` — depth: streams metadata, checks geometry, writes the report.
- `survey_archive.py` — breadth: CMR enumeration + archive-wide look-direction census.
- `reports/` — committed evidence (tracked).
- `.venv/`, `data/` — untracked.
