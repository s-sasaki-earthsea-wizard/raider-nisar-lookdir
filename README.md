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
That is now testable. This harness shows, against products live in the archive:

1. Every NISAR GUNW carries `science/LSAR/identification/lookDirection`
   (`'Left'` / `'Right'`), so RAiDER never has to guess.
2. **Left-looking NISAR GUNWs exist today.** The ASF beta stack over central
   Spain (track 159, ascending) is left-looking on all 7 granules.
3. The metadata is not taken on faith: the look side is re-derived from the
   state vectors and footprint embedded in the same product, and agrees on
   8/8 granules.
4. With the hardcoded `right`, RAiDER buffers the weather model on the wrong
   side of the scene for those granules.

The generated report is [`reports/lookdir_verification.md`](reports/lookdir_verification.md).

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
make verify   # stream metadata, re-derive geometry, regenerate the report
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
- `Makefile` — `setup` / `verify` / `verify-local` / `clean`.
- `verify_lookdir.py` — the harness; streams metadata, checks geometry, writes the report.
- `reports/lookdir_verification.md` — committed evidence (tracked).
- `.venv/`, `data/` — untracked.
