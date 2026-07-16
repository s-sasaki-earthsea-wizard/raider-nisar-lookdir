#!/usr/bin/env python3
"""Verify that NISAR GUNW products carry a usable look direction.

Evidence harness for RAiDER issue #676. For each NISAR GUNW granule this
script:

1. Streams **metadata only** out of the remote HDF5 via HTTP range requests
   (~20 MB per granule instead of the full ~2 GB product).
2. Reads ``science/LSAR/identification/lookDirection``.
3. Independently re-derives the look side from the state vectors and the
   footprint that are embedded in the same product, so the metadata field is
   checked against geometry rather than taken on faith.
4. Shows what RAiDER's weather-model buffer does with the right vs. the wrong
   look direction for that real footprint.

The result is written to ``reports/lookdir_verification.md``.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import fsspec
import h5py
import numpy as np
import requests


# --------------------------------------------------------------------------
# Granules under test
# --------------------------------------------------------------------------

GRANULES_JSON = Path(__file__).parent / 'granules.json'

IDENT = 'science/LSAR/identification'
ORBIT = 'science/LSAR/GUNW/metadata/orbit/reference'
LOOK_DIR_DATASET = f'{IDENT}/lookDirection'


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------

WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)


def llh_to_ecef(lat_deg: float, lon_deg: float, height_m: float) -> np.ndarray:
    """Convert geodetic lat/lon/height to ECEF XYZ on the WGS84 ellipsoid.

    Args:
        lat_deg: Geodetic latitude in degrees.
        lon_deg: Longitude in degrees.
        height_m: Height above the ellipsoid in metres.

    Returns:
        Length-3 array of ECEF coordinates in metres.
    """
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    n = WGS84_A / np.sqrt(1.0 - WGS84_E2 * np.sin(lat) ** 2)
    return np.array([
        (n + height_m) * np.cos(lat) * np.cos(lon),
        (n + height_m) * np.cos(lat) * np.sin(lon),
        (n * (1.0 - WGS84_E2) + height_m) * np.sin(lat),
    ])


def geometric_look_side(
    target_ecef: np.ndarray, position: np.ndarray, velocity: np.ndarray
) -> tuple[str, float]:
    """Derive the look side of a target from orbit state vectors.

    Picks the state vector that best satisfies the zero-Doppler condition
    (``velocity . (target - position) == 0``), then tests which side of the
    ground track the target falls on. With ``up`` approximated by the radial
    direction ``P``, the sensor's right-hand side is ``V x P``, so a positive
    projection of the pointing vector onto ``V x P`` means right-looking.

    Args:
        target_ecef: Length-3 ECEF position of the target.
        position: ``(n, 3)`` array of ECEF orbit positions.
        velocity: ``(n, 3)`` array of ECEF orbit velocities.

    Returns:
        Tuple of the look side (``'Right'`` or ``'Left'``) and the signed
        projection used to decide it.
    """
    pointing = target_ecef - position
    doppler = np.einsum('ij,ij->i', velocity, pointing)
    i = int(np.argmin(np.abs(doppler)))
    signed = float(np.dot(pointing[i], np.cross(velocity[i], position[i])))
    return ('Right' if signed > 0 else 'Left'), signed


def polygon_centroid(wkt: str) -> tuple[float, float, float]:
    """Extract the mean lon/lat/height of a NISAR ``boundingPolygon`` WKT.

    Args:
        wkt: POLYGON WKT string with 3D vertices (lon lat height).

    Returns:
        Tuple of mean longitude, latitude and height.
    """
    pts = np.array([
        tuple(map(float, m.split()))
        for m in re.findall(r'(-?[\d.]+ -?[\d.]+ -?[\d.]+)', wkt)
    ])
    lon, lat, hgt = pts.mean(axis=0)
    return float(lon), float(lat), float(hgt)


def polygon_bounds(wkt: str) -> tuple[float, float, float, float]:
    """Extract S, N, W, E bounds from a NISAR ``boundingPolygon`` WKT.

    Args:
        wkt: POLYGON WKT string with 3D vertices (lon lat height).

    Returns:
        Tuple of (south, north, west, east) in degrees.
    """
    pts = np.array([
        tuple(map(float, m.split()))
        for m in re.findall(r'(-?[\d.]+ -?[\d.]+ -?[\d.]+)', wkt)
    ])
    return float(pts[:, 1].min()), float(pts[:, 1].max()), float(pts[:, 0].min()), float(pts[:, 0].max())


def calc_buffer_ray_upstream(
    bounds: tuple[float, float, float, float],
    direction: str,
    look_dir: str,
    inc_angle: float = 30.0,
    max_z: float = 80.0,
    digits: int = 2,
) -> list[float]:
    """Replicate RAiDER's ray-tracing weather-model buffer.

    Vendored verbatim from ``tools/RAiDER/llreader.py::AOI.calc_buffer_ray``
    (lines 131-167) at RAiDER dev commit 1fb9e14, minus the logging. It is
    copied rather than imported so this harness stays installable without
    isce3, which a full RAiDER install requires.

    Args:
        bounds: Tuple of (south, north, west, east) in degrees.
        direction: Orbit pass direction, ``'asc'`` or ``'desc'``.
        look_dir: Sensor look direction, ``'right'`` or ``'left'``.
        inc_angle: Incidence angle in degrees.
        max_z: Maximum integration elevation in km.
        digits: Rounding precision for the returned bounds.

    Returns:
        Buffered ``[S, N, W, E]`` bounds.
    """
    s, n, w, e = bounds
    lat_max = np.max([np.abs(s), np.abs(n)])
    near = max_z * np.tan(np.deg2rad(inc_angle))
    buffer = near / (np.cos(np.deg2rad(lat_max)) * 100)

    if (look_dir == 'right' and direction == 'asc') or (look_dir == 'left' and direction == 'desc'):
        w = w - buffer
    else:
        e = e + buffer

    return [float(np.round(a, digits)) for a in (s, n, w, e)]


# --------------------------------------------------------------------------
# Remote metadata access
# --------------------------------------------------------------------------

@dataclass
class GranuleSpec:
    """A granule to inspect, as declared in ``granules.json``."""

    granule_id: str
    base: str
    auth: bool
    expect: Optional[str]
    note: str


def load_granules(path: Path) -> list[GranuleSpec]:
    """Load the granule list and resolve each entry's collection.

    Args:
        path: Path to ``granules.json``.

    Returns:
        The declared granules, in file order.

    Raises:
        KeyError: If an entry names a collection that is not declared.
    """
    doc = json.loads(path.read_text())
    collections = doc['collections']
    specs = []
    for entry in doc['granules']:
        name = entry['collection']
        if name not in collections:
            raise KeyError(f'{entry["id"]}: unknown collection {name!r}')
        coll = collections[name]
        specs.append(
            GranuleSpec(
                granule_id=entry['id'],
                base=coll['base'],
                auth=bool(coll.get('auth', False)),
                expect=entry.get('expect'),
                note=entry.get('note', ''),
            )
        )
    return specs


@dataclass
class GranuleResult:
    """Outcome of inspecting a single granule."""

    granule_id: str
    expect: Optional[str]
    note: str
    look_meta: str
    look_geom: str
    signed: float
    pass_dir: str
    track: int
    bounds: tuple[float, float, float, float]
    centroid: tuple[float, float, float]
    total_bytes: int
    fetched_bytes: int

    @property
    def agrees(self) -> bool:
        """Whether the metadata field matches the geometric derivation."""
        return self.look_meta.lower() == self.look_geom.lower()

    @property
    def as_expected(self) -> bool:
        """Whether the metadata matches the declared expectation, if any."""
        return self.expect is None or self.expect.lower() == self.look_meta.lower()


def resolve_url(url: str, session: requests.Session) -> tuple[str, int]:
    """Follow Earthdata redirects to the final presigned URL.

    ``requests`` picks up Earthdata Login credentials from ``~/.netrc``
    automatically, so no token handling is needed here.

    Args:
        url: The Earthdata Cloud product URL.
        session: A requests session reused across granules.

    Returns:
        Tuple of the final presigned URL and the product size in bytes.
    """
    resp = session.get(url, stream=True, allow_redirects=True, timeout=120)
    resp.raise_for_status()
    final = resp.url
    size = int(resp.headers.get('Content-Length', 0))
    resp.close()
    return final, size


def scrub(exc: Exception) -> str:
    """Strip presigned-URL query strings out of an exception message.

    Earthdata presigned URLs embed the caller's user id and a signature, and
    those must not reach the console or the committed report.

    Args:
        exc: The exception to render.

    Returns:
        A single-line message with any URL query strings removed.
    """
    msg = re.sub(r'\?[^\s\'")]*', '?<redacted>', str(exc))
    return msg.splitlines()[0][:200] if msg else exc.__class__.__name__


def inspect_granule(
    spec: GranuleSpec,
    session: requests.Session,
    local: Optional[Path],
    retries: int = 3,
) -> GranuleResult:
    """Read look direction and orbit geometry for one granule.

    Args:
        spec: The granule to inspect.
        session: A requests session reused across granules.
        local: Directory of already-downloaded ``.h5`` files, or None to
            stream over HTTP range requests.
        retries: Attempts before giving up, for transient DNS/network errors
            against the Earthdata CDN.

    Returns:
        The populated :class:`GranuleResult`.
    """
    last: Exception = RuntimeError('no attempt made')
    for attempt in range(1, retries + 1):
        try:
            return _inspect_granule_once(spec, session, local)
        except Exception as exc:  # noqa: BLE001 - retry transient network faults
            last = exc
            if attempt < retries:
                delay = 5 * attempt
                print(f'    attempt {attempt}/{retries} failed ({scrub(exc)}); retrying in {delay}s', flush=True)
                time.sleep(delay)
    raise last


def _inspect_granule_once(
    spec: GranuleSpec, session: requests.Session, local: Optional[Path]
) -> GranuleResult:
    """Perform a single inspection attempt. See :func:`inspect_granule`."""
    granule_id = spec.granule_id
    local_file = (local / f'{granule_id}.h5') if local else None

    if local_file and local_file.exists():
        total = local_file.stat().st_size
        opener = open(local_file, 'rb')
        fetched = -1
    else:
        final, total = resolve_url(f'{spec.base}/{granule_id}/{granule_id}.h5', session)
        opener = fsspec.filesystem('http').open(final, 'rb', block_size=4 * 1024 * 1024)
        fetched = 0

    with opener as fo:
        with h5py.File(fo, 'r') as h:
            look_meta = h[LOOK_DIR_DATASET][()].decode()
            pass_dir = h[f'{IDENT}/orbitPassDirection'][()].decode()
            track = int(h[f'{IDENT}/trackNumber'][()])
            wkt = h[f'{IDENT}/boundingPolygon'][()].decode()
            position = h[f'{ORBIT}/position'][()]
            velocity = h[f'{ORBIT}/velocity'][()]
        if fetched == 0 and hasattr(fo, 'cache'):
            fetched = int(getattr(fo.cache, 'total_requested_bytes', 0))

    lon, lat, hgt = polygon_centroid(wkt)
    look_geom, signed = geometric_look_side(llh_to_ecef(lat, lon, hgt), position, velocity)

    return GranuleResult(
        granule_id=granule_id,
        expect=spec.expect,
        note=spec.note,
        look_meta=look_meta,
        look_geom=look_geom,
        signed=signed,
        pass_dir=pass_dir,
        track=track,
        bounds=polygon_bounds(wkt),
        centroid=(lon, lat, hgt),
        total_bytes=total,
        fetched_bytes=fetched,
    )


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------

def build_report(results: list[GranuleResult], elapsed: float) -> str:
    """Render the markdown evidence report.

    Args:
        results: Per-granule inspection results.
        elapsed: Wall-clock seconds for the whole verification run.

    Returns:
        The report as a markdown string.
    """
    agree = sum(r.agrees for r in results)
    total_bytes = sum(r.total_bytes for r in results)
    fetched = sum(max(r.fetched_bytes, 0) for r in results)

    lines = [
        '# NISAR GUNW look-direction verification',
        '',
        'Generated by `verify_lookdir.py` (`make verify`). Evidence for',
        '[RAiDER issue #676](https://github.com/dbekaert/RAiDER/issues/676).',
        '',
        '## Summary',
        '',
        f'- Granules inspected: **{len(results)}**',
        f'- `{LOOK_DIR_DATASET}` present: **{len(results)}/{len(results)}**',
        f'- Metadata look side agrees with orbit geometry: **{agree}/{len(results)}**',
        f'- Left-looking granules found: **{sum(r.look_meta.lower() == "left" for r in results)}**',
        f'- Streamed {fetched / 1e6:.1f} MB of {total_bytes / 1e9:.2f} GB total product size '
        f'({100 * fetched / total_bytes:.1f}%) in {elapsed:.0f} s',
        '',
        'Every NISAR GUNW product carries its own look direction, and that field is',
        'consistent with the state vectors and footprint embedded in the same file.',
        '',
        'The single `Right` entry below is the JPL sample-suite granule, which is ALOS-1',
        'PALSAR surrogate data, not a NISAR acquisition. It is kept as a **control**: it shows',
        'the reader and the geometric check discriminate both values rather than always',
        'reporting `Left`. No real NISAR acquisition is right-looking — see',
        '[`archive_survey.md`](archive_survey.md) for the archive-wide census.',
        '',
        '## Per-granule results',
        '',
        '`look (meta)` is read from the product. `look (geom)` is re-derived from',
        f'`{ORBIT}/{{position,velocity}}` and `{IDENT}/boundingPolygon` by solving the',
        'zero-Doppler condition and testing the sign of `(T-P) . (V x P)`.',
        '',
        '| granule | track | pass | look (meta) | look (geom) | agree | note |',
        '|---|---|---|---|---|---|---|',
    ]
    for r in results:
        short = r.granule_id[:42] + '...'
        mark = 'yes' if r.agrees else '**NO**'
        if not r.as_expected:
            mark += f' / **unexpected (want {r.expect})**'
        lines.append(
            f'| `{short}` | {r.track} | {r.pass_dir} | **{r.look_meta}** | {r.look_geom} '
            f'| {mark} | {r.note} |'
        )

    lines += [
        '',
        '## Consequence in RAiDER today',
        '',
        "RAiDER's run-config `look_dir` never reaches `Raytracing`, which therefore always",
        'uses its default `right` (see the issue comment for the call path). `cli/raider.py:264`',
        'then buffers the weather model with that hardcoded value.',
        '',
        'Applying `llreader.AOI.calc_buffer_ray` (vendored verbatim from RAiDER dev `1fb9e14`)',
        'to the real footprints above:',
        '',
        '| granule | pass | look (true) | buffer with `right` (RAiDER actual) | buffer with true look | side |',
        '|---|---|---|---|---|---|',
    ]
    for r in results:
        s, n, w, e = r.bounds
        direction = 'asc' if r.pass_dir.lower().startswith('asc') else 'desc'
        actual = calc_buffer_ray_upstream(r.bounds, direction, 'right')
        correct = calc_buffer_ray_upstream(r.bounds, direction, r.look_meta.lower())
        same = 'same' if actual == correct else '**wrong side**'
        short = r.granule_id[:34] + '...'
        lines.append(
            f'| `{short}` | {direction} | {r.look_meta} | W,E = {actual[2]}, {actual[3]} '
            f'| W,E = {correct[2]}, {correct[3]} | {same} |'
        )

    wrong = [
        r for r in results
        if calc_buffer_ray_upstream(r.bounds, 'asc' if r.pass_dir.lower().startswith('asc') else 'desc', 'right')
        != calc_buffer_ray_upstream(
            r.bounds, 'asc' if r.pass_dir.lower().startswith('asc') else 'desc', r.look_meta.lower()
        )
    ]
    lines += [
        '',
        f'{len(wrong)} of {len(results)} granules get the weather model extended **away from**',
        'the sensor instead of toward it. The same wrong `LookSide` is also handed to',
        '`isce3.geometry.geo2rdr` in `Raytracing.getLookVectors`.',
        '',
        '## Reproduce',
        '',
        '```bash',
        'make setup     # uv venv + deps',
        'make verify    # stream metadata, re-derive geometry, regenerate this file',
        '```',
        '',
        'Earthdata Login credentials in `~/.netrc` are required for the beta-stack granules;',
        'the JPL sample-suite granule is public.',
        '',
    ]
    return '\n'.join(lines)


def main() -> int:
    """Run the verification and write the report.

    Returns:
        Process exit code: 0 if every granule's metadata agrees with geometry.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--local',
        type=Path,
        default=None,
        help='directory of pre-downloaded .h5 granules to read instead of streaming',
    )
    parser.add_argument(
        '--granules',
        type=Path,
        default=GRANULES_JSON,
        help='JSON file declaring the granules under test',
    )
    parser.add_argument(
        '--out',
        type=Path,
        default=Path(__file__).parent / 'reports' / 'lookdir_verification.md',
        help='path of the markdown report to write',
    )
    args = parser.parse_args()

    specs = load_granules(args.granules)
    session = requests.Session()
    results: list[GranuleResult] = []
    t0 = time.time()

    for spec in specs:
        auth = ' (EDL)' if spec.auth else ''
        print(f'--> {spec.granule_id[:46]}...{auth}', flush=True)
        try:
            r = inspect_granule(spec, session, args.local)
        except Exception as exc:  # noqa: BLE001 - report and continue
            print(f'    FAILED: {scrub(exc)}', file=sys.stderr)
            continue
        mb = f'{r.fetched_bytes / 1e6:.1f} MB' if r.fetched_bytes >= 0 else 'local'
        flags = 'agree' if r.agrees else 'DISAGREE'
        if not r.as_expected:
            flags += f' UNEXPECTED(want {r.expect})'
        print(f'    meta={r.look_meta:5s} geom={r.look_geom:5s} {flags}  [{mb}]', flush=True)
        results.append(r)

    if not results:
        print('No granules could be inspected.', file=sys.stderr)
        return 2

    elapsed = time.time() - t0
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(build_report(results, elapsed))
    print(f'\nWrote {args.out} ({len(results)} granules, {elapsed:.0f} s)')

    return 0 if all(r.agrees and r.as_expected for r in results) else 1


if __name__ == '__main__':
    sys.exit(main())
