"""
maketopo_source_selatan.py
===========================

Generates topography files (basal, mass fraction, surface/eta) for the
"Source Selatan_1" lahar scenario on the DEM domain (DEMPredict.tt3), and
produces verification plots (source area, hillshade, surface topo, and
initial lahar depth).

Pipeline: read DEM -> read shapefiles -> define basal/mfrac/eta ->
write *.tt3 files -> generate verification plots.
"""

# ---------------------------------------------------------------------------
# Standard / third-party imports
# ---------------------------------------------------------------------------
import os
import struct

import numpy as np
from pylab import *  # noqa: F401,F403  (kept as in the original script)
from scipy.interpolate import RegularGridInterpolator
from shapely.geometry import Polygon

# Clawpack imports
from clawpack.geoclaw import topotools


# ---------------------------------------------------------------------------
# Resolve path relative to the script directory
# ---------------------------------------------------------------------------
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def _abspath(filename):
    """Return an absolute path relative to this script's directory."""
    if os.path.isabs(filename):
        return filename
    return os.path.join(_SCRIPT_DIR, filename)


# ---------------------------------------------------------------------------
# FILE PATHS
# ---------------------------------------------------------------------------
DEM_PATH         = _abspath("DEMPredict.tt3")
SHAPEFILE_PATH   = _abspath("Source Selatan_1/Source Selatan_1.shp")   # lahar SOURCE (required)
OVERLAY_SHP_PATH = _abspath("Source Selatan/Source Selatan.shp")       # overlay layer (dashed blue)


# ---------------------------------------------------------------------------
# Lahar parameters
# ---------------------------------------------------------------------------
m0 = 0.59   # mass fraction inside the source area
h0 = 4.11   # initial lahar depth (meters)


# ---------------------------------------------------------------------------
# DEM metadata — DEMPredict.tt3
# ---------------------------------------------------------------------------
_NODATA_DEM        = -9999.0
_CELLSIZE_EXPECTED = 8.29   # nominal value ~= 8.289958903244 m

_PLOT_DPI = 200
_MAX_INCH = 20.0


# ---------------------------------------------------------------------------
# 1. Read DEM ASCII (tt3 format, 6-line header)
# ---------------------------------------------------------------------------
def _read_dem_ascii(path):
    """Read an ASCII tt3-format DEM file (6-line header) and return the grid."""
    with open(path, "r") as fh:
        ncols     = int(fh.readline().split()[0])
        nrows     = int(fh.readline().split()[0])
        xllcorner = float(fh.readline().split()[0])
        yllcorner = float(fh.readline().split()[0])
        cellsize  = float(fh.readline().split()[0])
        nodata    = float(fh.readline().split()[0])
        data      = np.loadtxt(fh)

    if abs(cellsize - _CELLSIZE_EXPECTED) > 0.1:
        print(f"  [WARNING] cellsize={cellsize:.6f} m, expected ~{_CELLSIZE_EXPECTED} m")

    Z = data.astype(float)
    Z[Z == nodata]      = np.nan
    Z[Z == _NODATA_DEM] = np.nan

    # FIX #16: flip Z so that row 0 = SOUTH
    Z = Z[::-1, :]

    # FIX #17: convert pixel edge -> pixel center
    half = cellsize / 2.0
    x1d  = (xllcorner + half) + np.arange(ncols) * cellsize
    y1d  = (yllcorner + half) + np.arange(nrows) * cellsize

    return x1d, y1d, Z, ncols, nrows, xllcorner, yllcorner, cellsize


(_x1d_dem, _y1d_dem, _Z_dem,
 _ncols, _nrows, _xll, _yll, _cellsize) = _read_dem_ascii(DEM_PATH)

xlower = float(_x1d_dem[0])
xupper = float(_x1d_dem[-1])
ylower = float(_y1d_dem[0])
yupper = float(_y1d_dem[-1])

nxpoints = int(round((xupper - xlower) / _cellsize)) + 1
nypoints = int(round((yupper - ylower) / _cellsize)) + 1

if nxpoints != _ncols:
    print(f"  [WARNING] nxpoints={nxpoints} != ncols={_ncols}")
if nypoints != _nrows:
    print(f"  [WARNING] nypoints={nypoints} != nrows={_nrows}")

# Domain verification
_xlower_ref = 427610.433
_xupper_ref = 444870.127
_ylower_ref = 9153774.650
_yupper_ref = 9171266.463
_tol = 0.5

for _name, _val, _ref in [
    ("xlower", xlower, _xlower_ref),
    ("xupper", xupper, _xupper_ref),
    ("ylower", ylower, _ylower_ref),
    ("yupper", yupper, _yupper_ref),
]:
    diff = abs(_val - _ref)
    if diff > _tol:
        print(f"  [WARNING] {_name}={_val:.3f} vs reference {_ref:.3f} "
              f"(diff={_val - _ref:.3f} m > tol={_tol:.3f} m)")
    else:
        print(f"  [OK] {_name}={_val:.3f} ~= {_ref:.3f}  (diff={_val - _ref:.4f} m) OK")

print("=" * 60)
print(f"[DEM] File      : {DEM_PATH}")
print(f"[DEM] Dimension : {_ncols} x {_nrows} pixels")
print(f"[DEM] Cellsize  : {_cellsize:.6f} m  (Pixel is Area -> pixel center)")
print(f"[DEM] X center  : {xlower:.3f} -> {xupper:.3f}")
print(f"[DEM] Y center  : {ylower:.3f} -> {yupper:.3f}")
print(f"[DEM] Elev min  : {np.nanmin(_Z_dem):.1f} m")
print(f"[DEM] Elev max  : {np.nanmax(_Z_dem):.1f} m")
nan_pct = np.sum(np.isnan(_Z_dem)) / _Z_dem.size * 100
print(f"[DEM] NaN       : {nan_pct:.2f}%")
print(f"[Grid output]   : {nxpoints} x {nypoints} points")
dx_check = (xupper - xlower) / (nxpoints - 1)
dy_check = (yupper - ylower) / (nypoints - 1)
print(f"[Grid output]   : dx={dx_check:.6f} m  dy={dy_check:.6f} m")
print(f"[FIX #16]       : Z flipped -> row 0 = SOUTH <-> y1d[0] = {ylower:.3f} m OK")
print(f"[FIX #17]       : xllcorner {_xll:.3f} + half {_cellsize / 2:.4f} = {xlower:.3f} m OK")
print("=" * 60)

_dem_interp = RegularGridInterpolator(
    (_y1d_dem, _x1d_dem),
    _Z_dem,
    method       = "nearest",
    bounds_error = False,
    fill_value   = np.nan,
)


# ---------------------------------------------------------------------------
# 2. Read shapefile
# ---------------------------------------------------------------------------
def _read_shp_polygon(path):
    """Read a polygon shapefile -> list of shapely.Polygon. Supports type 5 & 15."""
    polygons = []
    with open(path, 'rb') as f:
        f.read(4); f.read(20); f.read(4); f.read(4)
        shape_type = struct.unpack('<i', f.read(4))[0]
        f.read(64)

        if shape_type not in (5, 15):
            raise ValueError(
                f"Shape type {shape_type} is not Polygon (5) or PolygonZ (15).")

        while True:
            hdr = f.read(8)
            if len(hdr) < 8:
                break
            content_len = struct.unpack('>i', hdr[4:])[0] * 2
            rec_start   = f.tell()
            rec_stype   = struct.unpack('<i', f.read(4))[0]
            if rec_stype == 0:
                f.seek(rec_start + content_len)
                continue
            f.read(32)
            num_parts  = struct.unpack('<i', f.read(4))[0]
            num_points = struct.unpack('<i', f.read(4))[0]
            parts = [struct.unpack('<i', f.read(4))[0] for _ in range(num_parts)]
            pts   = [struct.unpack('<2d', f.read(16)) for _ in range(num_points)]
            remaining = content_len - (4 + 32 + 4 + 4 + num_parts * 4 + num_points * 16)
            if remaining > 0:
                f.read(remaining)
            parts.append(num_points)

            rings = []
            for i in range(num_parts):
                ring = pts[parts[i]:parts[i + 1]]
                if len(ring) < 3:
                    continue
                if ring[0] != ring[-1]:
                    ring = ring + [ring[0]]
                rings.append(ring)

            if rings:
                try:
                    poly = Polygon(
                        shell=rings[0],
                        holes=rings[1:] if len(rings) > 1 else []
                    )
                    if not poly.is_empty and poly.is_valid:
                        polygons.append(poly)
                    elif not poly.is_empty:
                        polygons.append(poly.buffer(0))
                except Exception as e:
                    print(f"  [WARNING] Failed to build Polygon: {e}")
    return polygons


def _load_shapefile(path, label=""):
    """Read a shapefile with geopandas (if available), fallback to internal parser."""
    if not os.path.isfile(path):
        print(f"  [WARNING] {label}: file not found -> '{path}'")
        return [], False

    try:
        import geopandas as gpd
        gdf   = gpd.read_file(path)
        polys = [g for g in gdf.geometry if g is not None and not g.is_empty]
        print(f"[{label}] Read using   : geopandas  ({len(polys)} features)")
        print(f"[{label}] CRS          : {gdf.crs}")
        return polys, True
    except ImportError:
        pass
    except Exception as e:
        print(f"  [WARNING] {label}: geopandas failed: {e}  -> trying internal parser")

    try:
        polys = _read_shp_polygon(path)
        print(f"[{label}] Read using   : internal SHP parser  ({len(polys)} features)")
        prj_path = path.replace(".shp", ".prj")
        try:
            with open(prj_path) as _pf:
                print(f"[{label}] CRS (prj)    : {_pf.read().strip()[:120]}...")
        except FileNotFoundError:
            print(f"[{label}] CRS (prj)    : not found")
        return polys, True
    except Exception as e:
        print(f"  [WARNING] {label}: internal SHP parser failed: {e}")
        return [], False


# ---------------------------------------------------------------------------
# lahar SOURCE — Source Selatan_1 (REQUIRED)
# ---------------------------------------------------------------------------
if not os.path.isfile(SHAPEFILE_PATH):
    raise FileNotFoundError(
        f"\n[ERROR] Lahar SOURCE shapefile not found: '{SHAPEFILE_PATH}'\n"
        f"  Make sure the .shp file (along with .dbf, .shx, .prj) exists\n"
        f"  inside the 'Source Selatan_1/' folder next to this script,\n"
        f"  or update SHAPEFILE_PATH.\n"
        f"  Expected file: 'Source Selatan_1.shp'"
    )

_polys_source, _source_ok = _load_shapefile(SHAPEFILE_PATH, label="Source")
_poly_source = (
    _polys_source[0]
    if len(_polys_source) == 1
    else _polys_source[0].union(*_polys_source[1:])
)

b = _poly_source.bounds
print(f"[Source] Shapefile : {SHAPEFILE_PATH}")
print(f"[Source] Bounds X  : {b[0]:.2f} -> {b[2]:.2f} m")
print(f"[Source] Bounds Y  : {b[1]:.2f} -> {b[3]:.2f} m")
print(f"[Source] Area      : {_poly_source.area:.2f} m^2")
print("=" * 60)

if b[0] < xlower or b[2] > xupper or b[1] < ylower or b[3] > yupper:
    print("  [WARNING] Part or all of the SOURCE shapefile is outside the DEM domain!")
    print(f"  DEM domain X : {xlower:.3f} -> {xupper:.3f}")
    print(f"  DEM domain Y : {ylower:.3f} -> {yupper:.3f}")
else:
    print("  [OK] SOURCE shapefile is inside the DEM domain OK")
print("=" * 60)

# ---------------------------------------------------------------------------
# OVERLAY shapefile — Source Selatan (visualization only, dashed blue)
# ---------------------------------------------------------------------------
_polys_overlay, _overlay_ok = _load_shapefile(OVERLAY_SHP_PATH, label="Overlay")

if _overlay_ok and _polys_overlay:
    bo = _polys_overlay[0].bounds
    for p in _polys_overlay[1:]:
        bx = p.bounds
        bo = (min(bo[0], bx[0]), min(bo[1], bx[1]),
              max(bo[2], bx[2]), max(bo[3], bx[3]))
    print(f"[Overlay] Shapefile : {OVERLAY_SHP_PATH}")
    print(f"[Overlay] Count     : {len(_polys_overlay)} features")
    print(f"[Overlay] Bounds X  : {bo[0]:.2f} -> {bo[2]:.2f} m")
    print(f"[Overlay] Bounds Y  : {bo[1]:.2f} -> {bo[3]:.2f} m")
    print(f"[Overlay] Note      : visualization layer only, NOT a separate lahar source")
else:
    print("  [INFO] Overlay shapefile not loaded — plotting continues without overlay.")
print("=" * 60)


def _plot_overlay(ax, linewidth=1.5, label="Source Selatan (overlay)"):
    """
    Draw all overlay shapefile features onto axes ax as a dashed blue
    line. No fill is applied. Safe to call even if the overlay is
    unavailable.
    """
    from shapely.geometry import MultiPolygon

    if not _overlay_ok or not _polys_overlay:
        return

    first = True
    for poly in _polys_overlay:
        geoms = list(poly.geoms) if isinstance(poly, MultiPolygon) else [poly]
        for geom in geoms:
            xo, yo = geom.exterior.xy
            ax.plot(
                xo, yo,
                color     = "blue",
                linewidth = linewidth,
                linestyle = "--",
                label     = label if first else "_nolegend_",
                zorder    = 9,
            )
            first = False
            for interior in geom.interiors:
                xi, yi = interior.xy
                ax.plot(xi, yi, color="blue", linewidth=linewidth,
                        linestyle="--", label="_nolegend_", zorder=9)


def _in_source(x_arr, y_arr):
    """Boolean array: True if the point (x, y) lies inside the SOURCE lahar polygon."""
    try:
        from shapely import contains_xy
        return contains_xy(_poly_source, x_arr, y_arr)
    except ImportError:
        from shapely.vectorized import contains
        return contains(_poly_source, x_arr, y_arr)


# ---------------------------------------------------------------------------
# 3. Basal topography function
# ---------------------------------------------------------------------------
def basal(x, y):
    """Basal topography elevation (meters). NaN filled with 0.0 so the solver doesn't crash."""
    pts = np.column_stack([np.ravel(y), np.ravel(x)])
    z   = _dem_interp(pts)
    z   = np.where(np.isnan(z), 0.0, z)
    return z.reshape(np.shape(x))


# ---------------------------------------------------------------------------
# 4. Mass fraction function
# ---------------------------------------------------------------------------
def mfrac(x, y):
    """Mass fraction: m0 inside the SOURCE lahar, 0 outside."""
    x_flat = np.ravel(x).astype(float)
    y_flat = np.ravel(y).astype(float)
    mask   = _in_source(x_flat, y_flat)
    return np.where(mask, m0, 0.0).reshape(np.shape(x))


# ---------------------------------------------------------------------------
# 5. Surface topo / eta function
# ---------------------------------------------------------------------------
def eta(x, y):
    """Upper material surface: basal + h0 inside the source, basal only outside."""
    B  = basal(x, y)
    mf = mfrac(x, y)
    return B + np.where(mf > 0, h0, 0.0)


# ---------------------------------------------------------------------------
# 6. Info
# ---------------------------------------------------------------------------
print("Maximum landslide depth: %.2f m" % h0)


# ---------------------------------------------------------------------------
# 7. Helper: write topography to a tt3 file
# ---------------------------------------------------------------------------
def _write_tt3(outfile, func, xlower, xupper, ylower, yupper, nx, ny):
    """Thin wrapper around topotools.topo3writer."""
    topotools.topo3writer(
        outfile, func,
        xlower, xupper,
        ylower, yupper,
        nx, ny,
    )


# ---------------------------------------------------------------------------
# Helper: compute figsize proportional to the UTM domain
# ---------------------------------------------------------------------------
def _domain_figsize(extra_ratio=1.0):
    """Compute figure size (inches) proportional to the DEM domain dimensions."""
    domain_w_px = nxpoints
    domain_h_px = nypoints
    scale = (_MAX_INCH * _PLOT_DPI) / max(domain_w_px, domain_h_px)
    fig_w = (domain_w_px * scale / _PLOT_DPI) * extra_ratio
    fig_h =  domain_h_px * scale / _PLOT_DPI
    return fig_w, fig_h


# ---------------------------------------------------------------------------
# 8. Write array directly to tt3
# ---------------------------------------------------------------------------
def _write_tt3_from_array(outfile, Z_southup, x1d, y1d, cellsize, nodata=-9999.0):
    """Write a Z array directly to a tt3 file (FIX #17)."""
    nrows, ncols = Z_southup.shape
    assert ncols == len(x1d), f"ncols mismatch: {ncols} vs {len(x1d)}"
    assert nrows == len(y1d), f"nrows mismatch: {nrows} vs {len(y1d)}"

    Z_write = Z_southup.copy()
    Z_write[np.isnan(Z_write)] = 0.0

    Z_file = Z_write[::-1, :]

    half = cellsize / 2.0
    xllcorner_out = x1d[0] - half
    yllcorner_out = y1d[0] - half

    with open(outfile, "w") as fh:
        fh.write(f"{ncols}\t\t\tncols\n")
        fh.write(f"{nrows}\t\t\tnrows\n")
        fh.write(f"{xllcorner_out:.10f}\t\txllcorner\n")
        fh.write(f"{yllcorner_out:.10f}\t\tyllcorner\n")
        fh.write(f"{cellsize:.10f}\t\tcellsize\n")
        fh.write(f"{nodata:.1f}\t\t\tnodata_value\n")
        for row in Z_file:
            fh.write("  ".join(f"{v:.4f}" for v in row) + "\n")

    print(f"  [OK] {outfile}")
    print(f"       ncols={ncols}  nrows={nrows}  cs={cellsize:.6f} m")
    print(f"       xllcorner (edge) = {xllcorner_out:.3f}")
    print(f"       yllcorner (edge) = {yllcorner_out:.3f}")
    print(f"       X center: {x1d[0]:.3f} -> {x1d[-1]:.3f}  ({x1d[-1] - x1d[0]:.1f} m)")
    print(f"       Y center: {y1d[0]:.3f} -> {y1d[-1]:.3f}  ({y1d[-1] - y1d[0]:.1f} m)")
    print(f"       Z: [{np.nanmin(Z_southup):.1f}, {np.nanmax(Z_southup):.1f}] m")


# ---------------------------------------------------------------------------
# 9. Main functions
# ---------------------------------------------------------------------------
def maketopo():
    """Write basal_topo.tt3 and mass_frac.tt3."""
    print("=" * 60)
    print("[maketopo] Writing basal_topo.tt3 (direct copy from DEMPredict)")
    _write_tt3_from_array("basal_topo.tt3", _Z_dem, _x1d_dem, _y1d_dem, _cellsize)

    print("[maketopo] Writing mass_frac.tt3")
    _write_tt3("mass_frac.tt3", mfrac, xlower, xupper, ylower, yupper, nxpoints, nypoints)
    print(f"  [OK] mass_frac.tt3")
    print("=" * 60)


def make_surface():
    """Write surface_topo.tt3 (eta)."""
    print("[make_surface] Writing surface_topo.tt3")
    _write_tt3("surface_topo.tt3", eta, xlower, xupper, ylower, yupper, nxpoints, nypoints)
    print("  [OK] surface_topo.tt3")


def make_plots():
    """
    Verification plots for the DEM, hillshade, and lahar depth.

    The 'Source Selatan/Source Selatan.shp' overlay is shown as a dashed
    blue line in every plot — a visualization layer only, it does NOT
    affect the lahar source.
    """
    fw, fh = _domain_figsize()
    fw_cb  = _domain_figsize(extra_ratio=1.15)[0]

    # Prepare hillshade
    basal_topo  = topotools.Topography("basal_topo.tt3", 3)
    dem_arr     = basal_topo.Z.copy()
    dem_display = np.where(np.isnan(dem_arr), np.nanmin(dem_arr), dem_arr)

    from matplotlib.colors import LightSource
    ls  = LightSource(azdeg=315, altdeg=45)
    rgb = ls.shade(
        dem_display,
        cmap       = plt.cm.terrain,
        blend_mode = "overlay",
        vert_exag  = 2.0,
        dx         = _cellsize,
        dy         = _cellsize,
    )

    half = _cellsize / 2.0
    hs_extent = [xlower - half, xupper + half, ylower - half, yupper + half]

    # -- Plot 1: Source polygon + overlay (diagnostic) ----------------------
    fig, ax = plt.subplots(figsize=(fw, fh), dpi=_PLOT_DPI)
    x_poly, y_poly = _poly_source.exterior.xy
    ax.fill(x_poly, y_poly, color="brown", alpha=0.4)
    ax.plot(x_poly, y_poly, color="brown", linewidth=2,
            label="Source Selatan_1 (lahar source)")
    _plot_overlay(ax)
    ax.set_xlabel("Easting (m UTM)")
    ax.set_ylabel("Northing (m UTM)")
    ax.set_title(
        "Lahar Source - Source Selatan_1.shp (source, brown)\n"
        "Overlay - Source Selatan (dashed blue, visualization)\n"
        "DEMPredict.tt3"
    )
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True)
    savefig("source_area.png", dpi=_PLOT_DPI, bbox_inches="tight")
    print("Created source_area.png")
    plt.close()

    # -- Plot 2: Hillshade + source outline + overlay ------------------------
    fig, ax = plt.subplots(figsize=(fw, fh), dpi=_PLOT_DPI)
    ax.imshow(rgb, interpolation="nearest", origin="lower",
              extent=hs_extent, aspect="equal")
    x_poly, y_poly = _poly_source.exterior.xy
    ax.plot(x_poly, y_poly, color="brown", linewidth=1.5,
            linestyle="-", label="Source Selatan_1 (source)", zorder=8)
    _plot_overlay(ax)
    ax.set_xlim(xlower - half, xupper + half)
    ax.set_ylim(ylower - half, yupper + half)
    ax.set_title(
        f"Basal Topo - DEMPredict  {_cellsize:.2f} m\n"
        f"Hillshade + Source (brown) + Overlay Selatan (dashed blue)"
    )
    ax.set_xlabel("Easting (m UTM)")
    ax.set_ylabel("Northing (m UTM)")
    ax.legend(loc="upper right", fontsize=8)
    savefig("basal_topo_hillshade.png", dpi=_PLOT_DPI, bbox_inches="tight")
    print("Created basal_topo_hillshade.png")
    plt.close()

    # -- Plot 3: Surface topo (eta) + overlay --------------------------------
    eta_topo = topotools.Topography("surface_topo.tt3", 3)
    fig, ax  = plt.subplots(figsize=(fw_cb, fh), dpi=_PLOT_DPI)
    eta_arr  = eta_topo.Z.copy()
    im2 = ax.imshow(eta_arr, origin="lower",
                    extent=hs_extent,
                    cmap="terrain", interpolation="nearest", aspect="equal")
    plt.colorbar(im2, ax=ax, label="Eta elevation (m)")
    _plot_overlay(ax)
    ax.set_title("Surface topo (eta) - DEMPredict / Source Selatan_1\n"
                 "Overlay: Source Selatan (dashed blue)")
    ax.set_xlabel("Easting (m UTM)")
    ax.set_ylabel("Northing (m UTM)")
    ax.legend(loc="upper right", fontsize=8)
    savefig("surface_topo.png", dpi=_PLOT_DPI, bbox_inches="tight")
    print("Created surface_topo.png")
    plt.close()

    # -- Plot 4: Lahar depth + overlay ----------------------------------------
    h = eta_topo.Z - basal_topo.Z
    fig, ax = plt.subplots(figsize=(fw_cb, fh), dpi=_PLOT_DPI)
    cm = ax.pcolormesh(eta_topo.X, eta_topo.Y, h, cmap="Reds")
    ax.set_aspect("equal", adjustable="datalim")
    plt.colorbar(cm, ax=ax, label="Lahar depth (m)")
    _plot_overlay(ax)
    ax.set_title("Landslide / Lahar Depth - DEMPredict / Source Selatan_1\n"
                 "Overlay: Source Selatan (dashed blue)")
    ax.set_xlabel("Easting (m UTM)")
    ax.set_ylabel("Northing (m UTM)")
    ax.legend(loc="upper right", fontsize=8)
    savefig("landslide_depth.png", dpi=_PLOT_DPI, bbox_inches="tight")
    print("Created landslide_depth.png")
    plt.close()

    # -- Statistics ------------------------------------------------------------
    print("\n-- DEM Statistics (basal_topo.tt3) -------------")
    print(f"  Elevation min  : {np.nanmin(basal_topo.Z):.1f} m")
    print(f"  Elevation max  : {np.nanmax(basal_topo.Z):.1f} m")
    print(f"  Elevation mean : {np.nanmean(basal_topo.Z):.1f} m")
    print(f"  Std dev        : {np.nanstd(basal_topo.Z):.1f} m")
    grad_y, grad_x = np.gradient(
        np.where(np.isnan(basal_topo.Z), 0, basal_topo.Z),
        _cellsize, _cellsize)
    slope_mean = np.nanmean(np.sqrt(grad_x**2 + grad_y**2))
    print(f"  Mean slope     : {slope_mean:.4f} m/m")
    nan_pct2 = np.sum(np.isnan(basal_topo.Z)) / basal_topo.Z.size * 100
    print(f"  NaN pixels     : {nan_pct2:.2f}%")

    h_max  = np.nanmax(h)
    h_mean = np.nanmean(h[h > 0]) if np.any(h > 0) else 0.0
    area_src = _poly_source.area
    print(f"\n-- Lahar Source Statistics (Source Selatan_1) --------")
    print(f"  Polygon area   : {area_src:.2f} m^2  ({area_src / 1e6:.4f} km^2)")
    print(f"  h0 (input)     : {h0:.2f} m")
    print(f"  h max grid     : {h_max:.2f} m")
    print(f"  h mean (>0)    : {h_mean:.2f} m")
    if _overlay_ok and _polys_overlay:
        total_area = sum(p.area for p in _polys_overlay)
        print(f"\n-- Overlay Statistics (Source Selatan) -----------")
        print(f"  Feature count  : {len(_polys_overlay)}")
        print(f"  Total area     : {total_area:.2f} m^2  ({total_area / 1e6:.4f} km^2)")
    print("--------------------------------------------------\n")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    maketopo()
    make_surface()
    make_plots()
