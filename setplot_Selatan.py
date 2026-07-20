import os
import subprocess
import glob
import matplotlib
matplotlib.rcParams['animation.embed_limit'] = 200

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.lines
import matplotlib.patheffects
import numpy as np
from clawpack.visclaw import colormaps, geoplot, gridtools


# ===========================================================================
# TFINAL SYNCHRONIZATION — from setrun_Selatan.py
# ===========================================================================
try:
    from setrun_Selatan import _TFINAL as TFINAL
    print(f"[setplot] TFINAL = {TFINAL} s ({TFINAL/60:.1f} min)")
except ImportError:
    TFINAL = 1620.0
    print(f"[setplot] WARNING: setrun_Selatan.py not found, fallback TFINAL={TFINAL}")


NUM_OUTPUT_FRAMES = 162   # output interval 10s
MOVIE_FPS         = 6
MOVIE_SLOWDOWN    = 0.5


# ===========================================================================
# MP4 CONFIGURATION
# ===========================================================================
MP4_FIGNOS = [1, 2, 6]

MP4_FIGNAMES = {
    1: "lahar_utama",
    2: "mass_fraction",
    6: "depth",
}

MP4_CRF = 20


# ===========================================================================
# DISPLAY TOGGLES
# ===========================================================================
# FIX (point 2): Kali Adem is now SHOWN, including in the transect panel.
SHOW_KALI_ADEM = True


# ===========================================================================
# DEM METADATA — identical to maketopo_Selatan.py (DEMPredict.tt3)
# This domain is ALREADY consistent with the Selatan Source (verified against
# the reference bounds in maketopo: xlower=427610.433 -> xupper=444870.127,
# ylower=9153774.650 -> yupper=9171266.463), so it is NOT changed.
# ===========================================================================
_ncols     = 2083
_nrows     = 2111
_xllcorner = 427606.288
_yllcorner = 9153770.505
_cellsize  = 8.289958903244
_half      = _cellsize / 2.0

xlower = _xllcorner + _half
xupper = xlower + (_ncols - 1) * _cellsize
ylower = _yllcorner + _half
yupper = ylower + (_nrows - 1) * _cellsize

_XLIM = [_xllcorner, _xllcorner + _ncols * _cellsize]
_YLIM = [_yllcorner, _yllcorner + _nrows * _cellsize]

print(f"[setplot] Domain X: {xlower:.3f} -> {xupper:.3f}")
print(f"[setplot] Domain Y: {ylower:.3f} -> {yupper:.3f}")
print(f"[setplot] View   X: {_XLIM[0]:.3f} -> {_XLIM[1]:.3f}")
print(f"[setplot] View   Y: {_YLIM[0]:.3f} -> {_YLIM[1]:.3f}")

# ---------------------------------------------------------------------------
# Figure layout
# ---------------------------------------------------------------------------
_fig_height = 12.0
_fig_width  = 12.8

_AX_MAP = [0.07, 0.34, 0.80, 0.59]
_AX_TR  = [0.07, 0.06, 0.80, 0.22]

sea_level = 0.0
dem_zmin  = 350.0
dem_zmax  = 2900.0

bouss = False
if bouss:
    i_eta = 9; i_hm = 5; i_pb = 6
else:
    i_eta = 7; i_hm = 3; i_pb = 4

M_ENCER = 0.3
M_PEKAT = 0.6

# ---------------------------------------------------------------------------
# FIX (residual red band at final minutes of the transect)
# ---------------------------------------------------------------------------
H_MIN_WET   = 1e-3
H_MIN_PEKAT = 0.05

# ---------------------------------------------------------------------------
# SOURCE coordinates — Source Selatan
# ---------------------------------------------------------------------------
# Original bounding box from "Source Selatan_1.shp" (lahar source, used to
# compute m0/h0 in maketopo). Verified: the center of this bbox
# (439164.988, 9165667.463) matches EXACTLY the SOURCE_POINT already used
# in the sensitivity scripts (in sync with setinput/maketopo).
_src_xmin = 438885.574
_src_xmax = 439444.402
_src_ymin = 9165049.500
_src_ymax = 9166285.426

src_cx = (_src_xmin + _src_xmax) / 2.0
src_cy = (_src_ymin + _src_ymax) / 2.0
print(f"[setplot] src_cx, src_cy = {src_cx:.3f}, {src_cy:.3f}  "
      f"(from Source Selatan_1.shp bounds)")

KALI_ADEM_X     = 439540.56
KALI_ADEM_Y     = 9161173.13
KALI_ADEM_LABEL = "Kali Adem"

# ---------------------------------------------------------------------------
# Kali Gendol transect (upstream -> downstream)
# ---------------------------------------------------------------------------
# Original coordinates from "Kali Gendol.shp" (1 polyline feature, 40 vertices,
# CRS WGS_1984_UTM_Zone_49S / EPSG:32749 -- same as this project's CRS,
# so used directly without reprojection). The vertex order in the original
# file is already upstream (near the Selatan source area, higher/northern Y)
# -> downstream (lower/southern Y), so it is NOT reversed.
# Hardcoded (not read at runtime) following the same pattern as the
# source polygon overlay -- so it is always available without depending on
# the shapefile when make plots is run.
_TRANSECT_WAYPOINTS = np.array([
    [438969.765, 9166235.287],   # upstream (B)
    [439057.078, 9165989.224],
    [439152.328, 9165766.973],
    [439199.953, 9165663.786],
    [439210.748, 9165523.186],
    [439268.957, 9165380.310],
    [439358.915, 9165210.977],
    [439443.582, 9165031.060],
    [439528.249, 9164941.101],
    [439612.916, 9164787.643],
    [439660.541, 9164692.392],
    [439549.416, 9164226.725],
    [439403.365, 9163945.208],
    [439276.365, 9163665.807],
    [439092.215, 9163380.056],
    [439035.065, 9163246.706],
    [439047.765, 9163056.206],
    [438971.565, 9162776.805],
    [439066.815, 9162599.005],
    [439111.265, 9162287.854],
    [439092.215, 9162154.504],
    [439250.965, 9161608.403],
    [439466.866, 9161259.152],
    [439555.766, 9161144.852],
    [439517.666, 9160986.102],
    [439689.116, 9160744.801],
    [439651.016, 9160516.201],
    [439828.816, 9160249.500],
    [439771.666, 9159970.100],
    [440063.767, 9159538.299],
    [439981.217, 9159322.398],
    [440108.217, 9159081.098],
    [440146.317, 9158935.048],
    [440101.867, 9158623.897],
    [440247.917, 9158141.296],
    [440406.667, 9157957.146],
    [440381.267, 9157576.145],
    [440514.618, 9157303.094],
    [440470.168, 9157144.344],
    [440470.168, 9157131.644],   # downstream (B')
])


def _build_transect_samples(n_pts=500):
    wp = _TRANSECT_WAYPOINTS
    d  = np.zeros(len(wp))
    for k in range(1, len(wp)):
        d[k] = d[k-1] + np.hypot(wp[k, 0]-wp[k-1, 0], wp[k, 1]-wp[k-1, 1])
    du   = np.linspace(0, d[-1], n_pts)
    xout = np.interp(du, d, wp[:, 0])
    yout = np.interp(du, d, wp[:, 1])
    return xout, yout, du, d[-1]


_TR_X, _TR_Y, _TR_DIST, _TR_LEN = _build_transect_samples(500)
R_MAX = _TR_LEN


# ---------------------------------------------------------------------------
# Kali Adem projection onto the transect (for the blue dashed line in the
# B-B' panel)
# ---------------------------------------------------------------------------
def _project_point_to_transect(px, py):
    """
    Project point (px, py) onto the transect polyline (_TR_X, _TR_Y) and
    return the distance-along-transect (arc length) of the nearest point.
    Used to place the Kali Adem vertical line in the transect panel at a
    position consistent with its real location, rather than a straight-line
    distance from the source.
    """
    dx = _TR_X - px
    dy = _TR_Y - py
    dist2 = dx * dx + dy * dy
    idx = int(np.argmin(dist2))
    return _TR_DIST[idx], float(np.sqrt(dist2[idx]))


_KALI_ADEM_TR_DIST, _KALI_ADEM_TR_OFFSET = _project_point_to_transect(KALI_ADEM_X, KALI_ADEM_Y)
print(f"[setplot] Kali Adem -> projected onto transect: distance={_KALI_ADEM_TR_DIST:.1f} m "
      f"(perpendicular offset from the transect line: {_KALI_ADEM_TR_OFFSET:.1f} m)")

# ---------------------------------------------------------------------------
# Lahar colormap (yellow -> dark red)
# ---------------------------------------------------------------------------
_lc = [
    (0.000, (1.00, 0.97, 0.40, 0.22)), (0.010, (1.00, 0.90, 0.10, 0.70)),
    (0.030, (1.00, 0.85, 0.00, 0.78)), (0.070, (1.00, 0.70, 0.00, 0.83)),
    (0.130, (1.00, 0.55, 0.00, 0.87)), (0.200, (0.95, 0.35, 0.00, 0.91)),
    (0.300, (0.85, 0.15, 0.00, 0.94)), (0.450, (0.65, 0.00, 0.00, 0.96)),
    (0.650, (0.40, 0.00, 0.00, 0.98)), (1.000, (0.10, 0.00, 0.00, 1.00)),
]
_n   = 512
_lcd = np.zeros((_n, 4))
_fr  = [c[0] for c in _lc]
_rg  = [c[1] for c in _lc]
for i in range(_n):
    f = i / (_n - 1)
    j = np.clip(np.searchsorted(_fr, f, side="right") - 1, 0, len(_fr) - 2)
    t = (f - _fr[j]) / max(_fr[j+1] - _fr[j], 1e-12)
    _lcd[i] = np.array(_rg[j]) + t * (np.array(_rg[j+1]) - np.array(_rg[j]))
lahar_cmap = mcolors.ListedColormap(_lcd, name="lahar_depth")
depth_cmap = lahar_cmap

_KALI_ADEM_COLOR      = "none"
_KALI_ADEM_EDGE_COLOR = "#0055CC"
_KALI_ADEM_MARKER     = "o"
_KALI_ADEM_MARKERSIZE = 6
_KALI_ADEM_ZORDER     = 25

# ---------------------------------------------------------------------------
# Styling for the transect line on the map (letters B / B')
# ---------------------------------------------------------------------------
_TRANSECT_LINE_COLOR = "white"
_TRANSECT_LINE_WIDTH = 1.2
_TRANSECT_LABEL_START = "B"
_TRANSECT_LABEL_END   = "B'"


# ===========================================================================
# OVERLAY — Source Selatan polygon coordinates (blue dashed line)
# ===========================================================================
# ORIGINAL coordinates from "Source Selatan.shp" (143 vertices, 1 feature,
# CRS EPSG:32749 -- same as this project). This is a VISUALIZATION-ONLY
# overlay layer (not the lahar source -- the lahar source still uses
# "Source Selatan_1.shp" via the src_cx/src_cy bounds above), following the
# exact same source/overlay separation pattern used in maketopo_Selatan.py.
# Hardcoded so it is always shown without depending on the shapefile when
# make plots is run.
_OVERLAY_COORDS_X = np.array([
    439040.489, 439052.598, 439062.797, 439068.370, 439081.195, 439089.952,
    439100.218, 439104.853, 439120.424, 439123.100, 439135.369, 439156.391,
    439166.534, 439155.610, 439169.468, 439191.779, 439204.057, 439202.205,
    439192.262, 439179.001, 439169.794, 439172.055, 439175.809, 439173.479,
    439173.011, 439175.974, 439180.811, 439187.495, 439201.861, 439216.617,
    439241.691, 439258.694, 439279.718, 439300.507, 439323.090, 439334.404,
    439347.573, 439349.235, 439364.662, 439376.190, 439379.211, 439382.301,
    439388.062, 439397.775, 439408.024, 439442.998, 439438.092, 439452.143,
    439465.854, 439471.757, 439481.872, 439480.396, 439488.266, 439493.339,
    439487.266, 439501.516, 439498.788, 439493.174, 439464.441, 439453.913,
    439445.743, 439431.083, 439430.304, 439427.792, 439417.573, 439382.718,
    439352.763, 439346.692, 439334.608, 439322.955, 439314.783, 439296.220,
    439279.799, 439264.089, 439241.109, 439225.800, 439202.513, 439188.852,
    439167.050, 439143.320, 439129.038, 439114.771, 439096.149, 439080.440,
    439084.509, 439079.397, 439081.697, 439096.474, 439104.388, 439106.139,
    439109.278, 439102.259, 439096.272, 439084.470, 439068.778, 439051.757,
    439040.323, 439021.516, 439009.113, 438998.182, 438978.784, 438960.925,
    438945.116, 438933.999, 438923.466, 438907.702, 438897.550, 438893.501,
    438891.390, 438889.292, 438882.384, 438878.605, 438873.082, 438859.983,
    438852.928, 438840.841, 438828.684, 438821.372, 438812.280, 438811.765,
    438812.420, 438828.195, 438839.366, 438835.104, 438843.902, 438838.848,
    438833.253, 438825.285, 438826.966, 438821.901, 438833.433, 438858.173,
    438873.392, 438887.771, 438905.140, 438930.073, 438949.165, 438968.390,
    438983.283, 438996.948, 439009.372, 439026.608, 439040.489,
])
_OVERLAY_COORDS_Y = np.array([
    9166370.166, 9166359.206, 9166357.488, 9166370.931, 9166356.612, 9166338.498,
    9166310.959, 9166299.283, 9166268.372, 9166251.041, 9166227.324, 9166200.450,
    9166187.697, 9166172.292, 9166146.750, 9166111.542, 9166104.604, 9166079.442,
    9166059.097, 9166030.313, 9166000.800, 9165979.399, 9165960.794, 9165953.203,
    9165936.529, 9165914.774, 9165898.286, 9165875.801, 9165856.093, 9165840.762,
    9165832.148, 9165835.807, 9165838.771, 9165839.802, 9165840.398, 9165833.787,
    9165819.616, 9165798.171, 9165786.814, 9165774.299, 9165751.167, 9165726.375,
    9165701.198, 9165674.460, 9165654.876, 9165620.735, 9165595.296, 9165550.196,
    9165500.193, 9165465.507, 9165425.879, 9165399.323, 9165354.301, 9165323.100,
    9165292.747, 9165261.589, 9165240.283, 9165200.632, 9165177.716, 9165148.516,
    9165119.021, 9165086.234, 9165067.382, 9165044.069, 9164998.036, 9164989.791,
    9164996.949, 9165014.607, 9165033.252, 9165062.694, 9165098.678, 9165139.406,
    9165170.767, 9165221.250, 9165232.501, 9165253.916, 9165280.849, 9165302.232,
    9165324.401, 9165356.143, 9165383.392, 9165416.075, 9165438.489, 9165463.565,
    9165487.491, 9165529.056, 9165564.107, 9165599.738, 9165623.501, 9165660.815,
    9165695.711, 9165716.350, 9165729.475, 9165761.882, 9165779.246, 9165786.502,
    9165792.346, 9165806.715, 9165818.244, 9165831.835, 9165844.582, 9165857.038,
    9165870.734, 9165878.811, 9165888.599, 9165905.305, 9165915.537, 9165925.173,
    9165936.939, 9165960.012, 9165982.420, 9165983.244, 9165977.938, 9165977.137,
    9165983.447, 9165992.597, 9165998.532, 9166007.275, 9166019.264, 9166031.332,
    9166054.392, 9166069.887, 9166090.504, 9166119.670, 9166146.159, 9166176.648,
    9166189.568, 9166204.272, 9166218.757, 9166269.027, 9166273.601, 9166280.678,
    9166285.405, 9166299.469, 9166309.243, 9166316.427, 9166325.062, 9166333.445,
    9166354.040, 9166368.249, 9166379.254, 9166375.439, 9166370.166,
])

print(f"[setplot] Source Selatan overlay: {len(_OVERLAY_COORDS_X)} vertices "
      f"(from Source Selatan.shp)")


def _plot_overlay(ax, linewidth=1.5, label="Source Selatan"):
    """Draw the source polygon outline as a blue dashed line."""
    ax.plot(
        _OVERLAY_COORDS_X, _OVERLAY_COORDS_Y,
        color     = "blue",
        linewidth = linewidth,
        linestyle = "--",
        label     = label,
        zorder    = 30,
    )


def _plot_transect_line(ax):
    """Draw the Kali Gendol transect line (upstream B -> downstream B') on the map."""
    ax.plot(
        _TR_X, _TR_Y,
        color     = _TRANSECT_LINE_COLOR,
        linewidth = _TRANSECT_LINE_WIDTH,
        linestyle = "-",
        zorder    = 29,
        solid_capstyle = "round",
    )
    ax.annotate(
        _TRANSECT_LABEL_START,
        xy=(_TR_X[0], _TR_Y[0]),
        xytext=(4, 4), textcoords="offset points",
        fontsize=9, fontweight="bold", color="white",
        path_effects=[matplotlib.patheffects.withStroke(linewidth=2.0, foreground="black")],
        zorder=31,
    )
    ax.annotate(
        _TRANSECT_LABEL_END,
        xy=(_TR_X[-1], _TR_Y[-1]),
        xytext=(4, 4), textcoords="offset points",
        fontsize=9, fontweight="bold", color="white",
        path_effects=[matplotlib.patheffects.withStroke(linewidth=2.0, foreground="black")],
        zorder=31,
    )

    # FIX (point 2): also show the blue marker + dashed line for Kali Adem
    # ON THE MAP at its real position (not projected), so its actual
    # location remains visible relative to the map.
    if SHOW_KALI_ADEM:
        ax.plot(KALI_ADEM_X, KALI_ADEM_Y,
                marker=_KALI_ADEM_MARKER, color=_KALI_ADEM_COLOR,
                markeredgecolor=_KALI_ADEM_EDGE_COLOR, markeredgewidth=1.5,
                markersize=_KALI_ADEM_MARKERSIZE, linestyle="None",
                zorder=_KALI_ADEM_ZORDER)


# ===========================================================================
# Load DEM & build hillshade (once, cached)
# ===========================================================================
_DEM_PATH = "DEMPredict.tt3" if os.path.isfile("DEMPredict.tt3") else "basal_topo.tt3"
print(f"[setplot] Hillshade DEM: {_DEM_PATH}")


def _read_tt3(path):
    try:
        with open(path) as fh:
            nc  = int(fh.readline().split()[0])
            nr  = int(fh.readline().split()[0])
            xll = float(fh.readline().split()[0])
            yll = float(fh.readline().split()[0])
            cs  = float(fh.readline().split()[0])
            nd  = float(fh.readline().split()[0])
            Z   = np.loadtxt(fh).astype(float)
        Z[Z == nd] = np.nan
        Z   = Z[::-1, :]
        h2  = cs / 2.0
        x1d = (xll + h2) + np.arange(nc) * cs
        y1d = (yll + h2) + np.arange(nr) * cs
        print(f"[setplot] DEM {nc}x{nr}  cs={cs:.6f} m")
        print(f"[setplot] DEM X: {x1d[0]:.3f} -> {x1d[-1]:.3f}")
        print(f"[setplot] DEM Y: {y1d[0]:.3f} -> {y1d[-1]:.3f}")
        print(f"[setplot] DEM Z: [{np.nanmin(Z):.0f}, {np.nanmax(Z):.0f}] m")
        return x1d, y1d, Z, cs
    except FileNotFoundError:
        print(f"[setplot] WARNING: {path} not found -- hillshade disabled")
        return None, None, None, None
    except Exception as e:
        print(f"[setplot] WARNING DEM read: {e}")
        return None, None, None, None


_dem_x, _dem_y, _dem_Z, _dem_cs = _read_tt3(_DEM_PATH)
_hillshade_rgba   = None
_hillshade_extent = None


def _build_hillshade():
    global _hillshade_rgba, _hillshade_extent
    if _dem_Z is None:
        return None, None
    if _hillshade_rgba is not None:
        return _hillshade_rgba, _hillshade_extent
    cs = _dem_cs
    Z  = _dem_Z.copy()
    nm = np.isnan(Z)
    if nm.any():
        from scipy.ndimage import distance_transform_edt
        _, idx = distance_transform_edt(nm, return_indices=True)
        Z[nm]  = Z[idx[0][nm], idx[1][nm]]
    dzdx = np.gradient(Z, cs, axis=1)
    dzdy = np.gradient(Z, cs, axis=0)
    az, alt = np.radians(360 - 315), np.radians(45)
    lx = np.cos(alt) * np.cos(az)
    ly = np.cos(alt) * np.sin(az)
    lz = np.sin(alt)
    nx = -dzdx; ny = -dzdy; nz = np.ones_like(Z)
    mag = np.sqrt(nx**2 + ny**2 + nz**2)
    nx /= mag; ny /= mag; nz /= mag
    hs = np.clip(nx*lx + ny*ly + nz*lz, 0, 1)
    hs = 0.30 + 0.65 * np.power(hs, 0.85)
    rgba = np.zeros((*hs.shape, 4), dtype=np.float32)
    rgba[..., :3] = hs[..., None]
    rgba[...,  3] = 1.0
    _hillshade_rgba = rgba
    h2 = cs / 2.0
    _hillshade_extent = [
        _dem_x[0]  - h2,
        _dem_x[-1] + h2,
        _dem_y[0]  - h2,
        _dem_y[-1] + h2,
    ]
    print(f"[setplot] Hillshade OK  shape={rgba.shape}")
    return _hillshade_rgba, _hillshade_extent


_build_hillshade()


def _draw_hillshade(ax):
    rgba, extent = _build_hillshade()
    if rgba is None:
        return
    for im in list(ax.images):
        if getattr(im, '_is_hs', False):
            im.remove()
    xlim_save = ax.get_xlim()
    ylim_save = ax.get_ylim()
    im = ax.imshow(rgba, extent=extent, origin="lower",
                   aspect="auto", zorder=0, interpolation="bilinear")
    im._is_hs = True
    ax.set_xlim(xlim_save)
    ax.set_ylim(ylim_save)


def _interp_topo_along_transect(xout, yout):
    if _dem_Z is None:
        return np.zeros_like(xout)
    from scipy.interpolate import RegularGridInterpolator
    f = RegularGridInterpolator((_dem_y, _dem_x), _dem_Z,
                                method="linear", bounds_error=False,
                                fill_value=np.nan)
    return f(np.column_stack([yout, xout]))


def _extract_all_patches(framesoln, var_indices):
    xs, ys = [], []
    vo = {v: [] for v in var_indices}
    for state in framesoln.states:
        g = state.patch
        XX, YY = np.meshgrid(g.x.centers, g.y.centers)
        xs.append(XX.ravel()); ys.append(YY.ravel())
        for v in var_indices:
            vo[v].append(state.q[v].ravel())
    res = {"x": np.concatenate(xs), "y": np.concatenate(ys)}
    for v in var_indices:
        res[v] = np.concatenate(vo[v])
    return res


def _draw_overlays(ax):
    """Draw the source polygon overlay + B-B' transect line (+ Kali Adem) on the map."""
    if SHOW_KALI_ADEM:
        ax.plot(KALI_ADEM_X, KALI_ADEM_Y,
                marker=_KALI_ADEM_MARKER, color=_KALI_ADEM_COLOR,
                markeredgecolor=_KALI_ADEM_EDGE_COLOR, markeredgewidth=1.5,
                markersize=_KALI_ADEM_MARKERSIZE, linestyle="None",
                zorder=_KALI_ADEM_ZORDER)

    _plot_overlay(ax)
    _plot_transect_line(ax)

    _transect_legend_line = matplotlib.lines.Line2D(
        [0], [0], color=_TRANSECT_LINE_COLOR, linewidth=_TRANSECT_LINE_WIDTH,
        linestyle="-"
    )
    _transect_legend_line.set_path_effects([
        matplotlib.patheffects.withStroke(linewidth=_TRANSECT_LINE_WIDTH + 2.0,
                                           foreground="black"),
    ])

    legend_handles = []
    legend_labels  = []

    if SHOW_KALI_ADEM:
        legend_handles.append(matplotlib.lines.Line2D(
            [0], [0], marker=_KALI_ADEM_MARKER, color="w",
            markerfacecolor=_KALI_ADEM_COLOR,
            markeredgecolor=_KALI_ADEM_EDGE_COLOR,
            markeredgewidth=1.5, markersize=_KALI_ADEM_MARKERSIZE,
            linestyle="None"
        ))
        legend_labels.append(KALI_ADEM_LABEL)

    legend_handles.append(matplotlib.lines.Line2D(
        [0], [0], color="blue", linewidth=1.5, linestyle="--"
    ))
    legend_labels.append("Source Selatan")

    legend_handles.append(_transect_legend_line)
    legend_labels.append("Transect B-B' (Kali Gendol)")

    ax.legend(legend_handles, legend_labels, fontsize=8, loc="lower right",
              framealpha=0.85, facecolor="white", edgecolor="gray")


# ===========================================================================
# MP4 EXPORT — stitch visclaw PNG frames with ffmpeg
# ===========================================================================

def _ffmpeg_available():
    try:
        subprocess.run(["ffmpeg", "-version"], capture_output=True, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def buat_mp4(plotdir):
    """
    Stitch visclaw PNG frames into an MP4 file per figure.
    Output MP4:
        <plotdir>/lahar_utama.mp4   (map + B-B' cross-section, Kali Gendol)
        <plotdir>/mass_fraction.mp4
        <plotdir>/depth.mp4
    """
    if not _ffmpeg_available():
        print("[mp4] ERROR: ffmpeg not found.")
        print("      Install with: sudo apt install ffmpeg")
        return

    if not isinstance(plotdir, str):
        plotdir = str(plotdir)

    plotdir = os.path.abspath(plotdir)
    os.makedirs(plotdir, exist_ok=True)

    if not os.path.isdir(plotdir):
        print(f"[mp4] ERROR: directory not found: {plotdir}")
        return

    print(f"\n[mp4] -- Exporting MP4 from PNG frames in {plotdir} --")

    for figno in MP4_FIGNOS:
        figname  = MP4_FIGNAMES.get(figno, f"fig{figno}")
        out_mp4  = os.path.join(plotdir, f"{figname}.mp4")
        listfile = os.path.join(plotdir, f"_framelist_fig{figno}.txt")

        pattern = os.path.join(plotdir, f"frame*fig{figno}.png")
        frames  = sorted(glob.glob(pattern))

        if not frames:
            print(f"[mp4] fig{figno}: no PNG frames found ({pattern})")
            continue

        print(f"[mp4] fig{figno} ({figname}): {len(frames)} frames -> {out_mp4}")

        if len(frames) < NUM_OUTPUT_FRAMES:
            print(f"  [WARNING] Only {len(frames)}/{NUM_OUTPUT_FRAMES} PNG frames "
                  f"found for fig{figno}! The resulting video may be "
                  f"TRUNCATED (not reaching t={TFINAL:.0f}s / "
                  f"{TFINAL/60:.1f} min). Make sure the entire plotting "
                  f"process (make plots / --parallel) has fully finished "
                  f"before running the MP4 export, then rerun "
                  f"buat_mp4() if needed.")

        try:
            with open(listfile, "w") as lf:
                for fp in frames:
                    lf.write(f"file '{os.path.abspath(fp)}'\n")
                    lf.write(f"duration {1.0 / MOVIE_FPS:.6f}\n")
                lf.write(f"file '{os.path.abspath(frames[-1])}'\n")

            cmd = [
                "ffmpeg", "-y",
                "-f", "concat",
                "-safe", "0",
                "-i", listfile,
                "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                "-c:v", "libx264",
                "-crf",    str(MP4_CRF),
                "-preset", "slow",
                "-pix_fmt", "yuv420p",
                out_mp4,
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode == 0 and os.path.isfile(out_mp4) and os.path.getsize(out_mp4) > 0:
                size_mb = os.path.getsize(out_mp4) / 1024 / 1024
                print(f"  OK  {out_mp4}  ({size_mb:.1f} MB)")
            else:
                print(f"  FAILED  ffmpeg error (fig{figno}):")
                for ln in result.stderr.strip().splitlines()[-6:]:
                    print(f"    {ln}")

        except Exception as e:
            print(f"  FAILED  Exception fig{figno}: {e}")

        finally:
            if os.path.isfile(listfile):
                try:
                    os.remove(listfile)
                except OSError:
                    pass

    print("[mp4] -- Done --\n")


# ===========================================================================
# setplot — main function called by visclaw
# ===========================================================================
def setplot(plotdata=None):
    import clawpack.dclaw.plot as dplot
    from clawpack.pyclaw import Solution

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()

    plotdata.clearfigures()
    plotdata.format = "ascii"

    def timeformat(t):
        from numpy import mod
        return f"{int(t/3600)}:{int(mod(t,3600)/60):02d}:{int(mod(t,60)):02d}"

    def ba_hillshade(current_data):
        pass

    def aa_map(current_data):
        from pylab import gca, ticklabel_format, xticks
        ax = gca()
        _draw_hillshade(ax)
        ax.set_xlim(_XLIM[0], _XLIM[1])
        ax.set_ylim(_YLIM[0], _YLIM[1])
        ax.set_title(f"t = {timeformat(current_data.t)}", fontsize=9, pad=3)
        ticklabel_format(useOffset=False)
        xticks(rotation=30, fontsize=7)
        ax.tick_params(axis='x', labelsize=7)
        ax.tick_params(axis='y', labelsize=7)
        ax.set_xlabel("Easting UTM (m)", fontsize=8, labelpad=2)
        ax.set_ylabel("Northing UTM (m)", fontsize=8)
        _draw_overlays(ax)

    # -----------------------------------------------------------------------
    # plot_var functions
    # -----------------------------------------------------------------------
    def land_nan(current_data):
        return np.full(current_data.q[0].shape, np.nan)

    def pure_water(current_data):
        q   = current_data.q
        h   = q[0]; hm = q[i_hm]
        eta = dplot.eta(current_data)
        with np.errstate(divide="ignore", invalid="ignore"):
            m = hm / h
        return np.where((h > 1e-3) & (m < M_ENCER), eta, np.nan)

    def lahar_depth(current_data):
        q  = current_data.q; h = q[0]; hm = q[i_hm]
        with np.errstate(divide="ignore", invalid="ignore"):
            m = hm / h
        return np.where((h > 1e-3) & (m >= M_ENCER), h, np.nan)

    def depth_wet(current_data):
        h = current_data.q[0]
        return np.where(h > 1e-3, h, np.nan)

    def mass_frac(current_data):
        q  = current_data.q; h = q[0]; hm = q[i_hm]
        with np.errstate(divide="ignore", invalid="ignore"):
            m = hm / h
        return np.where(h > 0.01, m, np.nan)

    # -----------------------------------------------------------------------
    # Figure 1 — Main lahar map + Kali Gendol transect
    # -----------------------------------------------------------------------
    plotfigure = plotdata.new_plotfigure(name="Banjir Lahar Dingin Merapi", figno=1)
    plotfigure.figsize   = (_fig_width, _fig_height)
    plotfigure.facecolor = "w"

    plotaxes = plotfigure.new_plotaxes("pcolor")
    plotaxes.axescmd    = (f"gcf().add_axes([{_AX_MAP[0]},{_AX_MAP[1]},"
                           f"{_AX_MAP[2]},{_AX_MAP[3]}])")
    plotaxes.title      = ""
    plotaxes.scaled     = False
    plotaxes.xlimits    = _XLIM
    plotaxes.ylimits    = _YLIM
    plotaxes.beforeaxes = ba_hillshade
    plotaxes.afteraxes  = aa_map

    plotitem = plotaxes.new_plotitem(plot_type="2d_pcolor")
    plotitem.plot_var    = land_nan
    plotitem.pcolor_cmap = plt.cm.gray
    plotitem.pcolor_cmin = dem_zmin; plotitem.pcolor_cmax = dem_zmax
    plotitem.add_colorbar       = False
    plotitem.amr_celledges_show = [0]
    plotitem.patchedges_show    = 0

    plotitem = plotaxes.new_plotitem(plot_type="2d_pcolor")
    plotitem.plot_var    = pure_water
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = sea_level - 5.0; plotitem.pcolor_cmax = sea_level + 20.0
    plotitem.add_colorbar       = False
    plotitem.amr_celledges_show = [0, 0, 0]
    plotitem.patchedges_show    = 0

    plotitem = plotaxes.new_plotitem(plot_type="2d_pcolor")
    plotitem.plot_var    = lahar_depth
    plotitem.pcolor_cmap = lahar_cmap
    plotitem.pcolor_cmin = 0.0; plotitem.pcolor_cmax = 15.0
    plotitem.add_colorbar       = True
    plotitem.colorbar_label     = "Lahar depth (m)"
    plotitem.amr_celledges_show = [0, 0, 0]
    plotitem.patchedges_show    = 0

    # Bottom panel: Kali Gendol transect
    plotaxes_tr = plotfigure.new_plotaxes("transect")
    plotaxes_tr.axescmd = (f"gcf().add_axes([{_AX_TR[0]},{_AX_TR[1]},"
                           f"{_AX_TR[2]},{_AX_TR[3]}])")
    plotaxes_tr.title   = ""

    def plot_transect(current_data):
        from clawpack.pyclaw import Solution
        from pylab import gca, grid, legend, nan
        import numpy as np

        pd  = current_data.plotdata
        sol = Solution(current_data.frameno, path=pd.outdir, file_format=pd.format)
        xo, yo, dist = _TR_X, _TR_Y, _TR_DIST

        eta_val = gridtools.grid_output_2d(sol, -1,   xo, yo)
        h       = gridtools.grid_output_2d(sol,  0,   xo, yo)
        hm      = gridtools.grid_output_2d(sol, i_hm, xo, yo)
        tp      = _interp_topo_along_transect(xo, yo)

        with np.errstate(divide="ignore", invalid="ignore"):
            m = np.where(h > H_MIN_PEKAT, hm / h, nan)
        wet = h > H_MIN_WET

        topo_min = float(np.nanmin(tp))
        topo_max = float(np.nanmax(tp))

        B = np.where(wet, eta_val - h, tp)

        pekat_mask_raw = (h > H_MIN_PEKAT) & (m >= M_PEKAT)

        ax = gca()
        ax.fill_between(dist, dem_zmin, tp, color=[.55, .55, .55], zorder=1)

        ax.fill_between(dist, B, np.where(wet & (m < M_ENCER), eta_val, nan),
                         color=[.30, .60, 1.], alpha=.85,
                         label=f"Mud flow (m<{M_ENCER})", zorder=3)
        ax.fill_between(dist, B,
                         np.where(wet & (m >= M_ENCER) & (m < M_PEKAT), eta_val, nan),
                         color=[1., .55, 0.], alpha=.85,
                         label=f"Hyperconcentrated ({M_ENCER}<=m<{M_PEKAT})", zorder=4)
        ax.fill_between(dist, B, np.where(pekat_mask_raw, eta_val, nan),
                         color=[.65, 0., 0.], alpha=.90,
                         label=f"Lahar (m>={M_PEKAT})", zorder=5)

        ax.plot(dist, np.where(wet, eta_val, nan), "k-", lw=.8, zorder=6)

        # FIX (point 2): Kali Adem now appears IN THE TRANSECT PANEL as a
        # blue dashed vertical line at its projected position along the
        # B-B' transect (computed once in _project_point_to_transect),
        # plus a round marker at the intersection with the topography profile.
        if SHOW_KALI_ADEM and 0.0 <= _KALI_ADEM_TR_DIST <= _TR_LEN:
            ax.axvline(x=_KALI_ADEM_TR_DIST, color=_KALI_ADEM_EDGE_COLOR,
                       linewidth=1.4, linestyle="--", zorder=7, alpha=.9,
                       label=f"{KALI_ADEM_LABEL} ({_KALI_ADEM_TR_DIST:.0f} m)")
            _tp_at_ka = _interp_topo_along_transect(
                np.array([KALI_ADEM_X]), np.array([KALI_ADEM_Y]))[0]
            ax.plot(_KALI_ADEM_TR_DIST, _tp_at_ka,
                    marker=_KALI_ADEM_MARKER, color=_KALI_ADEM_COLOR,
                    markeredgecolor=_KALI_ADEM_EDGE_COLOR, markeredgewidth=1.5,
                    markersize=_KALI_ADEM_MARKERSIZE, linestyle="None", zorder=8)

        ax.annotate(_TRANSECT_LABEL_START, xy=(0.0, 1.0), xycoords="axes fraction",
                    xytext=(2, 2), textcoords="offset points",
                    fontsize=10, fontweight="bold", ha="left", va="bottom")
        ax.annotate(_TRANSECT_LABEL_END, xy=(1.0, 1.0), xycoords="axes fraction",
                    xytext=(-2, 2), textcoords="offset points",
                    fontsize=10, fontweight="bold", ha="right", va="bottom")

        ax.set_title(
            f"Kali Gendol B-B' transect  t={current_data.t:.0f}s  "
            f"hmax={h.max():.2f} m  mmean={np.nanmean(m):.2f}",
            fontsize=12, fontweight="bold", pad=6)
        ax.set_xlabel("Distance from upstream (m)", fontsize=8, labelpad=2)
        ax.set_ylabel("Elevation (m)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_xlim([0, _TR_LEN])
        ax.set_ylim(topo_min - 50, topo_max + 120)
        ax.grid(True)
        ax.legend(fontsize=7, loc="upper right")

    plotaxes_tr.afteraxes = plot_transect
    plotitem = plotaxes_tr.new_plotitem(plot_type="2d_pcolor")
    plotitem.plot_var    = land_nan
    plotitem.pcolor_cmap = plt.cm.gray
    plotitem.pcolor_cmin = 0.0; plotitem.pcolor_cmax = 1.0
    plotitem.add_colorbar       = False
    plotitem.amr_celledges_show = [0]
    plotitem.patchedges_show    = 0

    # -----------------------------------------------------------------------
    # Figure 2 — Sediment mass fraction
    # -----------------------------------------------------------------------
    plotfigure = plotdata.new_plotfigure(name="Mass Fraction", figno=2)
    plotfigure.kwargs    = {"figsize": (_fig_width, _fig_height)}
    plotfigure.facecolor = "w"

    plotaxes = plotfigure.new_plotaxes("pcolor")
    plotaxes.title      = "Sediment mass fraction (m = hm/h)"
    plotaxes.scaled     = False
    plotaxes.xlimits    = _XLIM; plotaxes.ylimits = _YLIM
    plotaxes.beforeaxes = ba_hillshade; plotaxes.afteraxes = aa_map

    plotitem = plotaxes.new_plotitem(plot_type="2d_pcolor")
    plotitem.plot_var    = land_nan; plotitem.pcolor_cmap = plt.cm.gray
    plotitem.pcolor_cmin = dem_zmin; plotitem.pcolor_cmax = dem_zmax
    plotitem.add_colorbar       = False
    plotitem.amr_celledges_show = [0]; plotitem.patchedges_show = 0

    plotitem = plotaxes.new_plotitem(plot_type="2d_pcolor")
    plotitem.plot_var    = mass_frac
    plotitem.pcolor_cmap = colormaps.blue_yellow_red
    plotitem.pcolor_cmin = 0.0; plotitem.pcolor_cmax = 0.65
    plotitem.add_colorbar       = True
    plotitem.colorbar_label     = "Sediment mass fraction m (-)"
    plotitem.amr_celledges_show = [0, 0, 0]; plotitem.patchedges_show = 0

    # -----------------------------------------------------------------------
    # Figure 6 — Total flow depth
    # -----------------------------------------------------------------------
    plotfigure = plotdata.new_plotfigure(name="Depth", figno=6)
    plotfigure.kwargs    = {"figsize": (_fig_width, _fig_height)}
    plotfigure.facecolor = "w"

    plotaxes = plotfigure.new_plotaxes("pcolor")
    plotaxes.title      = "Lahar flow depth h (m)"
    plotaxes.scaled     = False
    plotaxes.xlimits    = _XLIM; plotaxes.ylimits = _YLIM
    plotaxes.beforeaxes = ba_hillshade; plotaxes.afteraxes = aa_map

    plotitem = plotaxes.new_plotitem(plot_type="2d_pcolor")
    plotitem.plot_var    = land_nan; plotitem.pcolor_cmap = plt.cm.gray
    plotitem.pcolor_cmin = dem_zmin; plotitem.pcolor_cmax = dem_zmax
    plotitem.add_colorbar       = False
    plotitem.amr_celledges_show = [0]; plotitem.patchedges_show = 0

    plotitem = plotaxes.new_plotitem(plot_type="2d_pcolor")
    plotitem.plot_var    = depth_wet
    plotitem.pcolor_cmap = depth_cmap
    plotitem.pcolor_cmin = 0.0; plotitem.pcolor_cmax = 15.0
    plotitem.add_colorbar       = True
    plotitem.colorbar_label     = "Depth h (m)"
    plotitem.amr_celledges_show = [0, 0, 0]; plotitem.patchedges_show = 0

    # -----------------------------------------------------------------------
    # Figure 3 — Depth vs distance-from-source scatter
    # -----------------------------------------------------------------------
    plotfigure = plotdata.new_plotfigure(name="scatter_depth", figno=3)
    plotfigure.kwargs = {"figsize": (10, 5)}

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title    = "Lahar depth vs distance"
    plotaxes.xlimits  = [0, R_MAX]
    plotaxes.ylimits  = "auto"
    plotaxes.grid     = True

    def scatter_depth(current_data):
        from pylab import gca, grid, title, xlabel, ylabel
        ax = gca()
        pd  = current_data.plotdata
        sol = Solution(current_data.frameno, path=pd.outdir, file_format=pd.format)
        d   = _extract_all_patches(sol, [0])
        x, y, h = d["x"], d["y"], d[0]
        wet = h > 1e-3
        if wet.sum() == 0:
            ax.text(0.5, 0.5, "No flow", ha="center", va="center",
                    transform=ax.transAxes, color="gray")
        else:
            r = np.sqrt((x[wet] - src_cx)**2 + (y[wet] - src_cy)**2)
            ax.plot(r, h[wet], ".", color=[.65, 0, 0], ms=1, alpha=.6)
        if SHOW_KALI_ADEM:
            rka = np.sqrt((KALI_ADEM_X - src_cx)**2 + (KALI_ADEM_Y - src_cy)**2)
            ax.axvline(x=rka, color=_KALI_ADEM_EDGE_COLOR, lw=1.2, ls=":",
                       zorder=5, alpha=.85, label=f"{KALI_ADEM_LABEL} ({rka:.0f} m)")
            ax.legend(fontsize=7)
        ax.set_xlim(0, R_MAX)
        xlabel("Distance from source (m)"); ylabel("h (m)")
        title(f"Depth vs distance  t={current_data.t:.0f}s"); grid(True)

    plotaxes.afteraxes = scatter_depth
    plotitem = plotaxes.new_plotitem(plot_type="2d_pcolor")
    plotitem.plot_var    = land_nan; plotitem.pcolor_cmap = plt.cm.gray
    plotitem.pcolor_cmin = 0.0; plotitem.pcolor_cmax = 1.0
    plotitem.add_colorbar       = False
    plotitem.amr_celledges_show = [0]; plotitem.patchedges_show = 0

    # -----------------------------------------------------------------------
    # Figure 9 — Speed vs distance-from-source scatter
    # -----------------------------------------------------------------------
    plotfigure = plotdata.new_plotfigure(name="scatter_speed", figno=9)
    plotfigure.kwargs = {"figsize": (10, 5)}

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title   = "Lahar speed vs distance"
    plotaxes.xlimits = [0, R_MAX]
    plotaxes.ylimits = "auto"
    plotaxes.grid    = True

    def scatter_speed(current_data):
        from pylab import gca, grid, title, xlabel, ylabel
        ax = gca()
        pd  = current_data.plotdata
        sol = Solution(current_data.frameno, path=pd.outdir, file_format=pd.format)
        d   = _extract_all_patches(sol, [0, 1, 2])
        x, y, h, hu, hv = d["x"], d["y"], d[0], d[1], d[2]
        wet = h > 1e-2
        if wet.sum() == 0:
            ax.text(0.5, 0.5, "No flow", ha="center", va="center",
                    transform=ax.transAxes, color="gray")
        else:
            r   = np.sqrt((x[wet] - src_cx)**2 + (y[wet] - src_cy)**2)
            spd = np.sqrt(hu[wet]**2 + hv[wet]**2) / h[wet]
            ax.plot(r, spd, "+", color="r", ms=2, alpha=.6)
        if SHOW_KALI_ADEM:
            rka = np.sqrt((KALI_ADEM_X - src_cx)**2 + (KALI_ADEM_Y - src_cy)**2)
            ax.axvline(x=rka, color=_KALI_ADEM_EDGE_COLOR, lw=1.2, ls=":",
                       zorder=5, alpha=.85, label=f"{KALI_ADEM_LABEL} ({rka:.0f} m)")
            ax.legend(fontsize=7)
        ax.set_xlim(0, R_MAX)
        xlabel("Distance from source (m)"); ylabel("Speed (m/s)")
        title(f"Speed vs distance  t={current_data.t:.0f}s"); grid(True)

    plotaxes.afteraxes = scatter_speed
    plotitem = plotaxes.new_plotitem(plot_type="2d_pcolor")
    plotitem.plot_var    = land_nan; plotitem.pcolor_cmap = plt.cm.gray
    plotitem.pcolor_cmin = 0.0; plotitem.pcolor_cmax = 1.0
    plotitem.add_colorbar       = False
    plotitem.amr_celledges_show = [0]; plotitem.patchedges_show = 0

    # -----------------------------------------------------------------------
    # Timing plots
    # -----------------------------------------------------------------------
    def make_timing_plots(plotdata):
        from clawpack.visclaw import plot_timing_stats
        try:
            tdir = plotdata.plotdir + "/_timing_figures"
            os.system(f"mkdir -p {tdir}")
            plot_timing_stats.make_plots(
                outdir=plotdata.outdir, make_pngs=True, plotdir=tdir,
                units={"comptime": "minutes", "simtime": "minutes",
                       "cell": "millions"})
            os.system(f"cp {plotdata.outdir}/timing.* {tdir}")
        except Exception:
            print("*** Error making timing plots")

    otherfigure = plotdata.new_otherfigure(name="timing",
                                           fname="_timing_figures/timing.html")
    otherfigure.makefig = make_timing_plots

    # -----------------------------------------------------------------------
    # MP4 export
    # -----------------------------------------------------------------------
    def _make_mp4_otherfigure(plotdata):
        buat_mp4(plotdata.plotdir)

    otherfigure_mp4 = plotdata.new_otherfigure(
        name="mp4_export",
        fname="lahar_utama.mp4",
    )
    otherfigure_mp4.makefig = _make_mp4_otherfigure

    # -----------------------------------------------------------------------
    # Output settings
    # -----------------------------------------------------------------------
    plotdata.printfigs           = True
    plotdata.print_format        = "png"
    plotdata.print_framenos      = "all"
    plotdata.print_gaugenos      = "all"
    plotdata.print_fignos        = "all"
    plotdata.html                = True
    plotdata.html_homelink       = "../README.html"
    plotdata.latex               = True
    plotdata.latex_figsperline   = 2
    plotdata.latex_framesperline = 1
    plotdata.latex_makepdf       = False
    plotdata.parallel            = True

    return plotdata


# ===========================================================================
# Manual entry point
# ===========================================================================
if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        _plotdir = sys.argv[1]
    else:
        _candidates = ["_plots", "_output/_plots", "../_plots"]
        _plotdir = None
        for _c in _candidates:
            if os.path.isdir(_c):
                _plotdir = _c
                break
        if _plotdir is None:
            print("Usage: python setplot.py [plotdir]")
            print("  plotdir  -- directory containing visclaw PNG frames (frame*fig*.png)")
            print("  Example  -- python setplot.py _plots")
            sys.exit(1)

    print(f"[mp4] Manual export from: {_plotdir}")
    buat_mp4(_plotdir)
