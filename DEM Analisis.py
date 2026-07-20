import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ============================================================
# CONFIGURATION
# ============================================================
OUTPUT_DIR  = "/home/venom/Downloads/Khavid/dclaw-main/examples/radial_slide/_output"
DEM_IN      = "/home/venom/Downloads/Khavid/dclaw-main/examples/radial_slide/DEMSkenario_3.tt3"
SHP_IN      = "/home/venom/Downloads/Khavid/dclaw-main/examples/radial_slide/Sebelum Lahar_2.shp"
DEM_OUT     = "/home/venom/Downloads/Khavid/dclaw-main/examples/radial_slide/DemFinal.tt3"
PLOT_DIR    = "/home/venom/Downloads/Khavid/dclaw-main/examples/radial_slide"

EROSI_DEPTH = 8.44
FGOUT_FGNO  = 1
T_TARGET    = 1620.0
H_THRESHOLD = 0.01
# ============================================================


# -----------------------------------------------------------
# Helper: truncate a colormap
# -----------------------------------------------------------
def truncate_colormap(cmap_name, minval=0.30, maxval=0.95, n=256):
    base = plt.cm.get_cmap(cmap_name)
    colors_array = base(np.linspace(minval, maxval, n))
    return mcolors.LinearSegmentedColormap.from_list(
        f'{cmap_name}_trunc', colors_array, N=n
    )

# ============================================================
# FIX 1a: CMAP_DEPOSIT — color composition from setplot.py (lahar_cmap)
# minimum alpha 0.22 (transparent yellow) -> dark red (alpha 1.00)
# DEPOSIT_VMAX = 15.0 m (consistent with setplot.py)
# ============================================================
DEPOSIT_VMAX = 15.0

_lahar_colors = [
    (0.000, (1.00, 0.97, 0.40, 0.22)),   # bright yellow, fairly transparent
    (0.010, (1.00, 0.90, 0.10, 0.70)),
    (0.030, (1.00, 0.85, 0.00, 0.78)),
    (0.070, (1.00, 0.70, 0.00, 0.83)),
    (0.130, (1.00, 0.55, 0.00, 0.87)),
    (0.200, (0.95, 0.35, 0.00, 0.91)),
    (0.300, (0.85, 0.15, 0.00, 0.94)),
    (0.450, (0.65, 0.00, 0.00, 0.96)),
    (0.650, (0.40, 0.00, 0.00, 0.98)),
    (1.000, (0.10, 0.00, 0.00, 1.00)),
]
_n = 512
_lahar_cmap_data = np.zeros((_n, 4))
_fracs = [c[0] for c in _lahar_colors]
_rgba  = [c[1] for c in _lahar_colors]
for i in range(_n):
    frac = i / (_n - 1)
    idx  = np.clip(np.searchsorted(_fracs, frac, side="right")-1, 0, len(_fracs)-2)
    t    = (frac - _fracs[idx]) / max(_fracs[idx+1] - _fracs[idx], 1e-12)
    _lahar_cmap_data[i] = (np.array(_rgba[idx])
                           + t * (np.array(_rgba[idx+1]) - np.array(_rgba[idx])))

CMAP_DEPOSIT = mcolors.ListedColormap(_lahar_cmap_data, name="lahar_depth")
CMAP_DEPOSIT.set_bad(alpha=0)
CMAP_DEPOSIT.set_under(alpha=0)

CMAP_EROSI = truncate_colormap('Blues_r', minval=0.25, maxval=0.95)
CMAP_EROSI.set_bad(alpha=0)


# -----------------------------------------------------------
# Read / write DEM  [UNCHANGED]
# -----------------------------------------------------------
def read_tt3(path):
    with open(path, 'r') as f:
        ncols     = int(f.readline().split()[0])
        nrows     = int(f.readline().split()[0])
        xllcorner = float(f.readline().split()[0])
        yllcorner = float(f.readline().split()[0])
        cellsize  = float(f.readline().split()[0])
        nodata    = float(f.readline().split()[0])
        data      = np.loadtxt(f)
    return {"ncols": ncols, "nrows": nrows, "xll": xllcorner,
            "yll": yllcorner, "cs": cellsize, "nodata": nodata}, data


def write_tt3(path, hdr, data):
    with open(path, 'w') as f:
        f.write(f"{hdr['ncols']}         ncols\n")
        f.write(f"{hdr['nrows']}         nrows\n")
        f.write(f"{hdr['xll']}      xllcorner\n")
        f.write(f"{hdr['yll']}      yllcorner\n")
        f.write(f"{hdr['cs']}          cellsize\n")
        f.write(f"{int(hdr['nodata'])}        nodata_value\n")
        np.savetxt(f, data, fmt='%.4f')


def make_hillshade(Z, cs):
    Z_fill = np.where(np.isnan(Z), np.nanmin(Z), Z) * 2.0
    dzdx = np.gradient(Z_fill, cs,  axis=1)
    dzdy = -np.gradient(Z_fill, cs, axis=0)
    alt  = np.radians(45)
    nx_  = -dzdx; ny_ = dzdy; nz_ = np.ones_like(Z_fill)
    mag  = np.sqrt(nx_**2 + ny_**2 + nz_**2)
    return 0.20 + 0.75 * np.clip(
        (-np.sin(alt)*nx_ - np.sin(alt)*ny_ + np.cos(alt)*nz_) / mag, 0, 1)


def build_shapefile_mask(hdr, shp_path):
    import geopandas as gpd
    from shapely.ops import unary_union

    cs  = hdr['cs']
    nx  = hdr['ncols']
    ny  = hdr['nrows']
    xll = hdr['xll']
    yll = hdr['yll']

    gdf  = gpd.read_file(shp_path)
    poly = unary_union(gdf.geometry)

    col_idx   = np.arange(nx)
    row_idx   = np.arange(ny)
    x_centers = xll + (col_idx + 0.5) * cs
    y_centers = yll + (ny - row_idx - 0.5) * cs

    XX, YY = np.meshgrid(x_centers, y_centers)

    try:
        from shapely import contains_xy
        mask = contains_xy(poly, XX.ravel(), YY.ravel()).reshape(ny, nx)
    except ImportError:
        from shapely.vectorized import contains
        mask = contains(poly, XX.ravel(), YY.ravel()).reshape(ny, nx)

    print(f"  Shapefile mask: {mask.sum():,} cells inside the polygon "
          f"({mask.sum()*cs**2/1e6:.4f} km²)")
    return mask


def read_fgout_full(output_dir, fgno, frame_num):
    t_file = os.path.join(output_dir, f"fgout{fgno:04d}.t{frame_num:04d}")
    q_file = os.path.join(output_dir, f"fgout{fgno:04d}.q{frame_num:04d}")

    with open(t_file, 'r') as f:
        lines_t = [l.strip() for l in f if l.strip()]
    t_val = float(lines_t[0].split()[0])
    nvar  = int(lines_t[1].split()[0])

    with open(q_file, 'r') as f:
        lines = f.readlines()

    patches = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1; continue
        parts = line.split()
        if len(parts) >= 2 and parts[1] == 'grid_number':
            i += 1; _    = int(lines[i].split()[0])
            i += 1; mx   = int(lines[i].split()[0])
            i += 1; my   = int(lines[i].split()[0])
            i += 1; xlow = float(lines[i].split()[0])
            i += 1; ylow = float(lines[i].split()[0])
            i += 1; dx   = float(lines[i].split()[0])
            i += 1; dy   = float(lines[i].split()[0])
            i += 1

            total_needed = nvar * mx * my
            collected = []
            while len(collected) < total_needed and i < len(lines):
                dl = lines[i].strip()
                i += 1
                if not dl:
                    continue
                p2 = dl.split()
                if len(p2) >= 2 and p2[1] == 'grid_number':
                    i -= 1; break
                for v in p2:
                    try:
                        collected.append(float(v))
                    except ValueError:
                        pass

            if len(collected) >= total_needed:
                arr = np.array(collected[:total_needed])
                q = arr.reshape((my, mx, nvar)).transpose(2, 0, 1)
                patches.append({
                    'q': q, 'mx': mx, 'my': my,
                    'xlow': xlow, 'ylow': ylow,
                    'dx': dx, 'dy': dy, 'nvar': nvar
                })
        else:
            i += 1

    return t_val, nvar, patches


def map_h_to_dem(patches, hdr, h_threshold=0.01):
    cs    = hdr['cs']
    ncols = hdr['ncols']
    nrows = hdr['nrows']
    xll   = hdr['xll']
    yll   = hdr['yll']

    h_on_dem = np.zeros((nrows, ncols), dtype=float)
    n_mapped  = 0

    for patch in patches:
        q    = patch['q']
        mx   = patch['mx']
        my   = patch['my']
        xlow = patch['xlow']
        ylow = patch['ylow']
        dx   = patch['dx']
        dy   = patch['dy']

        h_patch = q[0]

        col_arr = np.arange(mx, dtype=float)
        row_arr = np.arange(my, dtype=float)

        x_centers = xlow + (col_arr + 0.5) * dx
        y_centers = ylow + (row_arr + 0.5) * dy

        col_dem_arr = np.round((x_centers - xll) / cs).astype(int)
        row_dem_arr = (nrows - 1 - np.round((y_centers - yll) / cs)).astype(int)

        mask = h_patch > h_threshold
        rows_fg, cols_fg = np.where(mask)

        for r_fg, c_fg in zip(rows_fg, cols_fg):
            h_val   = h_patch[r_fg, c_fg]
            row_dem = row_dem_arr[r_fg]
            col_dem = col_dem_arr[c_fg]

            if 0 <= row_dem < nrows and 0 <= col_dem < ncols:
                if h_val > h_on_dem[row_dem, col_dem]:
                    h_on_dem[row_dem, col_dem] = h_val
                    n_mapped += 1

    print(f"  Total fgout cells mapped to DEM: {n_mapped:,}")
    print(f"  DEM cells with deposit > {h_threshold} m: "
          f"{(h_on_dem > h_threshold).sum():,}")
    if (h_on_dem > h_threshold).sum() > 0:
        hv = h_on_dem[h_on_dem > h_threshold]
        print(f"  h min={hv.min():.4f}  max={hv.max():.2f}  mean={hv.mean():.2f} m")
    return h_on_dem


# ===========================================================
# MAIN  [UNCHANGED]
# ===========================================================
print("=" * 60)
print("PART 1 — Read DEM & Erosion Shapefile")
print("=" * 60)

hdr, dem = read_tt3(DEM_IN)
valid_mask = dem != hdr['nodata']
print(f"  DEM grid  : {hdr['ncols']} x {hdr['nrows']}")
print(f"  Cellsize  : {hdr['cs']} m")
print(f"  xll={hdr['xll']:.1f}, yll={hdr['yll']:.1f}")
print(f"  Elevation : {dem[valid_mask].min():.1f} – {dem[valid_mask].max():.1f} m")
print(f"  Erosion   : {EROSI_DEPTH} m  <- h0 from maketopo.py")

shp_mask   = build_shapefile_mask(hdr, SHP_IN)
erosi_mask = valid_mask & shp_mask
dem_eroded = dem.copy().astype(float)
dem_eroded[erosi_mask] -= EROSI_DEPTH
print(f"\n  Eroded cells : {erosi_mask.sum():,}")
print(f"  Erosion volume : {erosi_mask.sum() * hdr['cs']**2 * EROSI_DEPTH / 1e6:.4f} million m³")

fgout_t_files = sorted([
    f for f in os.listdir(OUTPUT_DIR)
    if f.startswith(f"fgout{FGOUT_FGNO:04d}.t")
])
frame_times = {}
for fname in fgout_t_files:
    fn = int(fname.split('.t')[1])
    with open(os.path.join(OUTPUT_DIR, fname), 'r') as f:
        frame_times[fn] = float(f.readline().split()[0])

best_frame = min(frame_times, key=lambda fn: abs(frame_times[fn] - T_TARGET))
best_t     = frame_times[best_frame]

print("\n" + "=" * 60)
print("PART 2 — Read h from fgout & map onto DEM")
print("=" * 60)
print(f"\n  Frame: {best_frame:04d}  t={best_t:.1f}s ({best_t/60:.2f} min)")

t_final, nvar, patches = read_fgout_full(OUTPUT_DIR, FGOUT_FGNO, best_frame)
print(f"  NVAR={nvar}, number of patches={len(patches)}")
for ip, p in enumerate(patches):
    print(f"  Patch {ip}: mx={p['mx']}, my={p['my']}, "
          f"xlow={p['xlow']:.1f}, ylow={p['ylow']:.1f}, "
          f"dx={p['dx']:.4f}, dy={p['dy']:.4f}")
    h = p['q'][0]
    print(f"    h: min={h.min():.4f}  max={h.max():.4f}  "
          f"cells>{H_THRESHOLD}={(h>H_THRESHOLD).sum():,}")

print(f"\n  Mapping h -> DEM (explicit coordinates)...")
h_on_dem = map_h_to_dem(patches, hdr, h_threshold=H_THRESHOLD)

mask_deposit = valid_mask & (h_on_dem > H_THRESHOLD)
DemFinal     = dem_eroded.copy()
DemFinal[mask_deposit] += h_on_dem[mask_deposit]

print(f"\n  DEM updated:")
print(f"    Deposit cells : {mask_deposit.sum():,}")
if mask_deposit.sum() > 0:
    dv = h_on_dem[mask_deposit]
    print(f"    h min={dv.min():.4f}  max={dv.max():.2f}  mean={dv.mean():.2f} m")
    print(f"    Volume={( dv * hdr['cs']**2).sum()/1e6:.4f} million m³")

write_tt3(DEM_OUT, hdr, DemFinal)
print(f"\n  Saved: {DEM_OUT}")


# ===========================================================
# PART 2B — EXPORT DEPOSIT TO SHAPEFILE  [UNCHANGED]
# ===========================================================
print("\n" + "=" * 60)
print("PART 2B — Export Deposit to Shapefile")
print("=" * 60)

try:
    import geopandas as gpd
    import rasterio.features
    from rasterio.transform import from_origin
    from shapely.geometry import shape

    crs_ref = gpd.read_file(SHP_IN).crs
    deposit_array = np.where(mask_deposit, h_on_dem, 0).astype(np.float32)
    transform_raster = from_origin(
        west  = hdr['xll'],
        north = hdr['yll'] + hdr['nrows'] * hdr['cs'],
        xsize = hdr['cs'],
        ysize = hdr['cs']
    )
    shapes_gen = rasterio.features.shapes(
        deposit_array,
        mask      = mask_deposit.astype(np.uint8),
        transform = transform_raster
    )
    geometries = []; h_values = []
    for geom, value in shapes_gen:
        if value > H_THRESHOLD:
            geometries.append(shape(geom))
            h_values.append(float(value))

    if len(geometries) == 0:
        print("  No deposit geometry found -- shapefile not created.")
    else:
        cs2 = hdr['cs'] ** 2
        gdf_deposit = gpd.GeoDataFrame(
            {'h_dep_m'  : h_values,
             'vol_m3'   : [h * cs2 for h in h_values],
             'luas_m2'  : [geom.area for geom in geometries]},
            geometry=geometries, crs=crs_ref
        )
        gdf_dissolved = gdf_deposit.copy()
        gdf_dissolved['dummy'] = 1
        gdf_dissolved = (
            gdf_dissolved
            .dissolve(by='dummy', aggfunc={'h_dep_m':'mean','vol_m3':'sum','luas_m2':'sum'})
            .reset_index(drop=True)
        )
        gdf_dissolved['h_max_m'] = max(h_values)
        gdf_dissolved['h_min_m'] = min(h_values)

        SHP_DEPOSIT           = os.path.join(PLOT_DIR, "deposit_demfinal.shp")
        SHP_DEPOSIT_DISSOLVED = os.path.join(PLOT_DIR, "deposit_demfinal_dissolved.shp")
        gdf_deposit.to_file(SHP_DEPOSIT)
        gdf_dissolved.to_file(SHP_DEPOSIT_DISSOLVED)
        print(f"  {SHP_DEPOSIT}  ({len(gdf_deposit):,} polygons)")
        print(f"  {SHP_DEPOSIT_DISSOLVED}")
        print(f"    Volume : {gdf_deposit['vol_m3'].sum()/1e6:.4f} million m³")
        print(f"    Area   : {gdf_deposit['luas_m2'].sum()/1e6:.4f} km²")

except ImportError as e:
    print(f"  Library not found: {e}")
    print("    Run: pip install rasterio geopandas")
except Exception as e:
    print(f"  Failed to export shapefile: {e}")
    import traceback; traceback.print_exc()


# ===========================================================
# PART 3 — 4-Panel Visualization
# ===========================================================
print("\n" + "=" * 60)
print("PART 3 — Visualization")
print("=" * 60)

dem_old_f   = dem.astype(float);  dem_old_f[~valid_mask]   = np.nan
dem_ero_f   = dem_eroded.copy();  dem_ero_f[~valid_mask]   = np.nan
dem_final_f = DemFinal.copy();    dem_final_f[~valid_mask] = np.nan

valid_s3     = ~np.isnan(dem_old_f) & ~np.isnan(dem_final_f)
diff_total   = np.where(valid_s3, dem_final_f - dem_old_f, np.nan)
diff_erosi   = np.where(valid_s3 & (diff_total < -H_THRESHOLD), diff_total, np.nan)
diff_deposit = np.where(valid_s3 & (diff_total >  H_THRESHOLD), diff_total, np.nan)

xmin = hdr['xll'];  xmax = hdr['xll'] + hdr['ncols'] * hdr['cs']
ymin = hdr['yll'];  ymax = hdr['yll'] + hdr['nrows'] * hdr['cs']
extent = [xmin, xmax, ymin, ymax]

hs       = make_hillshade(dem_old_f, hdr['cs'])
vmin_dem = np.nanpercentile(dem_old_f, 2)
vmax_dem = np.nanpercentile(dem_old_f, 98)

hs_kw  = dict(cmap='gray',    extent=extent, origin='upper',
              aspect='equal', vmin=0.20, vmax=0.95)
dem_kw = dict(cmap='terrain', extent=extent, origin='upper',
              aspect='equal', alpha=0.45,
              vmin=vmin_dem, vmax=vmax_dem)
ov_kw  = dict(extent=extent,  origin='upper', aspect='equal')

dx_data = xmax - xmin
dy_data = ymax - ymin
data_ratio = dy_data / dx_data

PANEL_W_IN = 7.0
PANEL_H_IN = min(PANEL_W_IN * data_ratio, 22)
fig_w = PANEL_W_IN * 4 + 1.2 * 2 + 2.2 * 2 + 2.0
fig_h = PANEL_H_IN + 2.0

fig, axes = plt.subplots(
    1, 4,
    figsize=(fig_w, fig_h),
    gridspec_kw=dict(wspace=0.55, left=0.06, right=0.96, top=0.91, bottom=0.11)
)
ax1, ax2, ax3, ax4 = axes

fig.suptitle(
    f"DemFinal Verification — Merapi Lahar Erosion & Deposit (Scenario 3)\n"
    f"Erosion: {EROSI_DEPTH} m  |  Deposit fgout setrun_3  |  t = {t_final:.0f}s ({t_final/60:.0f} min)",
    fontsize=12, fontweight='bold', y=0.995
)

# ── FIX 2: font sizes increased for readability & to avoid overlap ──
FS_TITLE  = 13   # panel title
FS_LABEL  = 11   # xlabel/ylabel
FS_TICK   = 10   # tick labels
FS_CB     = 11   # colorbar label
FS_CBTICK = 10   # colorbar tick labels
FS_LEGEND = 11   # legend
FS_STATS  = 9    # stats text box

def add_colorbar(fig, ax, im, label, size="4%", pad=0.05, extend='neither'):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=size, pad=pad)
    cb  = fig.colorbar(im, cax=cax, extend=extend)
    cb.set_label(label, fontsize=FS_CB, labelpad=4)
    cb.ax.tick_params(labelsize=FS_CBTICK, pad=2)
    return cax, cb

def _style_ax(ax, title, show_ylabel=False):
    ax.set_title(title, fontsize=FS_TITLE, pad=5, fontweight='bold')
    ax.set_xlabel('Easting (m)', fontsize=FS_LABEL)
    ax.ticklabel_format(useOffset=False)
    ax.tick_params(axis='x', rotation=45, labelsize=FS_TICK, pad=2)
    ax.tick_params(axis='y', labelsize=FS_TICK, pad=2)
    if show_ylabel:
        ax.set_ylabel('Northing (m)', fontsize=FS_LABEL, labelpad=4)
    else:
        ax.set_ylabel('')
        ax.set_yticklabels([])

# Panel 1 — Original DEM
ax1.imshow(hs, **hs_kw)
im1 = ax1.imshow(dem_old_f, **dem_kw)
add_colorbar(fig, ax1, im1, 'Elevation (m)')
_style_ax(ax1, 'DEMSkenario_2', show_ylabel=True)

# Panel 2 — DEM after erosion
ax2.imshow(hs, **hs_kw)
im2 = ax2.imshow(dem_ero_f, **dem_kw)
add_colorbar(fig, ax2, im2, 'Elevation (m)')
_style_ax(ax2, f'DEMSkenario_2 Eroded\n(-{EROSI_DEPTH} m at source)')

import geopandas as gpd
gdf = gpd.read_file(SHP_IN)
gdf.boundary.plot(ax=ax2, color='cyan', linewidth=1.0, label='Source boundary')
ax2.legend(fontsize=FS_LEGEND, loc='lower right',
           framealpha=0.85, edgecolor='gray')

# Panel 3 — DemFinal + deposit overlay
# FIX 1: vmax=DEPOSIT_VMAX (15.0), alpha=0.88, extend='neither'
dep_ov = np.where(mask_deposit, h_on_dem, np.nan)

ax3.imshow(hs, **hs_kw)
im3a = ax3.imshow(dem_final_f, **dem_kw)
im3b = ax3.imshow(dep_ov, cmap=CMAP_DEPOSIT, alpha=0.88,
                   vmin=H_THRESHOLD, vmax=DEPOSIT_VMAX, **ov_kw)

divider3 = make_axes_locatable(ax3)
cax3a    = divider3.append_axes("right", size="4%", pad=0.05)
# FIX 2: pad 0.42 -> 0.70 to avoid overlap
cax3b    = divider3.append_axes("right", size="4%", pad=0.70)

cb3a = fig.colorbar(im3a, cax=cax3a, extend='neither')
cb3a.set_label('Elevation (m)', fontsize=FS_CB, labelpad=4)
cb3a.ax.tick_params(labelsize=FS_CBTICK, pad=2)

# FIX 1: extend='neither' (not 'max')
cb3b = fig.colorbar(im3b, cax=cax3b, extend='neither')
cb3b.set_label('Deposit h (m)', fontsize=FS_CB, labelpad=4)
cb3b.ax.tick_params(labelsize=FS_CBTICK, pad=2)

_style_ax(ax3, 'DemSkenario_2\n(Erosion + Scenario 2 deposit)')

# Panel 4 — Difference map
vmin_e = float(np.nanmin(diff_erosi)) if not np.all(np.isnan(diff_erosi)) else -1

ax4.imshow(hs, **hs_kw)
im_e = ax4.imshow(diff_erosi,   cmap=CMAP_EROSI,   alpha=0.85,
                   vmin=vmin_e, vmax=0, **ov_kw)
# FIX 1: vmax=DEPOSIT_VMAX (15.0), alpha=0.88
im_d = ax4.imshow(diff_deposit, cmap=CMAP_DEPOSIT, alpha=0.88,
                   vmin=H_THRESHOLD, vmax=DEPOSIT_VMAX, **ov_kw)

divider4 = make_axes_locatable(ax4)
cax4a    = divider4.append_axes("right", size="4%", pad=0.05)
# FIX 2: pad 0.42 -> 0.70 to avoid overlap
cax4b    = divider4.append_axes("right", size="4%", pad=0.70)

cb4a = fig.colorbar(im_e, cax=cax4a, extend='neither')
cb4a.set_label('Erosion (m)', fontsize=FS_CB, labelpad=4)
cb4a.ax.tick_params(labelsize=FS_CBTICK, pad=2)

# FIX 1: extend='neither' (not 'max')
cb4b = fig.colorbar(im_d, cax=cax4b, extend='neither')
cb4b.set_label('Deposit (m)', fontsize=FS_CB, labelpad=4)
cb4b.ax.tick_params(labelsize=FS_CBTICK, pad=2)

_style_ax(ax4, 'Total Difference Map\n(DEMSkenario_2 - DEMSkenario_1)')

if not np.all(np.isnan(diff_deposit)):
    dv    = diff_deposit[~np.isnan(diff_deposit)]
    stats = (f"Erosion: -{EROSI_DEPTH}m x {erosi_mask.sum():,} cells\n"
             f"Dep h : max {dv.max():.2f}m | avg {dv.mean():.2f}m\n"
             f"Vol   : {(dv*hdr['cs']**2).sum()/1e6:.4f} million m³\n"
             f"t     : {t_final:.0f}s")
else:
    stats = "Deposit: none"

ax4.text(0.02, 0.02, stats, transform=ax4.transAxes, fontsize=FS_STATS,
         va='bottom', fontfamily='monospace',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                   alpha=0.88, edgecolor='gray', linewidth=0.6))

out_map = os.path.join(PLOT_DIR, "DemSkenario_2_difference_map.png")
plt.savefig(out_map, dpi=200, bbox_inches='tight')
plt.show()
print(f"  Map saved: {out_map}")

# Deposit histogram
if not np.all(np.isnan(diff_deposit)):
    dv = diff_deposit[~np.isnan(diff_deposit)]
    fig2, ax = plt.subplots(figsize=(8, 5))
    n_bins = 60
    counts, edges, patches_bar = ax.hist(dv, bins=n_bins, edgecolor='none')
    norm_v = (edges[:-1] - edges[:-1].min()) / (edges[:-1].max() - edges[:-1].min() + 1e-12)
    for patch, nv in zip(patches_bar, norm_v):
        patch.set_facecolor(CMAP_DEPOSIT(nv))
    ax.axvline(dv.mean(),     color='navy',  ls='--', lw=1.5,
               label=f'Mean: {dv.mean():.2f} m')
    ax.axvline(np.median(dv), color='green', ls='-.', lw=1.5,
               label=f'Median: {np.median(dv):.2f} m')
    ax.set_xlabel('Deposit thickness (m)', fontsize=10)
    ax.set_ylabel('Cell count', fontsize=10)
    ax.set_title(f'DemFinal Deposit Thickness Distribution — t={t_final:.0f}s', fontsize=11)
    ax.legend(fontsize=9); ax.grid(True, alpha=0.4)
    plt.tight_layout()
    out_hist = os.path.join(PLOT_DIR, "DemSkenario_2_deposit_histogram.png")
    plt.savefig(out_hist, dpi=200, bbox_inches='tight')
    plt.show()
    print(f"  Histogram saved: {out_hist}")

print("\nDONE!")
