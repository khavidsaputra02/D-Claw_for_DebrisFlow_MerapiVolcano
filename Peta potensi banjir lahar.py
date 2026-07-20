import os
import sys
import warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.colors import LightSource
from matplotlib.ticker import FuncFormatter
from matplotlib.cm import ScalarMappable
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
from shapely.ops import unary_union
from pyproj import Transformer

warnings.filterwarnings("ignore")

# ═══════════════════════════════════════════════════════════════════════════
#  CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════

FILE_DEM     = "DEMPredict.tt3"
FILE_SHP_SEL = "Laharselatan.shp"
FILE_SHP_BAR = "Laharbarat-baratdaya.shp"   # <-- FIXED: name adjusted (has a hyphen)
FILE_SUNGAI  = "sungai.geojson"
FILE_JALAN   = "jalan.geojson"
FILE_PERMUKIMAN = "pemukiman_merapi_scale25ribu.gpkg"   # <-- NEW: settlement layer
LAYER_PERMUKIMAN = "permukiman_area"                     # layer name inside the .gpkg

WARNA_PERMUKIMAN = {"fill": "#5DADE2", "alpha": 0.85, "edge": "#1B4F72"}  # sky blue, solid

WORK_DIR = r"D:\Tugas Akhir\Peta Potensi Banjir Lahar Dingin"

_SEARCH_DIRS = [
    WORK_DIR,
    ".",
    os.path.dirname(os.path.abspath(__file__)),
]

def _resolve_path(filename: str) -> str:
    if os.path.isabs(filename) and os.path.exists(filename):
        return filename
    for d in _SEARCH_DIRS:
        p = os.path.join(d, os.path.basename(filename))
        if os.path.exists(p):
            return p
    return os.path.join(WORK_DIR, os.path.basename(filename))

CRS_UTM = "EPSG:32749"
CRS_WGS = "EPSG:4326"

NAMA_GUNUNG = "Gunung Merapi"
NAMA_PETA   = "COLD LAHAR FLOOD HAZARD POTENTIAL MAP"

WARNA_LAHAR = {"fill": "#FFE033", "alpha": 0.60, "edge": "#C8A800"}
RADIUS_LAHAR = 8564


# ═══════════════════════════════════════════════════════════════════════════
#  READ DEM ASCII GRID
# ═══════════════════════════════════════════════════════════════════════════

class DEMData:
    def __init__(self, fpath: str):
        self.fpath      = fpath
        self.header     = {}
        self.dem        = None
        self.hillshade  = None
        self.extent_utm = None
        self.extent_wgs = None
        self.peak_x_utm = None
        self.peak_y_utm = None
        self.peak_lon   = None
        self.peak_lat   = None
        self.peak_elev  = None
        self._tr_to_wgs = Transformer.from_crs(CRS_UTM, CRS_WGS, always_xy=True)
        self._load()
        self._compute_hillshade()
        self._detect_peak()

    def _load(self):
        print(f"[DEM] Reading {self.fpath} ...")
        with open(self.fpath, 'r') as f:
            for _ in range(6):
                line = f.readline().split()
                try:
                    self.header[line[1].lower()] = float(line[0])
                except (ValueError, IndexError):
                    self.header[line[0].lower()] = float(line[1])

        ncols  = int(self.header['ncols'])
        nrows  = int(self.header['nrows'])
        xll    = self.header['xllcorner']
        yll    = self.header['yllcorner']
        cell   = self.header['cellsize']
        nodata = self.header['nodata_value']
        print(f"  -> {nrows}x{ncols} pixels  |  cell={cell:.2f} m")

        self.dem = np.loadtxt(self.fpath, skiprows=6, dtype=np.float32)
        self.dem[self.dem == nodata] = np.nan

        xmax = xll + ncols * cell
        ymax = yll + nrows * cell
        self.extent_utm = (xll, xmax, yll, ymax)

        lon_ll, lat_ll = self._tr_to_wgs.transform(xll, yll)
        lon_ur, lat_ur = self._tr_to_wgs.transform(xmax, ymax)
        self.extent_wgs = (min(lon_ll, lon_ur), max(lon_ll, lon_ur),
                           min(lat_ll, lat_ur), max(lat_ll, lat_ur))

        valid = self.dem[~np.isnan(self.dem)]
        print(f"  -> Elevation: min={valid.min():.1f}m  max={valid.max():.1f}m")

    def _compute_hillshade(self):
        print("[DEM] Computing hillshade ...")
        dem_filled = self.dem.copy()
        mask_nan = np.isnan(dem_filled)
        if mask_nan.any():
            dem_filled[mask_nan] = np.nanmean(dem_filled)
        cell = self.header['cellsize']
        ls = LightSource(azdeg=315, altdeg=40)
        self.hillshade = ls.hillshade(dem_filled, vert_exag=2.5, dx=cell, dy=cell)

    def _detect_peak(self):
        iy, ix = np.unravel_index(np.nanargmax(self.dem), self.dem.shape)
        cell     = self.header['cellsize']
        xll      = self.header['xllcorner']
        yll      = self.header['yllcorner']
        nrows    = int(self.header['nrows'])
        ymax_utm = yll + nrows * cell
        px = xll + (ix + 0.5) * cell
        py = ymax_utm - (iy + 0.5) * cell
        self.peak_x_utm = px
        self.peak_y_utm = py
        lon, lat = self._tr_to_wgs.transform(px, py)
        self.peak_lon  = lon
        self.peak_lat  = lat
        self.peak_elev = float(np.nanmax(self.dem))
        print(f"  -> Peak: {self.peak_elev:.1f} m  @ UTM ({px:.0f}, {py:.0f})")


# ═══════════════════════════════════════════════════════════════════════════
#  LOAD & PROCESS LAHAR SHAPEFILES
# ═══════════════════════════════════════════════════════════════════════════

class LaharData:
    def __init__(self, shp_selatan: str, dem: DEMData, shp_barat=None):
        self.dem = dem
        self.gdf_per_source = {}
        self.gdf_gabungan_utm = None
        self.zona_gdf = {}

        if not os.path.exists(shp_selatan):
            raise FileNotFoundError(f"Shapefile not found: {shp_selatan}")
        self.gdf_per_source["selatan"] = self._load_shp(shp_selatan, "selatan")

        if shp_barat and os.path.exists(shp_barat):
            self.gdf_per_source["barat-baratdaya"] = self._load_shp(shp_barat, "barat-baratdaya")
        elif shp_barat:
            print(f"  [!] File '{shp_barat}' not found, skipped.")

        all_gdfs = list(self.gdf_per_source.values())
        if len(all_gdfs) == 1:
            self.gdf_gabungan_utm = all_gdfs[0].copy()
        else:
            self.gdf_gabungan_utm = gpd.GeoDataFrame(
                pd.concat(all_gdfs, ignore_index=True), crs=CRS_UTM)

        self._buat_zona_lahar()

    @staticmethod
    def _smooth_geom(geom, buffer_out=25, buffer_in=-22, quad_segs=32):
        smoothed = geom.buffer(buffer_out, quad_segs=quad_segs).buffer(buffer_in, quad_segs=quad_segs)
        return smoothed if not smoothed.is_empty else geom

    def _load_shp(self, fpath: str, label: str) -> gpd.GeoDataFrame:
        print(f"\n[lahar] Reading {os.path.basename(fpath)} ({label}) ...")
        gdf = gpd.read_file(fpath)
        if gdf.crs is None:
            gdf = gdf.set_crs(CRS_UTM)
        else:
            gdf = gdf.to_crs(CRS_UTM)
        n = len(gdf)
        if n > 1000:
            print(f"  -> {n} polygons, dissolving ...")
            union_geom = unary_union(gdf.geometry.values)
            merged = union_geom.buffer(6, quad_segs=16).buffer(-5, quad_segs=16)
            gdf = gpd.GeoDataFrame({"geometry": [merged], "sumber": [label]}, crs=CRS_UTM)
        gdf["geometry"] = gdf["geometry"].apply(
            lambda g: self._smooth_geom(g, buffer_out=25, buffer_in=-22, quad_segs=32))
        gdf["sumber"] = label
        return gdf

    def _buat_zona_lahar(self):
        print("\n[lahar] Building combined lahar zone ...")
        peak_utm = Point(self.dem.peak_x_utm, self.dem.peak_y_utm)
        union_all = unary_union(
            [unary_union(gdf.geometry.values) for gdf in self.gdf_per_source.values()])
        if union_all.is_empty or union_all.area < 1000:
            print("  [!] Lahar zone empty -> using full buffer")
            union_all = peak_utm.buffer(RADIUS_LAHAR, quad_segs=64)
        zona = self._smooth_geom(union_all, buffer_out=20, buffer_in=-18, quad_segs=48)
        gdf_utm = gpd.GeoDataFrame({"nama": ["Potensi Banjir Lahar Dingin"]},
                                    geometry=[zona], crs=CRS_UTM)
        self.zona_gdf = {"utm": gdf_utm, "wgs": gdf_utm.to_crs(CRS_WGS)}
        print(f"  -> Lahar zone: {zona.area/1e6:.2f} km2")


# ═══════════════════════════════════════════════════════════════════════════
#  LOAD ADDITIONAL VECTOR LAYERS
# ═══════════════════════════════════════════════════════════════════════════

def load_vektor_tambahan():
    sungai = None
    p = _resolve_path(FILE_SUNGAI)
    if os.path.exists(p):
        sungai = gpd.read_file(p).to_crs(CRS_UTM)
    jalan = None
    p = _resolve_path(FILE_JALAN)
    if os.path.exists(p):
        jalan = gpd.read_file(p).to_crs(CRS_UTM)
    return sungai, jalan


# NEW: dedicated function to load the settlement layer and clip it to the DEM extent
def load_pemukiman(dem: DEMData):
    """Reads the settlement layer from a GeoPackage, reprojects to CRS_UTM,
    then clips it to the DEM extent so it stays lightweight and relevant to the study area."""
    p = _resolve_path(FILE_PERMUKIMAN)
    if not os.path.exists(p):
        print(f"  [!] Settlement file '{FILE_PERMUKIMAN}' not found, skipped.")
        return None

    print(f"\n[pemukiman] Reading {os.path.basename(p)} (layer '{LAYER_PERMUKIMAN}') ...")
    try:
        gdf = gpd.read_file(p, layer=LAYER_PERMUKIMAN)
    except Exception as e:
        print(f"  [!] Failed to read layer '{LAYER_PERMUKIMAN}': {e}")
        print(f"      Trying to read the default layer instead ...")
        gdf = gpd.read_file(p)

    if gdf.crs is None:
        gdf = gdf.set_crs(CRS_UTM)
    else:
        gdf = gdf.to_crs(CRS_UTM)

    n_awal = len(gdf)
    xmin_utm, xmax_utm, ymin_utm, ymax_utm = dem.extent_utm
    bbox = (xmin_utm, ymin_utm, xmax_utm, ymax_utm)
    gdf = gdf.cx[bbox[0]:bbox[2], bbox[1]:bbox[3]]
    if len(gdf) == 0:
        print(f"  [!] No settlement polygons within the DEM extent (out of {n_awal} total).")
        return None

    print(f"  -> {len(gdf)} of {n_awal} settlement polygons fall within the DEM extent.")
    return gdf


# ═══════════════════════════════════════════════════════════════════════════
#  STATIC MAP
# ═══════════════════════════════════════════════════════════════════════════

def buat_peta_statis(dem: DEMData, lahar: LaharData, sungai, jalan, pemukiman=None, simpan_ke=None):
    if simpan_ke is None:
        simpan_ke = os.path.join(WORK_DIR, "peta_banjir_lahar_static.png")

    print("\n[peta] Creating static map ...")

    xmin_utm, xmax_utm, ymin_utm, ymax_utm = dem.extent_utm
    xmin_plot, xmax_plot = xmin_utm, xmax_utm
    ymin_plot, ymax_plot = ymin_utm, ymax_utm
    ext_dem = [xmin_utm, xmax_utm, ymin_utm, ymax_utm]

    L_in = 1.50
    R_in = 0.20
    T_in = 0.20
    B_in = 0.90
    CBAR_H_IN = 1.60

    data_w = xmax_plot - xmin_plot
    data_h = ymax_plot - ymin_plot
    ratio  = data_h / data_w

    ax_w_in = 16.0
    ax_h_in = ax_w_in * ratio

    fig_w = ax_w_in + L_in + R_in
    fig_h = ax_h_in + T_in + B_in + CBAR_H_IN

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=600)
    fig.patch.set_facecolor("white")

    ax = fig.add_axes([
        L_in / fig_w,
        (CBAR_H_IN + B_in) / fig_h,
        ax_w_in / fig_w,
        ax_h_in / fig_h,
    ])
    ax.set_facecolor("#c8c8c8")

    ax.imshow(dem.hillshade, extent=ext_dem, cmap="gray",
              vmin=0.0, vmax=1.0, aspect="auto", origin="upper",
              interpolation="bilinear", zorder=2)

    vmin_elev = float(np.nanpercentile(dem.dem, 2))
    vmax_elev = float(np.nanmax(dem.dem))
    dem_masked = np.ma.masked_invalid(dem.dem)
    ax.imshow(dem_masked, extent=ext_dem, cmap="gray",
              vmin=vmin_elev, vmax=vmax_elev, aspect="auto",
              origin="upper", interpolation="bilinear", zorder=3, alpha=0.30)
    ax.imshow(dem.hillshade, extent=ext_dem, cmap="gray",
              vmin=0.0, vmax=1.0, aspect="auto", origin="upper",
              interpolation="bilinear", zorder=4, alpha=0.35)

    ax.set_xlim(xmin_plot, xmax_plot)
    ax.set_ylim(ymin_plot, ymax_plot)

    style_per_src = {
        "selatan":         {"color": "#888888", "lw": 0.35, "ls": "-"},
        "barat-baratdaya": {"color": "#333333", "lw": 0.35, "ls": "-"},
    }
    for label, gdf in lahar.gdf_per_source.items():
        sty = style_per_src.get(label, {"color": "#888888", "lw": 0.25, "ls": "-"})
        gdf.plot(ax=ax, color="none", edgecolor=sty["color"],
                 linewidth=sty["lw"], linestyle=sty["ls"], alpha=0.55, zorder=5)

    if lahar.zona_gdf:
        lahar.zona_gdf["utm"].plot(ax=ax, color=WARNA_LAHAR["fill"],
                                   alpha=WARNA_LAHAR["alpha"],
                                   edgecolor=WARNA_LAHAR["edge"],
                                   linewidth=1.2, zorder=6)

    _buat_kontur(ax, dem)

    # NEW: plot the settlement layer - solid color (sky blue)
    if pemukiman is not None and len(pemukiman) > 0:
        pemukiman.plot(ax=ax, color=WARNA_PERMUKIMAN["fill"],
                        alpha=WARNA_PERMUKIMAN["alpha"],
                        edgecolor=WARNA_PERMUKIMAN["edge"],
                        linewidth=0.25, zorder=7.5)

    if jalan is not None:
        jalan.plot(ax=ax, color="#aaaaaa", linewidth=0.7, alpha=0.65, zorder=8)
    if sungai is not None:
        sungai.plot(ax=ax, color="#5599cc", linewidth=1.0, alpha=0.80, zorder=8)

    ax.plot(dem.peak_x_utm, dem.peak_y_utm,
            marker="^", ms=20, color="#FF2200", mec="white", mew=1.8, zorder=14)
    dy_label = (ymax_plot - ymin_plot) * 0.012
    ax.text(dem.peak_x_utm, dem.peak_y_utm + dy_label, NAMA_GUNUNG,
            fontsize=15, fontweight="bold", color="#111111",
            va="bottom", ha="center", linespacing=1.4,
            path_effects=[pe.withStroke(linewidth=2.5, foreground="white")],
            zorder=15)

    ax.set_xlim(xmin_plot, xmax_plot)
    ax.set_ylim(ymin_plot, ymax_plot)

    ax.set_xlabel("Easting UTM (m)", fontsize=16, color="#222222", labelpad=8)
    ax.set_ylabel("Northing UTM (m)", fontsize=16, color="#222222", labelpad=8)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{v:,.0f}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{v:,.0f}"))
    ax.tick_params(labelsize=14, colors="#222222", pad=4)
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center")
    for sp in ax.spines.values():
        sp.set_edgecolor("#333333")
    ax.grid(True, linestyle="--", linewidth=0.3, alpha=0.40, color="#888888", zorder=16)

    _buat_judul(ax)
    _buat_legenda(ax, ada_pemukiman=(pemukiman is not None and len(pemukiman) > 0))
    _buat_skala_dan_utara(ax, xmin_plot, ymin_plot, xmax_plot, ymax_plot)

    cb_l = (L_in + ax_w_in * 0.03) / fig_w
    cb_w = (ax_w_in * 0.70) / fig_w
    cb_h = 0.018
    cb_b = (CBAR_H_IN * 0.52) / fig_h

    cax = fig.add_axes([cb_l, cb_b, cb_w, cb_h])
    sm  = ScalarMappable(cmap="gray",
                         norm=mcolors.Normalize(vmin=vmin_elev, vmax=vmax_elev))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_label("Elevation (m asl)", fontsize=12, color="#222222", labelpad=3)
    cbar.ax.tick_params(labelsize=11, colors="#222222")
    cbar.outline.set_edgecolor("#444444")

    # FIXED: check FILE_SHP_BAR existence directly, instead of re-resolving the old name
    shp_bar_ada = os.path.exists(_resolve_path(FILE_SHP_BAR))
    shp_labels = " + ".join(
        os.path.basename(f) for f in
        ([FILE_SHP_SEL] + ([FILE_SHP_BAR] if shp_bar_ada else [])))
    fig.text(L_in / fig_w, 0.005,
             f"Source: PVMBG / Geological Agency, Ministry of Energy and Mineral Resources\n"
             f"DEM: {os.path.basename(FILE_DEM)}  |  Lahar: {shp_labels}  "
             f"|  Projection: UTM Zone 49S (EPSG:32749)",
             fontsize=9.5, color="#444444", va="bottom", ha="left", style="italic")

    plt.savefig(simpan_ke, dpi=600, bbox_inches=None,
                facecolor="white", edgecolor="none")
    plt.close()
    print(f"[peta] Saved -> {simpan_ke}")
    return simpan_ke


def _buat_kontur(ax, dem: DEMData):
    nrows = int(dem.header['nrows'])
    ncols = int(dem.header['ncols'])
    xmin_utm, xmax_utm, ymin_utm, ymax_utm = dem.extent_utm
    xs = np.linspace(xmin_utm, xmax_utm, ncols)
    ys = np.linspace(ymax_utm, ymin_utm, nrows)
    X, Y = np.meshgrid(xs, ys)
    vmax     = float(np.nanmax(dem.dem))
    vmin_pos = float(np.nanmin(dem.dem[dem.dem > 0])) if np.any(dem.dem > 0) else 0
    interval = max(50, round((vmax - vmin_pos) / 12 / 50) * 50)
    levels   = np.arange(round(vmin_pos / interval) * interval, vmax + interval, interval)
    cs = ax.contour(X, Y, dem.dem, levels=levels,
                    colors="#ffffff", linewidths=0.2, alpha=0.25, zorder=7)
    ax.clabel(cs, levels=levels[::3], fmt="%d m", fontsize=8,
              inline=True, inline_spacing=2, colors="#ffffff")


def _buat_judul(ax):
    props = dict(boxstyle="square,pad=0.6", facecolor="white",
                 edgecolor="#333333", linewidth=1.2, alpha=0.93)
    ax.text(0.99, 0.99, f"{NAMA_PETA}\n{NAMA_GUNUNG.upper()}",
            transform=ax.transAxes, fontsize=15, fontweight="bold",
            ha="right", va="top", color="#111111", bbox=props,
            zorder=20, linespacing=1.6)


def _buat_legenda(ax, ada_pemukiman=False):
    patches = [
        mpatches.Patch(facecolor=WARNA_LAHAR["fill"], edgecolor=WARNA_LAHAR["edge"],
                       alpha=0.8, label="Cold Lahar Flood Potential"),
    ]
    if ada_pemukiman:  # NEW
        patches.append(
            mpatches.Patch(facecolor=WARNA_PERMUKIMAN["fill"], edgecolor=WARNA_PERMUKIMAN["edge"],
                           alpha=WARNA_PERMUKIMAN["alpha"], label="Settlement Area"))
    patches.append(
        Line2D([0], [0], marker="^", color="none",
               markerfacecolor="#FF2200", markeredgecolor="#333333",
               markersize=14, label=f"{NAMA_GUNUNG} Peak"),
    )
    leg = ax.legend(handles=patches, loc="lower right", fontsize=13,
                    title="LEGEND", title_fontsize=15,
                    framealpha=0.95, edgecolor="#333333", frameon=True,
                    borderpad=0.8, handlelength=2.0,
                    facecolor="white", labelcolor="#111111")
    leg.get_title().set_color("#111111")
    leg.get_frame().set_linewidth(1.2)
    leg.set_zorder(30)  # FIXED: make sure the legend always stays above all map layers


def _buat_skala_dan_utara(ax, xmin, ymin, xmax, ymax):
    skala_m = 2000.0
    x0 = xmin + (xmax - xmin) * 0.04
    y0 = ymin + (ymax - ymin) * 0.04

    ax.plot([x0, x0 + skala_m], [y0, y0],
            color="#111111", linewidth=3.5, zorder=17, solid_capstyle="butt")
    tick_h = (ymax - ymin) * 0.007
    for xp in [x0, x0 + skala_m]:
        ax.plot([xp, xp], [y0 - tick_h/2, y0 + tick_h/2],
                color="#111111", lw=2, zorder=17)
    ax.text(x0 + skala_m/2, y0 + (ymax - ymin) * 0.010, "2 km",
            fontsize=12, ha="center", va="bottom", fontweight="bold",
            color="#111111", zorder=17,
            path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    ax.text(x0, y0 - (ymax - ymin) * 0.010, "0",
            fontsize=11, ha="center", va="top", color="#111111", zorder=17,
            path_effects=[pe.withStroke(linewidth=1.5, foreground="white")])

    gap     = (xmax - xmin) * 0.030
    arrow_h = (ymax - ymin) * 0.055
    ax_x    = x0 + skala_m + gap + (xmax - xmin) * 0.020
    ax_y0   = y0 - arrow_h * 0.25
    ax_y1   = ax_y0 + arrow_h
    ax.annotate("", xy=(ax_x, ax_y1), xytext=(ax_x, ax_y0),
                arrowprops=dict(arrowstyle="-|>", color="#111111", lw=2.2,
                                mutation_scale=18), zorder=21)
    ax.text(ax_x, ax_y1 + (ymax - ymin) * 0.010, "N",
            fontsize=18, fontweight="bold", ha="center", va="bottom",
            color="#111111", zorder=21,
            path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])


# ═══════════════════════════════════════════════════════════════════════════
#  INTERACTIVE MAP — Folium
# ═══════════════════════════════════════════════════════════════════════════

def buat_peta_interaktif(dem: DEMData, lahar: LaharData, sungai, jalan, pemukiman=None, simpan_ke=None):
    if simpan_ke is None:
        simpan_ke = os.path.join(WORK_DIR, "peta_banjir_lahar.html")
    try:
        import folium
        from PIL import Image
        import io, base64
    except ImportError as e:
        print(f"[folium] Missing dependency: {e}")
        return None

    print("\n[folium] Creating interactive map ...")
    lo, hi, la, ua = dem.extent_wgs
    m = folium.Map(location=[dem.peak_lat, dem.peak_lon], zoom_start=12, tiles=None)
    folium.TileLayer("CartoDB DarkMatter", name="Dark (default)").add_to(m)
    folium.TileLayer("CartoDB positron", name="Light").add_to(m)
    folium.TileLayer("OpenStreetMap", name="OpenStreetMap").add_to(m)

    hs_uri  = _hillshade_ke_png(dem)
    dem_uri = _dem_gray_ke_png(dem)

    folium.raster_layers.ImageOverlay(
        image=dem_uri, bounds=[[la, lo], [ua, hi]],
        opacity=0.60, name="DEM Elevation (grayscale)", show=True, interactive=False, zindex=1
    ).add_to(m)
    folium.raster_layers.ImageOverlay(
        image=hs_uri, bounds=[[la, lo], [ua, hi]],
        opacity=0.75, name="Hillshade Relief", show=True, interactive=False, zindex=2
    ).add_to(m)

    src_style = {
        "selatan":         {"color": "#888888", "weight": 1.2, "fillOpacity": 0},
        "barat-baratdaya": {"color": "#333333", "weight": 1.2, "fillOpacity": 0, "dashArray": "4 3"},
    }
    src_label = {
        "selatan": "Southern Lahar Boundary",
        "barat-baratdaya": "West-Southwest Lahar Boundary",
    }
    for label, gdf_utm in lahar.gdf_per_source.items():
        gdf_wgs = gdf_utm.to_crs(CRS_WGS)
        sty = src_style.get(label, {"color": "#aaaaaa", "weight": 1.0, "fillOpacity": 0})
        fg  = folium.FeatureGroup(name=src_label.get(label, label), show=True)
        folium.GeoJson(gdf_wgs.__geo_interface__,
                       style_function=lambda f, s=sty: s,
                       tooltip=src_label.get(label, label)).add_to(fg)
        fg.add_to(m)

    if lahar.zona_gdf:
        fg_zona = folium.FeatureGroup(name="Cold Lahar Flood Potential", show=True)
        folium.GeoJson(
            lahar.zona_gdf["wgs"].__geo_interface__,
            style_function=lambda f: {"fillColor": "#FFE033", "color": "#C8A800",
                                       "weight": 1.5, "fillOpacity": 0.55},
            highlight_function=lambda f: {"fillOpacity": 0.80, "weight": 3.0},
            popup=folium.Popup("<b>Cold Lahar Flood Potential</b>", max_width=250),
            tooltip="Cold Lahar Flood Potential",
        ).add_to(fg_zona)
        fg_zona.add_to(m)

    # NEW: settlement layer
    if pemukiman is not None and len(pemukiman) > 0:
        fg_p = folium.FeatureGroup(name="Settlement Area", show=True)
        folium.GeoJson(
            pemukiman.to_crs(CRS_WGS).__geo_interface__,
            style_function=lambda f: {"fillColor": "#5DADE2", "color": "#1B4F72",
                                       "weight": 0.8, "fillOpacity": 0.85},
            highlight_function=lambda f: {"fillOpacity": 1.0, "weight": 1.8},
            tooltip="Settlement Area",
        ).add_to(fg_p)
        fg_p.add_to(m)

    if sungai is not None:
        fg_s = folium.FeatureGroup(name="River Network", show=True)
        folium.GeoJson(sungai.to_crs(CRS_WGS).__geo_interface__,
                       style_function=lambda f: {"color": "#66AADD", "weight": 2.0, "opacity": 0.85}
                       ).add_to(fg_s)
        fg_s.add_to(m)

    if jalan is not None:
        fg_j = folium.FeatureGroup(name="Road Network", show=False)
        folium.GeoJson(jalan.to_crs(CRS_WGS).__geo_interface__,
                       style_function=lambda f: {"color": "#aaa", "weight": 1.2, "opacity": 0.7}
                       ).add_to(fg_j)
        fg_j.add_to(m)

    folium.Marker(
        location=[dem.peak_lat, dem.peak_lon],
        popup=folium.Popup(
            f"<b>{NAMA_GUNUNG}</b><br>Elevation: {dem.peak_elev:.0f} m asl<br>"
            f"Coordinates: {dem.peak_lat:.4f}S, {dem.peak_lon:.4f}E", max_width=220),
        tooltip=NAMA_GUNUNG,
        icon=folium.Icon(color="red", icon="fire", prefix="fa"),
    ).add_to(m)

    folium.LayerControl(position="topright", collapsed=False).add_to(m)
    m.get_root().html.add_child(folium.Element(_html_legenda()))
    m.get_root().html.add_child(folium.Element(_html_judul()))
    m.save(simpan_ke)
    print(f"[folium] Saved -> {simpan_ke}")
    return simpan_ke


def _hillshade_ke_png(dem: DEMData) -> str:
    from PIL import Image
    import io, base64
    hs   = (dem.hillshade * 255).astype(np.uint8)
    rgba = np.stack([hs, hs, hs,
                     (hs.astype(np.float32) * 0.85 + 40).clip(0, 220).astype(np.uint8)],
                    axis=-1)
    buf = io.BytesIO()
    Image.fromarray(rgba, mode="RGBA").save(buf, format="PNG")
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def _dem_gray_ke_png(dem: DEMData) -> str:
    from PIL import Image
    import io, base64
    dem_f = dem.dem.copy()
    dem_f[np.isnan(dem_f)] = np.nanmin(dem_f)
    vmin  = np.nanpercentile(dem_f, 2)
    vmax  = np.nanmax(dem_f)
    norm  = (dem_f - vmin) / max(vmax - vmin, 1)
    gray  = (norm * 180 + 30).clip(0, 255).astype(np.uint8)
    rgba  = np.stack([gray, gray, gray, np.full_like(gray, 200)], axis=-1)
    buf   = io.BytesIO()
    Image.fromarray(rgba, mode="RGBA").save(buf, format="PNG")
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def _html_legenda() -> str:
    return f"""
<div style="position:fixed;bottom:30px;left:20px;width:260px;
            background:#fff;border:1.5px solid #aaa;border-radius:8px;
            padding:12px 14px;font-family:Arial,sans-serif;font-size:12px;
            color:#222;z-index:9999;box-shadow:3px 3px 8px rgba(0,0,0,0.2);">
  <div style="font-weight:bold;font-size:13px;margin-bottom:8px;
              border-bottom:1px solid #ccc;padding-bottom:5px;">LEGEND</div>
  <div style="display:flex;align-items:center;margin-bottom:6px;">
    <div style="width:30px;height:14px;background:#FFE033;opacity:.9;
                border:1.5px solid #C8A800;margin-right:8px;flex-shrink:0;"></div>
    Cold Lahar Flood Potential</div>
  <div style="display:flex;align-items:center;margin-bottom:6px;">
    <div style="width:30px;height:14px;background:#5DADE2;opacity:.85;
                border:1.5px solid #1B4F72;margin-right:8px;flex-shrink:0;"></div>
    Settlement Area</div>
  <hr style="border:none;border-top:1px solid #ddd;margin:8px 0;">
  <div style="font-size:10px;color:#666;font-style:italic;">
    Source: PVMBG / Geological Agency<br>Ministry of Energy and Mineral Resources</div>
</div>"""


def _html_judul() -> str:
    return f"""
<div style="position:fixed;top:12px;left:50%;transform:translateX(-50%);
            background:#fff;border:1.5px solid #aaa;border-radius:8px;
            padding:8px 18px;font-family:Arial,sans-serif;font-size:13px;
            font-weight:bold;text-align:center;color:#111;z-index:9999;
            box-shadow:2px 2px 6px rgba(0,0,0,0.2);pointer-events:none;
            max-width:520px;">
  {NAMA_PETA}<br>
  <span style="font-size:12px;font-weight:normal;color:#444;">{NAMA_GUNUNG.upper()}</span>
</div>"""


# ═══════════════════════════════════════════════════════════════════════════
#  ENTRY POINT
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Cold Lahar Flood Hazard Potential Map")
    parser.add_argument("--dem",        default=FILE_DEM)
    parser.add_argument("--shp-sel",    default=FILE_SHP_SEL)
    parser.add_argument("--shp-bar",    default=FILE_SHP_BAR)
    parser.add_argument("--permukiman", default=FILE_PERMUKIMAN)   # NEW
    parser.add_argument("--nama",       default=NAMA_GUNUNG)
    parser.add_argument("--mode",       choices=["static", "interaktif", "keduanya"], default="keduanya")
    parser.add_argument("--out-static", default=None)
    parser.add_argument("--out-html",   default=None)
    args = parser.parse_args()

    FILE_DEM     = _resolve_path(args.dem)
    FILE_SHP_SEL = _resolve_path(args.shp_sel)
    FILE_SHP_BAR = _resolve_path(args.shp_bar)
    FILE_PERMUKIMAN = _resolve_path(args.permukiman)   # NEW
    NAMA_GUNUNG  = args.nama

    out_static = args.out_static or os.path.join(WORK_DIR, "peta_banjir_lahar_static.png")
    out_html   = args.out_html   or os.path.join(WORK_DIR, "peta_banjir_lahar.html")

    print("=" * 65)
    print("  COLD LAHAR FLOOD HAZARD POTENTIAL MAP")
    print("=" * 65)

    if not os.path.exists(FILE_DEM):
        print(f"[ERROR] DEM not found: {FILE_DEM}")
        sys.exit(1)
    if not os.path.exists(FILE_SHP_SEL):
        print(f"[ERROR] Southern shapefile not found: {FILE_SHP_SEL}")
        sys.exit(1)

    # FIXED: give an explicit warning if the west-southwest shapefile is not found
    if not os.path.exists(FILE_SHP_BAR):
        print(f"[WARNING] West-southwest shapefile NOT found at path: {FILE_SHP_BAR}")
        print(f"          Check the filename again (typo/different hyphen?) or the folder location.")

    shp_bar_resolved = FILE_SHP_BAR if os.path.exists(FILE_SHP_BAR) else None

    dem   = DEMData(FILE_DEM)
    lahar = LaharData(shp_selatan=FILE_SHP_SEL, dem=dem, shp_barat=shp_bar_resolved)
    sungai, jalan = load_vektor_tambahan()
    pemukiman = load_pemukiman(dem)   # NEW

    # FIXED: show a summary of which lahar sources were actually loaded
    print(f"\n[info] Lahar sources successfully loaded: {list(lahar.gdf_per_source.keys())}")

    if args.mode in ("static", "keduanya"):
        buat_peta_statis(dem, lahar, sungai, jalan, pemukiman=pemukiman, simpan_ke=out_static)

    if args.mode in ("interaktif", "keduanya"):
        buat_peta_interaktif(dem, lahar, sungai, jalan, pemukiman=pemukiman, simpan_ke=out_html)

    print("\n" + "=" * 65)
    print("  DONE")
    if args.mode in ("static", "keduanya"):
        print(f"  PNG  -> {out_static}")
    if args.mode in ("interaktif", "keduanya"):
        print(f"  HTML -> {out_html}")
    print("=" * 65)
