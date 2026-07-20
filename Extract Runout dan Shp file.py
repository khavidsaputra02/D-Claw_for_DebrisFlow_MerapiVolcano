import numpy as np, os, csv, re, sys
from scipy.spatial import KDTree

# ---------------------------------------------------------------------------
# Resolve path relative to the script directory — identical pattern to maketopo.py
# ---------------------------------------------------------------------------
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def _abspath(rel):
    return rel if os.path.isabs(rel) else os.path.join(_SCRIPT_DIR, rel)

OUTPUT_DIR   = _abspath("_output")
SETINPUT_DIR = _SCRIPT_DIR

# ---------------------------------------------------------------------------
# SOURCE_POINT — synchronized with src_cx / src_cy in setplot.py
# ---------------------------------------------------------------------------
SOURCE_POINT = (439164.988, 9165667.463)

H_THRESHOLD  = 0.01
CSV_OUTPUT   = os.path.join(SETINPUT_DIR, "hasil_sensitivitas.csv")

# UPDATED per latest Google Earth reading: Easting=439540.56, Northing=9161173.13
KALI_ADEM_COORD  = (439540.56, 9161173.13)
KALI_ADEM_RADIUS = 100.0

# Output shapefile CRS — UTM Zone 49S (matches simulation coordinates)
CRS_UTM = "EPSG:32749"

# ---------------------------------------------------------------------------
# MENIT_AKHIR (end time in minutes) — automatically synced with TFINAL from setrun.py
# ---------------------------------------------------------------------------
try:
    sys.path.insert(0, SETINPUT_DIR)
    from setrun import _TFINAL as _TFINAL_DETIK
    MENIT_AKHIR = _TFINAL_DETIK / 60.0
    print(f"[sensitivity] TFINAL = {_TFINAL_DETIK} s  "
          f"-> MENIT_AKHIR = {MENIT_AKHIR:.2f} min  (from setrun.py)")
except ImportError:
    MENIT_AKHIR = 27.0
    print(f"[sensitivity] WARNING: setrun.py not found or does not "
          f"export _TFINAL -- fallback MENIT_AKHIR = {MENIT_AKHIR} min")


# ===========================================================================
#  READ SIMULATION DATA
# ===========================================================================

def baca_dclaw_data(output_dir):
    params = {}
    f_path = os.path.join(output_dir, "dclaw.data")
    if not os.path.isfile(f_path):
        print("  [WARNING] dclaw.data not found")
        return params
    for line in open(f_path):
        line = line.strip()
        if not line or line.startswith("#") or "=:" not in line:
            continue
        v, n = line.split("=:", 1)
        params[n.strip().split("#")[0].strip()] = v.strip()
    return params


def baca_h0(setinput_dir):
    f_path = os.path.join(setinput_dir, "setinput.py")
    if not os.path.isfile(f_path):
        print("  [WARNING] setinput_1.py not found")
        return "?"
    for line in open(f_path):
        m = re.match(r'^h0\s*=\s*([0-9eE+\-\.]+)', line.strip())
        if m:
            return m.group(1)
    return "?"


def read_fgmax(f):
    data = np.loadtxt(f)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    print(f"  Cell count   : {len(data)}")
    print(f"  h_max range  : {data[:,2].min():.4f} - {data[:,2].max():.2f} m")
    return data[:,0], data[:,1], data[:,2]


def baca_fgout_frame(qfile):
    with open(qfile) as f:
        lines = f.readlines()
    mx   = int(lines[2].split()[0])
    my   = int(lines[3].split()[0])
    xlow = float(lines[4].split()[0])
    ylow = float(lines[5].split()[0])
    dx   = float(lines[6].split()[0])
    dy   = float(lines[7].split()[0])

    data = []
    for line in lines[8:]:
        vals = line.strip().split()
        if vals:
            data.append([float(v) for v in vals])
    data = np.array(data)
    h  = data[:, 0]
    xs = np.array([xlow + (i % mx + 0.5) * dx for i in range(mx * my)])
    ys = np.array([ylow + (i // mx + 0.5) * dy for i in range(mx * my)])
    return xs, ys, h[:mx * my], dx, dy


def cari_frame_terdekat(output_dir, target_detik):
    best_num  = None
    best_time = None
    best_diff = float('inf')
    for fname in sorted(os.listdir(output_dir)):
        if not fname.startswith('fgout0001.t'):
            continue
        suffix = fname[len('fgout0001.t'):]
        if not suffix.isdigit():
            continue
        try:
            t = float(open(os.path.join(output_dir, fname)).readline().split()[0])
        except Exception:
            continue
        diff = abs(t - target_detik)
        if diff < best_diff:
            best_diff = diff
            best_time = t
            best_num  = suffix
    return best_num, best_time, best_diff


def list_semua_frame(output_dir, t_min=0.0, t_max=None):
    """
    Collect ALL fgout0001.t* frames within the range [t_min, t_max] (seconds),
    sorted by ascending time. Used to aggregate (envelope) inundation across
    the whole simulation duration, not just a single last frame.

    Returns
    -------
    list of tuple (suffix, time_seconds), sorted ascending by time.
    """
    frames = []
    for fname in sorted(os.listdir(output_dir)):
        if not fname.startswith('fgout0001.t'):
            continue
        suffix = fname[len('fgout0001.t'):]
        if not suffix.isdigit():
            continue
        try:
            t = float(open(os.path.join(output_dir, fname)).readline().split()[0])
        except Exception:
            continue
        if t < t_min - 1e-6:
            continue
        if t_max is not None and t > t_max + 1e-6:
            continue
        frames.append((suffix, t))
    frames.sort(key=lambda item: item[1])
    return frames


# ===========================================================================
#  METRICS
# ===========================================================================

def hitung_luas(x, y, h):
    mask = h > H_THRESHOLD
    if not mask.any():
        print("  [WARNING] No inundated cells!")
        return 0.0
    xu = np.unique(np.round(x, 8))
    yu = np.unique(np.round(y, 8))
    dx = np.median(np.diff(xu)) if len(xu) > 1 else 1.0
    dy = np.median(np.diff(yu)) if len(yu) > 1 else 1.0
    if abs(x.mean()) < 360:
        lr = np.radians(y.mean())
        dx = dx * 111320 * np.cos(lr)
        dy = dy * 110540
        print("  Coordinates  : Geographic")
    else:
        print("  Coordinates  : UTM/Projected")
    print(f"  Resolution   : dx={dx:.2f} m  dy={dy:.2f} m")
    luas = mask.sum() * dx * dy
    print(f"  Inundated cells (h>{H_THRESHOLD}m) : {mask.sum()}")
    print(f"  Inundation area : {luas:.2f} m2  ({luas/1e6:.6f} km2)")
    return luas


def hitung_runout(x, y, h):
    """
    Straight-line (Euclidean) runout from SOURCE_POINT to the farthest
    inundated cell. SOURCE_POINT = src_cx/src_cy from setplot.py
    (439164.988, 9165667.463).
    """
    mask = h > H_THRESHOLD
    if not mask.any():
        print("  [WARNING] No inundated cells!")
        return 0.0
    xw, yw = x[mask], y[mask]
    sx, sy = SOURCE_POINT
    if abs(x.mean()) < 360:
        R  = 6371000
        p1 = np.radians(sy); p2 = np.radians(yw)
        dp = np.radians(yw - sy); dl = np.radians(xw - sx)
        a  = np.sin(dp/2)**2 + np.cos(p1)*np.cos(p2)*np.sin(dl/2)**2
        d  = R * 2 * np.arcsin(np.sqrt(a))
    else:
        d = np.sqrt((xw - sx)**2 + (yw - sy)**2)
    runout = d.max()
    print(f"  Runout distance : {runout:.2f} m  ({runout/1000:.4f} km)")
    return runout


def ekstrak_ketebalan_fgout(output_dir, target_menit,
                             target_coord=KALI_ADEM_COORD):
    target_detik = target_menit * 60.0
    frame_num, frame_time, selisih = cari_frame_terdekat(output_dir, target_detik)

    if frame_num is None:
        print("  [ERROR] No fgout file found!")
        return None

    qfile = os.path.join(output_dir, f'fgout0001.q{frame_num}')
    print(f"  Frame used      : fgout0001.q{frame_num}")
    print(f"  Frame time      : {frame_time:.1f} s  ({frame_time/60:.2f} min)")
    print(f"  Target diff     : {selisih:.1f} s")

    xs, ys, h, dx, dy = baca_fgout_frame(qfile)
    tx, ty = target_coord
    tree = KDTree(np.column_stack([xs, ys]))
    dist, idx = tree.query([tx, ty])

    # Statistics within radius
    jarak_semua  = np.sqrt((xs - tx)**2 + (ys - ty)**2)
    mask_radius  = jarak_semua <= KALI_ADEM_RADIUS
    mask_wet_rad = mask_radius & (h > H_THRESHOLD)

    hasil = {
        "frame"        : frame_num,
        "waktu_s"      : frame_time,
        "waktu_menit"  : frame_time / 60.0,
        "x_cell"       : float(xs[idx]),
        "y_cell"       : float(ys[idx]),
        "jarak_m"      : float(dist),
        "h"            : float(h[idx]),
        "tergenang"    : bool(h[idx] > H_THRESHOLD),
        "n_cell_radius": int(mask_radius.sum()),
        "n_cell_wet"   : int(mask_wet_rad.sum()),
        "h_mean_radius": float(h[mask_wet_rad].mean()) if mask_wet_rad.any() else 0.0,
        "h_max_radius" : float(h[mask_wet_rad].max())  if mask_wet_rad.any() else 0.0,
        "h_min_radius" : float(h[mask_wet_rad].min())  if mask_wet_rad.any() else 0.0,
    }

    print(f"\n  === KALI ADEM — t={frame_time:.0f}s ({frame_time/60:.2f} min) ===")
    print(f"  Target coordinate : ({tx:.2f}, {ty:.2f})")
    print(f"  Nearest cell      : ({hasil['x_cell']:.2f}, {hasil['y_cell']:.2f})")
    print(f"  Distance to target: {hasil['jarak_m']:.2f} m")
    print(f"  ─────────────────────────────────────────")
    if hasil['tergenang']:
        print(f"  Lahar thickness   : {hasil['h']:.4f} m  ({hasil['h']*100:.2f} cm)")
    else:
        print(f"  Lahar thickness   : {hasil['h']:.4f} m  (DRY / h < {H_THRESHOLD} m)")
    print(f"  ─────────────────────────────────────────")
    print(f"  Statistics within {KALI_ADEM_RADIUS:.0f} m radius:")
    print(f"    Total cells      : {hasil['n_cell_radius']}")
    print(f"    Inundated cells  : {hasil['n_cell_wet']}")
    if hasil['n_cell_wet'] > 0:
        print(f"    Mean h           : {hasil['h_mean_radius']:.4f} m")
        print(f"    Max h            : {hasil['h_max_radius']:.4f} m")
        print(f"    Min h            : {hasil['h_min_radius']:.4f} m")
    else:
        print("    No inundated cells within this radius.")
    print(f"  ─────────────────────────────────────────")
    return hasil


# ===========================================================================
#  DEPOSIT EXPORT (FULL-SIMULATION ENVELOPE) -> SHAPEFILE
# ===========================================================================

def ekspor_endapan_shp_envelope(output_dir, t_awal_menit=0.0, t_akhir_menit=None,
                                 shp_dir=None, shp_name="Laharselatan"):
    """
    Builds a MAXIMUM INUNDATION (envelope) shapefile by reading ALL fgout
    frames from t_awal_menit to t_akhir_menit (default: 0 to MENIT_AKHIR),
    then for each grid cell taking the MAXIMUM h that ever occurred across
    that time range — not just from a single last frame. This represents
    the entire area ever traversed/inundated by the lahar during the
    simulation from 0 to t_akhir_menit.

    Attributes per polygon (one polygon = one grid cell):
      - h_max_m    : MAXIMUM lahar thickness ever recorded in this cell
                     over the simulation (m)
      - h_max_cm   : same, in cm
      - t_max_s    : time (seconds) at which that maximum h occurred
      - t_max_min  : time (minutes) at which that maximum h occurred
      - t_wet_s    : time (seconds) this cell first became inundated
                     (h > H_THRESHOLD)
      - t_wet_min  : same, in minutes
      - x_ctr, y_ctr : cell center coordinates (m)
      - n_frame    : number of fgout frames processed for this envelope

    Parameters
    ----------
    output_dir    : simulation _output directory
    t_awal_menit  : start of time range (minutes), default 0
    t_akhir_menit : end of time range (minutes), default MENIT_AKHIR
    shp_dir       : directory to save the .shp (default: SETINPUT_DIR)
    shp_name      : filename without extension (default: Laharselatan)

    Returns
    -------
    path to the saved .shp file, or None on failure.
    """
    try:
        import geopandas as gpd
        from shapely.geometry import box
    except ImportError:
        print("  [ERROR] geopandas / shapely not installed.")
        print("          Run: pip install geopandas shapely")
        return None

    if shp_dir is None:
        shp_dir = SETINPUT_DIR
    if t_akhir_menit is None:
        t_akhir_menit = MENIT_AKHIR

    t_min_detik = t_awal_menit * 60.0
    t_max_detik = t_akhir_menit * 60.0

    frames = list_semua_frame(output_dir, t_min=t_min_detik, t_max=t_max_detik)
    if not frames:
        print(f"  [ERROR] No fgout frames in the range "
              f"{t_awal_menit:.1f}-{t_akhir_menit:.1f} min!")
        return None

    print(f"  Time range      : {t_awal_menit:.1f} - {t_akhir_menit:.1f} min "
          f"({t_min_detik:.0f} - {t_max_detik:.0f} s)")
    print(f"  Frame count     : {len(frames)}")
    print(f"  First frame     : t={frames[0][1]:.1f} s")
    print(f"  Last frame      : t={frames[-1][1]:.1f} s")

    xs = ys = dx = dy = None
    h_max      = None
    t_of_hmax  = None
    t_wet_first = None

    for i, (suffix, frame_time) in enumerate(frames):
        qfile = os.path.join(output_dir, f'fgout0001.q{suffix}')
        if not os.path.isfile(qfile):
            print(f"  [WARNING] File not found, skipped: {qfile}")
            continue

        xs_i, ys_i, h_i, dx_i, dy_i = baca_fgout_frame(qfile)

        if h_max is None:
            xs, ys, dx, dy = xs_i, ys_i, dx_i, dy_i
            h_max       = np.zeros_like(h_i)
            t_of_hmax   = np.zeros_like(h_i)
            t_wet_first = np.full_like(h_i, np.nan)
        elif len(h_i) != len(h_max):
            print(f"  [WARNING] Grid size for frame t={frame_time:.1f}s does "
                  f"not match, skipped (fgout grid may have changed).")
            continue

        # Update the maximum h + the time it occurred
        mask_baru = h_i > h_max
        h_max[mask_baru]     = h_i[mask_baru]
        t_of_hmax[mask_baru] = frame_time

        # Record the first time this cell became inundated (h > H_THRESHOLD)
        mask_wet_now  = h_i > H_THRESHOLD
        mask_belum_wet = np.isnan(t_wet_first)
        mask_set_wet  = mask_wet_now & mask_belum_wet
        t_wet_first[mask_set_wet] = frame_time

        if (i + 1) % 20 == 0 or (i + 1) == len(frames):
            print(f"    ... processed {i+1}/{len(frames)} frames "
                  f"(t={frame_time:.1f}s)")

    if h_max is None:
        print("  [ERROR] No frames were successfully read.")
        return None

    mask_wet = h_max > H_THRESHOLD
    n_wet    = mask_wet.sum()

    if n_wet == 0:
        print("  [WARNING] No inundated cells across this time range "
              "-- shapefile not created.")
        return None

    print(f"  Wet cells (h_max>{H_THRESHOLD} m) : {n_wet:,}  "
          f"(envelope of {len(frames)} frames)")

    half_dx = dx / 2.0
    half_dy = dy / 2.0

    xw          = xs[mask_wet]
    yw          = ys[mask_wet]
    hw          = h_max[mask_wet]
    tw_max      = t_of_hmax[mask_wet]
    tw_wet      = t_wet_first[mask_wet]

    geoms = [box(xi - half_dx, yi - half_dy, xi + half_dx, yi + half_dy)
             for xi, yi in zip(xw, yw)]

    gdf = gpd.GeoDataFrame(
        {
            "h_max_m"  : np.round(hw, 4),
            "h_max_cm" : np.round(hw * 100.0, 2),
            "t_max_s"  : np.round(tw_max, 1),
            "t_max_min": np.round(tw_max / 60.0, 3),
            "t_wet_s"  : np.round(tw_wet, 1),
            "t_wet_min": np.round(tw_wet / 60.0, 3),
            "x_ctr"    : np.round(xw, 3),
            "y_ctr"    : np.round(yw, 3),
            "n_frame"  : len(frames),
        },
        geometry=geoms,
        crs=CRS_UTM,
    )

    # Optionally dissolve if there are very many cells -- can be uncommented
    # to keep per-cell attributes intact
    # gdf = gdf.dissolve()

    os.makedirs(shp_dir, exist_ok=True)
    shp_path = os.path.join(shp_dir, f"{shp_name}.shp")
    gdf.to_file(shp_path, driver="ESRI Shapefile", encoding="utf-8")

    print(f"  Shapefile saved  : {shp_path}")
    print(f"  Polygon count    : {len(gdf):,}")
    print(f"  CRS              : {CRS_UTM}")
    print(f"  Attribute columns: {list(gdf.columns.drop('geometry'))}")
    return shp_path


# ===========================================================================
#  SAVE CSV
# ===========================================================================

def simpan_csv(params, luas, runout_lurus, h_fgout=None):
    kolom  = ["m0", "h0", "m_crit", "kref", "mu", "alpha_c", "rho_s", "rho_f", "phi"]
    header = kolom + [
        "luas_genangan_m2",
        "runout_lurus_m",
        "h_kali_adem_t27_m",
    ]
    ada = os.path.isfile(CSV_OUTPUT)
    with open(CSV_OUTPUT, "a", newline="") as f:
        w = csv.writer(f)
        if not ada:
            w.writerow(header)
        w.writerow([params.get(k, "") for k in kolom] + [
            round(luas,         2),
            round(runout_lurus, 2),
            round(h_fgout, 4) if h_fgout is not None else "",
        ])
    print(f"  CSV saved     : {CSV_OUTPUT}")


# ===========================================================================
#  MAIN
# ===========================================================================
if __name__ == "__main__":
    print("=" * 60)
    print("  D-Claw Sensitivity Metrics Extractor")
    print("  (+ Kali Adem extraction at t=27 min)")
    print("  (+ Deposit envelope export 0 - MENIT_AKHIR -> Laharselatan.shp)")
    print("=" * 60)
    print(f"  SOURCE_POINT : {SOURCE_POINT}  (synced with setplot.py)")

    print("\n[0] Reading simulation parameters ...")
    params = baca_dclaw_data(OUTPUT_DIR)
    params["h0"] = baca_h0(SETINPUT_DIR)

    info = [
        ("m0",      "Initial Solid Volume Fraction"),
        ("h0",      "Initial Lahar Depth (m)"),
        ("m_crit",  "Critical-State Volume Fraction"),
        ("kref",    "Hydraulic Permeability"),
        ("mu",      "Pore Fluid Viscosity (Pa.s)"),
        ("alpha_c", "Initial Compressibility"),
        ("rho_s",   "Solid Grain Density (kg/m3)"),
        ("rho_f",   "Pore Fluid Density (kg/m3)"),
        ("phi",     "Coulomb Friction Angle (deg)"),
    ]
    print(f"\n  {'Code':<10}  {'Value':<14}  Description")
    print(f"  {'-'*10}  {'-'*14}  {'-'*35}")
    for k, ket in info:
        print(f"  {k:<10}  {params.get(k,'?'):<14}  {ket}")

    print("\n[1] Reading fgmax ...")
    x, y, h = read_fgmax(os.path.join(OUTPUT_DIR, "fgmax0001.txt"))

    print("\n[2] Computing inundation area ...")
    luas = hitung_luas(x, y, h)

    print("\n[3] Computing runout ...")
    runout_lurus = hitung_runout(x, y, h)

    print(f"\n[4] Kali Adem thickness at MINUTE {MENIT_AKHIR:.0f} (from fgout) ...")
    hasil_fgout = ekstrak_ketebalan_fgout(OUTPUT_DIR, MENIT_AKHIR)

    print(f"\n[5] Exporting ENVELOPE deposit (0 to minute {MENIT_AKHIR:.0f}) "
          f"-> Laharselatan.shp ...")
    shp_path = ekspor_endapan_shp_envelope(
        output_dir    = OUTPUT_DIR,
        t_awal_menit  = 0.0,
        t_akhir_menit = MENIT_AKHIR,
        shp_dir       = SETINPUT_DIR,
        shp_name      = "Laharselatan",
    )

    print("\n[6] Saving to CSV ...")
    simpan_csv(params, luas, runout_lurus,
               h_fgout=hasil_fgout["h"] if hasil_fgout else None)

    print("\n" + "=" * 60)
    print("  FINAL RESULTS")
    print("=" * 60)
    for k, ket in info:
        print(f"  {k:<10} = {params.get(k,'?')}")
    print(f"  {'-'*55}")
    print(f"  Inundation area            : {luas:.2f} m2  ({luas/1e6:.6f} km2)")
    print(f"  Runout                     : {runout_lurus:.2f} m  "
          f"({runout_lurus/1000:.4f} km)")
    if hasil_fgout:
        print(f"  Kali Adem thickness t={MENIT_AKHIR:.0f}m  : "
              f"{hasil_fgout['h']:.4f} m  "
              f"({'INUNDATED' if hasil_fgout['tergenang'] else 'DRY'})")
    if shp_path:
        print(f"  Deposit shapefile (envelope 0-{MENIT_AKHIR:.0f} min) : {shp_path}")
    print("=" * 60)
