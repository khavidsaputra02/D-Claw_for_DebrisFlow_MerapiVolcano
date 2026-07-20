"""
setrun.py - Clawpack run-time parameters for the Merapi lahar simulation.

Pipeline: define domain from DEM metadata -> set time/timestepping/AMR
parameters -> define source flag region and gauges -> configure GeoClaw,
topo/init files, fgmax, fgout, and D-Claw material parameters.
"""

import os
import sys
import numpy as np

try:
    CLAW = os.environ["CLAW"]
except:
    raise Exception("*** Must first set CLAW environment variable")

from clawpack.amrclaw.data import FlagRegion
from clawpack.geoclaw import fgout_tools, fgmax_tools

bouss = False
if bouss:
    i_eta  = 10; i_hm = 6; num_eqn = 9
else:
    i_eta  = 8;  i_hm = 4; num_eqn = 7

# ---------------------------------------------------------------------------
# DEM metadata
# ---------------------------------------------------------------------------
_ncols    = 2083
_nrows    = 2111
_xllcorner = 427606.288
_yllcorner = 9153770.505
_cellsize  = 8.289958903244
_half      = _cellsize / 2.0

dom_xlower = _xllcorner + _half
dom_xupper = dom_xlower + (_ncols - 1) * _cellsize
dom_ylower = _yllcorner + _half
dom_yupper = dom_ylower + (_nrows - 1) * _cellsize

# ---------------------------------------------------------------------------
# SOURCE bounding box — Source Selatan_1.shp
# (taken directly from the shapefile bounds, EPSG:32749)
# ---------------------------------------------------------------------------
_src_xmin = 438885.574;  _src_xmax = 439444.402   # Source Selatan_1.shp
_src_ymin = 9165049.500; _src_ymax = 9166285.426
_src_cx   = (_src_xmin + _src_xmax) / 2.0         # 439164.988
_src_cy   = (_src_ymin + _src_ymax) / 2.0         # 9165667.463

_FLAG_BUF = 50.0
_TFINAL   = 1620.0
_NUM_OUTPUT_TIMES = int(round(_TFINAL / 10.0))   # 162

_AMR_FACTOR  = 2 * 2 * 2
_dx_level1   = _cellsize * _AMR_FACTOR
_num_cells_x = int(round((dom_xupper - dom_xlower) / _dx_level1))  # 260
_num_cells_y = int(round((dom_yupper - dom_ylower) / _dx_level1))  # 264
assert _num_cells_x > 0 and _num_cells_y > 0


def setrun(claw_pkg="dclaw"):
    from clawpack.clawutil import data
    assert claw_pkg.lower() == "dclaw"

    rundata = data.ClawRunData(claw_pkg, 2)
    clawdata = rundata.clawdata

    # -----------------------------------------------------------------------
    # Domain
    # -----------------------------------------------------------------------
    clawdata.num_dim  = 2
    clawdata.lower[0] = dom_xlower
    clawdata.upper[0] = dom_xupper
    clawdata.lower[1] = dom_ylower
    clawdata.upper[1] = dom_yupper
    clawdata.num_cells[0] = _num_cells_x
    clawdata.num_cells[1] = _num_cells_y
    clawdata.num_eqn    = num_eqn
    clawdata.num_aux    = 10
    clawdata.capa_index = 0

    # -----------------------------------------------------------------------
    # Time
    # -----------------------------------------------------------------------
    clawdata.t0      = 0.0
    clawdata.restart = False

    clawdata.output_style     = 1
    clawdata.num_output_times = _NUM_OUTPUT_TIMES   # 162 frames @ 10s interval
    clawdata.tfinal           = _TFINAL
    clawdata.output_t0        = True

    clawdata.output_format         = "ascii"
    clawdata.output_q_components   = "all"
    clawdata.output_aux_components = "none"
    clawdata.output_aux_onlyonce   = True

    # -----------------------------------------------------------------------
    # Timestepping — STABILITY FIX
    # -----------------------------------------------------------------------
    clawdata.verbosity   = 1
    clawdata.dt_variable = True
    clawdata.dt_initial  = 1.0e-6    # FIX: was 1e-4  -> start small
    clawdata.dt_max      = 30.0      # FIX: was 1e99  -> cap at 30 seconds
    clawdata.cfl_desired = 0.25      # FIX: was 0.40  -> conservative
    clawdata.cfl_max     = 0.40      # FIX: was 0.50
    clawdata.steps_max   = 100000    # FIX: was 5000  -> enough for 1620 s

    # -----------------------------------------------------------------------
    # Numerics
    # -----------------------------------------------------------------------
    clawdata.order             = 2
    clawdata.dimensional_split = "unsplit"
    clawdata.transverse_waves  = 2
    clawdata.num_waves         = 5
    clawdata.limiter           = [4, 4, 4, 4, 4]
    clawdata.use_fwaves        = True
    clawdata.source_split      = "godunov"
    clawdata.num_ghost         = 2
    clawdata.bc_lower[0] = "extrap"; clawdata.bc_upper[0] = "extrap"
    clawdata.bc_lower[1] = "extrap"; clawdata.bc_upper[1] = "extrap"
    clawdata.checkpt_style = 0

    # -----------------------------------------------------------------------
    # AMR — STABILITY FIX
    # -----------------------------------------------------------------------
    amrdata = rundata.amrdata
    amrdata.amr_levels_max      = 3
    amrdata.refinement_ratios_x = [2, 2, 2]
    amrdata.refinement_ratios_y = [2, 2, 2]
    # FIX: ratio_t [2,4,4] -> finer levels get more sub-steps
    amrdata.refinement_ratios_t = [2, 4, 4]   # was [2, 2, 2]

    amrdata.aux_type = ["center","center","yleft","center","center",
                        "center","center","center","center","center"]
    amrdata.flag_richardson     = False
    amrdata.flag2refine         = True
    amrdata.regrid_interval     = 1    # FIX: was 3 -> lahar moves fast
    amrdata.regrid_buffer_width = 2
    amrdata.clustering_cutoff   = 0.700000
    amrdata.verbosity_regrid    = 0

    # -----------------------------------------------------------------------
    # Flag regions — Source Selatan_1
    # -----------------------------------------------------------------------
    flagregions = rundata.flagregiondata.flagregions
    fr = FlagRegion(num_dim=2)
    fr.name     = "Source_Selatan_1"
    fr.minlevel = 3
    fr.maxlevel = 3
    fr.t1       = 0.0
    fr.t2       = 1e10
    fr.spatial_region_type = 1
    fr.spatial_region = [
        _src_xmin - _FLAG_BUF, _src_xmax + _FLAG_BUF,
        _src_ymin - _FLAG_BUF, _src_ymax + _FLAG_BUF,
    ]
    flagregions.append(fr)

    # -----------------------------------------------------------------------
    # Gauges — placed along the flow path from the source
    # -----------------------------------------------------------------------
    rundata.gaugedata.gauges = []
    for gid, dy in enumerate([-100, -1000, -3000, -5000], start=1):
        rundata.gaugedata.gauges.append([gid, _src_cx, _src_ymin + dy, 0.0, _TFINAL])

    # -----------------------------------------------------------------------
    # GeoClaw
    # -----------------------------------------------------------------------
    geo = rundata.geo_data
    geo.gravity            = 9.81
    geo.coordinate_system = 1
    geo.earth_radius       = 6367.5e3
    geo.coriolis_forcing   = False
    geo.sea_level          = 0.0
    geo.dry_tolerance      = 1.0e-3
    geo.friction_forcing   = True
    geo.manning_coefficient = 0.025
    geo.friction_depth     = 1e6

    rundata.refinement_data.variable_dt_refinement_ratios = True
    rundata.refinement_data.wave_tolerance = 0.01

    # -----------------------------------------------------------------------
    # Topo & init
    # -----------------------------------------------------------------------
    rundata.topo_data.topofiles.append([3, "basal_topo.tt3"])
    rundata.topo_data.topo_missing = 0.0
    rundata.dtopo_data

    rundata.qinitdclaw_data.qinitfiles.append([3, i_eta, "surface_topo.tt3"])
    rundata.qinitdclaw_data.qinitfiles.append([3, i_hm,  "mass_frac.tt3"])

    # -----------------------------------------------------------------------
    # fgmax
    # -----------------------------------------------------------------------
    rundata.fgmax_data.num_fgmax_val = 5
    fg = fgmax_tools.FGmaxGrid()
    fg.point_style = 2
    dx_fg = 25.0
    fg.x1 = dom_xlower + dx_fg/2; fg.x2 = dom_xupper - dx_fg/2
    fg.y1 = dom_ylower + dx_fg/2; fg.y2 = dom_yupper - dx_fg/2
    fg.dx = dx_fg
    fg.tstart_max = 0.0; fg.tend_max = 1e10
    fg.dt_check = 5.0; fg.min_level_check = 2
    fg.arrival_tol = 1e-2; fg.interp_method = 0
    rundata.fgmax_data.fgmax_grids.append(fg)

    # -----------------------------------------------------------------------
    # fgout
    # -----------------------------------------------------------------------
    fgout = fgout_tools.FGoutGrid()
    fgout.fgno = 1; fgout.point_style = 2; fgout.output_format = "ascii"
    fgout.nx = _ncols; fgout.ny = _nrows
    fgout.x1 = dom_xlower; fgout.x2 = dom_xupper
    fgout.y1 = dom_ylower; fgout.y2 = dom_yupper
    fgout.tstart = 0.0; fgout.tend = _TFINAL
    fgout.nout = _NUM_OUTPUT_TIMES + 1   # 163
    fgout.q_out_vars = [1, 4, 8]
    rundata.fgout_data.fgout_grids.append(fgout)

    # -----------------------------------------------------------------------
    # D-Claw material
    # -----------------------------------------------------------------------
    dc = rundata.dclaw_data
    dc.rho_s = 2500.0; dc.rho_f = 1150.0
    dc.m0    = 0.59;   dc.m_crit = 0.618
    dc.mref  = 0.6;    dc.kref   = 7.0e-11
    dc.phi   = 35;     dc.mu     = 0.005
    dc.alpha_c = 1.0e-7
    dc.delta   = 0.001; dc.c1 = 1
    dc.sigma_0 = 1.0e3
    dc.src2method  = 2; dc.alphamethod = 1
    dc.segregation = 0; dc.beta_seg = 0.0
    dc.chi0 = 0.5; dc.chie = 0.5
    dc.bed_normal = 0; dc.theta_input = 0.0
    dc.entrainment = 0; dc.entrainment_rate = 0.0
    dc.entrainment_method = 1; dc.me = 0.6

    rundata.pinitdclaw_data.init_ptype = 0

    rundata.flowgrades_data.flowgrades = [
        [1.0e-2, 1, 1, 2],
        [1.0e-2, 2, 1, 2],
    ]

    amrdata.dprint = False; amrdata.eprint = False; amrdata.edebug = False
    amrdata.gprint = False; amrdata.nprint = False; amrdata.pprint = False
    amrdata.rprint = False; amrdata.sprint = False; amrdata.tprint = False
    amrdata.uprint = False; amrdata.max1d = 300

    return rundata


if __name__ == "__main__":
    rundata = setrun(*sys.argv[1:])
    rundata.write()
