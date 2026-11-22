"""
Build the cranked-kite UAV in OpenVSP via the Python API.

Requires OpenVSP's Python bindings on the path:
    python3 -c "import openvsp"   ->  must not error

If that import fails, install the bindings that ship with your OpenVSP build:
    pip install <OpenVSP install dir>/python/openvsp
or add the OpenVSP python directory to PYTHONPATH.

Usage:
    python3 build_model.py                 # writes cranked_kite_v1.vsp3
    python3 build_model.py --no-masses     # geometry only
    python3 build_model.py --list-parms    # dump parameter names for debugging

Parameter names in the OpenVSP API have changed across releases. If a
SetParmVal call raises, run with --list-parms and match the names against
what your build reports.
"""

import argparse
import os
import sys

import geometry as G

try:
    import openvsp as vsp
except ImportError:
    sys.exit(
        "openvsp python module not found.\n"
        "Add the OpenVSP python directory to PYTHONPATH, or pip install it\n"
        "from your OpenVSP installation's python/ folder."
    )


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def setp(geom_id, group, name, value):
    """Set a parameter and fail loudly with a useful message."""
    pid = vsp.FindParm(geom_id, name, group)
    if not pid:
        raise RuntimeError(
            f"parameter '{name}' not found in group '{group}'. "
            f"Run with --list-parms to see what this OpenVSP build exposes."
        )
    vsp.SetParmVal(pid, value)
    return pid


def list_parms(geom_id):
    """Print every parameter on a geom, grouped. Debugging aid."""
    groups = {}
    for pid in vsp.GetGeomParmIDs(geom_id):
        g = vsp.GetParmGroupName(pid)
        groups.setdefault(g, []).append(vsp.GetParmName(pid))
    for g in sorted(groups):
        print(f"  [{g}]")
        for n in sorted(groups[g]):
            print(f"      {n}")


def set_airfoil(wing_id, xsec_index, thickness, airfoil_path):
    """
    Apply the reflex airfoil at one station.

    Tries the .dat file first. If it is missing, falls back to a NACA
    four-series with reflex-like camber so the model still builds - but a
    four-series is NOT a reflex section and will not trim a tailless wing.
    Get the real .dat file before you trust any VSPAERO pitching moment.
    """
    xsurf = vsp.GetXSecSurf(wing_id, 0)

    if airfoil_path and os.path.isfile(airfoil_path):
        vsp.ChangeXSecShape(xsurf, xsec_index, vsp.XS_FILE_AIRFOIL)
        vsp.Update()
        xsec = vsp.GetXSec(xsurf, xsec_index)
        vsp.ReadFileAirfoil(xsec, airfoil_path)
        vsp.Update()
        # scale the imported section to the requested thickness ratio
        tc = vsp.GetXSecParm(xsec, "ThickChord")
        if tc:
            vsp.SetParmVal(tc, thickness)
        return "file"

    vsp.ChangeXSecShape(xsurf, xsec_index, vsp.XS_FOUR_SERIES)
    vsp.Update()
    xsec = vsp.GetXSec(xsurf, xsec_index)
    for pname, pval in (("ThickChord", thickness),
                        ("Camber", 0.02),
                        ("CamberLoc", 0.40)):
        pid = vsp.GetXSecParm(xsec, pname)
        if pid:
            vsp.SetParmVal(pid, pval)
    vsp.Update()
    return "fallback"


# ---------------------------------------------------------------------------
# geometry construction
# ---------------------------------------------------------------------------

def build_wing(cfg, props, airfoil_path):
    """
    Three spanwise sections on one wing geom:
        XSec_1  inboard panel   root  -> crank
        XSec_2  outboard panel  crank -> tip
        XSec_3  winglet         tip   -> winglet tip
    OpenVSP mirrors about the XZ plane automatically, so everything below is
    defined for the right side only.
    """
    wid = vsp.AddGeom("WING")
    vsp.SetGeomName(wid, "MainWing")

    # A fresh WING geom has one section (XSec_1). Insert two more.
    vsp.InsertXSec(wid, 1, vsp.XS_FOUR_SERIES)
    vsp.InsertXSec(wid, 2, vsp.XS_FOUR_SERIES)
    vsp.Update()

    sections = [
        # (index, span, root chord, tip chord, LE sweep, dihedral)
        (1, cfg["p1_span"], cfg["p1_root_chord"], cfg["p1_tip_chord"],
         cfg["p1_sweep_le"], cfg["p1_dihedral"]),
        (2, cfg["p2_span"], cfg["p2_root_chord"], cfg["p2_tip_chord"],
         cfg["p2_sweep_le"], cfg["p2_dihedral"]),
        (3, cfg["wl_span"], cfg["wl_root_chord"], cfg["wl_tip_chord"],
         cfg["wl_sweep_le"], cfg["wl_dihedral"]),
    ]

    for idx, span, cr, ct, sweep, dihed in sections:
        grp = f"XSec_{idx}"
        # Drive the section by span and the two chords, so our numbers are
        # taken literally instead of being back-computed from area or AR.
        vsp.SetDriverGroup(wid, idx,
                           vsp.SPAN_WSECT_DRIVER,
                           vsp.ROOTC_WSECT_DRIVER,
                           vsp.TIPC_WSECT_DRIVER)
        setp(wid, grp, "Span", span)
        setp(wid, grp, "Root_Chord", cr)
        setp(wid, grp, "Tip_Chord", ct)
        setp(wid, grp, "Sweep", sweep)
        setp(wid, grp, "Sweep_Location", 0.0)   # 0.0 = leading edge
        setp(wid, grp, "Dihedral", dihed)
        setp(wid, grp, "SectTess_U", 12)
        vsp.Update()

    # ---- twist: washout builds up outboard -------------------------------
    # Twist is applied at the outboard end of each section, about the
    # quarter-chord so the leading edge does not swing fore and aft.
    setp(wid, "XSec_1", "Twist", cfg["p1_twist_tip"])
    setp(wid, "XSec_1", "Twist_Location", 0.25)
    setp(wid, "XSec_2", "Twist", cfg["p2_twist_tip"])
    setp(wid, "XSec_2", "Twist_Location", 0.25)
    setp(wid, "XSec_3", "Twist", 0.0)
    vsp.Update()

    # ---- airfoils at the four stations -----------------------------------
    modes = []
    for xsec_index, tc in ((0, cfg["tc_root"]),
                           (1, cfg["tc_crank"]),
                           (2, cfg["tc_tip"]),
                           (3, cfg["tc_tip"])):
        modes.append(set_airfoil(wid, xsec_index, tc, airfoil_path))
    vsp.Update()

    # ---- surface tessellation for a clean VLM mesh -----------------------
    for pname, pval in (("Tess_W", 33), ("LECluster", 0.25), ("TECluster", 0.8)):
        pid = vsp.FindParm(wid, pname, "WingGeom")
        if pid:
            vsp.SetParmVal(pid, pval)
    vsp.Update()

    return wid, modes


def add_control_surfaces(wid, cfg):
    """
    Elevons on the outer wing, drag rudders on the winglets.

    UStart and UEnd are normalised surface coordinates running root to tip
    across ALL sections, so the mapping from half-span fraction is not exactly
    linear. The values below are computed from cumulative section span, which
    is close enough for a conceptual model - check them in the GUI and nudge
    if the highlighted region does not match what you expect.
    """
    b1, b2, bw = cfg["p1_span"], cfg["p2_span"], cfg["wl_span"]
    total = b1 + b2 + bw

    def u_of_y(y):
        return y / total

    created = []

    # ---- elevons ---------------------------------------------------------
    y0 = cfg["elevon_span_start"] * (b1 + b2)
    y1 = cfg["elevon_span_end"] * (b1 + b2)
    ss = vsp.AddSubSurf(wid, vsp.SS_CONTROL)
    vsp.SetSubSurfName(wid, ss, "Elevon")
    for name, val in (("UStart", u_of_y(y0)),
                      ("UEnd", u_of_y(y1)),
                      ("Length_C_Start", 1.0 - cfg["elevon_hinge_xc"]),
                      ("Length_C_End", 1.0 - cfg["elevon_hinge_xc"]),
                      ("SE_Const_Flag", 1.0),
                      ("Surf_Type", vsp.SS_CONTROL_BOTH_SURF)):
        pid = vsp.FindParm(ss, name, "SS_Control")
        if pid:
            vsp.SetParmVal(pid, val)
    created.append("Elevon")

    # ---- drag rudders on the winglet -------------------------------------
    ss2 = vsp.AddSubSurf(wid, vsp.SS_CONTROL)
    vsp.SetSubSurfName(wid, ss2, "DragRudder")
    for name, val in (("UStart", u_of_y(b1 + b2 + 0.08 * bw)),
                      ("UEnd", u_of_y(total - 0.02 * bw)),
                      ("Length_C_Start", 1.0 - cfg["rudder_hinge_xc"]),
                      ("Length_C_End", 1.0 - cfg["rudder_hinge_xc"]),
                      ("SE_Const_Flag", 1.0),
                      ("Surf_Type", vsp.SS_CONTROL_BOTH_SURF)):
        pid = vsp.FindParm(ss2, name, "SS_Control")
        if pid:
            vsp.SetParmVal(pid, val)
    created.append("DragRudder")

    vsp.Update()
    return created


def add_point_masses(cfg, props):
    """
    Blank geoms carrying point masses, for the Mass Properties analysis.

    Structure mass is NOT added here - assign it through the wing's material
    density in the GUI, or accept the lumped value from geometry.py.
    """
    ids = []
    for name, mass, x in cfg["point_masses"]:
        bid = vsp.AddGeom("BLANK")
        vsp.SetGeomName(bid, name)
        for grp, pname, val in (("XForm", "X_Rel_Location", x),
                                ("XForm", "Y_Rel_Location", 0.0),
                                ("XForm", "Z_Rel_Location", 0.0),
                                ("Mass", "Point_Mass", 1.0),
                                ("Mass", "Mass", mass)):
            pid = vsp.FindParm(bid, pname, grp)
            if pid:
                vsp.SetParmVal(pid, val)
        vsp.Update()
        ids.append(bid)
    return ids


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None, help="output .vsp3 path")
    ap.add_argument("--no-masses", action="store_true")
    ap.add_argument("--no-controls", action="store_true")
    ap.add_argument("--list-parms", action="store_true")
    args = ap.parse_args()

    cfg = G.CONFIG
    props = G.properties(cfg)
    out = args.out or f"{cfg['name']}.vsp3"
    airfoil = cfg["airfoil_file"]

    vsp.VSPRenew()
    vsp.ClearVSPModel()

    wid, modes = build_wing(cfg, props, airfoil)

    if args.list_parms:
        print("Parameters on MainWing:")
        list_parms(wid)
        return

    if not args.no_controls:
        add_control_surfaces(wid, cfg)
    if not args.no_masses:
        add_point_masses(cfg, props)

    vsp.Update()

    errs = vsp.GetErrorMsgs() if hasattr(vsp, "GetErrorMsgs") else []
    vsp.WriteVSPFile(out, vsp.SET_ALL)

    # ---- verify OpenVSP's own area against our analytic value ------------
    vsp.SetAnalysisInputDefaults("CompGeom")
    rid = vsp.ExecAnalysis("CompGeom")
    try:
        vsp_area = sum(vsp.GetDoubleResults(rid, "Wet_Area"))
    except Exception:
        vsp_area = float("nan")

    print(f"wrote {out}")
    print(f"  airfoil mode per station : {modes}")
    if "fallback" in modes:
        print("  WARNING: reflex .dat file not found - a four-series stand-in was")
        print(f"           used. Put {airfoil} next to this script before trusting Cm.")
    print(f"  analytic S_ref           : {props['s_ref']:.5f} m^2")
    print(f"  OpenVSP wetted area      : {vsp_area:.5f} m^2  (~2x S_ref plus winglets)")
    print(f"  MAC                      : {props['mac']:.4f} m")
    print(f"  target CG                : {props['x_cg_target']*1000:.1f} mm from nose")
    if errs:
        print("  OpenVSP messages:")
        for e in errs:
            print(f"    {e}")
    print("\nNext: python3 run_analysis.py")


if __name__ == "__main__":
    main()
