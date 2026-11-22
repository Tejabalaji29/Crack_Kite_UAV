"""
Cranked-kite tailless UAV - geometry definition and derived properties.

Single source of truth. build_model.py and run_analysis.py both import from here,
so changing a number in CONFIG propagates through the whole toolchain.

All lengths in metres, angles in degrees, masses in kilograms.
Run directly to print the geometry report:  python3 geometry.py
"""

import math

# ----------------------------------------------------------------------------
# CONFIG - edit these for trade studies
# ----------------------------------------------------------------------------

CONFIG = {
    "name": "cranked_kite_v1",

    # Mass
    "mtow": 2.00,

    # Inboard panel: root -> crank
    "p1_span":        0.250,   # half-span of panel (per side)
    "p1_root_chord":  0.520,
    "p1_tip_chord":   0.300,
    "p1_sweep_le":    42.0,
    "p1_dihedral":     0.0,
    "p1_twist_root":   0.0,
    "p1_twist_tip":   -1.5,

    # Outboard panel: crank -> tip
    "p2_span":        0.200,
    "p2_root_chord":  0.300,
    "p2_tip_chord":   0.120,
    "p2_sweep_le":    28.0,
    "p2_dihedral":     2.0,
    "p2_twist_tip":   -4.0,

    # Winglet (not counted in reference area)
    "wl_span":        0.060,
    "wl_root_chord":  0.120,
    "wl_tip_chord":   0.070,
    "wl_sweep_le":    30.0,
    "wl_dihedral":    80.0,

    # Airfoil thickness ratios at the three main stations
    "tc_root":  0.12,
    "tc_crank": 0.10,
    "tc_tip":   0.09,

    # Airfoil data files - reflex sections. Put the .dat files next to this script.
    # MH 60 is the primary choice; EH 1.5/9 and MH 78 are valid substitutes.
    "airfoil_file": "mh60.dat",

    # Control surfaces, as fractions
    "elevon_span_start": 0.58,   # fraction of half-span
    "elevon_span_end":   0.95,
    "elevon_hinge_xc":   0.72,   # fraction of local chord
    "rudder_hinge_xc":   0.60,   # drag rudders on the winglets

    # Target stability
    "static_margin_target": 0.06,   # 6 % of MAC

    # Flight conditions
    "v_loiter": 17.0,   # m/s
    "v_dash":   28.0,   # m/s
    "rho":      1.225,  # kg/m^3 at sea level
    "nu":       1.46e-5,

    # Mass breakdown: (name, mass kg, x from nose m)
    "point_masses": [
        ("eo_sensor",        0.300, 0.120),
        ("battery_4s_3000",  0.320, 0.190),
        ("autopilot_servos", 0.150, 0.220),
        ("motor_esc",        0.180, 0.480),
    ],
    "structure_mass": 1.050,
}


# ----------------------------------------------------------------------------
# Panel maths
# ----------------------------------------------------------------------------

def _panel_area(b, cr, ct):
    """Planform area of one linearly tapered panel (one side only)."""
    return b * (cr + ct) / 2.0


def _panel_c2_integral(b, cr, ct):
    """Integral of c^2 dy over a linearly tapered panel."""
    return (b / 3.0) * (cr * cr + cr * ct + ct * ct)


def _panel_yc_integral(b, cr, ct):
    """
    Integral of y * c dy over a linearly tapered panel, y measured from the
    panel root.

        int_0^b y (cr + (ct-cr) y/b) dy = b^2 (cr + 2 ct) / 6

    Note the ordering: cr + 2ct, not 2cr + ct. Getting this backwards shifts
    y_MAC outboard and corrupts every downstream CG number.
    """
    return b * b * (cr + 2.0 * ct) / 6.0


def properties(cfg=None):
    """Compute every derived geometric and aerodynamic reference quantity."""
    c = dict(CONFIG if cfg is None else cfg)

    b1, cr1, ct1 = c["p1_span"], c["p1_root_chord"], c["p1_tip_chord"]
    b2, cr2, ct2 = c["p2_span"], c["p2_root_chord"], c["p2_tip_chord"]

    # --- reference area and span (wing only, winglets excluded) ---
    s_half = _panel_area(b1, cr1, ct1) + _panel_area(b2, cr2, ct2)
    s_ref = 2.0 * s_half
    b_ref = 2.0 * (b1 + b2)
    ar = b_ref ** 2 / s_ref

    # --- mean aerodynamic chord ---
    c2_int = _panel_c2_integral(b1, cr1, ct1) + _panel_c2_integral(b2, cr2, ct2)
    mac = c2_int / s_half

    # --- spanwise station of the MAC ---
    yc_int = _panel_yc_integral(b1, cr1, ct1)
    # panel 2 contribution, shifted outboard by b1
    yc_int += _panel_yc_integral(b2, cr2, ct2) + b1 * _panel_area(b2, cr2, ct2)
    y_mac = yc_int / s_half

    # --- leading-edge x at the MAC station ---
    tan1 = math.tan(math.radians(c["p1_sweep_le"]))
    tan2 = math.tan(math.radians(c["p2_sweep_le"]))
    if y_mac <= b1:
        x_le_mac = y_mac * tan1
    else:
        x_le_mac = b1 * tan1 + (y_mac - b1) * tan2

    x_mac_quarter = x_le_mac + 0.25 * mac

    # --- overall planform envelope ---
    x_tip_le = b1 * tan1 + b2 * tan2
    root_te = cr1
    tip_te = x_tip_le + ct2
    overall_length = max(root_te, tip_te)

    # --- taper and sweep summary ---
    taper_overall = ct2 / cr1
    # area-weighted quarter-chord sweep, a useful single number for reports
    sweep_c4_1 = math.degrees(math.atan(tan1 + (ct1 - cr1) * 0.25 / b1))
    sweep_c4_2 = math.degrees(math.atan(tan2 + (ct2 - cr2) * 0.25 / b2))
    a1 = _panel_area(b1, cr1, ct1)
    a2 = _panel_area(b2, cr2, ct2)
    sweep_c4 = (sweep_c4_1 * a1 + sweep_c4_2 * a2) / (a1 + a2)

    # --- loading and flight condition ---
    w = c["mtow"] * 9.80665
    wing_loading = c["mtow"] / s_ref

    def cl_at(v):
        return w / (0.5 * c["rho"] * v * v * s_ref)

    def reynolds(v):
        return v * mac / c["nu"]

    # --- mass and CG ---
    pm = c["point_masses"]
    pm_mass = sum(m for _, m, _ in pm)
    pm_moment = sum(m * x for _, m, x in pm)
    # structure assumed centred on the MAC quarter-chord as a first guess
    struct_x = x_mac_quarter
    total_mass = pm_mass + c["structure_mass"]
    x_cg = (pm_moment + c["structure_mass"] * struct_x) / total_mass

    # --- CG target window from the static-margin goal ---
    # VSPAERO supplies the real neutral point. For a swept tailless wing the
    # first-cut estimate puts it slightly aft of the MAC quarter-chord.
    x_np_estimate = x_le_mac + 0.28 * mac
    sm = c["static_margin_target"]
    x_cg_target = x_np_estimate - sm * mac
    x_cg_fwd = x_np_estimate - 0.08 * mac
    x_cg_aft = x_np_estimate - 0.05 * mac

    # --- control surface geometry ---
    def chord_at(y):
        if y <= b1:
            return cr1 + (ct1 - cr1) * (y / b1)
        return cr2 + (ct2 - cr2) * ((y - b1) / b2)

    y_el_1 = c["elevon_span_start"] * (b1 + b2)
    y_el_2 = c["elevon_span_end"] * (b1 + b2)
    el_c1 = chord_at(y_el_1) * (1.0 - c["elevon_hinge_xc"])
    el_c2 = chord_at(y_el_2) * (1.0 - c["elevon_hinge_xc"])
    elevon_area = 2.0 * (y_el_2 - y_el_1) * (el_c1 + el_c2) / 2.0

    return {
        "name": c["name"],
        "s_ref": s_ref,
        "b_ref": b_ref,
        "aspect_ratio": ar,
        "mac": mac,
        "y_mac": y_mac,
        "x_le_mac": x_le_mac,
        "x_mac_quarter": x_mac_quarter,
        "overall_length": overall_length,
        "x_tip_le": x_tip_le,
        "taper_overall": taper_overall,
        "sweep_c4": sweep_c4,
        "mtow": c["mtow"],
        "weight_n": w,
        "wing_loading": wing_loading,
        "cl_loiter": cl_at(c["v_loiter"]),
        "cl_dash": cl_at(c["v_dash"]),
        "re_loiter": reynolds(c["v_loiter"]),
        "re_dash": reynolds(c["v_dash"]),
        "mach_loiter": c["v_loiter"] / 343.0,
        "mach_dash": c["v_dash"] / 343.0,
        "mass_total": total_mass,
        "x_cg": x_cg,
        "x_np_estimate": x_np_estimate,
        "x_cg_target": x_cg_target,
        "x_cg_window": (x_cg_fwd, x_cg_aft),
        "elevon_area": elevon_area,
        "elevon_area_frac": elevon_area / s_ref,
        "config": c,
    }


def static_margin(x_cg, x_np, mac):
    """Static margin as a fraction of MAC. Positive means statically stable."""
    return (x_np - x_cg) / mac


def report(cfg=None):
    p = properties(cfg)
    c = p["config"]
    L = []
    A = L.append
    A("=" * 68)
    A(f"  {p['name']}  -  cranked-kite tailless UAV")
    A("=" * 68)
    A("")
    A("PLANFORM")
    A(f"  Reference area  S_ref      {p['s_ref']*1e4:8.1f} cm^2   ({p['s_ref']:.4f} m^2)")
    A(f"  Span            b_ref      {p['b_ref']*1000:8.1f} mm")
    A(f"  Aspect ratio    AR         {p['aspect_ratio']:8.3f}")
    A(f"  MAC                        {p['mac']*1000:8.1f} mm")
    A(f"  MAC station     y_mac      {p['y_mac']*1000:8.1f} mm from centreline")
    A(f"  MAC leading edge x         {p['x_le_mac']*1000:8.1f} mm from nose")
    A(f"  MAC quarter-chord x        {p['x_mac_quarter']*1000:8.1f} mm from nose")
    A(f"  Overall length             {p['overall_length']*1000:8.1f} mm")
    A(f"  Tip LE offset              {p['x_tip_le']*1000:8.1f} mm aft of nose")
    A(f"  Overall taper ratio        {p['taper_overall']:8.3f}")
    A(f"  Mean quarter-chord sweep   {p['sweep_c4']:8.2f} deg")
    A("")
    A("LOADING AND FLIGHT CONDITION")
    A(f"  MTOW                       {p['mtow']:8.3f} kg")
    A(f"  Wing loading               {p['wing_loading']:8.2f} kg/m^2")
    A(f"  CL at loiter {c['v_loiter']:.0f} m/s        {p['cl_loiter']:8.3f}")
    A(f"  CL at dash   {c['v_dash']:.0f} m/s        {p['cl_dash']:8.3f}")
    A(f"  Re at loiter               {p['re_loiter']:8.0f}")
    A(f"  Re at dash                 {p['re_dash']:8.0f}")
    A(f"  Mach loiter / dash         {p['mach_loiter']:8.4f} / {p['mach_dash']:.4f}")
    A("")
    A("MASS AND BALANCE")
    A(f"  Summed mass                {p['mass_total']:8.3f} kg")
    for n, m, x in c["point_masses"]:
        A(f"    {n:<20s}     {m:8.3f} kg at x = {x*1000:6.1f} mm")
    A(f"    {'structure':<20s}     {c['structure_mass']:8.3f} kg at MAC c/4")
    A(f"  Estimated CG               {p['x_cg']*1000:8.1f} mm from nose")
    A(f"  Estimated neutral point    {p['x_np_estimate']*1000:8.1f} mm  (VSPAERO gives the real value)")
    A(f"  CG target for {c['static_margin_target']*100:.0f} % SM      {p['x_cg_target']*1000:8.1f} mm from nose")
    A(f"  Acceptable CG window       {p['x_cg_window'][0]*1000:.1f} to {p['x_cg_window'][1]*1000:.1f} mm")
    sm_now = static_margin(p["x_cg"], p["x_np_estimate"], p["mac"])
    A(f"  Static margin as laid out  {sm_now*100:8.2f} %")
    if sm_now < 0.05:
        A("    -> CG too far aft. Move the battery forward.")
    elif sm_now > 0.08:
        A("    -> CG too far forward. Move the battery aft.")
    else:
        A("    -> within the 5-8 % target window.")
    A("")
    A("CONTROL SURFACES")
    A(f"  Elevon area (both)         {p['elevon_area']*1e4:8.1f} cm^2")
    A(f"  Elevon area / S_ref        {p['elevon_area_frac']*100:8.2f} %")
    if p["elevon_area_frac"] < 0.06:
        A("    -> below 6 %. Extend the elevons inboard toward the crank.")
    else:
        A("    -> adequate for a tailless wing of this sweep.")
    A("")
    A("VSPAERO REFERENCE VALUES - enter these in the analysis setup")
    A(f"  Sref = {p['s_ref']:.5f}    bref = {p['b_ref']:.4f}    cref = {p['mac']:.4f}")
    A(f"  X_cg = {p['x_cg']:.4f}")
    A("=" * 68)
    return "\n".join(L)


if __name__ == "__main__":
    print(report())
