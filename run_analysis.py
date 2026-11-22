"""
Run VSPAERO on the cranked-kite model and report the numbers that matter.

Usage:
    python3 run_analysis.py                     # loiter condition, alpha sweep
    python3 run_analysis.py --condition dash
    python3 run_analysis.py --alpha -4 14 19    # start, end, npts
    python3 run_analysis.py --no-stability

Outputs:
    <name>_<condition>_polar.csv   drag polar and L/D
    <name>_<condition>_summary.txt readable report
    prints CL_max-ish, L/D_max, neutral point and static margin
"""

import argparse
import csv
import math
import os
import sys

import geometry as G

try:
    import openvsp as vsp
except ImportError:
    sys.exit("openvsp python module not found - see build_model.py header.")


def _results_to_dict(rid):
    """Pull every double-valued field out of a VSPAERO result set."""
    out = {}
    for name in vsp.GetAllDataNames(rid):
        try:
            vals = vsp.GetDoubleResults(rid, name)
        except Exception:
            continue
        if vals:
            out[name] = vals
    return out


def _pick(d, *names):
    """First matching key, since VSPAERO field names vary by release."""
    for n in names:
        if n in d:
            return d[n]
    return None


def run_sweep(cfg, props, condition, alpha_spec, stability):
    vsp_file = f"{cfg['name']}.vsp3"
    if not os.path.isfile(vsp_file):
        sys.exit(f"{vsp_file} not found. Run build_model.py first.")

    vsp.VSPRenew()
    vsp.ClearVSPModel()
    vsp.ReadVSPFile(vsp_file)

    # Degenerate geometry is what the vortex lattice solver actually reads.
    vsp.SetAnalysisInputDefaults("VSPAERODegenGeom")
    vsp.ExecAnalysis("VSPAERODegenGeom")

    mach = props["mach_loiter"] if condition == "loiter" else props["mach_dash"]
    re = props["re_loiter"] if condition == "loiter" else props["re_dash"]
    v = cfg["v_loiter"] if condition == "loiter" else cfg["v_dash"]

    a_start, a_end, a_npts = alpha_spec

    A = "VSPAEROSweep"
    vsp.SetAnalysisInputDefaults(A)

    def seti(name, val):
        try:
            vsp.SetIntAnalysisInput(A, name, [int(val)])
        except Exception:
            pass

    def setd(name, val):
        try:
            vsp.SetDoubleAnalysisInput(A, name, [float(val)])
        except Exception:
            pass

    seti("AnalysisMethod", vsp.VORTEX_LATTICE)
    seti("RefFlag", vsp.MANUAL_REF_FLAG if hasattr(vsp, "MANUAL_REF_FLAG") else 0)
    setd("Sref", props["s_ref"])
    setd("bref", props["b_ref"])
    setd("cref", props["mac"])
    setd("Xcg", props["x_cg"])
    setd("Ycg", 0.0)
    setd("Zcg", 0.0)
    setd("MachStart", mach)
    setd("MachEnd", mach)
    seti("MachNpts", 1)
    setd("AlphaStart", a_start)
    setd("AlphaEnd", a_end)
    seti("AlphaNpts", a_npts)
    setd("BetaStart", 0.0)
    seti("BetaNpts", 1)
    setd("ReCref", re)
    seti("WakeNumIter", 5)
    seti("Symmetry", 1)
    seti("StabilityType", 1 if stability else 0)

    print(f"running VSPAERO: {condition}, M={mach:.4f}, Re={re:.0f}, "
          f"alpha {a_start} to {a_end} in {a_npts} points")
    rid = vsp.ExecAnalysis(A)
    return _results_to_dict(rid), v, mach, re


def summarise(cfg, props, data, v, mach, re, condition):
    alpha = _pick(data, "Alpha")
    cl = _pick(data, "CL", "CLtot")
    cd = _pick(data, "CDtot", "CD", "CDi")
    cdi = _pick(data, "CDi")
    cm = _pick(data, "CMy", "Cmy", "CM")

    if not alpha or not cl or not cd:
        print("could not find Alpha/CL/CD in the result set. Fields returned:")
        for k in sorted(data):
            print("   ", k)
        sys.exit(1)

    n = min(len(alpha), len(cl), len(cd))
    rows = []
    for i in range(n):
        ld = cl[i] / cd[i] if cd[i] else float("nan")
        rows.append({
            "alpha_deg": alpha[i],
            "CL": cl[i],
            "CD": cd[i],
            "CDi": cdi[i] if cdi and i < len(cdi) else "",
            "Cm": cm[i] if cm and i < len(cm) else "",
            "L_over_D": ld,
        })

    best = max(rows, key=lambda r: r["L_over_D"] if r["L_over_D"] == r["L_over_D"] else -1)

    # --- trimmed cruise: the alpha where CL matches the required CL --------
    cl_req = props["cl_loiter"] if condition == "loiter" else props["cl_dash"]
    trim = min(rows, key=lambda r: abs(r["CL"] - cl_req))

    # --- pitch stability --------------------------------------------------
    cm_alpha = None
    x_np = None
    if cm and len(cm) >= 2:
        # least-squares slope of Cm vs CL gives -(static margin)
        pairs = [(r["CL"], r["Cm"]) for r in rows if r["Cm"] != ""]
        if len(pairs) >= 2:
            xs = [p[0] for p in pairs]
            ys = [p[1] for p in pairs]
            mx = sum(xs) / len(xs)
            my = sum(ys) / len(ys)
            num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
            den = sum((x - mx) ** 2 for x in xs)
            if den:
                dcm_dcl = num / den
                sm = -dcm_dcl
                x_np = props["x_cg"] + sm * props["mac"]
                cm_alpha = dcm_dcl

    # --- power and endurance ---------------------------------------------
    # P_required = D * V / (eta_prop * eta_motor * eta_esc)
    eta = 0.55
    d_trim = 0.5 * cfg["rho"] * v * v * props["s_ref"] * trim["CD"]
    p_shaft = d_trim * v
    p_elec = p_shaft / eta
    # 4S 3000 mAh = 14.8 V * 3.0 Ah = 44.4 Wh, usable 85 %
    wh_usable = 14.8 * 3.0 * 0.85
    endurance_min = 60.0 * wh_usable / p_elec if p_elec > 0 else 0.0
    range_km = endurance_min * 60.0 * v / 1000.0

    L = []
    A_ = L.append
    A_("=" * 68)
    A_(f"  VSPAERO summary - {cfg['name']} - {condition}")
    A_("=" * 68)
    A_(f"  V = {v:.1f} m/s ({v*3.6:.0f} km/h)   M = {mach:.4f}   Re = {re:.0f}")
    A_(f"  Sref = {props['s_ref']:.5f} m^2   bref = {props['b_ref']:.4f} m   "
       f"cref = {props['mac']:.4f} m")
    A_(f"  Xcg  = {props['x_cg']*1000:.1f} mm from nose")
    A_("")
    A_("BEST L/D")
    A_(f"  alpha = {best['alpha_deg']:.2f} deg   CL = {best['CL']:.4f}   "
       f"CD = {best['CD']:.5f}   L/D = {best['L_over_D']:.2f}")
    A_("")
    A_(f"TRIMMED AT THE {condition.upper()} CONDITION  (CL required = {cl_req:.4f})")
    A_(f"  alpha = {trim['alpha_deg']:.2f} deg   CL = {trim['CL']:.4f}   "
       f"CD = {trim['CD']:.5f}   L/D = {trim['L_over_D']:.2f}")
    A_(f"  drag  = {d_trim:.2f} N")
    A_(f"  shaft power = {p_shaft:.1f} W    electrical power = {p_elec:.1f} W "
       f"(eta = {eta:.2f})")
    A_(f"  endurance   = {endurance_min:.1f} min on 44.4 Wh (85 % usable)")
    A_(f"  still-air range = {range_km:.1f} km")
    A_("")
    A_("PITCH STABILITY")
    if x_np is not None:
        sm = (x_np - props["x_cg"]) / props["mac"]
        A_(f"  dCm/dCL          = {cm_alpha:+.4f}")
        A_(f"  neutral point    = {x_np*1000:.1f} mm from nose")
        A_(f"  static margin    = {sm*100:+.2f} % MAC")
        if sm < 0:
            A_("    -> UNSTABLE. Move the CG forward immediately.")
        elif sm < 0.05:
            A_("    -> stable but twitchy. Move the CG forward a few mm.")
        elif sm > 0.08:
            A_("    -> over-stable. Sluggish pitch, wasted trim drag. CG aft.")
        else:
            A_("    -> inside the 5-8 % target window.")
        A_(f"  (analytic pre-estimate was {props['x_np_estimate']*1000:.1f} mm; "
           f"trust VSPAERO over that)")
    else:
        A_("  Cm not returned - rerun with stability enabled.")
    A_("")
    A_("CAVEATS TO STATE IN THE REPORT")
    A_("  Vortex lattice assumes attached inviscid flow: no stall prediction,")
    A_(f"  and at Re ~ {re:.0e} the real sections suffer laminar separation that")
    A_("  VLM cannot see. Treat CL beyond about 0.9 as fiction, and cross-check")
    A_("  the airfoil polars in XFLR5 at this Reynolds number.")
    A_("=" * 68)

    return rows, "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--condition", choices=["loiter", "dash"], default="loiter")
    ap.add_argument("--alpha", nargs=3, type=float, default=[-4.0, 14.0, 19.0],
                    metavar=("START", "END", "NPTS"))
    ap.add_argument("--no-stability", action="store_true")
    args = ap.parse_args()

    cfg = G.CONFIG
    props = G.properties(cfg)

    data, v, mach, re = run_sweep(
        cfg, props, args.condition,
        (args.alpha[0], args.alpha[1], int(args.alpha[2])),
        not args.no_stability,
    )
    rows, text = summarise(cfg, props, data, v, mach, re, args.condition)

    stem = f"{cfg['name']}_{args.condition}"
    with open(f"{stem}_polar.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(f"{stem}_summary.txt", "w") as f:
        f.write(text + "\n")

    print()
    print(text)
    print(f"\nwrote {stem}_polar.csv and {stem}_summary.txt")


if __name__ == "__main__":
    main()
