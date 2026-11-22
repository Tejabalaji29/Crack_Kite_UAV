"""
Trade studies on the cranked-kite planform.

Each study rebuilds the geometry, runs VSPAERO, and tabulates the result, so
you get the plots your report needs without clicking through the GUI 30 times.

Usage:
    python3 trade_study.py crank        # crank position 45..70 % half-span
    python3 trade_study.py sweep        # outboard sweep 20..40 deg
    python3 trade_study.py washout      # tip twist 0..-7 deg
    python3 trade_study.py --dry-run crank   # geometry only, no VSPAERO

Writes study_<name>.csv. --dry-run needs no OpenVSP install and is useful for
checking that the geometry stays sane before spending solver time.
"""

import argparse
import csv
import sys

import geometry as G

STUDIES = {
    # name: (label, list of (config overrides, x-value))
    "crank": "crank position (fraction of half-span)",
    "sweep": "outboard leading-edge sweep (deg)",
    "washout": "tip twist (deg)",
}


def variants(study):
    """Yield (xvalue, config) pairs, holding span and area as constant as possible."""
    base = G.CONFIG
    total_half = base["p1_span"] + base["p2_span"]

    if study == "crank":
        # Move the crank spanwise while holding span, root chord, tip chord and
        # reference AREA fixed. The chord at the crank is solved from the area
        # constraint, so each variant is a genuinely different break severity
        # rather than a relabelled straight-taper wing:
        #
        #   S_half = b1 (cr + cc)/2 + b2 (cc + ct)/2   with b1 + b2 = B
        #   =>  cc = (2 S_half - b1 cr - b2 ct) / B
        #
        s_half_target = G.properties(base)["s_ref"] / 2.0
        cr = base["p1_root_chord"]
        ct = base["p2_tip_chord"]
        for frac in (0.45, 0.50, 0.556, 0.60, 0.65, 0.70):
            c = dict(base)
            b1 = frac * total_half
            b2 = total_half - b1
            cc = (2.0 * s_half_target - b1 * cr - b2 * ct) / total_half
            c["p1_span"] = b1
            c["p2_span"] = b2
            c["p1_tip_chord"] = cc
            c["p2_root_chord"] = cc
            c["wl_root_chord"] = ct
            yield frac, c

    elif study == "sweep":
        for sw in (20.0, 24.0, 28.0, 32.0, 36.0, 40.0):
            c = dict(base)
            c["p2_sweep_le"] = sw
            yield sw, c

    elif study == "washout":
        for tw in (0.0, -1.5, -3.0, -4.0, -5.5, -7.0):
            c = dict(base)
            c["p2_twist_tip"] = tw
            c["p1_twist_tip"] = tw * 0.375   # keep the same distribution shape
            yield tw, c

    else:
        sys.exit(f"unknown study '{study}'. choose from: {', '.join(STUDIES)}")


def geometry_row(x, c):
    p = G.properties(c)
    return {
        "x": x,
        "crank_chord_mm": round(c["p1_tip_chord"] * 1000, 1),
        "crank_y_mm": round(c["p1_span"] * 1000, 1),
        "S_ref_m2": round(p["s_ref"], 5),
        "AR": round(p["aspect_ratio"], 4),
        "MAC_m": round(p["mac"], 5),
        "y_MAC_m": round(p["y_mac"], 5),
        "x_MAC_c4_mm": round(p["x_mac_quarter"] * 1000, 2),
        "sweep_c4_deg": round(p["sweep_c4"], 3),
        "taper": round(p["taper_overall"], 4),
        "elevon_frac_pct": round(p["elevon_area_frac"] * 100, 3),
        "CG_target_mm": round(p["x_cg_target"] * 1000, 2),
    }


def aero_row(x, c):
    """Build + solve one variant. Imported lazily so --dry-run needs no OpenVSP."""
    import build_model
    import run_analysis
    import openvsp as vsp

    p = G.properties(c)
    vsp.VSPRenew()
    vsp.ClearVSPModel()
    build_model.build_wing(c, p, c["airfoil_file"])
    vsp.Update()
    tmp = f"_study_{x}.vsp3".replace("-", "m")
    vsp.WriteVSPFile(tmp, vsp.SET_ALL)

    saved = dict(G.CONFIG)
    G.CONFIG.clear()
    G.CONFIG.update(c)
    G.CONFIG["name"] = tmp[:-5]
    try:
        data, v, mach, re = run_analysis.run_sweep(
            G.CONFIG, p, "loiter", (-2.0, 10.0, 13), True)
        rows, _ = run_analysis.summarise(G.CONFIG, p, data, v, mach, re, "loiter")
        best = max(rows, key=lambda r: r["L_over_D"])
        return {"LD_max": round(best["L_over_D"], 3),
                "alpha_LDmax": round(best["alpha_deg"], 2),
                "CL_LDmax": round(best["CL"], 4)}
    finally:
        G.CONFIG.clear()
        G.CONFIG.update(saved)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("study", choices=list(STUDIES))
    ap.add_argument("--dry-run", action="store_true",
                    help="geometry only, skip VSPAERO")
    args = ap.parse_args()

    print(f"study: {STUDIES[args.study]}")
    out = []
    for x, c in variants(args.study):
        row = geometry_row(x, c)
        if not args.dry_run:
            try:
                row.update(aero_row(x, c))
            except Exception as e:
                print(f"  VSPAERO failed at x={x}: {e}")
                row.update({"LD_max": "", "alpha_LDmax": "", "CL_LDmax": ""})
        out.append(row)
        print("  " + "  ".join(f"{k}={v}" for k, v in row.items()))

    path = f"study_{args.study}.csv"
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
