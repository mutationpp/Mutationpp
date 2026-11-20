"""Enthalpie de formation h_v(298) du composite liege/phenolique (P50), par
estimation elementaire.

Methode (bilan de combustion, loi de Hess) :

    hf(solide) = sum_i n_i hf(produit_i)  +  PCS

avec produits CO2(g) et H2O(l), et PCS le pouvoir calorifique superieur
estime par correlation sur l'analyse elementaire (C, H, O). La reference est
celle de Mutation++ (elements a l'etat standard a 298.15 K, carbone graphite
= 0), donc celle de la table h_g de cork_pyrolysis_gas.csv.

C'EST UNE ESTIMATION : les correlations de PCS sont precises a ~ +-1.5 MJ/kg
sur la biomasse ; aucune mesure calorimetrique du P50 n'existe dans l'article.
Les variantes (correlations, composition du liege, hf de la resine) donnent
l'incertitude.

Usage : python3 cork_hv.py
"""
import cork_pyrolysis_data as d

HF_CO2 = -393.51e3     # J/mol  CO2(g), NIST/CODATA
HF_H2O_L = -285.83e3   # J/mol  H2O(l)
M = d.M

# PCS [MJ/kg] en fonction des % masse (C, H, O)
CORRELATIONS = {
    "Channiwala-Parikh": lambda C, H, O: 0.3491*C + 1.1783*H - 0.1034*O,
    "Boie":              lambda C, H, O: 0.3517*C + 1.16225*H - 0.1109*O,
    "Dulong":            lambda C, H, O: 0.3383*C + 1.443*(H - O/8.0),
}


def hf_products(y):
    """sum n_i hf_i par kg de solide de fractions massiques y={C,H,O} [J/kg]."""
    n_co2 = 1000*y["C"]/M["C"]
    n_h2o = 1000*y["H"]/M["H"]/2
    return n_co2*HF_CO2 + n_h2o*HF_H2O_L


def hf_solid(y, corr):
    pcs = CORRELATIONS[corr](100*y["C"], 100*y["H"], 100*y.get("O", 0.0))*1e6
    return hf_products(y) + pcs, pcs


def resin_y():
    u = d.parse_formula(d.RESIN_FORMULA)
    mm = d.molar_mass(u)
    return {e: u.get(e, 0.0)*M[e]/mm for e in "CHO"}


def main():
    w_c, w_r = d.W_CORK, d.W_RESIN
    comps = {
        "litterature C62.4/H8.5/O28.4": d.normalize(d.CORK_MASS_PCT_LITERATURE),
        "reconstruit biochimique":      d.normalize(d.cork_elemental()),
    }
    yr = resin_y()
    HF_RESIN_BASE = -2143700.0   # J/kg, resine phenolique des bases TACOT/ZURAM
    print("composite = %.0f %% liege + %.0f %% resine (masse), char = C pur "
          "(hf_char = 0)\n" % (100*w_c, 100*w_r))
    rows = []
    for cn, yc in comps.items():
        for corr in CORRELATIONS:
            hf_c, pcs_c = hf_solid(yc, corr)
            hf_r, pcs_r = hf_solid(yr, corr)
            for rn, hfr in (("resine par correlation", hf_r),
                            ("resine = base TACOT/ZURAM", HF_RESIN_BASE)):
                hv = w_c*hf_c + w_r*hfr
                rows.append((cn, corr, rn, pcs_c/1e6, hf_c/1e6, hfr/1e6, hv/1e6))
    print("%-30s %-18s %-26s %6s %8s %8s %9s" % (
        "liege", "correlation PCS", "resine", "PCSliege", "hf_liege", "hf_res",
        "hv(298)"))
    for r in rows:
        print("%-30s %-18s %-26s %6.2f %8.3f %8.3f %9.3f" % r)
    ref = [r for r in rows if r[0].startswith("litt") and
           r[1] == "Channiwala-Parikh" and r[2].startswith("resine = base")][0]
    hv = [r[6] for r in rows]
    print("\nREFERENCE (litterature / Channiwala-Parikh / resine = base "
          "TACOT-ZURAM) : hv(298) = %.3f MJ/kg" % ref[6])
    print("plage sur les variantes : [%.3f ; %.3f] MJ/kg  (unites MJ/kg de "
          "composite)" % (min(hv), max(hv)))
    rv, rc = 480.0, 96.0
    print("hbar(298) = hv*rho_v/(rho_v-rho_c) = %.3f MJ/kg (reference)"
          % (ref[6]*rv/(rv-rc)))


if __name__ == "__main__":
    main()
