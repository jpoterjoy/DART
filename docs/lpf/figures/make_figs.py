"""Figures for docs/lpf (light and dark versions of each).
fig_localized_weights : localized weight 1 + rho (Ne w - 1) vs distance (illustration, exact formula)
fig_tempering         : Neff(beta) for one set of sharp weights and the beta that meets frac_neff (illustration)
fig_bgrid_ne          : bgrid prior RMSE vs ensemble size, LPF tempered / single pass (both with bins) vs EAKF
fig_bgrid_bins        : bgrid paired change in prior RMSE from merge-only binning (pf_bin_rho = 0.95)
bgrid numbers: lpf_build_test/bgrid ie*/ib* runs (ie_collect.py), installed LPF of 2026-10-27, 120 days, scored
days 60-120, 100 sites x U,V at level 3 every 6 h, cutoff 0.5, frac_neff 0.8, 3 seeds (Ne 1000: 1 seed).
usage: make_figs.py <out dir>"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out = sys.argv[1]
THEMES = {
    "light": dict(bg="#fcfcfb", ink="#0b0b0b", ink2="#52514e", grid="#e4e3df", s=["#2a78d6", "#eb6834", "#1baf7a"]),
    "dark":  dict(bg="#1a1a19", ink="#f0efec", ink2="#c3c2b7", grid="#34332f", s=["#3987e5", "#d95926", "#199e70"]),
}

def gc(z):
    z = np.abs(z); t = np.zeros_like(z)
    m1 = z <= 1; m2 = (z > 1) & (z < 2)
    a = z[m1]; t[m1] = -0.25*a**5 + 0.5*a**4 + 0.625*a**3 - 5/3*a**2 + 1
    a = z[m2]; t[m2] = a**5/12 - 0.5*a**4 + 0.625*a**3 + 5/3*a**2 - 5*a + 4 - 2/(3*a)
    return t

def style(ax, th):
    ax.set_facecolor(th["bg"])
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(th["grid"])
    ax.tick_params(colors=th["ink2"], labelsize=9, length=0)
    ax.grid(True, color=th["grid"], linewidth=1, linestyle="-"); ax.set_axisbelow(True)
    ax.xaxis.label.set_color(th["ink2"]); ax.yaxis.label.set_color(th["ink2"])

def fig(th, w=7.0, h=3.6):
    f, ax = plt.subplots(figsize=(w, h), dpi=150); f.patch.set_facecolor(th["bg"]); style(ax, th); return f, ax

def title(ax, th, t, sub):
    ax.set_title(t, loc="left", color=th["ink"], fontsize=11, fontweight="semibold", pad=22)
    ax.text(0, 1.03, sub, transform=ax.transAxes, color=th["ink2"], fontsize=8.5, va="bottom")

def save(f, name, mode):
    f.tight_layout(); f.savefig(os.path.join(out, "%s_%s.png" % (name, mode)), facecolor=f.get_facecolor()); plt.close(f)

for mode, th in THEMES.items():
    s = th["s"]
    # 1. localized weights
    f, ax = fig(th)
    d = np.linspace(0, 2.2, 400); rho = gc(d)
    for k, (a, lab) in enumerate(((3.0, "favoured particle (Ne·w = 3)"), (1.0, "neutral particle (Ne·w = 1)"),
                                  (0.1, "disfavoured particle (Ne·w = 0.1)"))):
        om = 1 + rho*(a - 1)
        ax.plot(d, om, color=s[k], lw=2, solid_capstyle="round", label=lab)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right", labelcolor=th["ink"])
    ax.set_xlim(0, 2.2); ax.set_xticks([0, 0.5, 1, 1.5, 2]); ax.set_ylim(0, 3.2)
    ax.set_xlabel("distance from the observation / cutoff"); ax.set_ylabel("localized weight ω")
    title(ax, th, "Localized weights fade to 1 away from the observation",
          "ω = 1 + ρ(d)·(Ne·w − 1), ρ = Gaspari–Cohn; beyond 2 × cutoff the observation has no effect")
    save(f, "fig_localized_weights", mode)

    # 2. tempering
    f, ax = fig(th)
    rng = np.random.default_rng(4); Ne = 80
    h = rng.normal(0, 1, Ne); lw = 0.5*(1.0 - h)**2/0.3; lw -= lw.min()
    b = np.linspace(0, 1, 401)
    neff = np.array([1/np.sum((np.exp(-bb*lw)/np.exp(-bb*lw).sum())**2) for bb in b])/Ne
    bs = b[np.argmin(np.abs(neff - 0.8))]
    ax.plot(b, neff, color=s[0], lw=2)
    ax.axhline(0.8, color=th["ink2"], lw=1)
    ax.plot([bs], [0.8], "o", ms=8, color=s[0], mec=th["bg"], mew=2)
    ax.annotate("β = %.2f meets frac_neff = 0.8\n(residual %.2f left for later passes)" % (bs, 1 - bs), (bs, 0.8),
                xytext=(bs + 0.08, 0.93), color=th["ink"], fontsize=8.5,
                arrowprops=dict(arrowstyle="-", color=th["ink2"], lw=1))
    ax.text(1.0, neff[-1] + 0.07, "untempered: Neff/Ne = %.2f" % neff[-1], color=th["ink"], fontsize=8.5, ha="right")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.05); ax.set_xlabel("tempering coefficient β"); ax.set_ylabel("Neff / Ne of exp(−β·lw)")
    title(ax, th, "Tempering keeps the effective ensemble size at the target",
          "One variable, Ne = 80, a sharp likelihood (illustration); each pass uses the largest β that meets frac_neff")
    save(f, "fig_tempering", mode)

    # 3. bgrid vs ensemble size
    f, ax = fig(th)
    ne = np.array([80, 160, 320, 1000])
    series = [("LPF tempered (pf_maxiter 20, KDDM, CF merge), bins", [1.049, 0.853, 0.772, 0.730], s[0], "o"),
              ("LPF single pass (QM merge), bins", [1.354, 1.120, 1.010, np.nan], s[1], "s"),
              ("EAKF, adaptive inflation", [0.688, 0.682, 0.689, 0.689], s[2], "D")]
    for lab, y, c, mk in series:
        y = np.array(y); ok = ~np.isnan(y)
        ax.plot(ne[ok], y[ok], color=c, lw=2, marker=mk, ms=6, mec=th["bg"], mew=1.5, label=lab)
        ax.text(ne[ok][-1]*1.08, y[ok][-1], "%.2f" % y[ok][-1], color=th["ink"], fontsize=8.5, va="center")
    ax.set_xscale("log", base=2); ax.set_xticks(ne); ax.set_xticklabels([str(n) for n in ne]); ax.set_xlim(70, 1500)
    ax.set_ylim(0.6, 1.45); ax.set_xlabel("ensemble size Ne"); ax.set_ylabel("prior RMSE (m/s, obs space)")
    lg = ax.legend(frameon=False, fontsize=8.5, loc="upper right", labelcolor=th["ink"])
    title(ax, th, "bgrid: the LPF approaches the EAKF as the ensemble grows",
          "bgrid, 100 sites × (u, v) every 6 h, days 60–120 scored, cutoff 0.5, frac_neff 0.8; 3 seeds (Ne 1000: 1)")
    save(f, "fig_bgrid_ne", mode)

    # 4. binning
    f, ax = fig(th)
    seeds = {("T", 80): [-6, -2, -2], ("T", 160): [-2, -6, -7], ("T", 320): [-5, -2, -3],
             ("S", 80): [-3, -3, -6], ("S", 160): [-2, -7, 1], ("S", 320): [-7, -4, -5]}
    means = {("T", 80): -3.4, ("T", 160): -4.9, ("T", 320): -3.4, ("S", 80): -3.9, ("S", 160): -2.6, ("S", 320): -5.2}
    x = np.arange(3); wbar = 0.22
    for k, (p, lab) in enumerate((("T", "tempered (pf_maxiter 20, KDDM, CF merge)"), ("S", "single pass (QM merge)"))):
        xs = x + (k - 0.5)*(wbar + 0.04)
        m = [means[(p, n)] for n in (80, 160, 320)]
        ax.bar(xs, m, width=wbar, color=s[k], label=lab, zorder=2)
        for xi, n in zip(xs, (80, 160, 320)):
            ax.plot([xi]*3, seeds[(p, n)], "o", ms=4, color=th["ink2"], zorder=3)
        for xi, mi in zip(xs, m):
            ax.text(xi, 1.6, "%.1f%%" % mi, ha="center", va="bottom", color=th["ink"], fontsize=8)
    ax.axhline(0, color=th["ink2"], lw=1)
    ax.set_xticks(x); ax.set_xticklabels(["Ne 80", "Ne 160", "Ne 320"]); ax.set_ylim(-9.5, 3.5)
    ax.set_ylabel("change in prior RMSE (%)")
    ax.legend(frameon=False, fontsize=8.5, loc="lower right", labelcolor=th["ink"])
    title(ax, th, "Merge-only binning (pf_bin_rho = 0.95) lowers RMSE and run time",
          "bgrid, paired with bins off; bars = mean of 3 seeds, dots = seeds; runs were also 15–27% faster")
    save(f, "fig_bgrid_bins", mode)
print("done")
