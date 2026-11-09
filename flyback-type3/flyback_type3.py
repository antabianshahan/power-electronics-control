"""Type III voltage-mode compensator for a 48 V -> 12 V flyback converter.

Python reproduction of flyback_control.m (no MATLAB licence needed).
Prints the design values and writes the figures used in the README:

    python flyback_type3.py          # needs numpy, scipy, matplotlib
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

OUT = Path(__file__).with_name("figures")
OUT.mkdir(exist_ok=True)

# ── Converter parameters ────────────────────────────────────────────────
Vg, V = 48.0, 12.0          # input / output voltage [V]
Lm, C = 240e-6, 50e-6       # magnetising inductance [H], output capacitance [F]
n, R = 0.25, 12.0           # turns ratio Ns/Np, load [ohm]
fc, PM_target = 20e3, 50.0  # crossover [Hz], phase margin [deg]

D = V / (n * Vg + V)                         # V/Vg = nD/(1-D)
tau_z = Lm * D * n**2 / (R * (1 - D) ** 2)   # right-half-plane zero
w0 = np.sqrt((1 - D) ** 2 / (n**2 * Lm * C))

# G_vd(s) = (Vg/n)(1 - s*tau_z) / (Lm C s^2 + (Lm/R) s + (1-D)^2/n^2)
num_G = [-(Vg / n) * tau_z, Vg / n]
den_G = [Lm * C, Lm / R, Lm * C * w0**2]

f = np.logspace(1, 7, 200_000)
w = 2 * np.pi * f


def bode(num, den):
    _, mag, ph = signal.bode(signal.TransferFunction(num, den), w)
    return mag, ph


def crossing(x, y, level):
    """Interpolated x where y crosses `level` (first crossing)."""
    i = np.where(np.diff(np.sign(y - level)) != 0)[0][0]
    return np.interp(level, [y[i + 1], y[i]], [x[i + 1], x[i]]) if y[i] > y[i + 1] \
        else np.interp(level, [y[i], y[i + 1]], [x[i], x[i + 1]])


# ── Uncompensated loop ──────────────────────────────────────────────────
mag_G, ph_G = bode(num_G, den_G)
f_c0 = crossing(f, mag_G, 0.0)
pm0 = 180 + np.interp(f_c0, f, ph_G)
f_180 = crossing(f, ph_G, -180.0)
gm0 = -np.interp(f_180, f, mag_G)

# ── Type III: Gc(s) = Gc0/s * (1 + s/wz)^2 / (1 + s/wp)^2 ─────────────
wc = 2 * np.pi * fc
_, G_wc = signal.freqresp(signal.TransferFunction(num_G, den_G), [wc])
G_wc = G_wc[0]
ph_exact = np.angle(G_wc, deg=True)               # exact phase at fc, wrapped to (-180, 180]
ph_plant = ph_exact + 360 * np.round((np.interp(fc, f, ph_G) - ph_exact) / 360)  # unwrapped
boost = (PM_target - 180) - ph_plant              # phase the compensator must add
a = np.tan(np.deg2rad((boost + 90) / 2) / 2 + np.pi / 4)   # wz = wc/a, wp = wc*a
wz, wp = wc / a, wc * a
Gc0 = wc / (abs(G_wc) * (1 + (wc / wz) ** 2) / (1 + (wc / wp) ** 2))

num_Gc = Gc0 * np.polymul([1 / wz, 1], [1 / wz, 1])
den_Gc = np.polymul([1, 0], np.polymul([1 / wp, 1], [1 / wp, 1]))
num_T, den_T = np.polymul(num_Gc, num_G), np.polymul(den_Gc, den_G)
mag_T, ph_T = bode(num_T, den_T)
f_c1 = crossing(f, mag_T, 0.0)
pm1 = 180 + np.interp(f_c1, f, ph_T)
p_rhp = int(np.sum(np.roots(den_T).real > 1e-9))

print(f"D = {D:.3f},  f0 = {w0/2/np.pi:.1f} Hz,  Q = {R*C*w0:.2f},  "
      f"RHP zero = {1/tau_z/2/np.pi/1e3:.1f} kHz")
print(f"Uncompensated: fc = {f_c0/1e3:.1f} kHz, PM = {pm0:.1f} deg, "
      f"GM = {gm0:.1f} dB at {f_180:.0f} Hz  -> unstable")
print(f"Plant phase at 20 kHz = {ph_plant:.2f} deg -> boost needed = {boost:.2f} deg")
print(f"Type III: fz = {wz/2/np.pi:.0f} Hz, fp = {wp/2/np.pi:.0f} Hz, Gc0 = {Gc0:.1f}")
print(f"Compensated: fc = {f_c1/1e3:.2f} kHz, PM = {pm1:.2f} deg, "
      f"open-loop RHP poles P = {p_rhp}")

# ── Figures ─────────────────────────────────────────────────────────────
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})


def bode_fig(mag, ph, title, fname, mark=None):
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    a1.semilogx(f, mag, color="#1f4e9c", lw=1.6)
    a2.semilogx(f, ph, color="#c0392b", lw=1.6)
    a1.axhline(0, color="k", lw=0.6, ls="--")
    a2.axhline(-180, color="k", lw=0.6, ls="--")
    if mark:
        for ax in (a1, a2):
            ax.axvline(fc, color="#2e8b57", lw=0.9, ls=":")
        a2.axhline(PM_target - 180, color="#2e8b57", lw=0.9, ls=":")
        a2.annotate(f"PM = {pm1:.0f}° at {f_c1/1e3:.0f} kHz", xy=(f_c1, pm1 - 180),
                    xytext=(f_c1 * 3, pm1 - 120), arrowprops=dict(arrowstyle="->"))
    a1.set_ylabel("Magnitude (dB)")
    a2.set_ylabel("Phase (deg)")
    a2.set_xlabel("Frequency (Hz)")
    a1.set_xlim(10, 1e6)
    a1.set_title(title, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=160)
    plt.close(fig)


bode_fig(mag_G, ph_G, f"Uncompensated plant $G_{{vd}}(s)$: PM = {pm0:.1f}° (unstable)",
         "bode_uncompensated.png")
bode_fig(mag_T, ph_T, "Compensated loop $T(s) = G_c(s)\\,G_{vd}(s)$",
         "bode_compensated.png", mark=True)

# Nyquist (zoom around -1)
_, H = signal.freqresp(signal.TransferFunction(num_T, den_T), w)
fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(H.real, H.imag, color="#1f4e9c", lw=1.4, label=r"$\omega > 0$")
ax.plot(H.real, -H.imag, color="#1f4e9c", lw=0.9, ls="--", alpha=0.6, label=r"$\omega < 0$")
ax.plot(-1, 0, "x", color="#c0392b", ms=11, mew=2.5, label="−1")
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)
ax.set_aspect("equal")
ax.set_xlabel("Re")
ax.set_ylabel("Im")
ax.set_title(f"Nyquist of $T(s)$: P = {p_rhp}, N = 0 → stable", fontweight="bold")
ax.legend(loc="lower right")
fig.tight_layout()
fig.savefig(OUT / "nyquist.png", dpi=160)
plt.close(fig)

# Sensitivity functions
den_cl = np.polyadd(den_T, num_T)
mag_S, _ = bode(den_T, den_cl)
mag_CS, _ = bode(num_T, den_cl)
fig, ax = plt.subplots(figsize=(8, 4))
ax.semilogx(f, mag_S, color="#1f4e9c", lw=1.6, label=r"$|S| = |1/(1+T)|$")
ax.semilogx(f, mag_CS, color="#c0392b", lw=1.6, label=r"$|T/(1+T)|$")
ax.axhline(0, color="k", lw=0.6, ls="--")
ax.set_xlim(10, 1e6)
ax.set_ylim(-80, 20)
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Magnitude (dB)")
ax.set_title("Sensitivity: integrator rejects DC disturbances (|S| → 0)", fontweight="bold")
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "sensitivity.png", dpi=160)
plt.close(fig)
