import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Load CSV exported by verlet_oscillator_example
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cmake-build-debug")
df = pd.read_csv(os.path.join(base, "verlet_comparison.csv"))

t = df["t"].to_numpy()

# Extract each method
xE = df["x_euler"].to_numpy();  vE = df["v_euler"].to_numpy();  EE = df["E_euler"].to_numpy()
xR = df["x_rk4"].to_numpy();    vR = df["v_rk4"].to_numpy();    ER = df["E_rk4"].to_numpy()
xV = df["x_verlet"].to_numpy(); vV = df["v_verlet"].to_numpy(); EV = df["E_verlet"].to_numpy()
xT = df["x_true"].to_numpy();   vT = df["v_true"].to_numpy();   ET = df["E_true"].to_numpy()

# ── Figure 1: Energy conservation ─────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

# Top panel: RK4 and Verlet zoomed in
ax1.plot(t, ER, label="RK4",             color="tab:blue",  linewidth=1.2)
ax1.plot(t, EV, label="Velocity Verlet", color="tab:green", linewidth=1.5)
ax1.axhline(ET[0], color="black", linestyle="--", linewidth=1.0, label="Exact energy")
ax1.set_ylabel("E(t)")
ax1.set_title("Energy Conservation: RK4 vs Velocity Verlet (zoomed)")
ax1.legend()
ax1.grid(True, alpha=0.4)

# Bottom panel: all three including Euler explosion
ax2.plot(t, EE, label="Euler",           color="tab:red",   linewidth=1.2)
ax2.plot(t, ER, label="RK4",             color="tab:blue",  linewidth=1.2)
ax2.plot(t, EV, label="Velocity Verlet", color="tab:green", linewidth=1.5)
ax2.axhline(ET[0], color="black", linestyle="--", linewidth=1.0, label="Exact energy")
ax2.set_xlabel("t")
ax2.set_ylabel("E(t)")
ax2.set_title("Energy Conservation: All Methods (Euler instability visible)")
ax2.legend()
ax2.grid(True, alpha=0.4)

plt.tight_layout()
plt.savefig("verlet_energy.png", dpi=150)
plt.show()

# ── Figure 2: Energy error (drift from exact) ─────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(t, ER - ET, label="RK4 drift",           color="tab:blue",  linewidth=1.2)
ax.plot(t, EV - ET, label="Velocity Verlet drift", color="tab:green", linewidth=1.5)
ax.axhline(0, color="black", linestyle="--", linewidth=0.8)
ax.set_xlabel("t")
ax.set_ylabel("E(t) − E_exact")
ax.set_title("Energy Drift Relative to Exact (Euler excluded — unstable)")
ax.legend()
ax.grid(True, alpha=0.4)
plt.tight_layout()
plt.savefig("verlet_energy_drift.png", dpi=150)
plt.show()

# ── Figure 3: Phase portrait ───────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(xR, vR, label="RK4",             color="tab:blue",  linewidth=0.8, alpha=0.8)
ax.plot(xV, vV, label="Velocity Verlet", color="tab:green", linewidth=1.2)
ax.plot(xT, vT, label="Exact",           color="black",     linewidth=1.0, linestyle="--")
ax.set_xlabel("x")
ax.set_ylabel("v")
ax.set_title("Phase Portrait (x, v) — Harmonic Oscillator")
ax.legend()
ax.grid(True, alpha=0.4)
ax.set_aspect("equal")
plt.tight_layout()
plt.savefig("verlet_phase.png", dpi=150)
plt.show()

# ── Figure 4: Position x(t) comparison ───────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(t, xR, label="RK4",             color="tab:blue",  linewidth=1.0, alpha=0.8)
ax.plot(t, xV, label="Velocity Verlet", color="tab:green", linewidth=1.2)
ax.plot(t, xT, label="Exact",           color="black",     linewidth=0.8, linestyle="--")
ax.set_xlabel("t")
ax.set_ylabel("x(t)")
ax.set_title("Position x(t): RK4 vs Velocity Verlet (Euler excluded — unstable)")
ax.legend()
ax.grid(True, alpha=0.4)
plt.tight_layout()
plt.savefig("verlet_position.png", dpi=150)
plt.show()