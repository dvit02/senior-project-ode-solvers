import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Load CSV files exported by dp45_example
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cmake-build-debug")
dp  = pd.read_csv(os.path.join(base, "dp45_vdp.csv"))
rk4 = pd.read_csv(os.path.join(base, "rk4_adaptive_vdp.csv"))

# ── Figure 1: x(t) comparison ─────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(rk4["t"], rk4["x"], label="Adaptive RK4 (step-doubling)", color="tab:blue",  linewidth=1.2)
ax.plot(dp["t"],  dp["x"],  label="Dormand-Prince RK45",          color="tab:orange", linewidth=1.2, linestyle="--")
ax.set_xlabel("t")
ax.set_ylabel("x(t)")
ax.set_title("Van der Pol x(t): Adaptive RK4 vs Dormand-Prince RK45 (mu=10)")
ax.legend()
ax.grid(True, alpha=0.4)
plt.tight_layout()
plt.savefig("dp45_xt.png", dpi=150)
plt.show()

# ── Figure 2: phase portrait ───────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(rk4["x"], rk4["v"], label="Adaptive RK4",        color="tab:blue",   linewidth=1.0, alpha=0.8)
ax.plot(dp["x"],  dp["v"],  label="Dormand-Prince RK45",  color="tab:orange", linewidth=1.2, linestyle="--")
ax.set_xlabel("x")
ax.set_ylabel("v")
ax.set_title("Phase Portrait (x, v) — Van der Pol (mu=10)")
ax.legend()
ax.grid(True, alpha=0.4)
plt.tight_layout()
plt.savefig("dp45_phase.png", dpi=150)
plt.show()

# ── Figure 3: step size distribution side by side ─────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

ax1.hist(np.diff(rk4["t"].to_numpy()), bins=40, color="tab:blue", alpha=0.7)
ax1.set_xlabel("Step size h")
ax1.set_ylabel("Count")
ax1.set_title("Adaptive RK4: step size distribution")
ax1.grid(True, alpha=0.4)

ax2.hist(np.diff(dp["t"].to_numpy()), bins=40, color="tab:orange", alpha=0.7)
ax2.set_xlabel("Step size h")
ax2.set_ylabel("Count")
ax2.set_title("Dormand-Prince RK45: step size distribution")
ax2.grid(True, alpha=0.4)

plt.tight_layout()
plt.savefig("dp45_stepsize.png", dpi=150)
plt.show()

# ── Figure 4: step size over time ─────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 4))
t_rk4 = rk4["t"].to_numpy()
t_dp  = dp["t"].to_numpy()
ax.plot(t_rk4[:-1], np.diff(t_rk4), label="Adaptive RK4",        color="tab:blue",   linewidth=1.0, alpha=0.8)
ax.plot(t_dp[:-1],  np.diff(t_dp),  label="Dormand-Prince RK45",  color="tab:orange", linewidth=1.0, alpha=0.8)
ax.set_xlabel("t")
ax.set_ylabel("h (step size)")
ax.set_title("Step Size h(t) over Integration: Adaptive RK4 vs Dormand-Prince RK45")
ax.legend()
ax.grid(True, alpha=0.4)
plt.tight_layout()
plt.savefig("dp45_h_over_time.png", dpi=150)
plt.show()
