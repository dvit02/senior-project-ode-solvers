import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

fixed = pd.read_csv("vdp_fixed.csv")
adap  = pd.read_csv("vdp_adaptive.csv")

# Adaptive step sizes from time grid (h_i = t_{i+1} - t_i)
tA = adap["t"].to_numpy()
hA = np.diff(tA)
tH = 0.5 * (tA[:-1] + tA[1:])  # midpoints for plotting h(t)

#
# FIGURE 1: x(t) comparison

plt.figure()
plt.plot(fixed["t"], fixed["x"], label="Fixed RK4: x(t)", linewidth=1.8)
plt.plot(adap["t"],  adap["x"],  label="Adaptive RK4: x(t)",
         linestyle="--", linewidth=2.0,
         marker="o", markersize=2.5, markevery=max(1, len(adap)//120))
plt.xlabel("t")
plt.ylabel("x")
plt.title("Van der Pol Oscillator: x(t) (Fixed vs Adaptive RK4)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("vdp_xt.png", dpi=200)
plt.show()

# FIGURE 2: phase portrait

plt.figure()
plt.plot(fixed["x"], fixed["v"], label="Fixed RK4 phase", linewidth=1.8)
plt.plot(adap["x"],  adap["v"],  label="Adaptive RK4 phase",
         linestyle="--", linewidth=2.0,
         marker="o", markersize=2.0, markevery=max(1, len(adap)//160))
plt.xlabel("x")
plt.ylabel("v")
plt.title("Van der Pol Oscillator: Phase Portrait (v vs x)")
plt.grid(True)
plt.legend()
plt.axis("equal")
plt.tight_layout()
plt.savefig("vdp_phase.png", dpi=200)
plt.show()


# FIGURE 3: adaptive step size h(t)

plt.figure()
plt.plot(tH, hA, linewidth=2.0)
plt.xlabel("t")
plt.ylabel("h (adaptive step size)")
plt.title("Adaptive RK4 Step Size Over Time (Van der Pol)")
plt.grid(True)
plt.tight_layout()
plt.savefig("vdp_h.png", dpi=200)
plt.show()

# FIGURE 4: density of steps along trajectory (nice visual)
# color = step size (smaller h => more points)


plt.figure()
plt.scatter(adap["x"][:-1], adap["v"][:-1], s=10, c=hA)  # no manual colormap for you
plt.xlabel("x")
plt.ylabel("v")
plt.title("Adaptive RK4: Step Density Along Limit Cycle (colored by h)")
plt.grid(True)
plt.axis("equal")
plt.tight_layout()
plt.savefig("vdp_phase_density.png", dpi=200)
plt.show()

print("Saved: vdp_xt.png, vdp_phase.png, vdp_h.png, vdp_phase_density.png")

