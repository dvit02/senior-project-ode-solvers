import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Read numerical solutions
import os; base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cmake-build-debug")
e = pd.read_csv(os.path.join(base, "oscillator_euler.csv"))
r = pd.read_csv(os.path.join(base, "oscillator_rk4.csv"))

# Parameters / IC (must match your C++ run)
omega = 1.0
x0 = 1.0
v0 = 0.0

# Use RK4 time grid for analytical (same grid anyway)
t = r["t"].to_numpy()

# Analytical solution
x_true = x0 * np.cos(omega * t) + (v0 / omega) * np.sin(omega * t)
v_true = -x0 * omega * np.sin(omega * t) + v0 * np.cos(omega * t)
E_true = 0.5 * (v_true**2 + (omega**2) * (x_true**2))

xe = e["x"].to_numpy()
ve = e["v"].to_numpy()
Ee = e["E"].to_numpy()

xr = r["x"].to_numpy()
vr = r["v"].to_numpy()
Er = r["E"].to_numpy()

# 1) Position x(t)

plt.figure()
plt.plot(e["t"], xe, label="Euler x(t)")
plt.plot(r["t"], xr, label="RK4 x(t)",
         linestyle="--", linewidth=2.0,
         marker="o", markersize=3, markevery=40, zorder=5)
plt.plot(t, x_true, label="Analytical x(t)", linewidth=1.5, alpha=0.7)
plt.xlabel("t")
plt.ylabel("x")
plt.legend()
plt.title("Harmonic Oscillator: Position")
plt.grid(True)
plt.show()

# 2) Velocity v(t)
plt.figure()
plt.plot(e["t"], ve, label="Euler v(t)")
plt.plot(r["t"], vr, label="RK4 v(t)",
         linestyle="--", linewidth=2.0,
         marker="o", markersize=3, markevery=40, zorder=5)
plt.plot(t, v_true, label="Analytical v(t)", linewidth=1.5, alpha=0.7)
plt.xlabel("t")
plt.ylabel("v")
plt.legend()
plt.title("Harmonic Oscillator: Velocity")
plt.grid(True)
plt.show()

# 3) Energy E(t)
plt.figure()
plt.plot(e["t"], Ee, label="Euler E(t)")
plt.plot(r["t"], Er, label="RK4 E(t)",
         linestyle="--", linewidth=2.0,
         marker="o", markersize=3, markevery=40, zorder=5)
plt.plot(t, E_true, label="Analytical E(t)", linewidth=1.5, alpha=0.7)
plt.xlabel("t")
plt.ylabel("E")
plt.legend()
plt.title("Harmonic Oscillator: Energy Drift")
plt.grid(True)
plt.show()

# 4) Error plots (log scale)

err_x_e = np.abs(xe - x_true)
err_x_r = np.abs(xr - x_true)

plt.figure()
plt.semilogy(t, err_x_e, label="|Euler - Analytical| (x)")
plt.semilogy(t, err_x_r, label="|RK4 - Analytical| (x)")
plt.xlabel("t")
plt.ylabel("absolute error (log scale)")
plt.legend()
plt.title("Position error vs Analytical")
plt.grid(True)
plt.show()

err_v_e = np.abs(ve - v_true)
err_v_r = np.abs(vr - v_true)

plt.figure()
plt.semilogy(t, err_v_e, label="|Euler - Analytical| (v)")
plt.semilogy(t, err_v_r, label="|RK4 - Analytical| (v)")
plt.xlabel("t")
plt.ylabel("absolute error (log scale)")
plt.legend()
plt.title("Velocity error vs Analytical")
plt.grid(True)
plt.show()

# 5) Energy drift relative to true energy
plt.figure()
plt.plot(t, Ee - E_true, label="Euler: E - E_true")
plt.plot(t, Er - E_true, label="RK4: E - E_true")
plt.xlabel("t")
plt.ylabel("Energy error")
plt.legend()
plt.title("Energy drift relative to analytical energy")
plt.grid(True)
plt.show()

# 6) Phase portrait (x vs v)

plt.figure()
plt.plot(xe, ve, label="Euler phase")
plt.plot(xr, vr, label="RK4 phase")
plt.plot(x_true, v_true, label="Analytical phase")
plt.xlabel("x")
plt.ylabel("v")
plt.legend()
plt.title("Phase portrait: (x, v)")
plt.grid(True)
plt.axis("equal")
plt.show()