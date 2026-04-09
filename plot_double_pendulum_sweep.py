from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def find_csv(root: Path) -> Path:
    candidates = [
        root / "outputs" / "double_pendulum" / "dpend_sweep.csv",
        root / "cmake-build-debug" / "dpend_sweep.csv",
        root / "cmake-build-debug-coverage" / "dpend_sweep.csv",
        root / "cmake-build-release" / "dpend_sweep.csv",
        root / "build" / "dpend_sweep.csv",
        root / "dpend_sweep.csv",
        ]

    for path in candidates:
        if path.exists():
            return path

    raise FileNotFoundError(
        "Could not find dpend_sweep.csv.\n"
        "Run double_pendulum_sweep first, or move the CSV into one of these locations:\n"
        + "\n".join(str(p) for p in candidates)
    )


def main() -> None:
    root = Path(__file__).resolve().parent
    csv_path = find_csv(root)
    df = pd.read_csv(csv_path)

    out_dir = root

    print(f"Loaded data from: {csv_path}")
    print(f"Saving plots to: {out_dir}")
    print("CSV columns:", list(df.columns))

    required = {"method", "tol_or_h", "drift", "rhs_calls", "steps"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in CSV: {sorted(missing)}")

    for col in ["tol_or_h", "drift", "rhs_calls", "steps"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["method", "tol_or_h", "drift", "rhs_calls", "steps"])

    plt.figure(figsize=(8, 5))
    for method, group in df.groupby("method"):
        group = group.sort_values("tol_or_h")
        plt.loglog(group["tol_or_h"], group["rhs_calls"], marker="o", label=method)

    plt.xlabel("tol_or_h")
    plt.ylabel("rhs_calls")
    plt.title("Double Pendulum Sweep: RHS Calls vs tol_or_h")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / "dpend_rhs_calls_vs_tol_or_h.png", dpi=150)

    plt.figure(figsize=(8, 5))
    for method, group in df.groupby("method"):
        group = group.sort_values("tol_or_h")
        plt.loglog(group["tol_or_h"], group["drift"].abs(), marker="o", label=method)

    plt.xlabel("tol_or_h")
    plt.ylabel("|drift|")
    plt.title("Double Pendulum Sweep: Energy Drift vs tol_or_h")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / "dpend_drift_vs_tol_or_h.png", dpi=150)

    plt.figure(figsize=(8, 5))
    for method, group in df.groupby("method"):
        group = group.sort_values("rhs_calls")
        plt.loglog(group["rhs_calls"], group["drift"].abs(), marker="o", label=method)

    plt.xlabel("rhs_calls")
    plt.ylabel("|drift|")
    plt.title("Double Pendulum Sweep: Work-Precision")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / "dpend_work_precision.png", dpi=150)

    plt.figure(figsize=(8, 5))
    for method, group in df.groupby("method"):
        group = group.sort_values("tol_or_h")
        plt.loglog(group["tol_or_h"], group["steps"], marker="o", label=method)

    plt.xlabel("tol_or_h")
    plt.ylabel("steps")
    plt.title("Double Pendulum Sweep: Steps vs tol_or_h")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / "dpend_steps_vs_tol_or_h.png", dpi=150)

    plt.show()


if __name__ == "__main__":
    main()