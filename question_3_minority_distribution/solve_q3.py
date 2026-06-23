from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid


BASE_DIR = Path(__file__).resolve().parent
IMAGE_DIR = BASE_DIR / "images"
DATA_DIR = BASE_DIR / "data"


def setup_plot_style() -> None:
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial Unicode MS", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.dpi"] = 140


def minority_profile(x: np.ndarray, impurity: np.ndarray, current_a: float, dn: float, q: float) -> np.ndarray:
    cumulative = cumulative_trapezoid(impurity, x, initial=0.0)
    integral_to_end = cumulative[-1] - cumulative
    return current_a / (dn * q) * integral_to_end / impurity


def main() -> None:
    setup_plot_style()
    IMAGE_DIR.mkdir(exist_ok=True)
    DATA_DIR.mkdir(exist_ok=True)

    q = 1.6e-19
    in_e = 0.01e-3
    dn_b = 2.0
    wb_cm = 0.05e-4
    na0 = 1.0e17
    eta_values = [1.0, 2.0, 4.0]

    u = np.linspace(0.0, 0.999, 600)
    x = u * wb_cm
    na_linear = na0 * (1.0 - u)
    n_linear = minority_profile(x, na_linear, in_e, dn_b, q)

    curves: list[tuple[str, np.ndarray, np.ndarray]] = [
        (r"$N_A(x)=N_A(0)(1-x/W_B)$", na_linear, n_linear)
    ]
    for eta in eta_values:
        na_exp = na0 * np.exp(-eta * u)
        n_exp = minority_profile(x, na_exp, in_e, dn_b, q)
        curves.append((rf"$N_A(x)=N_A(0)e^{{-{eta:g}x/W_B}}$", na_exp, n_exp))

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for label, _, nb in curves:
        ax.plot(x / 1e-4, nb / np.nanmax(nb), label=label)

    csv_rows: list[list[float]] = []
    for curve_index, (label, na_values, nb_values) in enumerate(curves):
        eta = 0.0 if curve_index == 0 else eta_values[curve_index - 1]
        for xi, ui, nai, ni in zip(x, u, na_values, nb_values):
            csv_rows.append([float(xi / 1e-4), float(ui), float(nai), float(ni), float(ni / np.nanmax(nb_values)), float(eta)])

    ax.set_xlabel(r"$x$ / $\mu m$")
    ax.set_ylabel("normalized minority electron concentration")
    ax.set_title("题目3：缓变基区少子浓度分布")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(IMAGE_DIR / "minority_distribution.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for label, _, nb in curves:
        ax.plot(x / 1e-4, nb, label=label)
    ax.set_xlabel(r"$x$ / $\mu m$")
    ax.set_ylabel(r"$n_B(x)$ / cm$^{-3}$")
    ax.set_title("题目3：缓变基区少子浓度绝对值")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(IMAGE_DIR / "minority_distribution_absolute.png")
    plt.close(fig)

    with (DATA_DIR / "minority_distribution.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x_um", "x_over_WB", "NA_cm-3", "nB_cm-3", "nB_normalized", "eta_0_means_linear"])
        writer.writerows(csv_rows)

    report = [
        "题目3：双极型晶体管缓变基区少子浓度分布",
        "",
        "计算式：nB(x)=InE/(DnB*q)*Integral[x,WB](N(t)dt)/N(x)",
        "参数：InE=0.01 mA, DnB=2 cm^2/s, WB=0.05 um, q=1.6e-19 C",
        "线性分布在 x=WB 处 NA=0，程序绘制到 0.999*WB。",
        "归一化浓度图片：images/minority_distribution.png",
        "绝对浓度图片：images/minority_distribution_absolute.png",
        "数据：data/minority_distribution.csv",
    ]
    (BASE_DIR / "result_report.txt").write_text("\n".join(report), encoding="utf-8")
    print(BASE_DIR / "result_report.txt")


if __name__ == "__main__":
    main()
