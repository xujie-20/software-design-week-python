from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq


BASE_DIR = Path(__file__).resolve().parent
IMAGE_DIR = BASE_DIR / "images"
DATA_DIR = BASE_DIR / "data"


class SiliconConstants:
    q = 1.602176634e-19
    k = 1.380649e-23
    t = 300.0
    eps0 = 8.8541878128e-14
    eps_si = 11.7 * eps0
    ni = 1.45e10

    @property
    def vt(self) -> float:
        return self.k * self.t / self.q


SI = SiliconConstants()


def setup_plot_style() -> None:
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial Unicode MS", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.dpi"] = 140


def abrupt_vbi(na: float, nd: float) -> float:
    return SI.vt * math.log(na * nd / SI.ni**2)


def abrupt_widths(na: float, nd: float, vr: float) -> tuple[float, float, float, float]:
    vbi = abrupt_vbi(na, nd)
    v = vbi + vr
    xn = math.sqrt(2 * SI.eps_si / SI.q * na / (nd * (na + nd)) * v)
    xp = math.sqrt(2 * SI.eps_si / SI.q * nd / (na * (na + nd)) * v)
    n0 = na * nd / (na + nd)
    emax = math.sqrt(2 * SI.q * n0 * v / SI.eps_si)
    return xn, xp, emax, vbi


def abrupt_field_profile(na: float, nd: float, vr: float, n: int = 900) -> tuple[np.ndarray, np.ndarray, float, float, float]:
    xn, xp, emax, _ = abrupt_widths(na, nd, vr)
    x = np.linspace(-xn, xp, n)
    e = np.where(x <= 0.0, SI.q * nd / SI.eps_si * (x + xn), SI.q * na / SI.eps_si * (xp - x))
    return x, np.maximum(e, 0.0), xn, xp, emax


def graded_vbi(a_grad: float) -> float:
    def f(vbi: float) -> float:
        x0 = (12 * SI.eps_si * vbi / (SI.q * a_grad)) ** (1.0 / 3.0)
        return vbi - SI.vt * math.log((a_grad * x0 / (2 * SI.ni)) ** 2)

    return brentq(f, 0.01, 2.0)


def graded_profile(a_grad: float, vr: float, n: int = 900) -> tuple[np.ndarray, np.ndarray, float, float, float]:
    vbi = graded_vbi(a_grad)
    xj = 0.5 * (12 * SI.eps_si * (vbi + vr) / (SI.q * a_grad)) ** (1.0 / 3.0)
    x = np.linspace(-xj, xj, n)
    emax = SI.q * a_grad / (2 * SI.eps_si) * xj**2
    e = emax * (1.0 - (x / xj) ** 2)
    return x, np.maximum(e, 0.0), xj, emax, vbi


def ionization_rate(e: np.ndarray | float, a_const: float, b_const: float, m: float = 1.0) -> np.ndarray:
    e_arr = np.asarray(e, dtype=float)
    alpha = np.zeros_like(e_arr)
    mask = e_arr > 1.0
    alpha[mask] = a_const * np.exp(-((b_const / e_arr[mask]) ** m))
    return alpha


def main() -> None:
    setup_plot_style()
    IMAGE_DIR.mkdir(exist_ok=True)
    DATA_DIR.mkdir(exist_ok=True)

    na = 1.0e17
    nd = 1.0e16
    a_grad = 1.0e20
    vr_sample = 40.0
    vr_values = np.linspace(0.0, 120.0, 121)
    field_vr_values = [0.0, 20.0, 40.0, 80.0, 120.0]

    x_abrupt, e_abrupt, _, _, emax_abrupt = abrupt_field_profile(na, nd, vr_sample)
    x_graded, e_graded, _, emax_graded, vbi_graded = graded_profile(a_grad, vr_sample)
    vbi_abrupt = abrupt_vbi(na, nd)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(x_abrupt / 1e-4, e_abrupt / 1e5, label=rf"突变结，$V_R$={vr_sample:g} V")
    ax.plot(x_graded / 1e-4, e_graded / 1e5, label=rf"线性缓变结，$V_R$={vr_sample:g} V")
    ax.set_xlabel(r"$x$ / $\mu m$")
    ax.set_ylabel(r"$E$ / $10^5$ V cm$^{-1}$")
    ax.set_title("题目4：势垒区电场分布")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(IMAGE_DIR / "field_distribution.png")
    plt.close(fig)

    fig, (ax_abrupt, ax_graded) = plt.subplots(1, 2, figsize=(12.0, 4.8), sharey=True)
    for vr in field_vr_values:
        x_multi, e_multi, *_ = abrupt_field_profile(na, nd, vr)
        ax_abrupt.plot(x_multi / 1e-4, e_multi / 1e5, label=f"{vr:g} V")

        x_multi_g, e_multi_g, *_ = graded_profile(a_grad, vr)
        ax_graded.plot(x_multi_g / 1e-4, e_multi_g / 1e5, label=f"{vr:g} V")

    ax_abrupt.set_title("突变结")
    ax_abrupt.set_xlabel(r"$x$ / $\mu m$")
    ax_abrupt.set_ylabel(r"$E$ / $10^5$ V cm$^{-1}$")
    ax_abrupt.grid(True, alpha=0.3)
    ax_abrupt.legend(title=r"$V_R$", fontsize=8)

    ax_graded.set_title("线性缓变结")
    ax_graded.set_xlabel(r"$x$ / $\mu m$")
    ax_graded.grid(True, alpha=0.3)
    ax_graded.legend(title=r"$V_R$", fontsize=8)

    fig.suptitle("题目4：多个反偏电压下的势垒区电场分布")
    fig.tight_layout()
    fig.savefig(IMAGE_DIR / "field_distribution_multi_voltage.png")
    plt.close(fig)

    alpha_n_a = 7.03e5
    alpha_n_b = 1.23e6
    alpha_p_a = 1.58e6
    alpha_p_b = 2.03e6

    rows: list[list[float]] = []
    for vr in vr_values:
        _, _, emax, _ = abrupt_widths(na, nd, float(vr))
        _, _, _, emax_g, _ = graded_profile(a_grad, float(vr), n=101)
        rows.append(
            [
                float(vr),
                emax,
                ionization_rate(emax, alpha_n_a, alpha_n_b).item(),
                ionization_rate(emax, alpha_p_a, alpha_p_b).item(),
                emax_g,
                ionization_rate(emax_g, alpha_n_a, alpha_n_b).item(),
                ionization_rate(emax_g, alpha_p_a, alpha_p_b).item(),
            ]
        )

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.semilogy(vr_values, [row[2] for row in rows], label="电子，突变结")
    ax.semilogy(vr_values, [row[3] for row in rows], label="空穴，突变结")
    ax.semilogy(vr_values, [row[5] for row in rows], "--", label="电子，线性缓变结")
    ax.semilogy(vr_values, [row[6] for row in rows], "--", label="空穴，线性缓变结")
    ax.set_xlabel(r"$V_R$ / V")
    ax.set_ylabel(r"碰撞电离率 / cm$^{-1}$")
    ax.set_title("题目4：碰撞电离率随反偏电压变化")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(IMAGE_DIR / "ionization_vs_voltage.png")
    plt.close(fig)

    with (DATA_DIR / "ionization_vs_voltage.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["VR_V", "Emax_abrupt_V_per_cm", "alpha_n_abrupt_cm-1", "alpha_p_abrupt_cm-1", "Emax_graded_V_per_cm", "alpha_n_graded_cm-1", "alpha_p_graded_cm-1"])
        writer.writerows(rows)

    report = [
        "题目4：PN 结势垒区电场分布和碰撞电离率",
        "",
        "硅材料参数：T=300 K, ni=1.45e10 cm^-3, eps_si=11.7eps0",
        f"突变结参数：NA={na:.2e} cm^-3, ND={nd:.2e} cm^-3, Vbi={vbi_abrupt:.4f} V",
        f"线性缓变结参数：a={a_grad:.2e} cm^-4, Vbi={vbi_graded:.4f} V",
        f"VR={vr_sample:g} V 时：突变结 Emax={emax_abrupt:.3e} V/cm，线性缓变结 Emax={emax_graded:.3e} V/cm",
        f"多电压电场图使用 VR={', '.join(f'{vr:g}' for vr in field_vr_values)} V。",
        "图片：images/field_distribution.png, images/field_distribution_multi_voltage.png, images/ionization_vs_voltage.png",
        "数据：data/ionization_vs_voltage.csv",
    ]
    (BASE_DIR / "result_report.txt").write_text("\n".join(report), encoding="utf-8")
    print(BASE_DIR / "result_report.txt")


if __name__ == "__main__":
    main()
