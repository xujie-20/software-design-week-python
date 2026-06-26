from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid
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


def abrupt_field_profile(na: float, nd: float, vr: float, n: int = 1400) -> tuple[np.ndarray, np.ndarray]:
    xn, xp, _, _ = abrupt_widths(na, nd, vr)
    x = np.linspace(-xn, xp, n)
    e = np.where(x <= 0.0, SI.q * nd / SI.eps_si * (x + xn), SI.q * na / SI.eps_si * (xp - x))
    return x, np.maximum(e, 0.0)


def graded_vbi(a_grad: float) -> float:
    def f(vbi: float) -> float:
        x0 = (12 * SI.eps_si * vbi / (SI.q * a_grad)) ** (1.0 / 3.0)
        return vbi - SI.vt * math.log((a_grad * x0 / (2 * SI.ni)) ** 2)

    return brentq(f, 0.01, 2.0)


def graded_field_profile(a_grad: float, vr: float, n: int = 1400) -> tuple[np.ndarray, np.ndarray]:
    vbi = graded_vbi(a_grad)
    xj = 0.5 * (12 * SI.eps_si * (vbi + vr) / (SI.q * a_grad)) ** (1.0 / 3.0)
    x = np.linspace(-xj, xj, n)
    emax = SI.q * a_grad / (2 * SI.eps_si) * xj**2
    e = emax * (1.0 - (x / xj) ** 2)
    return x, np.maximum(e, 0.0)


def ionization_rate(e: np.ndarray | float, a_const: float, b_const: float, m: float = 1.0) -> np.ndarray:
    e_arr = np.asarray(e, dtype=float)
    alpha = np.zeros_like(e_arr)
    mask = e_arr > 1.0
    alpha[mask] = a_const * np.exp(-((b_const / e_arr[mask]) ** m))
    return alpha


def multiplication_factors_from_field(x: np.ndarray, e: np.ndarray) -> tuple[float, float]:
    alpha_n = ionization_rate(e, 7.03e5, 1.23e6)
    alpha_p = ionization_rate(e, 1.58e6, 2.03e6)

    diff_np = alpha_n - alpha_p
    inner_n = cumulative_trapezoid(diff_np, x, initial=0.0)
    term_n = alpha_n * np.exp(np.clip(-inner_n, -100.0, 100.0))
    integral_n = np.trapezoid(term_n, x)

    diff_pn = alpha_p - alpha_n
    cumulative_p = cumulative_trapezoid(diff_pn, x, initial=0.0)
    inner_p = cumulative_p[-1] - cumulative_p
    term_p = alpha_p * np.exp(np.clip(-inner_p, -100.0, 100.0))
    integral_p = np.trapezoid(term_p, x)

    mn = math.inf if 1.0 - integral_n <= 0.0 else 1.0 / (1.0 - integral_n)
    mp = math.inf if 1.0 - integral_p <= 0.0 else 1.0 / (1.0 - integral_p)
    return mn, mp


def main() -> None:
    setup_plot_style()
    IMAGE_DIR.mkdir(exist_ok=True)
    DATA_DIR.mkdir(exist_ok=True)

    na = 1.0e17
    nd = 1.0e16
    a_grad = 1.0e20
    vr_values = np.linspace(0.0, 120.0, 121)
    rows: list[list[float]] = []
    mn_abrupt_values: list[float] = []
    mp_abrupt_values: list[float] = []
    mn_graded_values: list[float] = []
    mp_graded_values: list[float] = []

    for vr in vr_values:
        x_abrupt, e_abrupt = abrupt_field_profile(na, nd, float(vr))
        mn_abrupt, mp_abrupt = multiplication_factors_from_field(x_abrupt, e_abrupt)
        mn_abrupt_values.append(mn_abrupt)
        mp_abrupt_values.append(mp_abrupt)

        x_graded, e_graded = graded_field_profile(a_grad, float(vr))
        mn_graded, mp_graded = multiplication_factors_from_field(x_graded, e_graded)
        mn_graded_values.append(mn_graded)
        mp_graded_values.append(mp_graded)

        rows.append(
            [
                float(vr),
                mn_abrupt if math.isfinite(mn_abrupt) else np.nan,
                mp_abrupt if math.isfinite(mp_abrupt) else np.nan,
                mn_graded if math.isfinite(mn_graded) else np.nan,
                mp_graded if math.isfinite(mp_graded) else np.nan,
            ]
        )

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    mn_abrupt_plot = np.array([value if math.isfinite(value) and value < 100 else np.nan for value in mn_abrupt_values])
    mp_abrupt_plot = np.array([value if math.isfinite(value) and value < 100 else np.nan for value in mp_abrupt_values])
    mn_graded_plot = np.array([value if math.isfinite(value) and value < 100 else np.nan for value in mn_graded_values])
    mp_graded_plot = np.array([value if math.isfinite(value) and value < 100 else np.nan for value in mp_graded_values])
    ax.plot(vr_values, mn_abrupt_plot, label=r"$M_n$，突变结")
    ax.plot(vr_values, mp_abrupt_plot, label=r"$M_p$，突变结")
    ax.plot(vr_values, mn_graded_plot, "--", label=r"$M_n$，线性缓变结")
    ax.plot(vr_values, mp_graded_plot, "--", label=r"$M_p$，线性缓变结")
    ax.set_xlabel(r"$V_R$ / V")
    ax.set_ylabel("雪崩倍增因子")
    ax.set_title("题目5：雪崩倍增因子随反偏电压变化")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(IMAGE_DIR / "avalanche_multiplication.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    mn_abrupt_semilog = np.array([value if math.isfinite(value) and value > 0 else np.nan for value in mn_abrupt_values])
    mp_abrupt_semilog = np.array([value if math.isfinite(value) and value > 0 else np.nan for value in mp_abrupt_values])
    mn_graded_semilog = np.array([value if math.isfinite(value) and value > 0 else np.nan for value in mn_graded_values])
    mp_graded_semilog = np.array([value if math.isfinite(value) and value > 0 else np.nan for value in mp_graded_values])
    ax.semilogy(vr_values, mn_abrupt_semilog, label=r"$M_n$，突变结")
    ax.semilogy(vr_values, mp_abrupt_semilog, label=r"$M_p$，突变结")
    ax.semilogy(vr_values, mn_graded_semilog, "--", label=r"$M_n$，线性缓变结")
    ax.semilogy(vr_values, mp_graded_semilog, "--", label=r"$M_p$，线性缓变结")
    ax.set_xlabel(r"$V_R$ / V")
    ax.set_ylabel("雪崩倍增因子")
    ax.set_title("题目5：雪崩倍增因子半对数图")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(IMAGE_DIR / "avalanche_multiplication_semilog.png")
    plt.close(fig)

    with (DATA_DIR / "avalanche_multiplication.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["VR_V", "Mn_abrupt", "Mp_abrupt", "Mn_graded", "Mp_graded"])
        writer.writerows(rows)

    finite_mn_abrupt = [(vr, mn) for vr, mn in zip(vr_values, mn_abrupt_values) if math.isfinite(mn)]
    finite_mp_abrupt = [(vr, mp) for vr, mp in zip(vr_values, mp_abrupt_values) if math.isfinite(mp)]
    finite_mn_graded = [(vr, mn) for vr, mn in zip(vr_values, mn_graded_values) if math.isfinite(mn)]
    finite_mp_graded = [(vr, mp) for vr, mp in zip(vr_values, mp_graded_values) if math.isfinite(mp)]
    linear_mn_abrupt = [(vr, mn) for vr, mn in finite_mn_abrupt if mn < 100]
    linear_mp_abrupt = [(vr, mp) for vr, mp in finite_mp_abrupt if mp < 100]
    linear_mn_graded = [(vr, mn) for vr, mn in finite_mn_graded if mn < 100]
    linear_mp_graded = [(vr, mp) for vr, mp in finite_mp_graded if mp < 100]

    report = [
        "题目5：雪崩倍增因子",
        "",
        "方法：按局域碰撞电离模型，分别对突变结和线性缓变结耗尽区电场进行数值积分。",
        f"突变结参数：NA={na:.2e} cm^-3, ND={nd:.2e} cm^-3",
        f"线性缓变结参数：a={a_grad:.2e} cm^-4",
        f"突变结在 0-{vr_values[-1]:.0f} V 范围内，可计算的最大 Mn≈{max(v for _, v in finite_mn_abrupt):.4g}，最大 Mp≈{max(v for _, v in finite_mp_abrupt):.4g}。",
        f"线性缓变结在 0-{vr_values[-1]:.0f} V 范围内，可计算的最大 Mn≈{max(v for _, v in finite_mn_graded):.4g}，最大 Mp≈{max(v for _, v in finite_mp_graded):.4g}。",
        f"线性图只显示 M<100 的有限区域：突变结 Mn 约显示到 {linear_mn_abrupt[-1][0]:.0f} V，Mp 约显示到 {linear_mp_abrupt[-1][0]:.0f} V；线性缓变结 Mn 约显示到 {linear_mn_graded[-1][0]:.0f} V，Mp 约显示到 {linear_mp_graded[-1][0]:.0f} V。",
        "图片：images/avalanche_multiplication.png, images/avalanche_multiplication_semilog.png",
        "数据：data/avalanche_multiplication.csv",
    ]
    (BASE_DIR / "result_report.txt").write_text("\n".join(report), encoding="utf-8")
    print(BASE_DIR / "result_report.txt")


if __name__ == "__main__":
    main()
