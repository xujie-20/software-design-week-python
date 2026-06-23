from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


BASE_DIR = Path(__file__).resolve().parent
IMAGE_DIR = BASE_DIR / "images"
DATA_DIR = BASE_DIR / "data"


def setup_plot_style() -> None:
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial Unicode MS", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.dpi"] = 140


def circle_points(radius: float = 1.0, center: tuple[float, float] = (0.0, 0.0), n: int = 500) -> tuple[np.ndarray, np.ndarray]:
    theta = np.linspace(0.0, 2.0 * np.pi, n)
    x = center[0] + radius * np.cos(theta)
    y = center[1] + radius * np.sin(theta)
    return x, y


def main() -> None:
    setup_plot_style()
    IMAGE_DIR.mkdir(exist_ok=True)
    DATA_DIR.mkdir(exist_ok=True)

    x, y = circle_points(radius=1.0)

    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    ax.plot(x, y, linewidth=2.2)
    ax.scatter([0], [0], s=20)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("题目2：函数绘制圆形")
    fig.tight_layout()
    fig.savefig(IMAGE_DIR / "circle.png")
    plt.close(fig)

    with (DATA_DIR / "circle_points.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x", "y"])
        writer.writerows([[float(a), float(b)] for a, b in zip(x, y)])

    report = [
        "题目2：绘制圆形",
        "",
        "本题使用 circle_points(radius, center, n) 函数生成圆上采样点。",
        "圆心：(0, 0)",
        "半径：1",
        "图片：images/circle.png",
        "数据：data/circle_points.csv",
    ]
    (BASE_DIR / "result_report.txt").write_text("\n".join(report), encoding="utf-8")
    print(BASE_DIR / "result_report.txt")


if __name__ == "__main__":
    main()
