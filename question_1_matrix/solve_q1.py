from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"


def matrix_to_text(matrix: np.ndarray) -> str:
    return np.array2string(matrix, precision=6, suppress_small=True)


def save_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    DATA_DIR.mkdir(exist_ok=True)

    a = np.array([[5, 23, 56], [43, 9, 11], [14, 37, 216]], dtype=float)
    b = np.array([[3, 7, 4], [13, 2, 99], [61, 76, 1]], dtype=float)
    identity = np.eye(3)

    result_6a_5b = 6 * a + 5 * b
    result_a_b_i = a - b + identity
    matrix_product = a @ b
    elementwise_product = a * b

    save_csv(
        DATA_DIR / "matrix_results.csv",
        ["name", "r1c1", "r1c2", "r1c3", "r2c1", "r2c2", "r2c3", "r3c1", "r3c2", "r3c3"],
        [
            ["A", *a.ravel().tolist()],
            ["B", *b.ravel().tolist()],
            ["6A+5B", *result_6a_5b.ravel().tolist()],
            ["A-B+I", *result_a_b_i.ravel().tolist()],
            ["A@B", *matrix_product.ravel().tolist()],
            ["A_elementwise_B", *elementwise_product.ravel().tolist()],
        ],
    )

    report = [
        "题目1：矩阵运算结果",
        "",
        "A =",
        matrix_to_text(a),
        "",
        "B =",
        matrix_to_text(b),
        "",
        "6*A + 5*B =",
        matrix_to_text(result_6a_5b),
        "",
        "A - B + I =",
        matrix_to_text(result_a_b_i),
        "",
        "A*B（Python 中写作 A @ B，矩阵乘法）=",
        matrix_to_text(matrix_product),
        "",
        "A.*B（Python 中写作 A * B，对应元素相乘）=",
        matrix_to_text(elementwise_product),
        "",
        "结论：A*B 与 A.*B 不相同。A*B 是线性代数矩阵乘法，A.*B 是对应元素逐项相乘。",
    ]
    (BASE_DIR / "result_report.txt").write_text("\n".join(report), encoding="utf-8")
    print(BASE_DIR / "result_report.txt")


if __name__ == "__main__":
    main()
