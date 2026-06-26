# 题目 3：缓变基区少子浓度分布说明

本文档解释 `solve_q3.py` 的计算目的、物理意义、公式含义、程序流程，以及生成图片和 CSV 数据的作用。

## 1. 题目要解决什么

题目 3 研究的是双极型晶体管在有源放大区工作时，基区中的少子电子浓度分布。

对 NPN 晶体管来说，基区通常是 P 型区，电子是基区中的少数载流子。发射结正偏时，电子从发射区注入基区；集电结反偏时，电子被集电区收集。于是基区中会形成一个随位置 $x$ 变化的少子浓度分布：

$$
n_B(x)
$$

本题要画出不同缓变杂质分布下的 $n_B(x)$ 曲线，并解释杂质分布参数对少子浓度分布的影响。

## 2. 程序使用的公式

程序采用文档中的少子浓度分布公式：

$$
n_B(x)=\frac{I_{nE}}{D_{nB}q}\frac{\int_x^{W_B}N(t)dt}{N(x)}
$$

其中：

- $n_B(x)$：基区位置 $x$ 处的少子电子浓度。
- $I_{nE}$：电子电流。
- $D_{nB}$：基区电子扩散系数。
- $q$：元电荷。
- $W_B$：基区宽度。
- $N(x)$：基区杂质浓度分布。

程序对应函数是：

```python
def minority_profile(x, impurity, current_a, dn, q):
    cumulative = cumulative_trapezoid(impurity, x, initial=0.0)
    integral_to_end = cumulative[-1] - cumulative
    return current_a / (dn * q) * integral_to_end / impurity
```

这里：

- `impurity` 对应 $N(x)$。
- `cumulative_trapezoid(impurity, x, initial=0.0)` 对 $N(t)$ 做累计积分。
- `integral_to_end` 对应 $\int_x^{W_B}N(t)dt$。
- 最后一行对应完整的 $n_B(x)$ 公式。

## 3. 为什么要做数值积分

公式中有积分：

$$
\int_x^{W_B}N(t)dt
$$

如果 $N(x)$ 很简单，可以手算。但程序同时考虑线性分布和指数分布，用数值积分更通用。

程序使用 SciPy 的：

```python
cumulative_trapezoid(...)
```

它的作用是用梯形法近似积分。对于一组离散点 $x_0,x_1,\ldots,x_n$，程序先算从左端到每个点的累计积分，再用：

```python
integral_to_end = cumulative[-1] - cumulative
```

得到从当前位置 $x$ 到基区末端 $W_B$ 的积分。

## 4. 程序使用的参数

程序中的参数为：

```python
q = 1.6e-19
in_e = 0.01e-3
dn_b = 2.0
wb_cm = 0.05e-4
na0 = 1.0e17
eta_values = [1.0, 2.0, 4.0]
```

对应物理量：

$$
q=1.6\times 10^{-19}\ \text{C}
$$

$$
I_{nE}=0.01\ \text{mA}=1.0\times 10^{-5}\ \text{A}
$$

$$
D_{nB}=2\ \text{cm}^2/\text{s}
$$

$$
W_B=0.05\ \mu m=0.05\times 10^{-4}\ \text{cm}
$$

$$
N_A(0)=1.0\times 10^{17}\ \text{cm}^{-3}
$$

## 5. 杂质分布模型

程序绘制了两类杂质分布。

### 5.1 线性缓变分布

线性分布为：

$$
N_A(x)=N_A(0)\left(1-\frac{x}{W_B}\right)
$$

当 $x=0$ 时，杂质浓度最大：

$$
N_A(0)
$$

当 $x=W_B$ 时，公式给出：

$$
N_A(W_B)=0
$$

由于分母中有 $N(x)$，在 $x=W_B$ 处会出现除以零的问题。因此程序没有取到正好 $W_B$，而是取到：

```python
u = np.linspace(0.0, 0.999, 600)
```

也就是最大位置为：

$$
x=0.999W_B
$$

### 5.2 指数缓变分布

指数分布为：

$$
N_A(x)=N_A(0)e^{-\eta x/W_B}
$$

程序中取：

$$
\eta=1,2,4
$$

$\eta$ 越大，杂质浓度沿基区下降得越快。由于 $n_B(x)$ 的计算公式中既有积分项，也有除以 $N(x)$ 的项，因此 $\eta$ 会明显改变少子浓度曲线形状。

## 6. 程序生成了哪些结果

程序生成两个图片和一个 CSV 文件。

### 6.1 归一化少子浓度图

文件：

```text
images/minority_distribution.png
```

纵坐标是：

$$
\frac{n_B(x)}{n_{B,\max}}
$$

也就是每条曲线除以自身最大值后的结果。这样做的优点是可以更清楚地比较不同杂质分布对曲线形状的影响。

### 6.2 绝对少子浓度图

文件：

```text
images/minority_distribution_absolute.png
```

纵坐标是：

$$
n_B(x)
$$

单位是：

$$
\text{cm}^{-3}
$$

这张图用于观察少子浓度的实际数量级。

### 6.3 CSV 数据

文件：

```text
data/minority_distribution.csv
```

列含义如下：

| 列名 | 含义 |
|---|---|
| `x_um` | 位置 $x$，单位为 $\mu m$ |
| `x_over_WB` | 归一化位置 $x/W_B$ |
| `NA_cm-3` | 杂质浓度 $N_A(x)$ |
| `nB_cm-3` | 少子电子浓度绝对值 $n_B(x)$ |
| `nB_normalized` | 归一化少子浓度 |
| `eta_0_means_linear` | 指数分布参数；`0` 表示线性分布 |

## 7. 图像应该怎么理解

从图像中可以看出：

1. 在 $x=W_B$ 附近，积分区间 $\int_x^{W_B}$ 变短，所以 $n_B(x)$ 逐渐趋近于 0。
2. 不同杂质分布会改变 $N(x)$ 和积分项，因此少子浓度曲线的弯曲程度不同。
3. 指数分布中 $\eta$ 越大，杂质浓度下降越快，少子浓度分布变化也越明显。
4. 归一化图适合比较形状，绝对值图适合报告数量级。

## 8. 程序整体流程

`solve_q3.py` 的整体流程为：

```text
1. 设置物理参数 q, InE, DnB, WB, NA0
2. 生成位置数组 x
3. 构造线性杂质分布 NA(x)
4. 构造多个指数杂质分布 NA(x)
5. 对每一种分布计算 nB(x)
6. 保存归一化图
7. 保存绝对浓度图
8. 保存 CSV 数据
9. 生成 result_report.txt
```

## 9. 可以写进实验报告的说明

可以把下面这段放进实验报告：

> 本题根据缓变基区少子浓度分布公式 $n_B(x)=\frac{I_{nE}}{D_{nB}q}\frac{\int_x^{W_B}N(t)dt}{N(x)}$ 进行计算。程序分别设置线性杂质分布 $N_A(x)=N_A(0)(1-x/W_B)$ 和指数杂质分布 $N_A(x)=N_A(0)e^{-\eta x/W_B}$，通过梯形积分方法求得不同位置的 $n_B(x)$。结果表明，杂质浓度分布会明显影响少子浓度曲线形状；在接近集电结边界 $x=W_B$ 时，由于积分区间趋近于零，少子浓度也趋近于零。

