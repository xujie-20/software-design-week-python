# 题目 4：PN 结电场分布和碰撞电离率说明

本文档解释 `solve_q4.py` 的作用、物理模型、公式含义、程序流程，以及生成图片和 CSV 数据的意义。

## 1. 题目要解决什么

题目 4 要求：

> 确定 PN 结势垒区内电场分布和碰撞电离率随反偏电压的变化关系。

程序做了两件事：

1. 计算 PN 结耗尽区中的电场分布 $E(x)$。
2. 根据电场计算电子和空穴的碰撞电离率 $\alpha_n$、$\alpha_p$，并观察它们随反偏电压 $V_R$ 的变化。

当前程序同时考虑两种结：

- 突变 PN 结。
- 线性缓变结。

## 2. 程序使用的硅材料参数

程序中定义了一个 `SiliconConstants` 类：

```python
class SiliconConstants:
    q = 1.602176634e-19
    k = 1.380649e-23
    t = 300.0
    eps0 = 8.8541878128e-14
    eps_si = 11.7 * eps0
    ni = 1.45e10
```

含义如下：

| 参数 | 含义 |
|---|---|
| $q$ | 元电荷 |
| $k$ | 玻尔兹曼常数 |
| $T$ | 温度，取 300 K |
| $\varepsilon_0$ | 真空介电常数 |
| $\varepsilon_s$ | 硅介电常数，取 $11.7\varepsilon_0$ |
| $n_i$ | 硅本征载流子浓度 |

热电压为：

$$
V_T=\frac{kT}{q}
$$

程序中由 `SI.vt` 返回。

## 3. 突变 PN 结模型

突变结假设 P 区和 N 区的掺杂浓度在冶金结处突然改变。

当前程序参数为：

$$
N_A=1.0\times 10^{17}\ \text{cm}^{-3}
$$

$$
N_D=1.0\times 10^{16}\ \text{cm}^{-3}
$$

### 3.1 内建电势

突变结内建电势为：

$$
V_{bi}=V_T\ln\left(\frac{N_A N_D}{n_i^2}\right)
$$

程序对应函数：

```python
abrupt_vbi(na, nd)
```

### 3.2 耗尽区宽度

在反偏电压 $V_R$ 下，总势垒电压为：

$$
V=V_{bi}+V_R
$$

N 区耗尽宽度：

$$
x_n =
\sqrt{
\frac{2\varepsilon_s}{q}
\frac{N_A}{N_D(N_A+N_D)}
V
}
$$

P 区耗尽宽度：

$$
x_p =
\sqrt{
\frac{2\varepsilon_s}{q}
\frac{N_D}{N_A(N_A+N_D)}
V
}
$$

程序对应函数：

```python
abrupt_widths(na, nd, vr)
```

### 3.3 电场分布

突变结耗尽区中的电场是分段线性的：

$$
E(x)=
\begin{cases}
\frac{qN_D}{\varepsilon_s}(x+x_n), & -x_n\le x\le 0 \\
\frac{qN_A}{\varepsilon_s}(x_p-x), & 0\le x\le x_p
\end{cases}
$$

程序对应函数：

```python
abrupt_field_profile(na, nd, vr)
```

电场最大值出现在结处 $x=0$。

## 4. 线性缓变结模型

线性缓变结假设结附近杂质浓度近似按位置线性变化，浓度梯度为 $a$。

当前程序取：

$$
a=1.0\times 10^{20}\ \text{cm}^{-4}
$$

### 4.1 内建电势

线性缓变结的内建电势满足：

$$
V_{bi}=V_T\ln\left[\left(\frac{a x_0}{2n_i}\right)^2\right]
$$

其中：

$$
x_0=\left(\frac{12\varepsilon_s V_{bi}}{qa}\right)^{1/3}
$$

因为 $V_{bi}$ 同时在等式左右两侧，程序使用 `brentq` 做数值求根。

程序对应函数：

```python
graded_vbi(a_grad)
```

### 4.2 耗尽区宽度和电场

反偏电压 $V_R$ 下，线性缓变结的半耗尽区宽度为：

$$
x_j=\frac{1}{2}
\left[
\frac{12\varepsilon_s(V_{bi}+V_R)}{qa}
\right]^{1/3}
$$

电场分布为抛物线：

$$
E(x)=E_{max}\left[1-\left(\frac{x}{x_j}\right)^2\right]
$$

其中：

$$
E_{max}=\frac{qa}{2\varepsilon_s}x_j^2
$$

程序对应函数：

```python
graded_profile(a_grad, vr)
```

## 5. 碰撞电离率模型

碰撞电离率表示载流子在单位距离内发生碰撞电离的概率强弱。电场越强，载流子获得的能量越大，碰撞电离率通常越高。

程序使用经验公式：

$$
\alpha=A\exp\left[-\left(\frac{B}{E}\right)^m\right]
$$

当前程序取：

$$
m=1
$$

电子参数：

$$
A_n=7.03\times 10^5,\quad B_n=1.23\times 10^6
$$

空穴参数：

$$
A_p=1.58\times 10^6,\quad B_p=2.03\times 10^6
$$

程序对应函数：

```python
ionization_rate(e, a_const, b_const, m=1.0)
```

当 $E$ 较小时，指数项很小，碰撞电离率接近 0；当 $E$ 增大时，$\alpha$ 会快速上升。

## 6. 程序生成了哪些图

### 6.1 单电压电场分布图

文件：

```text
images/field_distribution.png
```

这张图比较 $V_R=40V$ 时：

- 突变结的电场分布。
- 线性缓变结的电场分布。

突变结电场呈三角形或折线形，线性缓变结电场呈抛物线形。

### 6.2 多电压电场分布图

文件：

```text
images/field_distribution_multi_voltage.png
```

这张图展示多个反偏电压：

$$
V_R=0,20,40,80,120V
$$

下的电场分布。可以看到反偏电压越大，耗尽区越宽，电场峰值也越高。

### 6.3 碰撞电离率随电压变化图

文件：

```text
images/ionization_vs_voltage.png
```

横坐标为反偏电压：

$$
V_R
$$

纵坐标为碰撞电离率：

$$
\alpha
$$

程序使用半对数坐标绘图，因为碰撞电离率随电场变化非常快，普通线性坐标不容易看清趋势。

## 7. CSV 数据说明

文件：

```text
data/ionization_vs_voltage.csv
```

列含义如下：

| 列名 | 含义 |
|---|---|
| `VR_V` | 反偏电压 $V_R$，单位 V |
| `Emax_abrupt_V_per_cm` | 突变结最大电场 |
| `alpha_n_abrupt_cm-1` | 突变结电子碰撞电离率 |
| `alpha_p_abrupt_cm-1` | 突变结空穴碰撞电离率 |
| `Emax_graded_V_per_cm` | 线性缓变结最大电场 |
| `alpha_n_graded_cm-1` | 线性缓变结电子碰撞电离率 |
| `alpha_p_graded_cm-1` | 线性缓变结空穴碰撞电离率 |

## 8. 程序整体流程

`solve_q4.py` 的整体流程为：

```text
1. 设置硅材料常数
2. 设置突变结掺杂浓度 NA、ND
3. 设置线性缓变结浓度梯度 a
4. 计算 VR=40V 时的突变结和线性缓变结电场分布
5. 绘制单电压电场分布图
6. 绘制多反偏电压电场分布图
7. 对 VR=0 到 120V 逐点计算最大电场
8. 根据最大电场计算电子和空穴碰撞电离率
9. 绘制碰撞电离率随反偏电压变化图
10. 保存 CSV 数据
11. 生成 result_report.txt
```

## 9. 可以写进实验报告的说明

可以把下面这段放进实验报告：

> 本题分别建立突变 PN 结和线性缓变结的耗尽区电场模型。对于突变结，电场在耗尽区内呈分段线性分布，并在结处达到最大值；对于线性缓变结，电场呈抛物线分布。程序首先计算不同反偏电压下的耗尽区宽度和最大电场，然后利用经验公式 $\alpha=A\exp[-(B/E)^m]$ 计算电子和空穴的碰撞电离率。结果表明，随着反偏电压增大，耗尽区宽度和最大电场均增大，碰撞电离率也随之快速上升。

## 10. 注意事项

1. 题目原文没有给出唯一的掺杂浓度和浓度梯度，因此程序使用一组典型参数。
2. 碰撞电离率对电场非常敏感，参数变化会明显影响曲线。
3. 程序中碰撞电离率随电压变化图使用的是最大电场处的电离率，用来展示随反偏电压变化的主要趋势。

