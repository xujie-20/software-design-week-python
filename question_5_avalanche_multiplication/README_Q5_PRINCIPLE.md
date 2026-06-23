# 题目 5：雪崩倍增因子的思路与原理说明

本文档解释 `solve_q5.py` 的计算思路、物理含义、公式来源和程序实现方式。

## 1. Markdown 里怎么显示 LaTeX

在 Markdown 文档里，一般有两种写公式的方法。

行内公式用一对 `$` 包起来，例如：

```markdown
电子倍增因子记为 $M_n$，空穴倍增因子记为 $M_p$。
```

显示效果是：电子倍增因子记为 $M_n$，空穴倍增因子记为 $M_p$。

独立公式用两对 `$$` 包起来，例如：

```markdown
$$
M_n = \frac{1}{1-\int_{x_n}^{x_p}\alpha_n(x)
\exp\left[-\int_{x_n}^{x}(\alpha_n(t)-\alpha_p(t))dt\right]dx}
$$
```

显示效果是：

$$
M_n = \frac{1}{1-\int_{x_n}^{x_p}\alpha_n(x)
\exp\left[-\int_{x_n}^{x}(\alpha_n(t)-\alpha_p(t))dt\right]dx}
$$

VS Code 的 Markdown 预览、Typora、Obsidian、很多在线 Markdown 编辑器都支持这种写法。如果某个环境不显示公式，通常不是写法错了，而是那个 Markdown 渲染器没有启用数学公式支持。

## 2. 题目 5 要解决什么问题

题目 5 要求：

> 确定雪崩倍增因子随外加反偏电压的变化关系。

也就是说，我们要研究反偏电压 $V_R$ 改变时，电子和空穴在 PN 结耗尽区中发生碰撞电离后，载流子数量被放大的程度。

程序最终画的是：

- 突变结电子雪崩倍增因子 $M_n$ 随 $V_R$ 的变化。
- 突变结空穴雪崩倍增因子 $M_p$ 随 $V_R$ 的变化。
- 线性缓变结电子雪崩倍增因子 $M_n$ 随 $V_R$ 的变化。
- 线性缓变结空穴雪崩倍增因子 $M_p$ 随 $V_R$ 的变化。

输出图像：

- `images/avalanche_multiplication.png`
- `images/avalanche_multiplication_semilog.png`

输出数据：

- `data/avalanche_multiplication.csv`

CSV 中包含：

```text
VR_V, Mn_abrupt, Mp_abrupt, Mn_graded, Mp_graded
```

## 3. 物理直觉

PN 结反偏时，耗尽区内会形成很强的电场。反偏电压 $V_R$ 越大，耗尽区电场一般越强。

当电子或空穴穿过强电场区域时，会被电场加速。如果获得的能量足够大，它就可能撞击晶格中的价电子，使价电子跃迁到导带，产生新的电子-空穴对。这个过程叫做碰撞电离。

新产生的电子和空穴也会继续被电场加速，继续产生新的电子-空穴对，于是载流子数量被连锁放大。这种现象就是雪崩倍增。

因此，题目 5 的基本逻辑是：

```text
反偏电压 VR 增大
    -> 耗尽区电场 E(x) 增大
    -> 碰撞电离率 alpha_n(x), alpha_p(x) 增大
    -> 雪崩倍增因子 Mn, Mp 增大
```

当倍增因子趋向无穷大时，就接近雪崩击穿。

## 4. 第一步：先得到耗尽区电场

题目 5 需要用到题目 4 的电场分布。当前程序同时计算两种结：

- 突变 PN 结。
- 线性缓变结。

两种结的后续处理完全相同：先得到电场 $E(x)$，再由 $E(x)$ 计算碰撞电离率，最后积分得到 $M_n$ 和 $M_p$。

### 4.1 突变 PN 结电场

设：

- $N_A$：P 区受主浓度。
- $N_D$：N 区施主浓度。
- $V_R$：外加反偏电压。
- $V_{bi}$：PN 结内建电势。
- $\varepsilon_s$：硅介电常数。
- $q$：元电荷。

程序中内建电势计算为：

$$
V_{bi}=V_T\ln\left(\frac{N_A N_D}{n_i^2}\right)
$$

其中：

$$
V_T=\frac{kT}{q}
$$

突变结耗尽区在 N 区和 P 区的宽度分别为：

$$
x_n =
\sqrt{
\frac{2\varepsilon_s}{q}
\frac{N_A}{N_D(N_A+N_D)}
(V_{bi}+V_R)
}
$$

$$
x_p =
\sqrt{
\frac{2\varepsilon_s}{q}
\frac{N_D}{N_A(N_A+N_D)}
(V_{bi}+V_R)
}
$$

程序里对应函数是：

```python
abrupt_widths(na, nd, vr)
```

耗尽区位置范围取为：

$$
x \in [-x_n, x_p]
$$

突变结电场分布是分段线性的：

$$
E(x)=
\begin{cases}
\frac{qN_D}{\varepsilon_s}(x+x_n), & -x_n\le x \le 0 \\
\frac{qN_A}{\varepsilon_s}(x_p-x), & 0\le x \le x_p
\end{cases}
$$

程序里对应函数是：

```python
abrupt_field_profile(na, nd, vr)
```

它返回两个数组：

- `x`：耗尽区中的位置坐标。
- `e`：每个位置上的电场强度 $E(x)$。

### 4.2 线性缓变结电场

线性缓变结假设结附近杂质浓度近似随位置线性变化，浓度梯度为 $a$。当前程序取：

$$
a=1.0\times 10^{20}\ \text{cm}^{-4}
$$

线性缓变结的内建电势需要由下面的关系确定：

$$
V_{bi}=V_T\ln\left[\left(\frac{a x_0}{2n_i}\right)^2\right]
$$

其中：

$$
x_0=\left(\frac{12\varepsilon_s V_{bi}}{qa}\right)^{1/3}
$$

因为 $V_{bi}$ 同时出现在等式左右两边，程序使用 `brentq` 进行数值求解。对应函数是：

```python
graded_vbi(a_grad)
```

在反偏电压 $V_R$ 下，线性缓变结的半耗尽区宽度记为 $x_j$：

$$
x_j=\frac{1}{2}\left[\frac{12\varepsilon_s(V_{bi}+V_R)}{qa}\right]^{1/3}
$$

线性缓变结电场分布为抛物线：

$$
E(x)=E_{max}\left[1-\left(\frac{x}{x_j}\right)^2\right],
\quad -x_j\le x\le x_j
$$

其中：

$$
E_{max}=\frac{qa}{2\varepsilon_s}x_j^2
$$

程序里对应函数是：

```python
graded_field_profile(a_grad, vr)
```

它同样返回位置数组 `x` 和电场数组 `e`。

## 5. 第二步：由电场计算碰撞电离率

有了电场 $E(x)$ 之后，电子和空穴的碰撞电离率可以用经验公式计算。

程序采用：

$$
\alpha = A\exp\left[-\left(\frac{B}{E}\right)^m\right]
$$

其中：

- $\alpha$：碰撞电离率，单位通常为 $\text{cm}^{-1}$。
- $E$：电场强度，单位为 $\text{V/cm}$。
- $A,B,m$：材料相关经验参数。

程序中取 $m=1$。

电子参数：

$$
A_n=7.03\times 10^5,\quad B_n=1.23\times 10^6
$$

空穴参数：

$$
A_p=1.58\times 10^6,\quad B_p=2.03\times 10^6
$$

程序中对应函数是：

```python
ionization_rate(e, a_const, b_const, m=1.0)
```

在 `multiplication_factors_from_field()` 中：

```python
alpha_n = ionization_rate(e, 7.03e5, 1.23e6)
alpha_p = ionization_rate(e, 1.58e6, 2.03e6)
```

也就是分别得到：

$$
\alpha_n(x)
$$

和：

$$
\alpha_p(x)
$$

## 6. 第三步：计算雪崩倍增因子

题目给出的电子和空穴雪崩倍增因子本质上都是：

```text
1 / (1 - 某个碰撞电离积分)
```

如果积分值很小，则分母接近 1，倍增因子接近 1，说明几乎没有倍增。

如果积分值增大，分母变小，倍增因子变大。

如果积分值接近 1，分母接近 0，倍增因子会急剧增大，表示接近雪崩击穿。

### 6.1 电子倍增因子

程序采用的电子倍增因子形式是：

$$
M_n =
\frac{1}{
1 -
\int_{x_n}^{x_p}
\alpha_n(x)
\exp\left[
-\int_{x_n}^{x}(\alpha_n(t)-\alpha_p(t))dt
\right]dx
}
$$

这里要注意，程序中的左边界坐标实际是 $-x_n$，右边界是 $x_p$。为了写法简洁，公式中把积分区间记成从左边界到右边界。

代码对应部分：

```python
diff_np = alpha_n - alpha_p
inner_n = cumulative_trapezoid(diff_np, x, initial=0.0)
term_n = alpha_n * np.exp(np.clip(-inner_n, -100.0, 100.0))
integral_n = np.trapezoid(term_n, x)
mn = 1.0 / (1.0 - integral_n)
```

解释一下：

- `diff_np = alpha_n - alpha_p` 对应 $\alpha_n-\alpha_p$。
- `inner_n` 是从左边界到当前位置的累计积分。
- `term_n` 是外层积分的被积函数。
- `integral_n` 是外层积分结果。
- `mn` 就是电子倍增因子。

### 6.2 空穴倍增因子

空穴倍增因子形式类似，只是电子和空穴的角色交换：

$$
M_p =
\frac{1}{
1 -
\int_{x_n}^{x_p}
\alpha_p(x)
\exp\left[
-\int_{x}^{x_p}(\alpha_p(t)-\alpha_n(t))dt
\right]dx
}
$$

代码对应部分：

```python
diff_pn = alpha_p - alpha_n
cumulative_p = cumulative_trapezoid(diff_pn, x, initial=0.0)
inner_p = cumulative_p[-1] - cumulative_p
term_p = alpha_p * np.exp(np.clip(-inner_p, -100.0, 100.0))
integral_p = np.trapezoid(term_p, x)
mp = 1.0 / (1.0 - integral_p)
```

这里：

- `cumulative_p[-1] - cumulative_p` 表示从当前位置积分到右边界。
- `mp` 是空穴倍增因子。

## 7. 为什么要用数值积分

公式里有两层积分，而且电场 $E(x)$ 是分段函数，碰撞电离率又包含指数项：

$$
\alpha(x)=A\exp\left[-\frac{B}{E(x)}\right]
$$

如果手算，会非常繁琐。程序采用数值方法：

1. 把耗尽区 $[-x_n,x_p]$ 离散成很多点。
2. 在每个点上计算电场 $E(x)$。
3. 在每个点上计算 $\alpha_n(x)$ 和 $\alpha_p(x)$。
4. 用梯形积分近似内层积分和外层积分。
5. 得到 $M_n$ 和 $M_p$。

程序中用到两个积分函数：

```python
cumulative_trapezoid(...)
np.trapezoid(...)
```

`cumulative_trapezoid` 用于做“从边界到每个位置”的累计积分。

`np.trapezoid` 用于做整个区间上的普通积分。

## 8. 为什么程序里有 `np.clip`

代码中有：

```python
np.exp(np.clip(-inner_n, -100.0, 100.0))
```

这是为了避免指数函数溢出。

因为当电场很强时，碰撞电离率会变大，积分项也会变大，指数函数可能出现非常大的数。`np.clip` 把指数的输入限制在 `[-100, 100]` 范围内，防止程序因为数值溢出而报错。

这不会改变低电压区域的正常趋势，主要是为了让接近击穿时的计算更稳定。

## 9. 为什么线性图只画到某个电压附近

程序里线性图使用：

```python
mn_plot = np.array([value if math.isfinite(value) and value < 100 else np.nan for value in mn_values])
mp_plot = np.array([value if math.isfinite(value) and value < 100 else np.nan for value in mp_values])
```

也就是说，线性图只显示 $M<100$ 的部分。

原因是：接近雪崩击穿时，$M$ 会快速增大。如果把特别大的值也画在线性坐标上，前面大部分曲线会被压扁，看不清。

因此程序另外生成了半对数图：

```text
images/avalanche_multiplication_semilog.png
```

半对数图更适合观察快速上升的趋势。

## 10. 图像怎么看

### 10.1 线性图

文件：

```text
images/avalanche_multiplication.png
```

横坐标：

$$
V_R
$$

表示反偏电压。

纵坐标：

$$
M
$$

表示雪崩倍增因子。

如果 $M=1$，说明基本没有倍增。

如果 $M=10$，说明载流子数量大约放大了 10 倍。

如果 $M$ 快速变大，说明 PN 结正在接近雪崩击穿。

### 10.2 半对数图

文件：

```text
images/avalanche_multiplication_semilog.png
```

纵坐标用对数刻度，所以更适合看 $M$ 从 1 增大到几十、几百的过程。

## 11. 当前程序使用的参数

当前程序使用：

```python
na = 1.0e17
nd = 1.0e16
a_grad = 1.0e20
vr_values = np.linspace(0.0, 120.0, 121)
```

也就是：

$$
N_A=1.0\times 10^{17}\ \text{cm}^{-3}
$$

$$
N_D=1.0\times 10^{16}\ \text{cm}^{-3}
$$

$$
a=1.0\times 10^{20}\ \text{cm}^{-4}
$$

$$
V_R = 0,1,2,\cdots,120\ \text{V}
$$

这些是典型仿真参数。题目原文没有给唯一指定的掺杂浓度和电压范围，所以报告中要说明参数是假定值。如果老师给了指定参数，就改 `solve_q5.py` 中 `main()` 函数里的：

```python
na = 1.0e17
nd = 1.0e16
a_grad = 1.0e20
vr_values = np.linspace(0.0, 120.0, 121)
```

## 12. 程序整体流程

`solve_q5.py` 的整体流程可以概括为：

```text
1. 设置硅材料参数
2. 设置突变结掺杂浓度 NA、ND，以及线性缓变结浓度梯度 a
3. 生成反偏电压数组 VR
4. 对每一个 VR：
   4.1 计算突变结耗尽区宽度 xn、xp 和电场 E(x)
   4.2 计算线性缓变结半耗尽区宽度 xj 和电场 E(x)
   4.3 根据 E(x) 计算 alpha_n(x)、alpha_p(x)
   4.4 对 alpha 进行数值积分
   4.5 分别得到突变结和线性缓变结的 Mn、Mp
5. 保存 CSV 数据
6. 绘制线性图和半对数图
7. 生成 result_report.txt
```

## 13. 可以写进实验报告的简短说明

可以把下面这段放入实验报告中：

> 在反偏 PN 结中，耗尽区内强电场会使载流子获得较高动能，并通过碰撞电离产生新的电子-空穴对，从而形成雪崩倍增。本文分别建立突变结和线性缓变结的耗尽区电场模型，计算不同反偏电压下的电场分布，然后由经验公式 $\alpha=A\exp[-(B/E)^m]$ 分别计算电子和空穴的碰撞电离率，最后对倍增因子公式进行数值积分，得到电子倍增因子 $M_n$ 和空穴倍增因子 $M_p$ 随反偏电压 $V_R$ 的变化关系。计算结果表明，随着 $V_R$ 增大，电场增强，碰撞电离率迅速上升，雪崩倍增因子也快速增大，说明器件逐渐接近雪崩击穿状态。

## 14. 注意事项

1. 倍增因子对电场非常敏感，所以掺杂浓度、电压范围、碰撞电离率参数改变后，曲线会明显变化。
2. 当分母 $1-\int(\cdots)dx$ 接近 0 时，$M$ 会急剧增大，这是接近击穿的表现。
3. 当前程序已经同时计算突变结和线性缓变结。两者使用同一个倍增因子积分函数，区别只在于电场 $E(x)$ 的计算模型不同。
