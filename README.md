# 软件设计题目 Python 程序说明

本目录把 5 道题目拆成了 5 个独立文件夹。每个文件夹中都有一个 `solve_q*.py` 程序，运行后会在本题文件夹内生成结果、图片和数据。

## 环境


需要的 Python 库：

```powershell
python -m pip install -U numpy matplotlib scipy
```

## 运行方式

在当前工作空间根目录运行对应题目的脚本即可。例如：

```powershell
python question_1_matrix\solve_q1.py
python question_2_circle\solve_q2.py
python question_3_minority_distribution\solve_q3.py
python question_4_pn_field_ionization\solve_q4.py
python question_5_avalanche_multiplication\solve_q5.py
```

每个脚本都会在自己的文件夹里生成 `result_report.txt`。有图片的题目会生成 `images/`，有数据表的题目会生成 `data/`。

## 文件夹说明

### `question_1_matrix`

题目 1：矩阵运算。

程序：`solve_q1.py`

输出：

- `result_report.txt`：矩阵 A、B 以及运算结果说明。
- `data/matrix_results.csv`：矩阵运算结果表。

实现内容：

- 计算 `6*A + 5*B`。
- 计算 `A - B + I`。
- 比较矩阵乘法 `A @ B` 和对应元素相乘 `A * B`。

### `question_2_circle`

题目 2：用函数绘制圆形。

程序：`solve_q2.py`

输出：

- `result_report.txt`：圆形绘制说明。
- `images/circle.png`：圆形图像。
- `data/circle_points.csv`：圆上采样点。

实现内容：

- 使用 `circle_points(radius, center, n)` 函数生成圆上坐标。
- 使用 `matplotlib` 绘制圆形。

### `question_3_minority_distribution`

题目 3：双极型晶体管缓变基区少子浓度分布。

程序：`solve_q3.py`

输出：

- `result_report.txt`：公式、参数和说明。
- `images/minority_distribution.png`：归一化少子浓度分布曲线。
- `images/minority_distribution_absolute.png`：少子浓度绝对值 `nB(x)` 曲线。
- `data/minority_distribution.csv`：计算数据。

实现内容：

- 根据 `nB(x)=InE/(DnB*q)*Integral[x,WB](N(t)dt)/N(x)` 进行数值积分。
- 绘制线性杂质分布和指数杂质分布下的归一化少子浓度曲线。
- 额外绘制纵坐标为绝对浓度 `nB(x)` 的曲线。

### `question_4_pn_field_ionization`

题目 4：PN 结势垒区电场分布和碰撞电离率。

程序：`solve_q4.py`

输出：

- `result_report.txt`：仿真参数和关键结果。
- `images/field_distribution.png`：突变结和线性缓变结电场分布。
- `images/field_distribution_multi_voltage.png`：多个反偏电压下的突变结和线性缓变结电场分布。
- `images/ionization_vs_voltage.png`：碰撞电离率随反偏电压变化。
- `data/ionization_vs_voltage.csv`：电场最大值和电离率数据。

实现内容：

- 计算突变结内建电势、耗尽区宽度和电场分布。
- 计算线性缓变结内建电势、耗尽区宽度和电场分布。
- 绘制多个反偏电压下的电场分布图。
- 使用硅材料碰撞电离率经验式计算电子和空穴电离率。

### `question_5_avalanche_multiplication`

题目 5：雪崩倍增因子随反偏电压变化。

程序：`solve_q5.py`

输出：

- `result_report.txt`：仿真参数和结果说明。
- `images/avalanche_multiplication.png`：突变结和线性缓变结雪崩倍增因子线性图。
- `images/avalanche_multiplication_semilog.png`：突变结和线性缓变结雪崩倍增因子半对数图。
- `data/avalanche_multiplication.csv`：反偏电压、突变结倍增因子和线性缓变结倍增因子数据。

实现内容：

- 采用题目 4 中的突变 PN 结和线性缓变结电场。
- 根据电子和空穴碰撞电离率进行数值积分。
- 得到电子倍增因子 `Mn` 和空穴倍增因子 `Mp`。

## 参数说明

题目 4 和题目 5 的原文只给出了公式和物理背景，没有给出唯一的掺杂浓度、浓度梯度和电压范围。因此程序中采用了一组典型硅 PN 结参数：

- `NA = 1.0e17 cm^-3`
- `ND = 1.0e16 cm^-3`
- `a = 1.0e20 cm^-4`
- `VR = 0 到 120 V`

如果老师给了指定参数，直接修改对应脚本开头 `main()` 函数中的这些变量即可。
