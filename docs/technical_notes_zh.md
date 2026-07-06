# pyCFRAM 技术说明文档

## 1. 项目概述

pyCFRAM 是 CFRAM（Climate Feedback–Response Analysis Method，气候反馈-响应分析方法）的 Python 封装实现。该框架以 Fortran RRTMG 辐射传输模型为核心引擎，通过 Python 完成数据预处理、矩阵求解、并行调度和可视化，实现了对温度变化的定量归因分解。

当前版本（v0.01）已成功复现 Wu et al. (2025, *J. Climate*) 论文中 Fig.3（温度分解空间图）、Fig.4（PAP 柱状图）和 Fig.5（气溶胶分种类空间图）的结果。

### 1.1 CFRAM 方法原理

CFRAM 的核心思想是将观测到的温度变化 ΔT 分解为各物理过程的贡献：

```
ΔT = -(∂R/∂T)⁻¹ × [Σ ΔF_radiative + Σ ΔF_non-radiative]
```

其中：
- `(∂R/∂T)` 是 Planck 反馈矩阵（38×38，37 个大气层 + 地表），描述温度对辐射通量的偏导数
- `ΔF_radiative` 包括水汽、云、气溶胶、CO₂、O₃、太阳辐射、地表反照率、地表温度等 9 项辐射强迫
- `ΔF_non-radiative` 包括潜热通量、感热通量、地表过程等非辐射强迫
- 动力项（atmdyn/sfcdyn/ocndyn/dry）由暖态辐射不平衡 `frc_full` 经同一 Planck 逆矩阵求得（见 §2.1、algorithm_spec §8.1）。其中 `dry = atmdyn + sfcdyn`，且 `sfcdyn = ocndyn + lhflx + shflx`（机器精度严格成立——lhflx/shflx 强迫仅置于地表行）。仅在 OLD-CFRAM 风格的 closure 出图中，`atmdyn` 才被重定义为残差 `ΔT_obs − Σ其他项`

### 1.2 参考文献

**CFRAM 方法论**
- Lu, J., and M. Cai, 2009: A new framework for isolating individual feedback processes in coupled general circulation climate models. Part I. *Climate Dynamics*, 32, 873–885.
- Cai, M., and J. Lu, 2009: A new framework for isolating individual feedback processes in coupled general circulation climate models. Part II. *Climate Dynamics*, 32, 887–900.

**CFRAM-A 框架（首次将气溶胶纳入 CFRAM）**
- Zhang, T., Y. Deng, J. Chen, S. Yang, P. Gao, and H. Zhang, 2022: Disentangling physical and dynamical drivers of the 2016/17 record-breaking warm winter in China. *Environmental Research Letters*, 17(7), 074024. https://doi.org/10.1088/1748-9326/ac79c1
  首个把气溶胶作为独立辐射强迫项显式加入 CFRAM 分解的工作。论文原文："the effect of aerosols … has not been included in the previous CFRAM analysis"。给出 5 物种分解（BC/OC/Sulfate/Sea salt/Dust），用 MERRA-2 驱动 + 离线 RRTMG v5 做辐射传输。pyCFRAM 沿用这一 CFRAM-A 概念框架（含 RRTMG + MERRA-2 驱动），并用 Python 编排层重新实现：ERA5 + MERRA-2 自动预处理、Python 端全矩阵 Planck 求逆、multiprocessing 按格点并行。
  注意：Zhang 2022 的气溶胶光学用 **MAM4/CAM6**；pyCFRAM 用 **GOCART 查找表**（承自 Wu 等 2025 参考代码），两者对同一 MERRA-2 混合比给出不同的物种级 forcing。

**CFRAM 应用（极端事件归因）**
- Wu, Q., Q. Li, T. Zhang, X. Sun, S. Yang, and X. Hu, 2025: Quantitative attribution of 2013 and 2022 extremely hot events in China: insights from a climate feedback–response perspective. *Journal of Climate*, 38(17), 4331–4349.

---

## 2. 系统架构

### 2.1 Fortran-Python 分工

```
┌─────────────────────────────────────────────────────────┐
│ Fortran 辐射引擎 (RRTMG 或 Fu, 由 case.yaml.scheme 选)  │
│                                                         │
│  输入：单柱大气廓线 + 地表状态（base 和 perturbed）     │
│  计算：~11 次 rad_driver 辐射传输调用                   │
│        + nlayer+1 次 Planck 矩阵构建（逐层 ±0.5K 扰动） │
│        + 6 次 per-species 气溶胶扰动 (RRTMG only)       │
│        + Fu only: dual-MC sub-column overlap (base/warm)│
│  输出：9 项 bulk forcing + cloud_lw/sw 拆分 +           │
│        6 物种 forcing + Planck 逆矩阵                   │
└──────────────────────┬──────────────────────────────────┘
                       │ (临时二进制文件)
┌──────────────────────▼──────────────────────────────────┐
│ Python（分解求解 + 并行调度）                           │
│                                                         │
│  辐射项：  dT_i = -(∂R/∂T)⁻¹ × frc_i    (numpy 矩阵乘)│
│  非辐射项：dT_lhflx/shflx (强迫仅置地表行,同一 Planck) │
│  气溶胶分种类：dT_bc/oc/sulf/seas/dust   (同一 Planck)  │
│  动力项：dT_atmdyn/sfcdyn/ocndyn/dry  (frc_full,同Planck)│
│         dry=atmdyn+sfcdyn; sfcdyn=ocndyn+lhflx+shflx     │
│  并行：multiprocessing Pool, 每格点独立进程              │
└─────────────────────────────────────────────────────────┘
```

### 2.2 为什么不用 OpenMP？

RRTMG 内部的 Fortran module 全局变量在线程间共享，导致数据竞争。尝试 OpenMP 并行后出现静默错误（无输出但不崩溃）。

采用 Python multiprocessing 方案：每个 worker 进程独立运行一个 Fortran 可执行文件（cfram_rrtmg_1col），通过进程隔离天然避免了 RRTMG 的线程安全问题。在 256 核的 hqlx204 上，80 进程并行，每个 case（9801 格点）约 19 分钟完成。

### 2.3 关键设计决策

| 决策 | 选择 | 理由 |
|------|------|------|
| dT 求解位置 | Python（非 Fortran） | 新增分解项不需要改 Fortran，只需提供 forcing |
| 并行方式 | multiprocessing（非 OpenMP） | RRTMG 线程不安全 |
| 非辐射项 | 用户提供 forcing，Planck 矩阵求解 | 保持通用性 |
| 气溶胶分种类 | per-species forcing × 自计算 Planck 矩阵 | 半自计算 |
| 动力项 (atmdyn/sfcdyn/ocndyn/dry) | 由 frc_full 经 Planck 逆求得（非残差） | sfcdyn=ocndyn+lhflx+shflx 严格成立 |
| 输入格式 | 标准 NetCDF（base + perturbed 两个状态） | 不预设 event/climatology |

---

## 3. 目录结构

```
pyCFRAM/
├── run_case.py              # 一键入口：python3 run_case.py eh13
├── requirements.txt         # Python 依赖
├── configs/
│   └── defaults.yaml        # 全局配置（辐射参数、气压层、气溶胶映射、出图设置）
├── cases/                   # 每个 case 一个子目录
│   ├── eh13/
│   │   ├── case.yaml        # case 配置（输入路径、出图区域）
│   │   ├── input/           # 输入数据（NetCDF，符合 input_spec）
│   │   ├── output/          # CFRAM 计算结果
│   │   └── figures/         # 出图结果
│   └── eh22/
├── core/                    # Python 核心模块
│   ├── config.py            # 统一配置加载器（含 get_drdt_eval/probe/get_co2_handling/get_q_handling）
│   └── constants.py         # 物理常量 + RRTMG 固有参数
│                            #（气溶胶光学已 inline 进 run_parallel_python.py:compute_aerosol_column）
├── fortran/                 # Fortran 辐射引擎（RRTMG + Fu，由 case.yaml 选择）
│   ├── cfram_rrtmg.f90      # RRTMG 单柱主程序（make sed 出 cfram_rrtmg_1col.f90）
│   ├── cfram_fu_1col.f90    # Fu 单柱主程序（dual-MC sub-column 重叠）
│   ├── Fu/                  # Fu 源码（cas_fu_radiation.f, fu_helpers.f, para.file）
│   ├── rad_driver.f90       # RRTMG 调用驱动
│   ├── rad_driver_lw.f90    # LW 通量计算（calc_drdt 用）
│   ├── drdt.f90             # Planck 矩阵 ∂R/∂T（含 centered FD ±0.5K 分支）
│   ├── math.f90             # 矩阵求逆（LAPACK DGETRF/DGETRS）
│   ├── writeout.f90         # 输出工具
│   ├── makefile             # 双 toolchain：intel (ifort+MKL) / gnu (gfortran+conda LAPACK)
│   ├── data_prep/aerosol/   # GOCART 气溶胶查找表（opticsBands_*.nc）
│   ├── RRTMG_LW-master/     # 上游 RRTMG 长波库（含 rrtmg_lw.nc）
│   └── RRTMG_SW-master/     # 上游 RRTMG 短波库（含 rrtmg_sw.nc）
├── scripts/                 # 工作流脚本
│   ├── build_case_input.py         # ERA5 + MERRA-2 → 标准 NetCDF 输入（case.yaml 驱动）
│   ├── download_era5_flux.py       # ERA5 PL+SL 下载（CDS API）
│   ├── download_merra2_aerosol.py  # MERRA-2 M2I3NVAER 下载（NASA GES DISC）
│   ├── extract_full_field.py       # 输入 NetCDF → Fortran 二进制（含 GOCART 气溶胶光学）
│   ├── run_parallel_python.py      # 并行 CFRAM 计算（multiprocessing）
│   ├── plot_fig3_independent.py    # Fig.3 双列对比图（paper | independent）
│   ├── plot_fig3_self.py           # Fig.3 温度分解空间图（基于 paper_data 复跑）
│   ├── plot_fig3.py                # Fig.3 paper_data 直出（验证用）
│   ├── plot_fig4.py                # Fig.4 PAP 柱状图
│   ├── plot_fig5.py                # Fig.5 气溶胶分种类
│   └── validate_vs_paper.py        # 对比 surface dT 与 Wu et al. 结果
├── data/                    # 数据源模块
│   ├── source_base.py              # DataSource 抽象基类 + 工厂注册器
│   ├── era5_source.py              # ERA5DailySource：Wu et al. warm-period 方法全场版
│   └── merra2_aerosol.py           # MERRA-2 气溶胶加载器（垂直 log-p + 水平双线性插值）
├── plotting/                # 可视化模块
└── docs/                    # 文档
    ├── algorithm_spec.md    # 算法规格说明
    ├── input_spec.md        # 输入数据格式规范
    ├── technical_notes_en.md # 精简英文版技术笔记（社区版）
    └── technical_notes_zh.md # 本文档（详细中文版）
```

---

## 4. 安装与编译

### 4.1 环境要求

| 组件 | 版本 / 说明 |
|------|-------------|
| Python | 3.8+（推荐 conda/miniconda 环境） |
| gfortran | ≥ 9.0（编译 RRTMG） |
| LAPACK/BLAS | 桌面 Linux 通常自带；HPC 集群须 `module load lapack`（见 §13.9） |
| GLIBC | ≥ 2.25（netCDF4 Python 包依赖；老机器如 hqlx74 不满足，需换机） |
| 磁盘 | 代码 ~500 MB；ERA5+MERRA-2 数据约 **40 GB** |
| 账号 | Copernicus CDS（下载 ERA5）+ NASA Earthdata（下载 MERRA-2） |

### 4.2 Python 环境准备

```bash
# 推荐 conda
conda create -n pycfram python=3.10 -y
conda activate pycfram

# 安装依赖
git clone git@github.com:lzhenn/pyCFRAM.git
cd pyCFRAM
pip install -r requirements.txt
```

`requirements.txt` 列出了 numpy、netCDF4、matplotlib、cartopy、scipy、pyyaml、cdsapi 等。

### 4.3 账号与凭据配置

本文档假设组员**已具备** CDS 与 NASA Earthdata 账号。凭据文件需放置在以下位置：

**CDS（ERA5 下载）**：`~/.cdsapirc`
```
url: https://cds.climate.copernicus.eu/api
key: <YOUR_UID>:<YOUR_API_KEY>
```
申请地址：`https://cds.climate.copernicus.eu/` → Profile → API key

**NASA Earthdata（MERRA-2 下载）**：`~/.netrc`（权限 600）
```
machine urs.earthdata.nasa.gov login <USERNAME> password <PASSWORD>
```
申请地址：`https://urs.earthdata.nasa.gov/` → 注册后还需在 Profile → Applications 里批准 "NASA GESDISC DATA ARCHIVE"。

### 4.4 Fortran 编译

```bash
cd fortran
make                       # 默认 RRTMG: cfram_rrtmg_1col
make fu                    # Fu 引擎: cfram_fu_1col（runtime nlev）
make TOOLCHAIN=intel       # ifort + MKL（hqlx220/204 生产环境，默认）
make TOOLCHAIN=gnu         # Mac 本地 gfortran + conda LAPACK
cd ..
```

可执行文件（两者均通过 `data_prep/plev.dat` 文件大小推导运行时 nlev，单一 binary 处理任意垂直网格 17/19/30/37/…）：
- `cfram_rrtmg_1col`：RRTMG 单柱版
- `cfram_fu_1col`：Fu 单柱版
- `cfram_fu_1col_sp`：Fu 单精度调试版（与 OLD CFRAM 默认 ifort 行为对齐，validation 用）

每个 case 的 `case.yaml` 通过 `radiation.scheme` 字段选择引擎（推荐）：

```yaml
radiation:
  scheme: rrtmg      # 或 fu；将来扩展见 core/config.py:RADIATION_SCHEMES
```

`run.executable: <binary>` 仍作为 legacy 直传名 binary 的逃生口保留。

**Fu 引擎 dual MC sub-column overlap**（2026-05-10 修复）：
Fu 引擎需要两套 Monte Carlo 子列云重叠模式，分别从 `cc_base` 和
`cc_warm` 抽样：

- `base_no_cloud` —— 由 `cc_base` 抽样，用于 case 0/2/3/4/5/6/7 + drdt
- `warm_no_cloud` —— 由 `cc_warm` 抽样，用于 case 1/8/9（warm/cloud/full）

apple-to-apple 对齐 OLD CFRAM 源码（GW-base.f 写 base_no_cloud_out，
GW-warm.f 写 warm_no_cloud_out，GW-cloud.f 读 warm 版）。如果两个状态
共用一份 pattern（pyCFRAM 2026-05-10 之前的行为），warm 云态下 RT 看到
的 sub-col cloud fraction 与实际加载的 cc_warm 错配，会出现 +1.8 K 的
CLDL/CLDS 全场镜像翻转 bug。修复后单柱所有 13 个 dT 项与 collab
`partial_T_1.grd` 偏差 ≤ 0.003 K（详见 `session_log.md` 2026-05-10 条目）。

可选地把 OLD CFRAM 的精确 MC 输出（`base/warm_no_cloud_out_1.dat`）拷贝到
`data_prep/{base,warm}_no_cloud_seed.dat`，Fu binary 会读取这两份文件做
bit-perfect 单柱复现。

---

## 5. 运行流程

### 5.1 使用 ERA5 + MERRA-2 独立驱动（标准流程）

整个流程包含 5 步：**下载 → 构建输入 → 运行 CFRAM → 出图 → 验证**。每步给出可直接执行的命令；所有路径均相对项目根 `pyCFRAM/`。

#### Step 1：下载 ERA5 + MERRA-2 源数据

数据量估计：ERA5 daily ~5 GB，MERRA-2 aerosol ~35 GB（2003–2022 年 8 月全量）。

```bash
# ERA5 pressure-level + single-level（通过 CDS API）
python3 scripts/download_era5_flux.py

# MERRA-2 M2I3NVAER 3-hourly 气溶胶（通过 NASA GES DISC）
python3 scripts/download_merra2_aerosol.py
```

下载后目录结构：

```
era5_data/
├── daily/
│   ├── era5_pl_{var}_{YYYY}08.nc          # var ∈ {t,q,o3,cc,ciwc,clwc}
│   └── era5_sl_{YYYY}08/
│       ├── data_stream-oper_stepType-instant.nc   # skt, sp
│       ├── data_stream-oper_stepType-accum.nc     # ssrd, ssr, tisr, slhf, sshf
│       └── data_stream-oper_stepType-max.nc       # mx2t
└── merra2/
    └── M2I3NVAER_{YYYYMMDD}.nc4           # 每天一个文件
```

两个下载脚本均可断点续传（重复运行会跳过已存在文件）。

#### Step 2：构建标准 NetCDF 输入

```bash
python3 scripts/build_case_input.py --case eh13
python3 scripts/build_case_input.py --case eh22
```

该脚本从 `cases/<case>/case.yaml` 的 `source:` 段读取配置（warm_days、clim_years、CO2 值、气溶胶数据源），执行：

1. ERA5 PL 变量 → daily mean → warm-days 平均 → base / perturbed 两态
2. ERA5 SL accum 字段（solar, lhflx, shflx）经 `×6 ÷86400` 累积修正转 W/m²
3. MERRA-2 气溶胶 6 物种（bc/ocphi/ocpho/sulf/ss/dust）垂直 log-p + 水平双线性插值到 ERA5 网格
4. 写出 5 个 NetCDF 文件到 `cases/<case>/input/`

输出：

```
cases/eh13/input/
├── base_pres.nc          # (1, 37, 97, 121) T, q, o3, camt, cliq, cice, co2, 6 种气溶胶
├── base_surf.nc          # (1, 97, 121) ts, ps, solar, albedo
├── perturbed_pres.nc
├── perturbed_surf.nc
└── nonrad_forcing.nc     # ΔF lhflx, shflx (W/m²)
```

构建过程会打印关键指标（solar/lhflx/shflx domain mean 等），若数值与 Step 5 验证表相差 5× 以上，检查 ×6 累积修正（见 [第 13 节排错](#13-常见坑与排错)）。

#### Step 3：运行 CFRAM 计算

```bash
python3 run_case.py eh13
python3 run_case.py eh22
```

`run_case.py` 默认依次执行四步：`build → extract → run → plot`。若 Step 2 已完成，可跳过 build：

```bash
python3 run_case.py eh13 --step extract    # NetCDF → Fortran 二进制 + GOCART 光学属性
python3 run_case.py eh13 --step run --nproc 40   # multiprocessing 并行调 cfram_rrtmg_1col
python3 run_case.py eh13 --step plot       # 生成 fig3_independent.png
```

80 进程并行时单 case 约 **20 分钟**（9801 格点，~8.3 pts/s）。

#### Step 4：查看结果

```
cases/eh13/output/cfram_result.nc       # 分解结果（变量列表见 §8）
cases/eh13/figures/fig3_independent.png # 双列对比图：paper | independent
```

额外出图（可选）：

```bash
python3 scripts/plot_fig5.py --case eh13 eh22   # Fig.5 气溶胶分种类
```

#### Step 5：验证（对比 Wu et al. paper_data）

前提：`paper_data/` 目录已放置（仅验证需要，非标准流程必需）。

```bash
python3 scripts/validate_vs_paper.py eh22
python3 scripts/validate_vs_paper.py eh13
```

期望输出（关键行）：

```
=== EH22: CORRECTED comparison: surface dT ===
  term            paper      indep     corr
  cld            2.0986     2.5684    0.952
  wv             1.1597     1.2035    0.998
  ...
  lhflx         -3.5229    -3.6152    0.993
  shflx         -3.3168    -3.3560    0.993
```

合格标准（§9 详细表）：
- **wv、lhflx、shflx 的 corr > 0.99** — pipeline 核心正确
- **cld、aerosol 的 corr > 0.90** — 空间模式一致，数值上存在已知的路径 A vs B 互补偏差

### 5.2 case.yaml 标准格式

以 `cases/eh22/case.yaml` 为例：

```yaml
case_name: EH22
description: "East China extreme heat, Aug 2022 vs 2003–2022 climatology"

source:
  type: era5_daily
  data_dir: era5_data/daily
  temporal:
    event_year: 2022
    event_month: 8
    clim_years: [2003, 2022]
    warm_days: [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
    # 可选：若不指定 warm_days，会用以下自动检测
    warm_detect:
      threshold_pct: 90
      min_consecutive: 3
      detect_region:
        lat: [27, 35]
        lon: [103, 122]
  co2:
    source: constant
    base_ppmv: 394.23          # 气候态，来自 paper_data 反推
    perturbed_ppmv: 416.89      # 事件年
  aerosol:
    source: merra2
    data_dir: era5_data/merra2

input:
  base_pres: input/base_pres.nc
  base_surf: input/base_surf.nc
  perturbed_pres: input/perturbed_pres.nc
  perturbed_surf: input/perturbed_surf.nc
  nonrad_forcing: input/nonrad_forcing.nc

plot:
  key_region:
    lon: [103, 122]
    lat: [27, 35]
```

### 5.3 添加新 case

1. 复制现有 case 目录作为模板：

```bash
cp -r cases/eh22 cases/my_case
```

2. 编辑 `cases/my_case/case.yaml`：修改 `case_name`、`source.temporal`（event_year、warm_days）、`co2`、`plot.key_region`。

3. 运行：

```bash
python3 scripts/build_case_input.py --case my_case
python3 run_case.py my_case
```

如果是非 ERA5 数据源（如 CMIP6 模式输出），需在 `data/` 下新增 DataSource 子类（参考 `era5_source.py`）并用 `@register_source('your_type')` 注册，再在 `case.yaml` 里设 `source.type: your_type`。

---

## 6. 输入数据格式

详见 [input_spec.md](input_spec.md)。核心要求：

- **base_pres.nc / perturbed_pres.nc**：大气廓线（37层），包含温度、比湿、臭氧、云量/水含量、CO₂、6种气溶胶
- **base_surf.nc / perturbed_surf.nc**：地表状态（温度、气压、太阳辐射、反照率）
- **nonrad_forcing.nc**（可选）：非辐射 forcing（潜热、感热、地表过程），以及气溶胶分种类 forcing（bc, oc, sulf, seas, dust）

pyCFRAM 不限定两个状态的定义方式——可以是事件期 vs 气候态、未来 vs 历史、敏感性实验 vs 控制实验等。

---

## 7. 配置系统

### 7.1 全局配置 `configs/defaults.yaml`

包含所有与 case 无关的参数：

| 配置项 | 说明 |
|--------|------|
| `radiation.scon` | 太阳常数 (1360.98 W/m²) |
| `radiation.co2_ppmv` | CO₂ 默认浓度 (395 ppmv，输入数据有则覆盖) |
| `radiation.ch4_ppmv` | CH₄ 浓度 (1.6 ppmv) |
| `radiation.n2o_ppmv` | N₂O 浓度 (0.28 ppmv) |
| `radiation.cloud_re_ice/liq` | 云粒子有效半径 (5/20 microns) |
| `grid.pressure_levels` | 37 层标准气压层 (hPa) |
| `aerosol.species` | 气溶胶种类 → GOCART 查找表映射 |
| `run.nproc` | 并行进程数 (auto = 所有 CPU) |
| `plotting.levels_dt` | 温度分解图色阶 |

### 7.2 Case 配置 `cases/<name>/case.yaml`

每个 case 独立的配置：输入文件路径、出图参数等。

### 7.3 配置加载

所有脚本通过 `core/config.py` 统一加载配置：

```python
from core.config import load_case, defaults, get_plev, get_aerosol_map
```

---

## 8. 输出格式

### 8.1 cfram_result.nc

每个 case 的 CFRAM 计算结果保存为 NetCDF：

| 变量 | 维度 | 说明 |
|------|------|------|
| `dT_q` | (lev, lat, lon) | 水汽贡献的温度变化 |
| `dT_cloud` | (lev, lat, lon) | 云贡献 |
| `dT_aerosol` | (lev, lat, lon) | 气溶胶总贡献（RRTMG one-at-a-time） |
| `dT_co2` | (lev, lat, lon) | CO₂ 贡献 |
| `dT_o3` | (lev, lat, lon) | 臭氧贡献 |
| `dT_ts` | (lev, lat, lon) | 地表温度辐射效应 |
| `dT_solar` | (lev, lat, lon) | 太阳辐射贡献 |
| `dT_albedo` | (lev, lat, lon) | 地表反照率贡献 |
| `dT_warm` | (lev, lat, lon) | 总辐射变化 |
| `dT_lhflx` | (lev, lat, lon) | 潜热通量贡献（需 nonrad_forcing） |
| `dT_shflx` | (lev, lat, lon) | 感热通量贡献 |
| `dT_sfcdyn` | (lev, lat, lon) | 地表过程贡献（= ocndyn + lhflx + shflx，严格成立） |
| `dT_ocndyn` | (lev, lat, lon) | 海洋/次表层热输送（frc_full 地表行 − Δlh − Δsh） |
| `dT_bc/oc/sulf/seas/dust` | (lev, lat, lon) | 气溶胶分种类贡献 |
| `dT_observed` | (lev, lat, lon) | 观测温度变化 (T_perturbed - T_base) |
| `dT_atmdyn` | (lev, lat, lon) | 大气动力贡献（frc_full 大气层行，非残差） |
| `dT_dry` | (lev, lat, lon) | = atmdyn + sfcdyn（frc_full 整柱响应） |
| `dT_cloud_lw/sw` | (lev, lat, lon) | 云 LW/SW 拆分（精确可加 = dT_cloud） |
| `frc_*` | (lev, lat, lon) | 对应的辐射 forcing (W/m²) |

### 8.2 level 维度

38 层 = 37 个大气层（TOA→地表）+ 1 个地表层。地表层为最后一个 index。

---

## 9. 验证结果（独立 ERA5+MERRA-2 pipeline, v3）

以下结果基于 2026-04-12 的 v3 版本（×6 累积修正、CO2 对齐、warm_days 对齐完成后）。完整重跑命令：

```bash
python3 scripts/build_case_input.py --case eh22
python3 run_case.py eh22
python3 scripts/validate_vs_paper.py eh22
```

### 9.1 EH22（2022 年 8 月事件 vs 2003–2022 气候态）

Domain mean 地表 ΔT 与 Wu et al. paper_data 对比（空间裁剪至 20–40°N）：

| 分解项 | paper | v3 (indep) | 空间 corr | 说明 |
|---|---|---|---|---|
| Water Vapor | 1.16 K | 1.20 K | **0.998** | 核心项，近完美 |
| Cloud | 2.10 K | 2.57 K | 0.952 | 符号正确；数值差异见 §9.3 |
| Aerosol (total) | 1.11 K | 1.38 K | 0.829 | 同上 |
| CO₂ | 0.077 K | 0.081 K | 0.964 | ✓ |
| Albedo | 0.220 K | 0.196 K | 0.972 | ✓ |
| Solar | −0.032 K | −0.029 K | 0.757 | 值小，相关性对小扰动敏感 |
| LHFLX | −3.52 K | −3.62 K | **0.993** | 核心项 |
| SHFLX | −3.32 K | −3.36 K | **0.993** | 核心项 |

### 9.2 EH13（2013 年 8 月事件 vs 2003–2022 气候态）

| 分解项 | paper | v3 (indep) | 空间 corr | 说明 |
|---|---|---|---|---|
| Water Vapor | −1.35 K | −1.43 K | **0.997** | 核心项 |
| Cloud | 1.33 K | 0.60 K | 0.945 | 偏低但空间模式一致 |
| Aerosol (total) | 0.62 K | 1.01 K | 0.873 | 偏高但空间模式一致 |
| CO₂ | 0.003 K | 0.002 K | 0.156 | 值极小，corr 不具参考性 |
| Albedo | 0.082 K | 0.058 K | 0.962 | ✓ |
| LHFLX | −1.95 K | −1.80 K | **0.992** | 核心项 |
| SHFLX | −2.49 K | −2.39 K | **0.992** | 核心项 |

### 9.3 Cloud / Aerosol 互补偏差的来源

Paper（路径 B）采用**离线预计算气溶胶 forcing + Planck 矩阵求逆**，本质为一阶线性分解；本实现（路径 A）采用 **RRTMG 内直接注入气溶胶混合比 + CFRAM partial perturbation**，保留 aerosol-cloud 辐射耦合的非线性。两者在 cloud 与 aerosol 各自的分配上会产生系统性偏差，但**cloud + aerosol 合计**与 paper corr > 0.95：

- EH22: cloud+aer paper=3.21 K, indep=3.95 K
- EH13: cloud+aer paper=1.95 K, indep=1.62 K

这是**方法学差异**而非 bug，选择路径 A 的理由详见项目根目录 `README.md` 及开发笔记。物种数上，本实现比 paper 更精细（6 种 vs 5 种，OC 拆为亲水/疏水）。

---

## 10. 关键代码说明

### 10.1 Fortran 辐射引擎 (`fortran/cfram_rrtmg.f90`)

主程序流程：

```fortran
! 对每个格点 (ilat, ilon)：
! 1. 确定有效大气层数 nlayer (基于地表气压)
! 2. 调用 rad_driver × 11 次：
!    - base 态 + warm 态 → 总变化
!    - 逐项替换 (CO₂, q, ts, O₃, solar, albedo, cloud, aerosol) → 各项 forcing
! 3. 调用 calc_drdt → 构建 Planck 矩阵 (∂R/∂T)
! 4. 矩阵求逆 → drdt_inv
! 5. 输出 9 项 forcing + drdt_inv
```

### 10.2 Python 并行运行器 (`scripts/run_parallel_python.py`)

```python
def process_column(args):
    ilat, ilon = args
    # 1. 提取单柱数据
    # 2. 计算气溶胶光学属性
    # 3. 写临时二进制文件到 tmpdir
    # 4. 运行 cfram_rrtmg_1col
    # 5. 读取 forcing + drdt_inv
    # 6. dT = -drdt_inv @ frc  (所有项统一)
    # 7. 清理 tmpdir，返回结果

with Pool(nproc) as pool:
    for ilat, ilon, result in pool.imap_unordered(process_column, tasks):
        # 收集结果到全场数组
```

### 10.3 气溶胶光学属性计算 (`core/aerosol_optics.py`)

从 MERRA-2 气溶胶混合比 → RRTMG 波段分辨的 AOD/SSA/g：

```
mixing_ratio × air_density × layer_thickness × mass_extinction_coeff → AOD per band
```

查找表来自 GOCART（opticsBands_*.nc），包含 30 个波段（14 SW + 16 LW）的消光系数，依赖相对湿度。

当前限制：paper_data 提供 6 种聚合物种（bc, ocphi, ocpho, sulf, ss, dust），但 GOCART 有 15 个 bins（DU001-005, SS001-005 等）。每种聚合物种映射到单一代表性 bin，这是气溶胶精度的主要瓶颈。

---

## 11. 性能指标

| 指标 | 值 |
|------|------|
| 格点数 | 81 × 121 = 9,801 |
| 并行度 | 80 进程 (hqlx204, 256 核) |
| 单 case 耗时 | ~19 分钟 |
| 单格点耗时 | ~1 秒 |
| 吞吐率 | 8.6 pts/s |
| 内存 | ~2 GB/进程（主进程加载全场数据 ~500 MB） |

---

## 12. 后续开发方向

已完成项（移出本节）：
- ✅ 气溶胶 15-bin → 6 种独立物种映射（已经是每物种独立光学属性 + RRTMG `iaer=10`）
- ✅ ERA5+MERRA-2 自动预处理（`build_case_input.py` + `case.yaml` 驱动）
- ✅ 支持任意气压层数（runtime 从 `data_prep/plev.dat` 大小推导 nlev，单一 binary 处理任意网格）
- ✅ CMIP6 接口（`scripts/build_cesm2_official.py` + `cases/cesm2_4xco2*` 系列；hybrid σ → plev 重网格 + Phase A 1850 climatology O3 注入）
- ✅ 单柱 climlab 验证框架（见 §14）+ 高阶 CFRAM 选项（midstate Planck、centered FD、co2_handling=midstate）
- ✅ Fu 引擎接入 + dual-MC sub-column 重叠 bug 修复（apple-to-apple 对齐 OLD CFRAM）

待办：

1. **路径 B 备用模块**（低优先级）：实现离线气溶胶 forcing 预计算 + Planck 求逆，与 Wu et al. 一比一数值对齐，供对照实验用
2. **CO2 / O3 空间变化**：当前 CO2 用常数（paper_data 反推），O3 用 ERA5；升级为 MERRA-2 M2I3NVCHM CO2 或 CAMS 3D 场，O3 切到 MERRA-2 减少数据源异质性
3. **f2py 封装 RRTMG**：避免临时文件 I/O，预计单格点速度提升 2-3 倍
4. **2 阶 CFRAM**：单柱验证（§14）显示 1 阶 CFRAM 在 RCE 上有 ~0.2 K mid-trop floor，源于对流调整 mass redistribution + R_TTT + CC 二阶；要进一步逼近真值需显式 Hessian
5. **敏感性实验框架**：单柱 ablation（关闭 cloud / aerosol / WV 各一项）

---

## 13. 常见坑与排错

### 13.1 cloud dT 符号翻转 / lhflx 量级偏小 ~6×

**症状**：
- validate_vs_paper.py 显示 `cld` 为负值（paper 为正）
- `lhflx` forcing domain mean ≈ −1 W/m²（paper ≈ −7）

**根因**：ERA5 CDS 返回的 6-hourly accumulation 字段**每步只含 1 小时累积**，并非完整的 6 小时累积。`data/era5_source.py:_sixhourly_accum_to_daily_wm2()` 中必须保留 `×6` 外推因子：

```python
daily_total = reshaped.sum(axis=1) * 6    # ← 关键的 ×6
return daily_total / SECS_PER_DAY
```

**自检**：build 阶段日志中 solar ~415 W/m²、lhflx ~−7 W/m²、shflx ~−5 W/m² 属正常；若 solar ~70 W/m² 则说明 ×6 缺失。

### 13.2 CO2 偏差导致 `dT_co2` 数量级错误

**症状**：`dT_co2` 与 paper 偏差 > 50%。

**根因**：`case.yaml` 中 CO2 值使用了今日大气浓度（~415 ppmv），未对齐 Wu et al. 使用的气候态值。

**修复**：
- EH22: `base_ppmv: 394.23, perturbed_ppmv: 416.89`
- EH13: `base_ppmv: 395.54, perturbed_ppmv: 396.14`

（这些值从 `paper_data/input_check.nc` 的 co2 字段 domain mean 反推。自建 case 需另行确定气候态 CO2。）

### 13.3 Warm days 自动检测与 paper 不一致

**症状**：`build_case_input.py` 打印的 warm period 范围与 paper 描述差 1-2 天。

**修复**：在 `case.yaml` 的 `temporal` 段手动指定 `warm_days` 列表（0-indexed August day），覆盖自动检测：

```yaml
warm_days: [2,3,4,...,23]    # Aug 3-24
```

### 13.4 Surface 层读取 index 错误

**症状**：`dT_cloud` 等在地表的数值看起来只有几个百分之一，远小于预期。

**根因**：`cfram_result.nc` 中 `lev` 维度**surface 在最后（`lev[-1] = 1013 hPa`）**，不是首位。Paper 的 `partial_t.nc` 则是 surface 在首位（`lev[0] = 1013 hPa`）。

**修复**：读地表层时，独立结果用 `d[-1, :, :]`，paper 结果用 `d[0, :, :]`。`validate_vs_paper.py` 已区分处理。

### 13.5 netCDF4 / GLIBC 版本错误

**症状**：`ImportError: GLIBC_2.25 not found`。

**根因**：机器 GLIBC 过旧（如 RHEL 7 早期）。

**修复**：换一台新机器跑，或用 conda 环境并设置 `LD_LIBRARY_PATH` 指向 conda 的 libc。

### 13.6 Fortran OpenMP 并行崩溃 / 静默错误

**症状**：给 `cfram_rrtmg.f90` 加 `!$OMP PARALLEL DO` 后多线程运行无输出，或输出全 0。

**根因**：RRTMG 内部 module 有 SAVE 变量，线程不安全。

**修复**：不要用 OpenMP。使用 `run_parallel_python.py` 的 multiprocessing 方案——每个 worker 进程独立调用 `cfram_rrtmg_1col`，进程隔离规避线程竞争。

### 13.7 CDS / Earthdata 下载中断

**症状**：`download_*.py` 在下到一半时网络超时。

**修复**：两个下载脚本均支持断点续传（检查目标文件存在性跳过）。直接重跑即可：

```bash
python3 scripts/download_era5_flux.py
python3 scripts/download_merra2_aerosol.py
```

### 13.8 MERRA-2 气溶胶插值后出现负值

**症状**：`base_pres.nc` 里 bc/oc 等字段含极小负值（~ −1e-12）。

**根因**：log-p 插值 + 浮点舍入。

**修复**：`merra2_aerosol.py` 的 `load_merra2_aerosol()` 末尾已执行 `np.maximum(result, 0.0)` 强制非负。无需手动处理。

### 13.9 HPC 集群上 LAPACK/BLAS 未加载导致 Fortran 瞬间退出

**症状**：`run_case.py` 在"run"阶段几秒内完成（正常约 20 分钟），日志显示吞吐率异常高（>1000 pts/s），`cfram_result.nc` 中各项分解结果全为 NaN 或 0。

**根因**：HPC 集群上 LAPACK/BLAS 通常由 `module` 系统管理，不加载时动态库不在 `LD_LIBRARY_PATH` 中，`cfram_rrtmg_1col` 启动即退出（无报错输出）。Python worker 未检查 Fortran 子进程返回码，将失败计入"done"，导致结果全为空。

**排查**：

```bash
# 检查动态库是否可找到（有 "not found" 则说明缺库）
ldd fortran/cfram_rrtmg_1col | grep -i 'not found'

# 直接运行可执行文件查看错误信息
cd fortran && ./cfram_rrtmg_1col
```

**修复**：在编译和运行前加载对应的 LAPACK 模块（具体名称因集群而异）：

```bash
module load lapack       # 或 module load openblas / intel-mkl 等，视集群而定
module load blas
```

也可在 `fortran/makefile` 中将 LAPACK 改为静态链接（`-l:liblapack.a -l:libblas.a`），一次编译后无需每次加载模块。

**注意**：`§4.1` 环境要求表中"LAPACK/BLAS 系统自带即可"的说法适用于普通 Linux 桌面，在 HPC 集群上须显式加载模块。

---

## 14. 单柱 climlab 验证方法

### 14.1 动机

实地观测/再分析数据驱动的 CFRAM 分解（EH13/EH22、CESM2 4×CO2）受多种误差源混合：观测不确定性、气溶胶光学近似、Manabe 假设、对流参数化等。当 closure (`dT_obs = Σ dT_X + dT_dry`) 残差大时,难以独立判断错误来自数据还是 CFRAM 方法本身。

`cases/climlab_4xco2` 提供一个**完全可控的理想验证**：用 climlab 软件包跑两个 clear-sky RCE 平衡态(1×CO2 vs 4×CO2 Manabe RH-fixed),把这两个平衡作为 base/warm 喂给 pyCFRAM。这是一个 1×1 单柱、纯辐射 + 干对流调整的简化系统,没有云、气溶胶、动力,因此 CFRAM 的闭合性、Planck 矩阵正确性、forcing 路径选择都暴露得很清楚。

- ΔTs ≈ +4.59 K (climlab 收敛态)
- 30-level 默认 σ-grid (climlab),pyCFRAM runtime nlev 自适应,无需重编 Fortran
- 单柱 RRTMG 跑时 ~1.5 s,Fu ~2 s

### 14.2 执行

```bash
# 1. 跑 climlab 双态 RCE,写入 cases/climlab_4xco2/input/ 和 cases/climlab_4xco2_fu/input/
/path/to/conda/python experiments/climlab_validation/run_rce_4xco2.py

# 2. CFRAM 分解 (1×1 自动识别,无 mp.Pool 启动开销)
python3 run_case.py climlab_4xco2     --step run     # RRTMG
python3 run_case.py climlab_4xco2_fu  --step run     # Fu

# 3. 4-panel 垂直 profile (Σ closure + co2 + q + dry)
python3 scripts/plot_singlecol_profile.py climlab_4xco2                    # RRTMG only
python3 scripts/plot_singlecol_profile.py climlab_4xco2 climlab_4xco2_fu   # 双引擎 overlay
```

### 14.3 高阶 CFRAM 选项 (仅适用于 RCE clear-sky 单柱)

1 阶 CFRAM 在 clear-sky 单柱中线性化误差成为主要误差源。pyCFRAM 提供 4 个 opt-in `radiation:` flag 逐步消除 Taylor 展开中的特定项。**所有 flag 默认关闭**,且**绝不应在 ERA5/CESM2 真实大气 case 上启用** —— 它们假设辐射近平衡、无非辐射耦合。

| Flag | 默认 | 作用 |
|------|------|------|
| `drdt_eval: midstate` | `base` | Planck 矩阵在 (T_base+T_warm)/2 算,cancel 单变量 R_TT 项 (1-变量 2 阶 Taylor)。+1 rad_driver/cell |
| `drdt_probe: centered` | `onesided` | calc_drdt 用 T_j ± 0.5K centered FD 替代 +1K 单边,cancel J 矩阵构造内部的 R_TT 污染。calc_drdt 内 rad_driver_lw 翻倍 |
| `co2_handling: midstate` | `base` | `frc_co2 = R(T_mid, q_mid, co2_warm) − R(T_mid, q_mid, co2_base)` 在 mid-state 大气算,cancel ∂²R/∂T∂C cross-term。+2 rad_driver/cell |
| `q_handling: midstate` / `feedback` | `independent` | 仅诊断用。`midstate` 在 Manabe RH-fixed RCE 上让 upper-trop 闭合**变差 3×** —— baseline 的 R_Tq cross-term 恰好与约束 RH 解对齐,cancel 反而丢失"正确"贡献。`feedback` 翻转 R_Tq 符号(单边 warm 路径) |

数学背景:midstate 对称差展开只消单变量偶次项 (R_TT、R_qq、R_CC),混合 cross-term (R_TC、R_Tq) 仍存活。Cancel mixed term 需要把对应的 `frc_X` **和** baseline 都搬到 midstate (即 `co2_handling: midstate`)。详见 `temp.md`。

### 14.4 已验证闭合性

`climlab_4xco2` 上 `dT_obs − (dT_co2 + dT_q + dT_ts + dT_dry)` 按高度带 RMS:

| 配置 | mid-trop (450–650 hPa) | upper-trop (200–400 hPa) | strato (<200 hPa) |
|------|:----------------------:|:------------------------:|:-----------------:|
| baseline 1 阶 | 0.45 K | 0.07 K | 0.30 K |
| + `drdt_eval: midstate` | 0.45 K | 0.07 K | **0.28 K** (stratosphere 闭合) |
| + `co2_handling: midstate` | **0.24 K (-47%)** | 0.04 K | 0.28 K |
| + `drdt_probe: centered` | **0.21 K (-55%)** | **0.02 K (-68%)** | **0.25 K** |

剩余 ~0.2 K mid-trop 残差是 1 阶 CFRAM 在该 setup 上的**物理 floor**:climlab 干对流调整的 mass redistribution(非辐射、本质上在 J⁻¹·ΔR 框架之外)、R_TTT 高阶 T 项、Clausius–Clapeyron Δq 二阶项。要继续下降需 2 阶 CFRAM 显式 Hessian —— 算法结构需重新设计,工程量较大。

### 14.5 引擎交叉验证 (RRTMG vs Fu)

两个引擎用相同的 climlab 输入跑出来,`dT_co2` / `dT_q` 个体量级差 10–15 %(不同 RT solver),但总 ECS-equivalent 响应在数值噪声内一致。这是独立于 OLD CFRAM 参考的有用 cross-check —— 任一引擎单独的 bug 都不会同时影响另一个。

---

## 15. Phase 3：DAMIP 多模式支持

### 15.1 目标与架构

在 ERA5+MERRA-2、单一 CESM2 4×CO2 数据集之外，pyCFRAM 新增一个注册式数据源
`data/cmip6_damip_source.py`（`@register_source('cmip6_damip')`），使其能吃
**任意 CMIP6 模式**的**任意 DAMIP 单强迫实验**（`hist-aer`/`hist-GHG`/
`hist-nat`/`hist-stratO3`……）。架构上完全复用已有链路：

```
case.yaml (source.type: cmip6_damip)
  → run_case.py --step build（通用分发：非 cesm2_cmip6/None 都走这条路）
  → scripts/build_case_input.py（SOURCE_MODULES 注册表触发 import 拿到插件）
  → data/cmip6_damip_source.py :: build_states()
        ├─ data/cmip6_common.py       （模式无关数值/探测机械）
        ├─ configs/damip_models.d/<model>.yaml   （模式专属怪癖，可选）
        └─ configs/damip_experiments.yaml        （实验语义，可选）
  → validate_states() 写盘前拒 NaN → 4 个 pres/surf NC + nonrad_forcing.nc + provenance.json
  → run_case.py --step run（scripts/run_parallel_python.py + fortran/cfram_rrtmg_1col，零改动）
  → cfram_result.nc → scripts/write_run_summary.py → <case>.summary.txt
```

`core/` 与 `fortran/` 全程零改动——这是合同 M5 验收口径的硬性要求，由
`tests/test_damip_userguide_example.py::test_phase3_diff_never_touches_core_or_fortran`
对整条 Phase 3 分支做 `git diff --name-only` 机械化断言。完整架构图、IO 契约
（`build_states()` 返回值、`provenance.json` 字段）、`build_states()` 十步流程
逐步说明见 `docs/m4_damip_module.md`；每模式支持矩阵与已知问题库见
`docs/m5_multimodel_userguide.md`。

### 15.2 跨模式异构处理：真实踩坑记录

13 个候选 hist-aer 模式里只有 3 个发布全部所需变量——**缺云/缺 O₃/缺 rsdt 是主
路径，不是边缘情况**。以下每一条都是本阶段对真实下载的 ESGF 数据实测踩到的坑，
不是预想的假设场景（呼应 §13 的踩坑记录风格）：

- **日历**：`cmip6_common.decode_time()` 用 `cftime.num2date()` 替代早期
  CESM2-only 的 `days/365.0` 算法（后者对除 noleap 外的日历全部静默算错）。M5
  八模式跨 5 种日历：noleap（CESM2）、365_day（GISS/CanESM5）、gregorian
  （IPSL/MIROC6）、proleptic_gregorian（MRI）、360_day（HadGEM3，每月恒 30 天）。
- **IPSL-CM6A-LR 的 hybrid 垂直坐标，一次性踩了三个坑**：
  1. `cl`/`clw`/`cli` 文件**完全没有 `formula_terms` 属性**，自动探测无字符串
     可读——`configs/damip_models.d/IPSL-CM6A-LR.yaml` 显式给
     `vertical: {scheme: hybrid_ap_b, ap: ap, b: b}` 覆盖。
  2. `ap`/`b` 系数发布在**层界面**（`nlev+1` 个值）而非**层中点**（`nlev` 个值，
     与数据本身的层维度一致）——`_read_hybrid_coeffs()` 相邻界面值配对平均，
     这是**精确**而非近似的处理：因为 `p = ap·p0 + b·ps` 对 `ap`/`b` 是线性的，
     两个界面气压平均恰好等于该层的中点气压。
  3. IPSL 的混合层顺序是**地表→TOA**（`b` 从 ~1 降到 ~0），与 CESM2 原生
     `hyam`/`hybm`（TOA→地表，`hybrid_to_plev_mass_conserving` 唯一验证过的
     约定）相反。`_read_hybrid_coeffs()` 从 `b` 的单调方向自动探测并反转系数
     数组；调用方必须用**同样的方式**反转云场数据本身——只反转一侧会把云廓线
     悄悄打乱成一个"看起来合理但物理错误"的近零云场（本 bug 最初不是崩溃发现
     的，而是靠对照真实原始数据的云量数量级做 sanity check 才定位到）。
- **CNRM-CM6-1 的非标准属性名**：`lev` 坐标携带 `formula_term`（单数），不是
  标准 CF 的 `formula_terms`（复数）——CERFACS 的 CMOR 配置怪癖。自动探测找不到
  该字符串，静默落到错误的（当作已在 plev 上）分支。同样用 models.d 的
  `vertical:` 覆盖解决。
- **HadGEM3-GC31-LL 的 hybrid-height 坐标（真正无法支持的坐标类型）**：
  `formula_terms = "a: lev b: b orog: orog"`，即 `z = a + b·orog`——**高度**混合
  坐标而非**气压**混合坐标。`detect_vertical()` 对无法识别的 `formula_terms`
  正确抛出 `ValueError`；`cmip6_damip_source.py` 捕获这个异常，降级为
  `cloud=SKIPPED` 而不是让整个 build 崩溃。真正支持 hybrid-height 需要模式自身
  温度廓线做高度→气压反演，是另一种物理转换，不是系数命名的变体——列为
  known-issue（`docs/plan_ph3.md` §12 R4），不强行适配。
- **MRI-ESM2-0 的 `o3` 发布在比同模式其它变量更粗的原生网格上**（64×128 vs
  ta/hus/cl 等的 160×320——CMIP6 并不要求同一 table 内所有变量共享一套网格）。
  加了 `cmip6_common.regrid_horizontal_bilinear()`（双线性、经度周期）。这个
  regrid 还牵出**第二个坑**：CMIP6 标准的下层大气缺测掩膜（高地形区，如青藏
  高原上空 1000 hPa 层为 NaN）如果不先处理，双线性插值会把 NaN 扩散到相邻的
  **有效**格点——不只是原本缺测的点。`cmip6_common.
  fill_nan_hold_toward_surface()` 在 regrid 前把地表端的 NaN 段用最浅有效层
  的值填住，regrid 之后，最终、真正 ps-aware 的 `fill_subsurface()`（用**目标**
  网格自己的 ps）仍然照常在下游跑一遍——前者只是 regrid 这一步的 NaN 传播
  防护，不是最终的子地表处理。
- **单位换算**：`cl` 是 %（÷100 转 0–1 分数）；`o3` 是 mol/mol（×48/29 转
  kg/kg，`VMR_TO_MMR`，与既有 `scripts/inject_cesm_o3.py` 同一常数）。O₃ 单位
  换算漏做是 **5 个数量级**的错误，不是小误差——正是 2026-05-10/12
  session_log 记录的"huss-as-O3"数据损坏事故的近亲（协作者的 Fortran CFRAM
  曾把比湿误当 O3 喂进辐射方案）。`cmip6_damip_source.py` 的 O₃ 分支刻意让这
  类错误难以复现：换算只在一处、读入 mol/mol 字段后立即执行一次。

### 15.3 ESGF 下载客户端（`data/esgf_fetch.py` + `scripts/download_damip.py`）

纯 stdlib（只用 `urllib`，不依赖 `requests`/`intake-esgf`/`esgpull`），免认证
下载（DKRZ/CEDA 经典 Solr 搜索 + HTTP Range 请求）。三个真实踩坑：

- **未过滤 variant 导致下载体积膨胀**：IPSL-CM6A-LR 一个模式就发布 10 个
  `hist-aer` variant，不加过滤的 `--model IPSL-CM6A-LR` dry-run 报出 54.8 GB
  （应为约 20.4 GB/variant）。修复：`download_damip.py` 的 `MODEL_DEFAULTS` 表
  给每个模式一个默认 `(variant_label, grid_label)`——代码注释明确写"这不是可选
  的体积优化，遗漏它是正确性 bug"。
- **跨节点 replica 重复，且 Solr 的 `replica` 标志不可靠**：同一物理文件常被
  多个 ESGF 数据节点索引（如 CNRM 的家节点 + DKRZ 镜像）。直觉修法——按
  `replica=false` 过滤——**行不通**：部分模式（实测：MRI-ESM2-0）在 DKRZ/CEDA
  上全部被标记 `replica=true`（因为这两个节点都不是该模式的家节点），
  `replica=false` 会返回**零结果**，即使数据其实可下载。改用文件名去重
  （`search_damip_files()`），并**保留每个文件名对应的所有已知 URL**（而非只
  留第一个），为下面的镜像回退铺路。
- **`checksum`/`checksum_type` 是 Solr 多值字段**：返回的是列表（如
  `checksum_type=['SHA256']`）而非标量。这个 bug 只有在一次真实的完整下载
  （1.5 GB）走到校验步骤时才会暴露——dry-run 根本走不到这一步。
- **CNRM-CM6-1 家节点超时，镜像可用**：`esg1.umr-cnrm.fr` 对 `cl` 的 17 个
  分片文件中 7 个超时；`esgf3.dkrz.de` 完全可用。`download_files()` 现在按序
  重试每个已知镜像 URL，而非第一个 URL 失败就放弃。
- **`http.client.RemoteDisconnected` 不是 `urllib.error.URLError` 的子类**：
  在为 WP-M5.3 准备 NorESM2-LM dry-run 时，`list_files()` 的 `urlopen()`
  调用因此未捕获崩溃。`esgf_fetch.py` 里所有原本只捕获 `URLError` 的
  except 子句（`search_datasets` ×2、`list_files`、`fetch`）现在都同时捕获
  `OSError`（commit `15a2845`）。

### 15.4 如何接入新模式（NorESM2-LM 实测例子）

WP-M5.3 用 NorESM2-LM（不在 M4/M5 八模式清单内）验证了"接入新模式只需
`configs/damip_models.d/<model>.yaml` + `cases/<case>/case.yaml` 两个 yaml 文件，
零 Python 改动"这条路径。2026-07-06 hqlx210 实测：build 一次成功；run 87.5 秒
（13824 格点，NorESM2-LM 原生 96×144 网格）；`cfram_result.nc` 零 NaN；恒等式
`dT_sfcdyn = dT_ocndyn + dT_lhflx + dT_shflx` 残差 ~4.6e-15 K（机器精度）；全球
均值 `dT_observed[sfc] = -1.05 K`（净冷却，符号正确），北半球 `-1.72 K` 冷于
南半球 `-0.37 K`（符合人为气溶胶集中北半球的预期）。完整 yaml 内容 + 永久回归
门 `tests/test_damip_userguide_example.py` 说明见
`docs/m5_multimodel_userguide.md` §"如何接入新模式"一节。

### 15.5 覆盖率与测试

```bash
coverage run -m pytest tests/ -q
coverage report -m --include="*/data/cmip6_common.py,*/data/cmip6_damip_source.py,*/data/esgf_fetch.py"
```

2026-07-06 复核（新增 regrid/hybrid-height try-except/多镜像回退/OSError 捕获
等逻辑之后）：`cmip6_common.py` 98%、`cmip6_damip_source.py` 92%、
`esgf_fetch.py` 99%，均远高于合同 60% 硬指标，相对 WP-M4.6 首次测量的
98%/93%/98% 无实质回退。

