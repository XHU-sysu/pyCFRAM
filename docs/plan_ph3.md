# pyCFRAM Phase 3 执行计划 v2（M4 DAMIP 基础接入 + M5 多模式扩展与自定义接口）

> **版本**：v2（2026-07-05）。v1 由 Opus 尽调起草；v2 为督导审计修订版——**逐条实测核对**了仓库代码、合同原文、ESGF 三个索引节点（DKRZ/CEDA/LLNL）、13 个候选模式的 hist-aer 变量矩阵、免认证下载链路、本地与 hqlx210 环境。修正 **6 处事实错误 / 空白**，其中两处会直接导致 v1 按图施工失败（8 模式清单中 5 个缺云变量；7/13 模式缺 rsdt 而 v1 无 solar 兜底）。
> 交付物：本文件（`docs/plan_ph3.md`）。给 **执行者（lead session）+ sonnet / haiku 级 subagent** 的作业书：读完即可分阶段落地。
> 覆盖合同 `contract/contract_10.pdf` 的 **Phase 3**：M4（第 13–15 周，20%）+ M5（第 16–19 周，17%）。
> 上游：Phase 2 已完成并本地提交（commit `9401a7e`/`e51cf47`，分支 `feat/m2-m3-lapse-rate-kernel`，**PR 尚未开**——见 §8 WP-0 的 PR 叠链处理）。
> **已拍板（用户 2026-07-05 确认）**：主 experiment = **hist-aer**；ESGF **从零下载**；模式清单原则上用推荐清单——v2 依据实测数据可得性做了**必要调整**（见 §1.4，调整理由全部有 ESGF 检索实证）。

---

## v2 相对 v1 的修订摘要（督导审计结论，先读这个）

| # | v1 的问题 | v2 修订 | 证据 |
|---|---|---|---|
| 1 | **8 模式清单未做可得性核实**：v1 清单中 CESM2/CanESM5/NorESM2-LM/ACCESS-ESM1-5/GFDL-ESM4 五个模式的 hist-aer **没发布 cl/clw/cli**（云项全跳），其中 4 个连 **rsdt** 都没有 | 清单重排（§1.4）：M4 三模式 = IPSL-CM6A-LR + MRI-ESM2-0 + CESM2；M5 八模式按"6 个全变量 + 2 个缺变量演示 skip"配置 | 2026-07-05 DKRZ+CEDA 双索引联查，逐模式变量 facet 计数（§1.4 表） |
| 2 | **rsdt 缺失无兜底**：v1 §6 决策树只覆盖云/O3/气溶胶/通量/rsus+rsds，没有 solar 的缺省路径，而这是 M4 首批（CESM2）就会踩的主路径 | 新增 `solar=ANALYTIC` 兜底（解析式年均 TOA 入射）；hist-aer 太阳强迫冻结→两态 solar 相同→frc_solar≡0，兜底不伤分解 | CESM2/NorESM2/ACCESS/GFDL/FGOALS hist-aer Amon 均无 rsdt |
| 3 | **intake-esgf 引入不必要**：本地与 hqlx210 均未装；hqlx210 conda base 是 NFS 共享环境不宜随意装包；且 CMIP6 下载**免认证**（实测 HTTP 206 + HDF5 magic，无需甲方 ESGF 账号） | `esgf_fetch.py` 用**纯 stdlib（urllib）**实现 search+下载：DKRZ/CEDA 经典 Solr API 实测可用（LLNL 已下线重定向）；intake-esgf 降为可选。"等甲方账号"不再是前置 | 实测：`esgf-data.dkrz.de` / `esgf.ceda.ac.uk` 200 OK（hqlx210 直连 2s）；MRI ts 文件 Range 请求免认证 206 |
| 4 | **build 后 hook 设计与写盘器校验冲突**：v1 §6 允许"O₃/子地表掩膜作为 build 后可选 hook"，但路径 A 的 `build_case_input.validate_states` 在**写盘前**拒绝一切 non-finite 值——NaN 留到 hook 阶段根本到不了写盘 | 子地表 HOLD 填充、O₃ 注入、单位转换**全部收进 `build_states()` 内部**完成；不存在 build 后 hook | `scripts/build_case_input.py:61-77` |
| 5 | **仓库落地细节空白**（Sonnet 会卡住的坑）：`build_case_input.py:239` 硬编码 `import data.era5_source` 触发注册；case.yaml 必须含 `input:` 块；huss 不在路径 A 写盘清单且不能无脑加；nonrad 符号两条路径表述相反；CMIP6 `cl` 是 %、`o3` 是 mol/mol 需转换 | §2 逐条写明 + §4.1 给出**完整可抄** case.yaml；§6.3 nonrad 符号统一律 | 见 §2 各条的 file:line 引用 |
| 6 | **PR 节奏与合同不符**：v1 只在 M5 末开一个 PR；合同验收流程是"甲方收到**各检查点** PR 后 7 个工作日 Code Review"——M4 需要自己的 PR；且 Phase 2 的 PR 还没开，分支叠链需要先解 | 新增 WP-M4.7（M4 PR）；WP-0 增加 PR 叠链决策（先开 Ph2 PR，M4 PR 以其为 base 或等其合并） | 合同第三条"验收流程"+ 第十一条"甲方评估与回复时限" |

其余 v1 内容（DAMIP 协议解读、base↔perturbed 配对、路径 A 架构选型、cmip6_common 抽取 + 回归金标、models.d 防冲突、决策树 + summary.txt 思路）经核实**全部成立**，v2 原样保留并细化。

---

## Context（为什么做这件事）

pyCFRAM 目前只能吃 ERA5 再分析、MERRA-2 气溶胶、以及一套 CESM2 CMIP6 4×CO₂ 数据（`data/cesm2_cmip6_source.py` + `scripts/build_cesm2_official.py`，硬编码、未走插件注册）。合同 Phase 3 要求把工具升级为 **CMIP6 DAMIP 单强迫（single-forcing）多模式**可用：给定一个 DAMIP 实验（hist-aer / hist-GHG / hist-nat / hist-stratO3 …），自动配对 base↔perturbed 两态、跨模式吸收异构（不同日历 / 垂直坐标 / 压力层 / variant 标签 / 缺变量），端到端跑出 CFRAM 分解。M5 进一步要求"**用户只改 `case.yaml`（+ 可选 `data/<source>_source.py`）就能接新模式，不碰 `core/` 与 `fortran/`**"。

**期望结果**：新增一个干净的注册式数据源插件 `data/cmip6_damip_source.py`，把已有 CMIP6 机械（hybrid→plev、气候态、反照率、O₃ 注入、子地表填充）抽成模式无关的 `data/cmip6_common.py` 复用；下游 CFRAM 引擎（Fortran + `run_parallel_python.py`）**零改动**。M4 用 hist-aer 跑通 ≥3 模式，M5 扩到 ≥8 模式（≥6 通过、其余文档化 known-issue），并交付用户自定义接口 + 全套文档。

**合同边界（重要，减负条款）**：
- "生产规模的批量计算"**不在服务范围**——Phase 3 仅做**最小数据集功能验证**（每个 experiment 取若干年份子集气候态），不做多年多情景科学跑批。模式 native 网格（~1-2°，CESM2 为 192×288=55296 列）单气候态两态在 hqlx210 上 ≈7 分钟/case，完全可行。
- 每个检查点验收 = **PR + 甲方 7 个工作日 Code Review**（逾期未回复视为通过）。
- **检查点独立性**：M4 启动不依赖 M3 验收通过；M3 文档 latest draft 即可作为 M4 残差解释引用源。

---

## 0. TL;DR（30 秒版）

**要做的事**：给 pyCFRAM 加一个 **DAMIP 数据源插件**，让它能吃任意 CMIP6 模式的单强迫实验。

- **架构落点**：`data/cmip6_damip_source.py` 实现 `@register_source('cmip6_damip')` 的 `DataSource` 子类，`build_states()` 返回 `(base_state, perturbed_state, nonrad_forcing)` 三字典 → 走**已有的通用写盘器** `scripts/build_case_input.py`。**不写新 build 脚本、不改 `core/`、不改 `fortran/`**。
- **base↔perturbed 配对**：同一单强迫实验的**首个十年（1850–1859）气候态 vs 末个十年气候态**。自洽、契合合同数据口径。
- **跨模式异构吸收**（真正的工程量）：日历、垂直坐标、发布压力层、variant、grid_label、经纬序、单位、**缺变量（云/rsdt/o3 缺失是主路径不是边缘）**——全部 config 驱动 + 运行时探测。
- **可选量缺省处理**：缺云→云项跳过；缺 o3→注入 CESM 1850 气候态；缺 rsdt→解析式 solar；缺通量→无 nonrad。每次 build+run 落 `*.summary.txt`。

**M4 验收硬指标**（合同原文）：≥3 个 DAMIP 模式（**合同建议清单**：CESM2、IPSL-CM6A-LR、MRI-ESM2-0、CanESM5、NorESM2-LM 中选）× ≥1 个 single-forcing experiment 端到端运行；气溶胶/臭氧可选缺省逻辑 + 每次运行生成 `*.summary.txt`（引用 M3 文档结论）。**v2 选型：IPSL-CM6A-LR、MRI-ESM2-0、CESM2**（前两个全变量，CESM2 演示 skip 路径；三个都在合同建议清单内）。
**M5 验收硬指标**：≥8 模式尝试、≥6 端到端通过，失败模式 README 中文 known-issue；用户自定义接口（仅 yaml + 可选 source.py，core/fortran 零改动）+ ≥1 可运行示例；新模块 pytest 行覆盖率 ≥ 60%；文档覆盖 Phase 2/3 全部新增功能。

**进场先过 3 道门**：① **回归金标**（重构 CESM2 机械后 `cesm2_4xco2_official` build 输出不变——金标只能在 hqlx210 生成，见 WP-M4.1）；② ESGF 首个 hist-aer 子集在 hqlx210 落盘并能 netCDF4 打开（**不等甲方账号，免认证已实证**）；③ 首个 DAMIP case `--step build` 产出合规 4+1 个输入 NC。

---

## 1. 背景：DAMIP 是什么、CFRAM 拿它做什么

### 1.1 DAMIP 单强迫协议（Gillett et al. 2016, GMD）

DAMIP（Detection & Attribution MIP）在历史期（1850–2014/2020）跑一批**只让单个强迫因子随时间变化、其余全部冻结在 1850** 的实验：

| experiment_id | 时变的强迫 | 其余（冻结在 1850） | Tier | 备注 |
|---|---|---|---|---|
| **hist-aer**（★主验收） | 人为气溶胶 | GHG（含 CO₂）、O₃、太阳、火山、土地利用 | 1 | 发布模式最多、信号强；直接检验 pyCFRAM 气溶胶/云链路 |
| hist-GHG | well-mixed GHG（CO₂/CH₄/N₂O…） | 气溶胶、O₃、自然 | 1 | **CO₂ 时变** → config 给时变 CO₂ |
| hist-nat | 太阳 + 火山 | 全部人为 | 1 | 信号弱、火山年际噪声大 |
| hist-stratO3 | 平流层臭氧 | 其余 | 2 | **O₃ 时变** → 必须用模式自带 o3，不注入气候态 |
| hist-sol / hist-volc / hist-CO2 | 单一自然/CO₂ | 其余 | 2/3 | 发布少，M5 可选 |

**关键推论（决定各强迫项怎么处理）**：
- hist-aer：CO₂/O₃/太阳两态**相同**（1850 值）→ `co2_base == co2_warm`（284.7 ppm，与 `cesm2_4xco2_official` 的 `co2_pictrl_ppm` 同源）；O₃ 用模式值或注入气候态（两态同场，frc_o3≡0）；**solar 两态相同 → frc_solar≡0，因此 solar 缺失可以用解析兜底而不伤任何 dT_X**；气溶胶本身是信号，但 Amon 里没有 3D 气溶胶 mmr（在 AERmon 且非 GOCART 6 种）→ pyCFRAM 的气溶胶分解**跳过**（M4 交付物 2 的核心场景，runner 已有自动检测：`run_parallel_python.py:639` `aer_max < 1e-15 → skip_aerosol`）。
- hist-GHG：CO₂ 两态**不同**（1850 vs 末十年）→ config 给 `base_ppmv/perturbed_ppmv`。
- hist-stratO3：O₃ 两态**不同** → **必须用模式 o3**，不能注入 1850 气候态（否则 frc_o3≡0，抹掉信号）。

### 1.2 base↔perturbed 配对（CFRAM 两态法的落地）

CFRAM 比较两个大气态。为**隔离单个强迫因子的历史信号**，取**同一实验**的：
- **base** = 首十年（默认 1850–1859）月气候态；
- **perturbed（warm）** = 末十年气候态。**per-model 默认**：时间到 2020-12 的模式用 [2011, 2020]；只到 2014-12 的模式（CESM2、GISS-E2-1-G，实测确认）用 [2005, 2014]。写进各自 `configs/damip_models.d/<model>.yaml`。

则 `dT_observed = warm.ta − base.ta` ≈ 该单强迫因子在历史期造成的温度变化。**自洽**（同一模式、同一实验、同一 variant，无跨实验漂移/单位对齐问题），且直接契合合同数据口径"每个 experiment 取若干年份子集（例如 1850-1859 与 last decade 月平均）"。

> **已否决的替代方案（决策记录）**：① base=piControl → 需额外下载、有模式漂移 detrend 问题；② base=historical(all-forcing) → 两态混入多种强迫。均比"同实验首末十年"更复杂且更易错，按奥卡姆取最简。

### 1.3 CFRAM 拿它做什么

把 `dT_observed` 喂进 pyCFRAM **已有的** RRTMG/Fu 引擎与 `run_parallel_python.py`，输出逐过程分解 `dT_q / dT_cloud / dT_co2 / dT_albedo / dT_atmdyn / dT_sfcdyn / dT_ocndyn / …`。**Phase 3 不碰引擎**——只把"任意 DAMIP 模式的两态"翻译成 `docs/input_spec.md` 规定的输入 NC。Phase 2 的辐射核 LR 分解（`run_case.py --step lr` / `--step lr-attr`，已实装于 `run_case.py:77-81`）对 DAMIP case 同样可用：case.yaml 加 `lapse_rate:` 块即免费复用。

### 1.4 ★ 数据可得性实测矩阵（2026-07-05，DKRZ+CEDA 双索引；v2 新增，选型依据）

对 13 个候选模式逐一检索 `CMIP6 / DAMIP / hist-aer / Amon` 的数据集计数。pyCFRAM 需要 14 个变量：`ta hus ts ps cl clw cli rsdt rsds rsus hfls hfss huss o3`。

| 模式 | 缺失变量 | ta 时间跨度 | variant | grid | 日历 | 判定 |
|---|---|---|---|---|---|---|
| **IPSL-CM6A-LR** | **无**（14/14 全，含 o3） | 1850–2020 | r1i1p1f1 ×10 | **gr** | gregorian | ★M4 首选（调试基准） |
| **MRI-ESM2-0** | **无**（14/14 全，含 o3） | 1850–2020 | r1i1p1f1 ×5 | gn | proleptic_greg. | ★M4 第二 |
| **CNRM-CM6-1** | **无**（14/14 全，含 o3） | 1850–2020 | **r1i1p1f2** ×10 | **gr** | gregorian | M5（f2 变体覆盖） |
| MIROC6 | o3 | 1850–2020 | r1i1p1f1 ×9 | gn | gregorian | M5 |
| GISS-E2-1-G | o3 | 1850–**2014** | r1i1p1f1 ×11（有 p1/p3） | gn | 365_day | M5 |
| HadGEM3-GC31-LL | o3 | 1850–2020 | **r1i1p1f3** ×55 | gn | **360_day** | M5（日历+f3 压力测试） |
| CanESM5 | cl, clw, cli, o3 | 1850–2020 | r1i1p1f1 ×30（有 p1/p2） | gn | 365_day | M5（云 skip 演示） |
| **CESM2** | cl, clw, cli, **rsdt**, o3 | 1850–**2014** | r1i1p1f1, r3i1p1f1 | gn | noleap | ★M4 第三（合同建议 + skip 演示） |
| NorESM2-LM | cl, clw, cli, rsdt, o3 | 1850–2020 | r1i1p1f1 ×3 | gn | noleap | 备选/known-issue 候选 |
| ACCESS-ESM1-5 | cl, clw, cli, rsdt, o3 | 1850–2020 | r1i1p1f1 ×3 | gn | proleptic_greg. | 备选 |
| GFDL-ESM4 | cl, clw, cli, rsdt, o3 | 1850–2020 | r1i1p1f1 ×1 | **gr1** | noleap | 备选 |
| BCC-CSM2-MR | cl（有 clw/cli）, o3 | — | r1i1p1f1 ×3 | gn | 365_day | 备选（三缺一→云仍跳） |
| FGOALS-g3 | cl, clw, cli, rsdt, o3 | — | r1i1p1f1 ×3 | gn | 365_day | 备选 |

**由此拍板（v2 修订 v1 清单的理由）**：
- **M4 三模式 = IPSL-CM6A-LR + MRI-ESM2-0 + CESM2**。前两个全变量（云、o3、rsdt 齐活，调试时任何异常都能排除"缺数据"因素）；CESM2 在合同建议清单内且恰好是"缺云+缺rsdt+缺o3"的 skip 全家桶——**M4 交付物 2（可选缺省处理）的天然验收样本**。v1 的 CanESM5 移到 M5（它也缺云，M4 阶段两个 skip 样本冗余）。
- **M5 八模式 = IPSL-CM6A-LR、MRI-ESM2-0、CNRM-CM6-1、MIROC6、GISS-E2-1-G、HadGEM3-GC31-LL、CanESM5、CESM2**。前 6 个云项活跃（全变量或仅缺 o3），后 2 个演示 skip；变体覆盖 f1/f2/f3、日历覆盖 gregorian/proleptic/365_day/360_day/noleap、网格覆盖 gn/gr——异构矩阵全覆盖且 "≥6/8 通过"有 6 个全变量模式打底，**几乎不可能不达标**。
- 备选池（谁挂了换谁）：NorESM2-LM、ACCESS-ESM1-5、GFDL-ESM4、BCC-CSM2-MR、FGOALS-g3——注意它们全部缺云，端到端仍可通过（skip 是特性），但科学演示效果打折。
- **v1 清单的问题**：8 个里 5 个缺云、4 个缺 rsdt——若按 v1 执行，M5 大半案例云项全零，且 M4 首个模式 CESM2 会在"没有 solar 兜底"上直接卡死。

> 检索脚本已验证可用（纯 stdlib，见 §7），执行期还可随时重跑刷新矩阵。**注意**：以上是 index 元数据；个别文件可能存在"索引在、data node 死"的情况，下载期用 §7 的多副本回退。

---

## 2. pyCFRAM 现状：DAMIP 往哪接（已逐文件核对，含 file:line）

### 2.1 两条 build 路径（DAMIP 走干净的那条）

```
路径 A（注册式插件，干净，ERA5 用）★DAMIP 目标：
  case.yaml source.type ──► core.config.load_case
                            data.source_base.get_source(cfg)      # 工厂 (source_base.py:113)
                            └─► @register_source('x') 的 DataSource 子类
                                  .build_states() → (base, pert, nonrad) 三字典
                            scripts/build_case_input.py            # 通用写盘器
                                  validate_states → write_pres_nc / write_surf_nc / write_nonrad_nc
                                  → cases/<case>/input/{base,perturbed}_{pres,surf}.nc + nonrad_forcing.nc

路径 B（bespoke，硬编码，CESM2 4×CO₂ 现用）——不复用，仅抽取其机械：
  run_case.py:55-59: if src_type=='cesm2_cmip6': build_cesm2_official.py
                      + inject_cesm_o3.py + mask_subsurface_layers.py（build 后 hook）
```

**结论**：DAMIP 用**路径 A**。但注意路径 A 与路径 B 的一个关键差异（v1 遗漏）：**路径 A 在写盘前跑 `validate_states`（`build_case_input.py:61-77`），拒绝一切 non-finite 值并对气溶胶做量级检查**。路径 B 的"build 后 hook"模式（先写占位、再 inject O₃、再 mask 子地表）在路径 A 走不通——NaN 活不过 validate。因此 DAMIP 的 `build_states()` 必须在**返回前**完成：子地表 HOLD 填充、O₃ 注入/置零、全部单位转换。§6 决策树按此收口。

### 2.2 已有可复用机械（`data/cesm2_cmip6_source.py`，CESM2 专属，需泛化）

| 函数 | 作用 | CESM2 硬编码点（M4 要泛化） |
|---|---|---|
| `list_files(raw_dir, exp_subdir)` | glob `*.nc` 按 `<var>_` 前缀归类 | 文件名模板假定 CESM2 |
| `years_to_month_indices(...)` (L46) | 选年份月段 | **假定 noleap 365 天/年**（`days/365.0`），其它日历错位 |
| `annual_climo_from_monthly(...)` (L64) | 日加权年气候态、fillvalue→NaN | noleap 月天数硬编码 |
| `hybrid_to_plev_mass_conserving(...)` (L196) | hybrid σ-p → plev，柱质量守恒 | 系数名 `a/b/p0`（其它模式多为 `ap/b`） |
| `compute_albedo(rsus, rsds)` (L269) | 反照率、极夜安全 | 通用，直接抽 |
| `load_climo_pres(...)` (L100) | 读全部变量 → 气候态字典 | 变量名/表假定 CESM2 全变量都在；**`cl` 的 %→fraction `/100` 转换在 L132**（泛化时保留此语义） |

`compute_albedo`、`hybrid_to_plev_mass_conserving` 的数值内核模式无关；泛化主要在**日历（改 cftime）**、**hybrid 系数名探测（formula_terms）**、**缺变量容错**、**文件发现**四处。

### 2.3 通用写盘器（`scripts/build_case_input.py`）——可用，但有 3 个 v1 未提的硬点

1. **注册触发是硬编码 import**（L239：`import data.era5_source`）。新增源必须同步在这里加 `import data.cmip6_damip_source`，或者（推荐，一次到位）改成按 `cfg['source']['type']` 的注册表映射动态 import——这属于根目录脚本，不在 core/，改动合规。**忘了这一条，`get_source` 会报 "Unknown source type"，Sonnet 会白排查半小时。**
2. **`write_surf_nc` 只写 `SURF_2D_VARS = ['ts','ps','solar','albedo']`（L25），不含 `huss`**。runner 端 huss 是可选的（`run_parallel_python.py:585-591`：缺失→−999→Fortran 用 HOLD 兜底），所以 RRTMG 主路径没它也能跑。但 DAMIP 模式全都有 huss（实测 13/13），Fu 引擎表面行要用。**正确做法**：让 `write_surf_nc` 支持"可选变量——state 里有才写"；**禁止**把 `'huss'` 无脑加进 `SURF_2D_VARS`——那会让不产 huss 的 ERA5 源写出 `huss=0`，而 0 不会触发 runner 的 |x|>900 缺失判定，等于**悄悄改变现有 ERA5 case 的 Fu 行为**（回归金标之外的隐性破坏）。
3. **`validate_states` 的气溶胶量级检查**（kg/kg ≤1e-5）：DAMIP 源 aerosol=skip → 全零数组，检查自然通过 ✓。

### 2.4 需要泛化的调度点（root 脚本，非 core/fortran，改动合规）

`run_case.py:51-63` build 分发目前白名单限死：`era5_daily / era5_date_range / era5_merra2 / None`。要改成**通用**：任何已注册且非 `cesm2_cmip6` 的源 → 走 `build_case_input.py`。改完后"注册 + yaml"即可接新源——这正是 M5 验收口径。注意保留 `cesm2_cmip6` 的路径 B 分支原样（回归金标护栏）。

### 2.5 输入规范硬约束（`docs/input_spec.md` + 实码核对）

- `*_pres.nc`：`lev` 单位 hPa、序 **surface→TOA**（pyCFRAM 内部 `[::-1]` 翻转，`run_parallel_python.py:579`）；3D 变量 `(time=1, lev, lat, lon)`；变量名 `ta_lay,q,o3,camt,cliq,cice,co2,bc,ocphi,ocpho,sulf,ss,dust`。
- `*_surf.nc`：`ts[K], ps[Pa], solar[W/m²], albedo[0–1]`（+可选 `huss[kg/kg]`）。
- `nonrad_forcing.nc`：runner 读 `[0, 0, :, :]`（time 0、**lev index 0**，`run_parallel_python.py:673`）→ 写盘时 surface 值必须放 lev 轴第 0 位。路径 A 的 `write_nonrad_nc` 用 `lev_out=[1013, …]` 且 surface 在 index 0 ✓；路径 B 用 `[1000,…,1]` surface 也在 index 0 ✓——两种布局 runner 通吃，DAMIP 走路径 A 现成的即可。
- **nonrad 符号统一律（v1 两处表述打架，这里定死）**：runner 计算 `dT_X = −drdt⁻¹·frc_X`；恒等式 `dT_sfcdyn = dT_ocndyn + dT_lhflx + dT_shflx` 要求 **frc = Δ(向下地表能量通量)**。CMIP6 的 `hfls/hfss` 是**向上为正** → `frc_lhflx = −(warm.hfls − base.hfls)`（`build_cesm2_official.py:126-134` 先例）；ERA5 的 `slhf/sshf` 是向下为正 → `frc = +(event − clim)`（`era5_source.py:390` 先例）。两条先例物理上同一条规则，DAMIP 用 CMIP6 变量 → **取负号**。
- CO₂ 作 3D 常数场（mol/mol，`ppmv×1e-6`）；`camt` 0–1（CMIP6 `cl` 是 % → **÷100**）；`o3` 输入规范是 **kg/kg 质量混合比**，CMIP6 `o3` 是 mol/mol 摩尔分数 → **×48/29**（`inject_cesm_o3.py:40` `VMR_TO_MMR` 先例）。
- **case.yaml 必须包含 `input:` 块**（`core/config.py:40-47` 靠它解析路径；缺了 run step 直接 KeyError）。v1 §4.1 示例漏了，v2 §4.1 已补全。

---

## 3. 跨模式异构矩阵（逐维定策略；v2 依实测更新）

| 维度 | 实测变体（§1.4） | 探测/泛化策略 | 落点 |
|---|---|---|---|
| **日历** | noleap/365_day（CESM2/CanESM5/GISS/BCC/FGOALS…）、gregorian/proleptic（IPSL/MRI/CNRM/MIROC/ACCESS）、**360_day（HadGEM3）** | **弃自算天数，改用 cftime**（hqlx210 已装 1.6.2，Mac 1.6.4）：`cftime.num2date(time[:], units, calendar)` 得 `(year,month)`；日加权用该日历真实月长（`date.daysinmonth` 或按 calendar 查表；360_day 恒 30） | `cmip6_common.decode_time()` |
| **垂直坐标（cl/clw/cli）** | ① hybrid `a,b,p0`（CESM2 系，p=a·p0+b·ps）② hybrid `ap,b`（p=ap+b·ps，IPSL/MRI/CNRM/HadGEM 等 CMOR 主流）③ HadGEM 可能为 hybrid height（`lev` standard_name=atmosphere_hybrid_height_coordinate）→ 探测到即 known-issue 或气压近似 | 读 `cl` 变量 `formula_terms` 属性 + 坐标 standard_name 探测；统一成 `p(k)=ap_eff+b·ps`（`ap_eff=a·p0` 或 `ap`）。config 可显式覆盖 | `cmip6_common.detect_vertical()` + models.d yaml `vertical:` |
| **发布压力层（ta/hus/o3）** | 主流 plev19；个别 plev39 | 运行时读实际 `plev`（Pa），log-p 线性插值到 case 统一目标 plev；**目标=源时插值必须恒等**（单测断言） | `cmip6_common.interp_plev_to_target()` |
| **variant_label** | f1 主流；**CNRM=f2**、**HadGEM3=f3**；CanESM5/GISS 有 p1/p2/p3 之分 | models.d yaml 给 default；运行时 glob 兜底（取字典序首个并告警） | `discover_variant()` |
| **grid_label** | gn 主流；**IPSL/CNRM=gr**、GFDL=gr1 | models.d yaml 给 default；glob 兜底 | 同上 |
| **经度/纬度** | lon 0–360 vs −180–180；lat 两向 | 统一 lon 0–360 升序、lat S→N | `cmip6_common.normalize_grid()` |
| **单位** | `cl` %→0–1（÷100）；`o3` mol/mol→kg/kg（×48/29）；ps Pa；plev Pa→hPa | 按 `units` 属性归一，缺属性按 CMIP6 标准假定并告警 | 各 loader |
| **缺变量（主路径！）** | **cl/clw/cli 缺：7/13 模式；rsdt 缺：5/13；o3 缺：10/13**；hfls/hfss/rsus/rsds/huss 实测 13/13 全有 | §6 决策树（v2 新增 solar=ANALYTIC 兜底） | §6 |

> **设计原则**：`cmip6_common.py` 只放模式无关的数值/探测函数；模式专属 quirk 全部进 `configs/damip_models.d/<model>.yaml`（数据而非代码）——加新模式 = 加一个 yaml 文件，不碰 Python（M5 验收口径）。

---

## 4. 架构与文件布局（给子代理的"写哪里"地图）

```
data/
  cmip6_common.py         # [M4] 模式无关 CMIP6 机械：decode_time / detect_vertical /
                          #      hybrid_to_plev（复用现有质量守恒内核）/ interp_plev_to_target /
                          #      normalize_grid / compute_albedo / analytic_solar /
                          #      o3_climatology（抽自 inject_cesm_o3 的插值逻辑）/
                          #      fill_subsurface（抽自 mask_subsurface_layers 的 HOLD 逻辑）/
                          #      discover_files / discover_variant
  cmip6_damip_source.py   # [M4] @register_source('cmip6_damip')；build_states() 见 §5.1
  esgf_fetch.py           # [M4] 纯 stdlib ESGF 客户端：Solr search + file list + HTTP 下载
                          #      （断点续传 + checksum 校验 + 多副本回退），无第三方依赖
  cesm2_cmip6_source.py   # [M4] 重构成薄壳调 cmip6_common（4×CO₂ build 行为不变——回归金标）
configs/
  damip_models.d/         # [M4 起即用] 每模式一个 yaml（v1 放单文件、M5 才拆——v2 改为第一天就
    CESM2.yaml            #   用目录 glob，省掉后期迁移 + 天然免扇出冲突）
    IPSL-CM6A-LR.yaml
    MRI-ESM2-0.yaml
    ...
  damip_experiments.yaml  # [M4] 逐 experiment 语义（哪些强迫时变、CO₂/O₃ 默认）
scripts/
  download_damip.py       # [M4] 下载 CLI（调 data/esgf_fetch.py；--dry-run 报体积）
  make_damip_case.py      # [M5] 从 models.d 脚手架生成 cases/<name>/case.yaml
cases/
  damip_ipsl_histaer/     # [M4] 首个端到端 case（另两个：mri / cesm2）
  damip_<model>_histaer/  # [M5] 逐模式
docs/
  m4_damip_module.md      # [M4] 模块技术文档（架构/IO/异构策略/summary 规范）
  m5_multimodel_userguide.md  # [M5] 用户指南 + known-issues 表 + 自定义模式 howto
  plan_ph3.md             # ← 本文件
tests/
  test_cmip6_common.py    # [M4] 日历×3 / hybrid 质量守恒 / plev 插值恒等 / 归一 / solar 解析
  test_damip_source.py    # [M4] build_states 在合成 fixture 上（含缺变量路径）
  test_damip_regression.py# [M4] 金标：CESM2 4×CO₂ build 重构前后不变
  test_damip_userguide_example.py # [M5] 自定义接口示例 + core/fortran 零改动断言
  data/damip_smoke/       # 合成迷你多模式 NetCDF（~KB 级，不同日历/系数名/缺变量组合）
run_case.py               # [M4] build 分发泛化（§2.4）
scripts/build_case_input.py # [M4] ①注册 import 泛化 ②write_surf_nc 支持可选 huss
```

### 4.1 `case.yaml` 完整模板（v2 补全 input/lapse_rate 块；用户唯一要写的东西）

```yaml
case_name: DAMIP_IPSL_histaer
description: "IPSL-CM6A-LR hist-aer, first decade vs last decade climatology"

source:
  type: cmip6_damip
  model: IPSL-CM6A-LR          # → configs/damip_models.d/IPSL-CM6A-LR.yaml 查 quirk
  experiment: hist-aer         # → configs/damip_experiments.yaml 查语义
  variant: r1i1p1f1            # 可省，省则用 models.d default，再兜底 glob 自动发现
  grid_label: gr               # 可省，同上
  raw_dir: raw_data/cmip6_damip
  base_years: [1850, 1859]
  warm_years: [2011, 2020]     # CESM2/GISS 用 [2005, 2014]（数据止于 2014-12）
  co2:                         # hist-aer：CO₂ 冻结 1850 → 两态相同
    source: constant
    base_ppmv: 284.7
    perturbed_ppmv: 284.7      # 注意键名是 *_ppmv（source_base.get_co2 的实际键，非 *_ppm）
  o3: auto                     # auto|use_model|climatology|skip（语义见 §6.1）
  aerosol:
    source: zero               # source_base.get_aerosol 已支持 'zero'（全零→runner 自动 skip）

grid:
  pressure_levels: [1, 5, 10, 20, 30, 50, 70, 100, 150, 200, 250, 300,
                    400, 500, 600, 700, 850, 925, 1000]   # TOA→sfc（与 cesm2_official 同款
                    # 19 层 CMIP6 标准 plev19；源=目标时插值恒等，误差最小）

input:                          # ★必须有——core/config.load_case 靠它解析路径
  base_pres: input/base_pres.nc
  base_surf: input/base_surf.nc
  perturbed_pres: input/perturbed_pres.nc
  perturbed_surf: input/perturbed_surf.nc
  nonrad_forcing: input/nonrad_forcing.nc

radiation:
  scheme: rrtmg

lapse_rate:                     # 可选：加了就能 --step lr 免费复用 Phase 2 模块
  kernels: [CloudSat, GFDL]
  kernel_months: annual
  sky: all-sky
  tropopause: climkern

run:
  nproc: auto

plot:
  key_region: {lon: [0, 360], lat: [-90, 90]}
```

### 4.2 `configs/damip_models.d/<model>.yaml` 骨架（加模式 = 加一个文件）

```yaml
# configs/damip_models.d/IPSL-CM6A-LR.yaml
institution_id: IPSL
default_variant: r1i1p1f1
default_grid: gr
calendar: gregorian            # 仅作 sanity 对照；实际解码一律信 NC 的 calendar 属性
warm_years_default: [2011, 2020]
vertical: {scheme: auto}       # auto=运行时 formula_terms 探测；可显式 {scheme: hybrid_ap_b, ap: ap, b: b}
```

```yaml
# configs/damip_models.d/CESM2.yaml
institution_id: NCAR
default_variant: r1i1p1f1
default_grid: gn
calendar: noleap
warm_years_default: [2005, 2014]   # hist-aer 数据止于 2014-12（实测）
vertical: {scheme: hybrid_ab_p0, a: a, b: b, p0: p0}
missing_ok: [cl, clw, cli, rsdt, o3]   # 声明性标注：这些缺失是已知且可接受的（→ §6 决策树）
```

> loader 规则：`glob('configs/damip_models.d/*.yaml')` 合并成字典；`case.yaml source.model` 查不到时**不硬失败**——运行时探测兜底 + 告警（这样"接入全新模式"甚至可以不写 models.d 文件，M5 示例的最小路径）。

### 4.3 `configs/damip_experiments.yaml`（逐实验语义）

```yaml
hist-aer:     {varying: [aerosol], co2: fixed_1850, o3_default: auto,      note: "气溶胶不在 Amon → aerosol 恒 skip"}
hist-GHG:     {varying: [ghg],     co2: time_varying, o3_default: auto}
hist-nat:     {varying: [solar, volcanic], co2: fixed_1850, o3_default: auto}
hist-stratO3: {varying: [o3],      co2: fixed_1850, o3_default: use_model}  # 必须用模式 o3
historical:   {varying: [all],     co2: time_varying, o3_default: use_model}
```

---

## 5. base↔perturbed 配对与单强迫一致性检查

### 5.1 `build_states()` 流程（`cmip6_damip_source`；每步标注复用来源）

1. `discover_files` + `discover_variant`：定位 model/experiment/variant/grid 下每变量文件（含多时间分片归并）。
2. `decode_time`（cftime）：逐月 `(year, month)`；取 `base_years` / `warm_years` 月索引。**首末月完整性检查**：所选年段的月数必须 == 年数×12，缺月即 fail-fast 报哪个文件缺哪段。
3. 逐变量 `annual_climo`（该日历真实月长做日加权、fillvalue→NaN）。
4. 云三件套齐全 → `detect_vertical` + `hybrid_to_plev_mass_conserving`（若在 hybrid 层）→ camt=cl/100, cliq=clw, cice=cli；不齐 → 三者全零 + 标记 SKIPPED（§6）。
5. ta/hus/o3 若源 plev ≠ 目标 plev → `interp_plev_to_target`（log-p 线性）。
6. `normalize_grid` 统一经纬；`compute_albedo(rsus, rsds)`；solar：有 rsdt 用 rsdt，无 → `analytic_solar(lat)`（§6.2）。
7. CO₂ 常数场（`get_co2`，注意键名 `base_ppmv/perturbed_ppmv`）；O₃ 按 §6.1（use_model→×48/29；climatology→`o3_climatology()` 注入，两态同场）。
8. **`fill_subsurface`**（抽自 `mask_subsurface_layers.py` 的既有策略）：`lev > ps` 或 ta 为 NaN 的胞 → `ta_lay=ts`、`q/o3` HOLD（复制最低真实层——**RRTMG 不吃 0 值 H₂O/O₃**，`mask_subsurface_layers.py:33-37` 有明文警告）、云/气溶胶=0；顺带 clip 负 cliq/cice。**此步之后 state 内不允许残留任何 NaN**（validate_states 会拒）。
9. nonrad：`lhflx = −(warm.hfls − base.hfls)`，`shflx = −(warm.hfss − base.hfss)`（§2.5 符号统一律）。
10. 组装三字典返回；同时把每个变量的来源/处理决定记进 `self.provenance`（供 §6.3 summary 落盘）。

### 5.2 单强迫一致性 sanity（写进 summary，不作硬门）

- hist-aer：断言 `co2_base == co2_warm`（config 保证）；若用了模式 o3，检查两态域均相对差 < 1%（O₃ 冻结，应近равны）；`|solar_warm − solar_base|` 域均 < 0.5 W/m²（太阳冻结）。
- hist-GHG：反向断言 CO₂ 两态**确实不同**。
- 不满足只告警（提示数据/配置有误），不阻断。

---

## 6. 可选量缺省处理 + `summary.txt`（M4 交付物 2 的核心）

### 6.1 决策树（每变量族；**全部在 `build_states()` 内落地**，无 build 后 hook——理由见 §2.1）

```
云 (cl/clw/cli):  三者齐全且非全 NaN → 用（%→/100，hybrid→plev 质量守恒）
                  否则 → camt=cliq=cice=0，标记 cloud=SKIPPED
                  （BCC 只缺 cl 有 clw/cli → 仍判不齐全 → 跳；不做部分云）
O3:  o3==use_model            → 模式 o3 ×48/29 → kg/kg（hist-stratO3 必须走这条）
     o3==auto:   模式有 o3    → 同 use_model，标记 o3=MODEL
                 模式无 o3    → 注入 CESM 1850 气候态（raw_data/ozone_*_1850clim_c090420.nc，
                                复用 inject_cesm_o3 的 vmr→mmr + log-p/lat 插值，两态同场
                                → frc_o3≡0），标记 o3=CLIMATOLOGY
     o3==climatology          → 强制注入（同上）
     o3==skip 或注入源文件缺失 → o3 = HOLD 一个极小常数剖面? ——否：直接 0 会踩 RRTMG 零值坑。
                                实现：o3=1e-12 kg/kg 常数（辐射上≈0），标记 o3=SKIPPED
solar (rsdt):  有 rsdt → 年均气候态；缺 → analytic_solar(lat)（§6.2），标记 solar=ANALYTIC
               （hist-aer 太阳冻结 → 两态同场 → frc_solar≡0，兜底不污染任何 dT_X）
气溶胶:  hist-aer Amon 无 3D 气溶胶 → source: zero（source_base 现成）→ 全零
         → runner 自动 skip_aerosol（run_parallel_python.py:639 实测阈值 1e-15），标记 aerosol=SKIPPED
通量 (hfls/hfss):  齐全（实测 13/13 都有）→ nonrad 按 §2.5 符号；理论缺失路径 → 返回空 dict
                   → 写盘器跳过 nonrad_forcing.nc → runner 不出 sfcdyn 族，标记 UNAVAILABLE
反照率 (rsus/rsds):  齐全（实测 13/13）→ compute_albedo；理论缺失 → albedo=0.15 常数 + 告警，标记 FALLBACK
huss:  有（实测 13/13）→ 写入 surf nc（写盘器需支持可选变量，§2.3-2）；缺 → 不写（runner HOLD 兜底）
```

### 6.2 `analytic_solar(lat)`（v2 新增；CESM2/NorESM2/ACCESS/GFDL/FGOALS 的主路径）

年均 TOA 入射日均通量只依赖纬度与太阳常数（轨道参数取现代值）：对一年逐日积分 `S0·cos(SZA)` 的日均（标准天文公式：赤纬 δ(day)、时角积分 `H0=arccos(−tanφ·tanδ)`，`Q̄day = S0/π·(H0·sinφ·sinδ + cosφ·cosδ·sinH0)`，年均取 365 日平均；S0=1361 W/m²）。~30 行纯 numpy，单测对照两个理论锚点：全球均值 ≈ S0/4 ≈ 340 W/m²、赤道年均 ≈ 417 W/m²。**两态用同一场**（hist-aer 太阳冻结），对分解的唯一影响是背景 SW 态的微小偏差，且有 rsdt 的模式可作交叉验证（用 IPSL 的 rsdt 对照解析式，纬向曲线应几乎重合——这本身就是一条单测）。

### 6.3 `cases/<case>/output/<case>.summary.txt` 规范（每次 run 生成）

由独立脚本 `scripts/write_run_summary.py` 产出（`run_case.py` 在 run step 末尾追加调用；**不**耦合进 `run_parallel_python.py` 内部——runner 是 Phase 2 已验收代码，少碰）：
- 头：case / model / experiment / variant / grid / base_years / warm_years / 目标 plev / 网格尺寸 / 引擎。
- **过程活跃表**：每个物理过程 ACTIVE / SKIPPED / CLIMATOLOGY / ANALYTIC / FALLBACK / UNAVAILABLE + 一句话原因（数据从 build 阶段的 provenance JSON 读取，`cases/<case>/input/provenance.json`，由 build_states 顺手落盘）。
- **可加性残差**：`dT_observed[sfc] − Σ dT_X[sfc]` 的域均与 max，并引用 `docs/m3_methodology_comparison.md` 既有结论（CFRAM 一阶展开固有非线性，非 bug；合同明文要求引用 M3 文档）。
- **单强迫一致性**（§5.2 结果）。
- 尾：复现命令 + 输入文件 md5。

> summary 是 M4 交付物 2 的**可验收落点**——审阅者一眼看到这次 run 哪些过程被跳过、为什么。跳过是**特性不是 bug**。

---

## 7. ESGF 数据获取（v2 重写：免账号、纯 stdlib、hqlx210 直下）

### 7.1 实测事实（2026-07-05，全部有复现命令）

- **经典 Solr Search API 可用节点**：`esgf-data.dkrz.de`、`esgf.ceda.ac.uk`（含跨节点联邦视图）。**LLNL 已下线**（301 → esgf-1-5-bridge），不要用。hqlx210 直连 DKRZ/CEDA 实测 2s 内 200。
- **CMIP6 文件下载免认证**：对 CEDA THREDDS fileServer 的 Range 请求直接 206 + HDF5 magic。**无需等甲方 ESGF 账号**（合同写"甲方提供 ESGF 账号或机构访问"是保险条款，实测用不上；若某镜像今后要求登录，再启用账号即可）。
- 文件级时间过滤：CMIP6 文件名自带 `_YYYYMM-YYYYMM.nc` 后缀 → 从文件名解析时间段，只下与 base/warm 年段相交的分片。注意 CanESM5/IPSL 等把 1850–2020 打成单文件（IPSL ta 实测 1.65 GB/文件）——单文件的没法省，照单全收。
- **体积精算**：全变量模式（IPSL 档）≈ 6 个 3D 变量×~1.6GB + 8 个 2D×~0.2GB ≈ **11 GB/模式**；缺云模式 ≈ 5–6 GB。M5 八模式合计 **≈ 65–75 GB**，压合同 80 GB 预算上限——`--dry-run` 体积门必须先跑（超了就对单文件大户换分片变体或用 [2005,2014] 减一个十年分片）。
- **落盘地点**：hqlx210 `raw_data/cmip6_damip/<model>/<experiment>/`（`/home/lzhenn/work` 挂在 hqlx74:/export/r074，实测剩 1.5T ✓）。本地 Mac 只留 tests/data 迷你 fixture。下载在 hqlx210 上跑（数据在哪算在哪；ERA5 先例同此）。

### 7.2 `data/esgf_fetch.py` 设计（纯 stdlib，无第三方依赖）

```python
# 核心三函数（探测脚本已在本次尽调实测跑通，直接演化成模块即可）
def search_datasets(model, experiment, variables, variant=None, grid=None,
                    node="https://esgf-data.dkrz.de/esg-search/search") -> list[dict]
    # type=Dataset, facets 见附录 B；返回含 variant/grid/version 的数据集清单
def list_files(dataset_id, node=...) -> list[dict]
    # type=File, fields=title,size,url,checksum,checksum_type；解析 HTTPServer url
def fetch(url, dest, checksum=None, resume=True)
    # urllib + Range 断点续传；sha256 校验；单文件多副本 url 逐个回退
```

- 年段过滤：`filename_time_overlap(title, base_years, warm_years)`（解析 `_YYYYMM-YYYYMM.nc`）。
- 节点回退顺序：DKRZ → CEDA →（文件级 url 列表中的任意 data node）。
- `manifest.json` 落每文件 {url, size, checksum, mtime}，重跑幂等。
- **intake-esgf / esgpull 均降为"可选替代"**：不装、不依赖；README 提一句给未来用户。

### 7.3 变量清单（每模式每实验）

`ta hus ts ps cl clw cli rsdt rsds rsus hfls hfss huss` + 条件 `o3`（o3=auto 且模式有才下）。全部 `table_id=Amon`。

### 7.4 下载 CLI

```bash
# 全部在 hqlx210 上执行（嵌套 ssh，见 persistent_context.md）
python3 scripts/download_damip.py --model IPSL-CM6A-LR --experiment hist-aer --dry-run   # 只报文件清单+体积
python3 scripts/download_damip.py --model IPSL-CM6A-LR --experiment hist-aer             # 实下
python3 scripts/download_damip.py --all-m5 --dry-run                                     # 八模式总体积（80GB 门）
```

---

## 8. 工作分解（WP）· 子代理指派 · 检查点 · 回退

> 每 WP 标注 **Owner 档位 / 产出 / 前置 / Done 判据 / 回退**。依赖批次见 §9.2。

### 批次 0（打底，阻塞门）

**WP-0 · 分支 + 计划落盘 + PR 叠链决策 〔lead〕**
- 产出：从 `feat/m2-m3-lapse-rate-kernel` 切 `feat/m4-m5-damip`；本文件已落 `docs/plan_ph3.md`。
- **PR 叠链（v2 新增）**：Phase 2 的 PR 还没开。合同按检查点 PR 验收 → 正确顺序是**先开 M2/M3 的 PR**（`feat/m2-m3-lapse-rate-kernel` → main），M4 分支叠其上；M4 PR 用 base=`feat/m2-m3-lapse-rate-kernel`（GitHub stacked PR），Ph2 PR 合并后 retarget main。此事需用户点头（涉及对甲方的交付节奏），lead 在批次 0 向用户确认一次。
- Done：`git branch --show-current` = `feat/m4-m5-damip`；PR 策略有用户结论。
- 回退：无。

**WP-M4.0a · ESGF 探测与 dry-run 〔haiku；可立即开始，无任何前置〕**
- 产出：`data/esgf_fetch.py` + `scripts/download_damip.py`（含 --dry-run）。
- Done：对 M4 三模式 dry-run 各报出文件清单+体积；三模式合计 <35 GB；对 §1.4 矩阵抽查 2 个模式的变量 facet 与本计划一致。
- 回退：DKRZ 挂 → CEDA；两个都挂 → 文件级 url 直连任一 data node。

**WP-M4.0b · 首批 3 模式数据落盘 hqlx210 〔lead 触发，后台跑〕**
- 前置：M4.0a。**不等甲方账号**（免认证已实证）。
- Done（门②）：IPSL/MRI/CESM2 hist-aer 全变量文件落 `raw_data/cmip6_damip/`，逐个 netCDF4 打开成功；`ta` 有 `plev`、IPSL/MRI 的 `cl` 有 formula_terms；manifest checksum 全过。
- 回退：§7.2 多副本回退；个别文件死链 → 换 variant（IPSL 有 r1..r10 十个）。
- **可与 WP-M4.1 完全并行**（互不依赖）。

### 批次 A（模式无关内核 + 回归金标）

**WP-M4.1 · 抽 `cmip6_common.py` + CESM2 回归金标 〔sonnet〕**
- 产出：`data/cmip6_common.py`（§4 函数集：新写 `decode_time`/`detect_vertical`/`interp_plev_to_target`/`normalize_grid`/`analytic_solar`/`discover_*`；**平移复用** `hybrid_to_plev_mass_conserving`/`compute_albedo` 数值内核、`inject_cesm_o3` 的 `o3_climatology` 插值、`mask_subsurface_layers` 的 `fill_subsurface` HOLD 策略）；`cesm2_cmip6_source.py` 重构成薄壳。
- 前置：无（纯本地重构）。
- Done（门①，最高优先）：`tests/test_damip_regression.py` 断言重构后 CESM2 4×CO₂ build 产出不变。**金标生成的现实约束（v1 未写）**：cesm2 raw 数据（6.9+4.1+21 GB）只在 hqlx210——两层金标：
  - **本地层**：从真实 CESM2 raw 裁一个 4×4×全层迷你子集进 `tests/data/damip_smoke/cesm2_mini/`（lead 在 hqlx210 上裁好拉回，~几 MB），金标 = 迷你 build 输出的 md5/数值，重构前先跑一次锁定，本地 pytest 可复跑。
  - **远端层**：重构合入后在 hqlx210 全量 `run_case.py cesm2_4xco2_official --step build` 一次，input NC 与重构前 `md5` 对比（重构前先备份 `.md5sum` 清单）。
- 回退：金标红 → 回滚，`cmip6_common` 与 legacy 并存不动（DRY 让步于安全）。

### 批次 B（DAMIP 插件，串行于 A）

**WP-M4.2 · DAMIP 源插件 + 配置 〔sonnet〕**
- 产出：`data/cmip6_damip_source.py`（§5.1 十步流程）；`configs/damip_models.d/{IPSL-CM6A-LR,MRI-ESM2-0,CESM2}.yaml` + glob loader；`configs/damip_experiments.yaml`；provenance.json 落盘。
- 前置：WP-M4.1（用 common）+ WP-M4.0b（真数据在 hqlx210；本地开发用 smoke fixture）。
- Done：在 IPSL hist-aer 真数据上 `build_states()` 返回三字典，断言：变量名/维度=input_spec、lev 序 sfc→TOA、无 NaN、camt∈[0,1]、o3 量级 1e-8~1e-5 kg/kg（防 vmr/mmr 转换漏做——**huss-as-O3 事故（session_log 2026-05-12）的教训：O₃ 错 5 个量级会污染整套分解**）、nonrad 符号（全球均 Δ(−hfls) 与文献符号一致）。
- 回退：hybrid 探测失败 → models.d yaml 显式给系数名；plev 异常 → 打印实际 plev 人工确认。

**WP-M4.3 · 分发泛化 + 写盘器两处小改 〔sonnet〕**
- 产出：`run_case.py` build 分支通用化（§2.4）；`build_case_input.py` ①注册 import 泛化（§2.3-1）②`write_surf_nc` 可选 huss（§2.3-2，**只写 state 里有的**）。
- 前置：WP-M4.2。
- Done（门③）：`python3 run_case.py damip_ipsl_histaer --step build` 产出 4+1 个合规 NC + provenance.json；**回归护栏**：`eh13 --step build --dry-run` 与 `cesm2_4xco2_official` 分发路径均不受影响（跑一遍确认）。
- 回退：分发泛化波及 era5/cesm2 → 保留显式分支 + 新增通用兜底。

**WP-M4.4 · summary.txt 〔sonnet〕**
- 产出：`scripts/write_run_summary.py`（§6.3 规范，读 provenance.json + cfram_result.nc）；`run_case.py` run step 末尾追加调用。
- 前置：WP-M4.2（provenance 格式）。
- Done：IPSL run 后 summary 生成：aerosol=SKIPPED、o3=MODEL、cloud=ACTIVE、solar=ACTIVE；CESM2 run 后：cloud=SKIPPED、solar=ANALYTIC、o3=CLIMATOLOGY；含可加性残差 + M3 文档引用。
- 回退：无风险（独立后处理脚本）。

### 批次 C（端到端 + 测试，M4 收口）

**WP-M4.5 · 三模式 × hist-aer 端到端 〔lead 编排 + sonnet 执行〕**
- 产出：`cases/damip_{ipsl,mri,cesm2}_histaer/`；hqlx210 各跑出 `cfram_result.nc` + summary。
- 前置：M4.1–M4.4 + M4.0b。
- Done（=M4 端到端门）：3 case 均出 `cfram_result.nc`（无 NaN 泄漏、`dT_sfcdyn=dT_ocndyn+dT_lhflx+dT_shflx` 残差 <1e-10 K）+ summary 正确标记；**物理 sanity**：hist-aer 为净冷却强迫 → 全球均 `dT_observed[sfc]` ∈ (−3, 0) K，且 NH 冷于 SH（人为气溶胶集中北半球）——三模式方向一致即过，不设定量阈值（合同数值阈值表仅适用"对比参考结果"型验收，Phase 3 无外部参考）。
- 回退：某模式跑挂 → 数据/坐标问题记 known-issue 换备选池（§1.4）；代码问题修 common/source。
- **运行纪律**：rsync Mac→mini→hqlx210 排除 `fortran/cfram_*_1col*`，部署后远端 `make` 重编 + verify size/mtime（persistent_context.md 踩坑记录）。

**WP-M4.6 · 测试 + 覆盖率 〔haiku 起草，sonnet 补〕**
- 产出：`tests/test_cmip6_common.py`（noleap/gregorian/360_day 三日历合成 time 轴、hybrid 柱质量守恒、plev 插值恒等、normalize 幂等、analytic_solar 两锚点+对照 IPSL rsdt）、`test_damip_source.py`（smoke fixture 上全决策树分支：全变量/缺云/缺rsdt/缺o3）、`tests/data/damip_smoke/`（合成迷你 NC，不同日历×系数名×缺变量组合，~KB 级）。
- Done：`pytest tests/ -q` 全绿（**含 Phase 2 既有 20 个测试不回归**）；`pytest --cov=data.cmip6_common --cov=data.cmip6_damip_source --cov=data.esgf_fetch` 行覆盖 ≥60%（合同硬指标；esgf_fetch 的网络函数用 monkeypatch 假响应测解析逻辑）。
- 回退：合成 fixture 太繁 → 从真实数据裁 4×4 迷你 NC（cesm2_mini 先例）。

**WP-M4.7 · M4 PR 〔lead〕（v2 新增）**
- 产出：M4 全部 commit → PR（base 按 WP-0 决策）；描述含验收对照表（3 模式端到端实测、skip 逻辑演示、覆盖率数字、`docs/m4_damip_module.md`）。
- Done：PR 开出、CI/pytest 绿；甲方 7 个工作日 Code Review 时钟启动。
- **M5 不等 M4 验收回复**（合同只约束甲方回复时限，未禁止乙方继续）——批次 D 立即开工。

**───────── M4 验收门（≥3 模式端到端 + 可选跳过 + summary + PR）─────────**

### 批次 D（M5 多模式扩展）

**WP-M5.1 · 扩到 8 模式 〔下载 lead 串行；case 配置+运行 sonnet/haiku 扇出〕**
- 产出：`configs/damip_models.d/` 补齐 §1.4 八模式；各 `cases/damip_<model>_histaer/`；各自跑通。
- **执行序（v2 修订）**：①lead 先跑 `download_damip.py --all-m5 --dry-run` 过 80GB 门 → 在 hqlx210 **串行**批量下载其余 5 模式（单机带宽，扇出下载无意义且互相挤占）；②数据就位后按模式扇出 agent（每 agent：写 models.d yaml + case.yaml → build → run → 验证 → 回报），每 agent 只碰自己模式的文件，天然无冲突（models.d 目录制从 M4 就生效，无 v1 的"M5 才拆文件"迁移）。
- 前置：M4 门过。
- Done（=M5 端到端门）：8 模式尝试、**≥6 出 `cfram_result.nc`**；失败模式定位到具体报错与根因。
- 回退：某模式挂 → 备选池替换或记 known-issue（§1.4 已预判：HadGEM3 若 cl 是 hybrid-height 坐标则可能是第一个真 known-issue；CanESM5/CESM2 云 skip 是设计内行为不算失败）。

**WP-M5.2 · known-issues 文档化 〔haiku〕**
- 产出：`docs/m5_multimodel_userguide.md` known-issues 表 + README 挂链；逐条含**实测报错原文**（合同要求"中文档化具体 known-issue：数据缺值、坐标异常、单位不一致等"）。
- Done：每个未通过/降级模式在 README 有一行可核对记录（云 skip 的模式也列出，注明"设计内缺省处理"）。

**WP-M5.3 · 用户自定义接口 + 可运行示例 〔sonnet〕**
- 产出：`scripts/make_damip_case.py` 脚手架；docs "如何接入新模式/新 experiment" 章节；**一个不在八模式清单内的全新模式**（从备选池选，如 ACCESS-ESM1-5 或 NorESM2-LM——故意选缺云的，同时演示 skip 与自定义接入）仅通过 models.d yaml + case.yaml 接入并 build 通过。
- Done（=M5 可扩展门）：`tests/test_damip_userguide_example.py` 断言接入全程 `git diff --name-only` 中 `core/` 与 `fortran/` 零出现；示例 case build 出 4+1 NC。
- 回退：新模式需要新 quirk 类型（如 hybrid-height）→ 允许在 `cmip6_common` 加**通用**分支，文档说明这是"框架增强"而非"每模式改代码"。

**WP-M5.4 · 全套文档 + 覆盖率 〔sonnet 正文，haiku 排版〕**
- 产出：`docs/m4_damip_module.md` + `docs/m5_multimodel_userguide.md` 定稿；README Phase 3 段；`docs/technical_notes_{en,zh}.md` 更新；核对文档覆盖 **Phase 2 + Phase 3 全部新增功能**（合同 M5 交付物 3 原文）。
- Done：覆盖率 ≥60% 复核（加了新代码后重跑）；文档命令逐条可复现。

**WP-M5.5 · M5 PR 〔lead〕**
- 产出：M5 增量 → PR；描述含 M5 验收表（八模式清单、≥6 通过实测记录、known-issues、覆盖率、自定义接口示例）。
- Done：PR 开出、pytest 绿。合同 10 天免费缺陷修复期自 M5 终验起算——PR 后保留分支现场。

---

## 9. 子代理编排 playbook

### 9.1 档位分工
- **haiku**：esgf_fetch/download CLI 起草、逐模式 yaml、boilerplate 测试、known-issue 表、文档排版。
- **sonnet**：cmip6_common 重构、DAMIP 源插件、分发泛化、决策树落地、单模式端到端调试、文档正文、根因定位。
- **lead**：架构与回归金标把关、PR 叠链与用户沟通、hqlx210 下载/部署/批量运行、模式清单微调拍板、两个 PR。

### 9.2 依赖图与并行批次

```
WP-0（分支+PR叠链决策）
  ├─ WP-M4.0a（esgf_fetch+dry-run，haiku）──► WP-M4.0b（3模式下载@hqlx210，lead 后台）──┐
  └─ WP-M4.1（cmip6_common+金标，sonnet）★门①                                      │(数据)
                        ▼                                                          │
              WP-M4.2（DAMIP 源，sonnet）◄─────────────────────────────────────────┘
                        ▼
        ┌─ WP-M4.3（分发泛化+写盘器，sonnet）★门③
        └─ WP-M4.4（summary，sonnet）
                        ▼
              WP-M4.5（3 模式端到端@hqlx210，lead+sonnet）★M4 门
              WP-M4.6（测试，haiku→sonnet）贯穿批次 B/C
              WP-M4.7（M4 PR，lead）
──────────── M4 验收（不阻塞 D 开工）────────────
      WP-M5.1（8 模式：下载 lead 串行 → case 扇出）★M5 端到端门
              WP-M5.2（known-issues，haiku）
              WP-M5.3（自定义接口+示例，sonnet）★M5 可扩展门
              WP-M5.4（文档+覆盖率，sonnet/haiku）
              WP-M5.5（M5 PR，lead）
```

### 9.3 文件所有权边界（防并行冲突）
- M4.1 独占 `data/cmip6_common.py` + `data/cesm2_cmip6_source.py`；M4.2 独占 `data/cmip6_damip_source.py` + `configs/`；M4.3 独占 `run_case.py` + `scripts/build_case_input.py`；M4.4 独占 `scripts/write_run_summary.py`；M4.0a 独占 `data/esgf_fetch.py` + `scripts/download_damip.py`。
- M5.1 扇出：每模式 agent 只写 `configs/damip_models.d/<model>.yaml` + `cases/damip_<model>_histaer/`，零交集。
- 并行 agent 用 `Agent(isolation="worktree")`；`run_case.py` 等共享文件全程 lead/单 agent 串行改。

### 9.4 派活 prompt 模板（lead 逐 WP 填空）
```
subagent_type: general-purpose   model: <sonnet|haiku>   isolation: worktree
prompt: |
  你在 pyCFRAM（/Users/zhenningli/work/ust-jumper/pyCFRAM，分支 feat/m4-m5-damip）执行 Phase 3 的 <WP 编号>。
  先读 docs/plan_ph3.md 的 §<节> + §2（现状与 file:line 硬点）+ §6（决策树）。
  只允许改：<文件白名单>。不碰其它文件（尤其 core/ 与 fortran/——M5 验收断言这两目录零改动）。
  实现要点：<抄该 WP 的产出/Done/回退>。
  完成后跑 <该 WP Done 命令>，把实测结果（build 出几个 NC / 覆盖率数字 / e2e 是否出 cfram_result.nc + summary 关键行）原文贴回。
  远程纪律：hqlx210 事务走 persistent_context.md 的嵌套 ssh + rsync 链（排除 fortran/cfram_*_1col*，部署后 verify size+mtime）。
  不 git push、不开 PR（lead 统一做）。
```

### 9.5 集成协议
- 每 WP 一原子 commit，末尾 `Co-Authored-By: Claude ...`。
- lead 收编 worktree → 本地 `pytest` 全绿（**含 Phase 2 的 20 个既有测试**）+ 金标过 → rsync hqlx210 验证 → 推进下一批。
- 远程纪律同 §8 WP-M4.5。

---

## 10. 测试策略

- **单元**（pytest，本地跑，不碰网络/大数据）：
  ① `decode_time`：noleap/gregorian/360_day 各造合成 time 轴（`cftime.date2num` 正反构造），断言 (year,month) 与月长；
  ② `hybrid_to_plev` 柱质量守恒（沿用现有实现性质，`a,b,p0` 与 `ap,b` 两种系数路径都测）；
  ③ `interp_plev_to_target` 目标=源时恒等；
  ④ `normalize_grid` 经度 wrap/纬度翻转幂等；
  ⑤ `analytic_solar` 两锚点（340/417 W/m²）；
  ⑥ 决策树四分支（全变量/缺云/缺rsdt/缺o3）在 smoke fixture 上产出正确标记与合规 state；
  ⑦ `esgf_fetch` 解析逻辑（monkeypatch 假 Solr JSON + 假 HTTP 响应，不真联网）。
- **回归金标**（最高优先）：`test_damip_regression.py` 双层（本地 cesm2_mini fixture + hqlx210 全量 md5 清单，见 WP-M4.1）。
- **E2E smoke**：迷你 DAMIP case（smoke fixture 8×8）`--step build` → 本地 1×1 或 8×8 `run`（Mac gfortran toolchain 已支持，persistent_context "gnu toolchain" 节）→ 出 cfram_result.nc + summary。
- **覆盖率**：`pytest --cov=data.cmip6_common --cov=data.cmip6_damip_source --cov=data.esgf_fetch --cov-report=term-missing` ≥ 60%（合同硬指标）。
- **M5 断言**：`test_damip_userguide_example.py` 用 `git diff --name-only <base>` 断言 `core/`、`fortran/` 零改动。
- **既有测试不回归**：Phase 2 的 `test_kernels/test_lr_*/test_e2e_smoke`（20 passed + 1 skipped 基线）每批次全绿。

---

## 11. 交付物 ↔ 合同验收映射（对照 contract_10.pdf 原文）

### 11.1 M4 清单（合同：DAMIP-aware 预处理模块 + 可选缺省逻辑；第 13–15 周，20%）
- [ ] **回归金标 PASS**（CESM2 4×CO₂ build 不变；双层金标）
- [ ] `data/cmip6_common.py` + `data/cmip6_damip_source.py`（识别 experiment_id、自动 base↔single-forcing 配对、**动态检测各模式可用压力层并 Python 端插值到统一网格**——合同原文逐项落地）
- [ ] `data/esgf_fetch.py` + `scripts/download_damip.py`（从零下载链路，免认证实测）
- [ ] `run_case.py` + `build_case_input.py` 分发泛化
- [ ] **≥3 模式 × hist-aer 端到端**（IPSL/MRI/CESM2，全部在合同建议清单内；`cfram_result.nc` + summary）
- [ ] 气溶胶/O₃（+云/solar，超合同范围的同类缺省）自动跳过 + `*.summary.txt` 含 additivity 残差说明、引用 M3 文档结论（合同原文）
- [ ] `tests/`（新模块覆盖 ≥60%）+ smoke 数据
- [ ] `docs/m4_damip_module.md`
- [ ] **M4 PR**（甲方 7 工作日 Code Review 时钟）

### 11.2 M5 清单（第 16–19 周，17%）
- [ ] **≥8 模式尝试、≥6 端到端通过**；未通过模式 README **中文** known-issue（数据缺值、坐标异常、单位不一致等，含实测报错）
- [ ] 用户自定义接口：配置文档 + ≥1 可运行示例；**仅 case.yaml + 可选 data/<source>_source.py，core/ 与 fortran/ 零修改**（测试断言）
- [ ] 文档覆盖 Phase 2/3 全部新增功能；新模块 pytest 行覆盖率 ≥60%
- [ ] **M5 PR + Code Review**；此后 10 天免费缺陷修复期

---

## 12. 风险登记 + 回退矩阵（v2 重排）

| # | 风险 | 触发信号 | 等级 | 回退 |
|---|---|---|---|---|
| R1 | **重构破坏 CESM2 现有 build** | 金标红 | 高 | 门①拦截；回滚重构，common 与 legacy 并存 |
| R2 | ~~ESGF 凭据故障~~ → 降级为**节点可用性波动** | 下载 4xx/5xx/超时 | 低（免认证+双索引+多副本实测） | DKRZ↔CEDA 互备；file url 逐副本回退；真要凭据再启用甲方账号 |
| R3 | 模式缺云/o3/rsdt | 变量检索为空 | **已消化为主路径**（§1.4 预判 + §6 决策树） | 不是风险是特性；summary 标记 + known-issue |
| R4 | 异常日历/垂直坐标（360_day/hybrid-height） | cftime/formula_terms 解析异常 | 中 | cftime 通吃日历（单测锁）；HadGEM3 hybrid-height 若真不支持 → known-issue（计入允许失败 2/8） |
| R5 | 80GB 预算超支 | --all-m5 dry-run 报 >75GB | 中 | 单文件大户换 [2005,2014] 窗口少下分片；o3 只给有 o3 且 auto 命中的模式下 |
| R6 | 写盘器/分发改动波及 ERA5/cesm2 既有 case | eh13 dry-run 或 cesm2 分发行为变化 | 中 | WP-M4.3 Done 内建回归护栏；huss 只写"state 有" |
| R7 | 单强迫语义配错（hist-GHG 固定 CO₂ 之类） | §5.2 sanity 告警 | 低 | experiments.yaml 单一事实源 + sanity |
| R8 | O₃ 单位转换漏做（vmr↔mmr 5 个量级） | build 断言 o3 量级越界 | 中 | WP-M4.2 Done 显式断言（huss-as-O3 事故教训） |
| R9 | nonrad 符号搞反 | 恒等式残差正常但 dT_lhflx 符号与文献反 | 中 | §2.5 统一律 + M4.5 物理 sanity（NH 冷于 SH） |
| R10 | rsync 覆盖 Linux 二进制 | hqlx210 `Exec format error` | 中 | 排除 `fortran/cfram_*_1col*`、远端 make 重编（persistent_context 踩坑） |
| R11 | PR 叠链混乱（Ph2 未合并先动 Ph3） | review 范围含 Ph2 diff | 低 | WP-0 与用户定死 stacked PR 策略 |

**诊断阶梯**（某模式端到端结果可疑时按序查）：① summary 的过程标记是否符合 §1.4 预期 → ② provenance 单位转换记录（cl/100、o3×48/29）→ ③ 子地表填充统计（高地形区比例是否合理）→ ④ nonrad 符号（恒等式+半球对比）→ ⑤ hybrid 系数探测结果 vs models.d 声明 → ⑥ 换 IPSL 同流程对照（全变量基准模式）。

---

## 13. 时间线（对齐合同第 13–19 周）

| 周 | 里程碑 | 批次 | 关键子检查点 |
|---|---|---|---|
| 13 | M4 启动 | 0+A | 分支+PR 叠链决策；esgf_fetch dry-run（门②前半）；3 模式下载后台跑；cmip6_common 重构 + **金标（门①）** |
| 14 | M4 核心 | B | build_states 三字典（含 O₃ 量级断言）；分发泛化 + 写盘器；`damip_ipsl_histaer --step build` 出 4+1 NC（门③） |
| 15 | M4 收口 | C | summary；**3 模式端到端（M4 门）**；测试 ≥60%；`m4_damip_module.md`；**M4 PR** |
| 16 | M5 扩展 | D | --all-m5 dry-run 过 80GB 门；串行下载其余 5 模式；先扇出数据齐的模式 |
| 17 | M5 扩展 | D | **≥6/8 通过**；失败根因定位 |
| 18 | M5 接口 | D | 自定义接口 + 新模式仅 yaml 接入示例（可扩展门）；文档正文 |
| 19 | M5 收口 | D | 文档覆盖 Ph2/3 复核；覆盖率复跑；**M5 PR** |

---

## 附录 A · DAMIP experiment_id 速查
见 §1.1。Tier-1（hist-GHG/hist-aer/hist-nat）发布最广；hist-stratO3 必须 `o3: use_model`。

## 附录 B · ESGF Solr API 速查（实测可用形态）
```
GET https://esgf-data.dkrz.de/esg-search/search?
    project=CMIP6&activity_id=DAMIP&experiment_id=hist-aer&source_id=<model>
    &table_id=Amon&variable_id=<var>&type=Dataset|File
    &format=application/solr+json&limit=<n>&facets=variable_id,variant_label,grid_label
```
- 变量覆盖矩阵：`limit=0&facets=variable_id` 读 facet 计数（本计划 §1.4 即此法产出）。
- 文件清单：`type=File&fields=title,size,url,checksum,checksum_type`；url 取 `HTTPServer` 项 `|` 前段。
- CEDA 同构：`https://esgf.ceda.ac.uk/esg-search/search`。**LLNL 已死，勿用。**

## 附录 C · 常用命令
```bash
# 分支
git checkout feat/m2-m3-lapse-rate-kernel && git checkout -b feat/m4-m5-damip
# 下载（hqlx210 上；免认证）
ssh mini 'ssh lzhenn@hqlx210 "cd /home/lzhenn/work/ust-jumper/pyCFRAM && \
  nohup python3 -u scripts/download_damip.py --model IPSL-CM6A-LR --experiment hist-aer \
  > /tmp/damip_dl.log 2>&1 &"'
# 金标（重构前生成，重构后复跑）
pytest tests/test_damip_regression.py -q
# 端到端（hqlx210）
python3 run_case.py damip_ipsl_histaer --step build     # 本地或 hqlx210 均可（数据在哪跑哪）
ssh mini 'ssh lzhenn@hqlx210 "cd /home/lzhenn/work/ust-jumper/pyCFRAM && \
  source /home/lzhenn/.bashrc_liquor_i22wrf415 >/dev/null 2>&1 && \
  nohup python3 -u run_case.py damip_ipsl_histaer --step run --nproc 200 > /tmp/damip.log 2>&1 &"'
# 覆盖率
pytest tests/ -q --cov=data.cmip6_common --cov=data.cmip6_damip_source --cov=data.esgf_fetch --cov-report=term-missing
```

## 附录 D · 关键决策记录（v2 增补）
| 决策 | 取值 | 理由 |
|---|---|---|
| 主 experiment | hist-aer | 用户确认；发布最广、信号强 |
| base↔warm 配对 | 同实验首十年 vs 末十年气候态（warm 窗 per-model） | 自洽；CESM2/GISS 数据止于 2014 的实测约束 |
| M4 三模式 | **IPSL-CM6A-LR / MRI-ESM2-0 / CESM2** | 实测变量矩阵（§1.4）：两个全变量基准 + 一个 skip 全家桶演示；均在合同建议清单 |
| M5 八模式 | IPSL/MRI/CNRM/MIROC6/GISS/HadGEM3/CanESM5/CESM2 | 6 全变量打底 + 2 skip 演示；f1/f2/f3 + 五种日历 + gn/gr 全覆盖 |
| 数据源形态 | 注册式 DataSource 插件 | 走通用写盘器；满足 M5 口径 |
| ESGF 工具 | **纯 stdlib urllib（自写 esgf_fetch）** | 免认证实测；不给共享 NFS conda 装包；intake-esgf/esgpull 降为可选 |
| 下载地点 | **hqlx210 直下**（DKRZ/CEDA 实测 2s 可达） | 数据在哪算在哪；Mac 只留迷你 fixture |
| 缺省处理落点 | **全部在 build_states() 内**（无 build 后 hook） | validate_states 写盘前拒 NaN（§2.1） |
| solar 兜底 | analytic_solar(lat)（hist-aer 两态同场 frc_solar≡0） | 5/13 模式缺 rsdt 的实测现实 |
| O₃ 兜底 | CESM 1850 气候态注入（vmr×48/29） | inject_cesm_o3 先例；两态同场 frc_o3≡0 |
| nonrad 符号 | frc = Δ(向下地表通量) = −Δhfls（CMIP6） | build_cesm2_official/era5_source 两先例统一律 |
| 模式 quirk | configs/damip_models.d/<m>.yaml 目录 glob（M4 首日即用） | 加模式=加文件；免 M5 迁移 + 免扇出冲突 |
| 日历 | cftime 通吃（NC calendar 属性为准） | 360_day/noleap/gregorian 实测并存 |
| PR 节奏 | M4、M5 各一个 PR；先解 Ph2 PR 叠链 | 合同按检查点 PR 验收（7 工作日时钟） |

---

> **给执行者的最后一句**：进场先钉三道门——①CESM2 金标（别把 4×CO₂ build 跑坏；金标先于重构生成）、②首个 hist-aer 子集在 hqlx210 落盘可开（**今天就能下，不等任何账号**）、③`damip_ipsl_histaer --step build` 出合规 4+1 NC。三个最容易翻车的坑按出现顺序：`build_case_input.py:239` 的注册 import、O₃ 的 mol/mol→kg/kg（差 5 个量级，huss-as-O3 事故的近亲）、nonrad 的 −Δhfls 符号。缺云/缺 rsdt/缺 o3 不是异常是**主路径**（13 个候选里只有 3 个全变量），决策树 + summary 把"跳过"讲成可审计的特性。`core/` 与 `fortran/` 全程零改动——这既是奥卡姆最省，也正好是 M5 的验收口径。
