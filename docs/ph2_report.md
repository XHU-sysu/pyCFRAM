# Phase 2（M2/M3）完成汇报

**日期**：2026-07-02
**对应合同**：`contract/contract_10.pdf` Phase 2 — M2（辐射核对比，第3-8周，28%）+ M3（物理过程分解与对比文档，第9-12周，25%）
**执行依据**：`docs/plan.md`
**分支/提交**：`feat/m2-m3-lapse-rate-kernel`，基于 `9401a7e` + 后续文档修订

---

## 一、背景与目标

M2 要求给 pyCFRAM 增加一个"辐射核 Lapse-Rate 分解"能力：从 `cfram_result.nc` 的逐层温度响应中，抽出"垂直非均匀分量"（偏离地表均匀增温的部分），用 Kramer 核（ClimKern 里的 `CloudSat`）为主、GFDL 核为参考，算出对应的 TOA 辐射扰动 ΔR_LR，并与官方 ClimKern 库的 `calc_T_feedbacks` 做交叉验证。合同硬指标：全场空间相关性 ≥0.85，域均值相对误差 ≤15%。

M3 要求在此基础上，把"垂直非均匀增温是被哪些物理过程造出来的"归因到 CFRAM 已有的分解项（`dT_q`、`dT_atmdyn`、`dT_sfcdyn` 等），并产出方法论对比技术文档（≥4 张图）。

---

## 二、进场发现的阻塞项及处理

正式开工前的尽调发现两个必须先处理的问题：

1. **数据陈旧**：`cases/cesm2_4xco2_official/output/cfram_result.nc` 无论本地还是 hqlx210 规范路径，都是 2026-05-11 sfcdyn 修复（commit `fedcf4d`）之前的旧输出，恒等式检验 `dT_ocndyn+dT_lhflx+dT_shflx == dT_sfcdyn` 残差实测高达 190K（应为机器精度）。`dT_observed` 本身不受影响（M2 可用），但 M3 依赖的 `dT_atmdyn/dT_sfcdyn/dT_lhflx/dT_shflx` 全部是坏的。
   → 在 hqlx210 上重跑该 case（200 进程，394.7 秒，55296 格点），重跑后残差降到 **3.55×10⁻¹⁴ K**（机器精度），新鲜度门通过。旧文件已改名 `cfram_result.nc.bak_prefix_20260511` 留档。

2. **环境依赖**：`climkern` 库依赖 `esmpy`/`xesmf`，此前一直装不上。经 dry-run 确认 conda-forge 在 Apple Silicon 上可以正常解出这两个包，于是新建独立环境 `pycfram-kern`（不污染 conda base），成功装齐 esmpy 8.9.1 + xesmf 0.9.2 + climkern 1.2.1，`import climkern` 完全正常。

### 意外踩坑：rsync 不认 `.gitignore`

重跑前的代码同步（Mac → mini → hqlx210）过程中，Mac 本地用 gfortran 编译的 ARM64 版 `cfram_rrtmg_1col`（虽然在 `.gitignore` 里，但 rsync 是纯文件同步，不认 git 规则）被同步过去，**覆盖了 hqlx210 上原本正常的 Linux ifort 二进制**，导致首次重跑直接报 `Exec format error`。定位后在 hqlx210 上 `make` 重新编译修复，并把这条踩坑记入 `.claude/persistent_context.md`，避免下次重犯。

### 意外阻塞：ClimKern 核数据下载

ClimKern 自带的 `python -m climkern download` 会拉取一个 **5.3GB** 的完整压缩包（内含大量本次用不到的其它模式核/地表核数据），而 Zenodo 官方在 2026 年因爬虫流量普遍限速（非本项目专属问题，已查证 Zenodo 官方博客确认），实测 Mac 直连速度约 76KB/s，若走默认下载器预计要 10+ 小时。
→ 改用 `remotezip`（基于 HTTP Range 请求），只从远程 zip 里精确抽取所需的两个文件（`TOA_CloudSat_Kerns.nc` + `TOA_GFDL_Kerns.nc`，各约130MB，合计252MB，体积降到官方下载器的 1/20）。经测试从 mini 中转下载速度提升到约350KB/s（约4.6倍），最终252MB在几分钟内下载完成。

> 说明：抽取用的 Zenodo record（`18565513`）是通过官方 `doi.org` 解析 climkern 库自带的 concept DOI（`10.5281/zenodo.10223376`）得到的当前版本记录，不是另找的替代数据源。

---

## 三、代码改动清单

### 新增核心模块（`core/`，不依赖 xesmf/climkern，任何机器可跑）

| 文件 | 内容 |
|---|---|
| `core/kernels.py` | `KernelSet` 类：读取 `TOA_<name>_Kerns.nc`，自实现 bilinear 重网格（scipy，经度周期wrap + 纬度边界clip近似climkern的nearest外插） |
| `core/lr_kernel.py` | 对流层顶公式、层厚计算、ΔR_LR/ΔR_PL 核心积分公式、垂直插值（严格复刻 ClimKern `calc_T_feedbacks` 的 NaN 传播顺序） |
| `core/lr_attribution.py` | M3 route i：把 `core/lr_kernel.py` 的机制套用到任意 CFRAM 分解项 `dT_X` 上，得到逐过程 ΔR_LR |

### 新增脚本（`scripts/`）

| 文件 | 内容 |
|---|---|
| `scripts/compute_lr_kernel.py` | CLI：case → `lr_kernel.nc`（原生实现，无需 xesmf） |
| `scripts/validate_lr_vs_climkern.py` | 唯一依赖 xesmf/climkern 的脚本：跑 ClimKern 参照 + 交叉验证指标 + 报告 |
| `scripts/plot_lr_comparison.py` | native/ClimKern/差异三联图 + Kramer vs GFDL 差异图 |
| `scripts/compute_lr_attribution.py` | CLI：case → `lr_attribution.nc`（逐过程归因） |
| `scripts/plot_lr_attribution.py` | 逐过程归因地图 + 纬向廓线图 |
| `data/kernel_source.py` | 从已安装的 climkern 包定位/搬运核数据到 `data/kernels/`（gitignored，带 md5 manifest） |

### 修改的既有文件

- `run_case.py`：新增 `--step lr` / `--step lr-attr` 分发
- `cases/cesm2_4xco2_official/case.yaml`：新增 `lapse_rate` 配置块（选核/选月份/选晴空模式）
- `README.md`：新增"Lapse-Rate Kernel Module"章节 + Key Scripts 表格条目
- `.gitignore`：新增 `data/kernels/`、`.coverage`、`.pytest_cache/`
- `.claude/persistent_context.md`（本地未纳入版本控制）：追加 rsync 踩坑记录

### 新增测试（`tests/`）

| 文件 | 覆盖内容 |
|---|---|
| `tests/test_kernels.py` | KernelSet 读取/重网格/月份选择；含一条仅在 `pycfram-kern` 环境跑的 xesmf 交叉检验 |
| `tests/test_lr_kernel.py` | 对流层层厚求和恒等式、均匀增温→ΔR_LR≡0、常数解析核手算对拍、NaN传播顺序 |
| `tests/test_lr_attribution.py` | 逐过程归因函数、可加性残差计算 |
| `tests/test_e2e_smoke.py` | 端到端：8×8裁剪的迷你 case + 粗化核（各约0.5MB）驱动完整 CLI 命令 |
| `tests/data/smoke/` | 上述 e2e 测试用的小规模固定数据（cfram_result裁剪、perturbed_surf裁剪、两个粗化核文件） |

### 新增文档（`docs/`）

- `docs/m2_kernel_module.md` — M2 模块技术说明（算法、架构、环境搭建、验证结果、已知限制、实现过程中修复的bug）
- `docs/m3_methodology_comparison.md` — M3 方法论对比正文（合同硬性交付物）
- `docs/m3_route_decision_memo.md` — M3 路线选择备忘录（route i vs route ii）
- `docs/ph2_report.md` — 本文档

---

## 四、关键 bug：netCDF masked array 静默失真

实现过程中最严重的一处 bug：`netCDF4` 库对带 `_FillValue` 属性的变量默认返回 **masked array**，而 `np.array(masked_array)` 在**不经过 `np.ma.filled()`** 的情况下，会把被 mask 掉的元素（例如核数据里地下无效层）静默替换成**原始填充值**（约 1×10³⁶）而不是 NaN。

首次跑通 `compute_lr_kernel.py` 时，CloudSat 核算出的 ΔR_LR 域均值达到 **-4.39×10³⁷ W/m²**（正常量级应为个位数 W/m²），排查后定位到这个问题，修复为全程使用 `np.ma.filled(arr, np.nan)`。修复后数值立即恢复正常量级（-1.90 W/m²，简单算术平均；面积加权平均见下节）。这个坑已写入 `docs/m2_kernel_module.md` 的"Gotcha"章节，供后续维护者参考。

---

## 五、数值验证结果

### M2：交叉验证（全球 `cesm2_4xco2_official`，192×288）

原生模块与 ClimKern 官方库喂**完全同一份**温度响应数据（`dT_observed`，即 `perturbed.ta − base.ta`），分别过核后做面积加权（cos纬度）空间相关性和域均值相对误差对比：

| 核 | 空间相关性 corr | 域均值相对误差 | 合同门槛（0.85 / 15%） |
|---|---|---|---|
| CloudSat（Kramer，主） | **0.9997** | **1.32%** | ✅ PASS |
| GFDL（参考） | **0.9999** | **0.02%** | ✅ PASS |

两者都大幅超过合同门槛，也超过计划里给自己定的"内部绿线"（corr≥0.98，rel_diff≤5%）。三联图（native / ClimKern / 差异）目视上两个 panel 几乎无法区分，差异图接近全场机器噪声水平，仅在南北极附近（约2°范围内）因核数据原生网格（2°×2.5°）在极点附近外推产生局部色标饱和，影响面积可忽略。

### M3：逐过程归因（route i）

对 `q, co2, o3, solar, albedo, cloud, aerosol, lhflx, shflx, atmdyn, ocndyn` 共11个 CFRAM 分解项逐一过核，得到各自的 ΔR_LR 贡献：

```
CloudSat   可加性残差(逐点绝对值均值) = 1.40 W/m² (占总量均值的21.8%)
GFDL       可加性残差(逐点绝对值均值) = 1.75 W/m² (占总量均值的23.1%)
```

**重要澄清**：这22%左右的残差是**逐点**统计量，不代表域均值系统偏差。实测面积加权**域均值**层面，Σ各过程ΔR_LR（CloudSat：−4.356 W/m²）与M2总量（−4.335 W/m²）相差**不到0.5%**——也就是说逐过程归因在全球平均意义上几乎精确，22%的残差是空间上相互抵消的散点（集中在相邻过程局部互相"过冲/欠冲"的地方），不是系统性偏差。这个非线性残差符合预期：CFRAM分解本身是能量收支的一阶（线性化）展开，各过程温度响应的辐射效应加总不严格等于总响应的辐射效应，这与session_log.md里此前记录的气溶胶/云可加性残差是同一性质，不是bug，文档里明确说明不强行闭合到0。

**各过程域均值贡献**（CloudSat核，W/m²）：水汽`q`=+9.09（最大正贡献，热带对流层上层增暖放大的经典信号）、潜热通量`lhflx`=−10.43（最大负贡献）、大气动力`atmdyn`=−4.80、海洋动力`ocndyn`=−2.99、反照率`albedo`=+2.82（南北极海冰/积雪消融驱动）、感热通量`shflx`=+1.69、CO2=+1.42、云`cloud`=−1.16，O3/太阳辐射/气溶胶≈0（本case未配置气溶胶强迫，O3在base/warm两态均固定为CESM 1850气候态）。

---

## 六、测试与覆盖率

```
pytest tests/ -v --cov=core --cov-report=term-missing
20 passed, 1 skipped in 1.35s
```

跳过的1项是仅能在 `pycfram-kern` 环境跑的 xesmf 交叉检验（重网格实现 vs xesmf 官方实现），已在该环境下单独验证通过（corr>0.999，最大绝对误差<核值域1%）。

新增模块行覆盖率：

| 模块 | 覆盖率 |
|---|---|
| `core/kernels.py` | 94% |
| `core/lr_kernel.py` | 100% |
| `core/lr_attribution.py` | 100% |

均远超合同M5口径的60%最低要求。

---

## 七、检查清单（建议核对顺序）

1. **验证报告**（数字最直接）：`cases/cesm2_4xco2_official/output/lr_validation_report.txt`
2. **6张图**（`cases/cesm2_4xco2_official/figures/`）：
   - `fig_lr_comparison_CloudSat.png` / `fig_lr_comparison_GFDL.png` — native vs ClimKern 三联图
   - `fig_lr_kramer_vs_gfdl.png` — 两核差异
   - `fig_lr_attribution_CloudSat.png` / `fig_lr_attribution_GFDL.png` — 逐过程归因地图
   - `fig_lr_zonal_profile.png` — 纬向廓线（对标论文Fig 3a风格）
3. **技术文档**：`docs/m2_kernel_module.md`（模块说明+已知限制+踩坑记录）、`docs/m3_methodology_comparison.md`（方法论对比正文）、`docs/m3_route_decision_memo.md`（route选择依据）
4. **测试**：`pytest tests/ -v --cov=core`（应为20 passed 1 skipped）
5. **PR**：`feat/m2-m3-lapse-rate-kernel` 分支
6. `scripts/plot_alaska_profile.py`（已修改）、`scripts/setup_paperdata_input.py`（未跟踪）是会话前遗留改动，本次未涉及，未纳入本次提交范围。

---

## 八、待办与后续

- **Route ii（ω动力代理）**未实施：`docs/m3_route_decision_memo.md`已记录决策依据（全仓库grep确认raw_data无任何`wap`/omega数据，需先补充CMIP6变量下载才可能开展），列为未来可选拓展项。
- M2/M3 均已通过合同数值门槛，PR待Code Review后合并。
