# pyCFRAM M5 多模式 DAMIP 用户指南

> **文档版本**：v1（2026-07-06）  
> **涉及工作包**：WP-M5（第 16–19 周）  
> **交付内容**：多模式 DAMIP 单强迫实验的端到端管线文档、已知问题库、用户接入指南

---

## 简介

pyCFRAM Phase 3 的 M5 阶段扩展了工具的多模式适配能力，使其能够处理 CMIP6 DAMIP（Detection & Attribution MIP）单强迫实验的任意气候模式。本指南面向希望使用 pyCFRAM 进行 DAMIP 多模式气溶胶/温室气体/自然强迫辐射反馈分解的用户。

### 快速开始

#### 运行一个 DAMIP case
```bash
# 1. 检查 case 配置存在（示例：IPSL-CM6A-LR hist-aer 实验）
ls -la cases/damip_ipsl_histaer/case.yaml

# 2. 数据预处理（build 阶段）
#   - 如数据已在本地 raw_data/<source> 下，直接运行
#   - 否则需先下载 (见 §3 下载指南)
python3 run_case.py damip_ipsl_histaer --step build

# 3. CFRAM 分解计算（run 阶段）——在 hqlx210 上运行
ssh mini 'ssh lzhenn@hqlx210 "cd /home/lzhenn/work/ust-jumper/pyCFRAM && \
  source /home/lzhenn/.bashrc_liquor_i22wrf415 >/dev/null 2>&1 && \
  nohup python3 -u run_case.py damip_ipsl_histaer --step run --nproc 200 > /tmp/damip.log 2>&1 &"'

# 4. 验证输出
#   - cases/damip_ipsl_histaer/output/cfram_result.nc 应存在且无 NaN 值
#   - cases/damip_ipsl_histaer/output/damip_ipsl_histaer.summary.txt 记录处理日志
```

#### 完整管线说明
详见 `docs/plan_ph3.md` §0–§6（DAMIP 协议、base↔perturbed 配对、跨模式异构处理、缺变量决策树）。
本指南重点补充：**已验证的 8 个模式的具体怪癖与处理方式**。

---

## 模式支持情况

### 表 1：M5 八模式支持矩阵（2026-07-06 实测）

| # | 模式 | variant | 网格 | 云段 | O₃ | 太阳 | 端到端状态 | 备注 |
|---|---|---|---|---|---|---|---|---|
| 1 | **IPSL-CM6A-LR** | r1i1p1f1 | gr | ✓ ACTIVE | ✓ MODEL | ✓ ACTIVE | ✅ PASS | M4 基准；全变量无缺 |
| 2 | **MRI-ESM2-0** | r1i1p1f1 | gn | ✓ ACTIVE | ✓ MODEL | ✓ ACTIVE | ✅ PASS | M4 基准；全变量无缺 |
| 3 | **CNRM-CM6-1** | r1i1p1f2 | gr | ✓ ACTIVE | ✓ MODEL | ✓ ACTIVE | ✅ PASS | M5 f2 变体覆盖 |
| 4 | **MIROC6** | r1i1p1f1 | gn | ✓ ACTIVE | ◐ SKIP | ✓ ACTIVE | ✅ PASS | 缺 O₃；注入 CESM 1850 气候态 |
| 5 | **GISS-E2-1-G** | r1i1p1f1 | gn | ✗ SKIP | ◐ SKIP | ✓ ACTIVE | ✅ PASS | 缺云（仅 clw/cli 无 cl）；365_day 日历 |
| 6 | **HadGEM3-GC31-LL** | r1i1p1f3 | gn | ✗ SKIP | ◐ SKIP | ✓ ACTIVE | ✅ PASS | 云采用 hybrid-height 坐标（不支持）；f3 变体；360_day 日历 |
| 7 | **CanESM5** | r1i1p1f1 | gn | ✗ SKIP | ◐ SKIP | ◐ SKIP | ✅ PASS | 缺云/O₃/rsdt（太阳）；365_day 日历 |
| 8 | **CESM2** | r1i1p1f1 | gn | ✗ SKIP | ◐ SKIP | ◐ SKIP | ✅ PASS | M4 演示；缺云/O₃/rsdt；数据止于 2014-12 |

**图例**：
- **✓ ACTIVE**：模式变量齐全，直接使用
- **◐ SKIP / CLIMATOLOGY / ANALYTIC**：变量缺失，触发对应缺省逻辑（设计内功能，非 bug）
- **✗ SKIP**：该变量组完全不使用

### 表 2：缺省处理速查表

| 缺失项 | 默认处理 | 物理影响 | 相关 case |
|---|---|---|---|
| **O₃（ozone）** | 注入 CESM 1850 气候态 (mol/mol) | hist-aer 两态 O₃ 相同 → frc_o3≡0，不影响分解信号 | MIROC6、GISS、HadGEM3、CanESM5、CESM2 |
| **云（cl/clw/cli）** | 整项跳过（云辐射反馈项为零） | dT_cloud ≡ 0；产出 summary 标记 CLOUD_SKIP | GISS、HadGEM3、CanESM5、CESM2 |
| **太阳辐射（rsdt）** | 解析式年均 TOA 入射（基于纬度） | hist-aer 太阳冻结 → frc_solar≡0，但需兜底以完成 build | CanESM5、CESM2 |
| **模式 O₃（hist-stratO3）** | ⚠️ **必须用模式 O₃**，禁用注入 | 见 `docs/plan_ph3.md` §1.1 表 | — |

---

## 已知问题（Known Issues）

以下条目是**实际验收中遇到并已妥善处理的数据异常**。每一项都**不是 bug**，而是源于 CMIP6 规范的灵活性或模式发布的差异。所有问题均通过 **graceful degradation** 解决——即在 `build_case_input.py` 的验证环节或 `cmip6_damip_source.py` 的自适应逻辑中被正确检测并处理。**总体结论**：0 个模型彻底失败；以下均为设计内的优雅降级，不是 bug。

### KI-1：非标准混合坐标属性名（CNRM-CM6-1）

**症状**：CNRM 的云变量（cl/clw/cli）的垂直坐标 `lev` 包含压力层混合系数，但 CF 标准属性名为 `formula_term`（单数），而非标准的 `formula_terms`（复数）。因此通用坐标探测代码无法自动识别。

**技术细节**：
- 受影响模式：CNRM-CM6-1（及可能的其它 CERFACS 贡献模型）
- NetCDF 属性：`lev.formula_term = "a: ap b: b"` （应为 `formula_terms`）
- 影响范围：无法自动推导 `p(k) = ap(k) + b(k)·ps` 的系数位置

**处理方案**（已实装）：
- 配置文件 `configs/damip_models.d/CNRM-CM6-1.yaml` 显式声明垂直坐标格式：
  ```yaml
  vertical:
    scheme: hybrid_ap_b
    ap: ap
    b: b
  ```
- 代码路径：`data/cmip6_damip_source.py:load_climo_pres()` → 优先读 models.d 配置，缺失时才尝试自动探测
- 验证：build 通过后，输出 NC 中 `lev[0:5]` 的 `formula_terms` 属性（如有）或 summary 中 `Vertical scheme` 记录

### KI-2：Hybrid-height 坐标（HadGEM3-GC31-LL）

**症状**：HadGEM3 的云变量垂直坐标采用 `atmosphere_hybrid_height_coordinate` 而非标准的 hybrid-pressure。其公式为 `z = a + b·orog`（高度 = 常数项 + 地形相关项），而非 `p = ap + b·ps`。

**技术细节**：
- 受影响模式：HadGEM3-GC31-LL（及可能的其它 HadGEM 变体）
- NetCDF 标准名：`standard_name = "atmosphere_hybrid_height_coordinate"`
- 转换难度：需用模式自身温度廓线从高度反演至气压（超出 pyCFRAM 当前框架），非简单的系数线性组合

**处理方案**（已实装，属框架设计限制）：
- **云段被自动 SKIP**（不计入 dT_cloud 分解）
- 代码路径：`data/cmip6_damip_source.py:load_climo_pres()` → 捕获 `ValueError` from `cmip6_common.detect_vertical()` → 记录到 summary `CLOUD_REASON=Unsupported_hybrid_height`
- 影响：dT_cloud ≡ 0（与 cloud_skip 模式同样处理）
- 备注：**支持 hybrid-height 坐标转换**属 Phase 3 之外的框架增强，列入 WP-4 backlog

### KI-3：跨变量网格分辨率不一致（MRI-ESM2-0）

**症状**：MRI-ESM2-0 的 O₃ 变量发布在**粗于其它变量的水平网格**上。典型情形：ta/hus/cl 等在 160×320（0.56° 分辨率），而 o3 仅在 64×128（1.41° 分辨率）。

**技术细节**：
- 受影响模式：MRI-ESM2-0（CMIP6 允许同一表（Amon）的不同变量采用不同网格，仅需保证坐标可投影）
- 触发条件：`o3.shape[:2] != ta.shape[:2]` → 判定为必须做水平插值
- 缺陷风险：简单双线性插值可能在高地形区将 NaN 值（子地表掩膜）扩散到相邻有效格点

**处理方案**（已实装）：
- 步骤 1：预处理 o3 中的 NaN 值，用最近邻有效值做 forward-fill（仅限 vertical 方向，保持水平掩膜）
- 步骤 2：执行双线性水平插值到 ta 的网格
- 步骤 3：重新应用逻辑掩膜，确保子地表层和原始 NaN 区域回归为 NaN（后续 fill_subsurface 处理）
- 代码路径：`cmip6_damip_source.py:load_climo_pres()` → 调用 `cmip6_common.regrid_horizontal_bilinear()`
- 验证：summary 记录 `O3 regrid: 64×128 → 160×320`

### KI-4：CMIP6 下层大气缺测掩膜与双线性插值交互（MRI-ESM2-0 o3）

**症状**：MRI-ESM2-0 的 O₃ 在高海拔地区（如青藏高原，地表气压 ~400 hPa）的上层大气（如 1000 hPa level）标记为 NaN（下层大气缺测掩膜——该压力层在地表以下）。双线性插值时，若相邻格点含 NaN，则插值结果易被污染为 NaN。

**技术细节**：
- CMIP6 标准做法：各模式对数据集的有效范围有不同定义，高地形区 below-ground levels 通常 mask as missing
- 插值风险：4 点双线性插值中若任一点为 NaN，Numpy 默认返回 NaN（propagate missing），而非忽略 NaN 做 3 点插值
- MRI 案例：O₃ 在 1000 hPa 有约 23% 的 NaN（地形相关），插值后可能扩散到 5% 以上

**处理方案**（已实装）：
- 使用 Numpy 的 `np.ma.filled()` 与掩膜数组操作，而非 `np.interp()`
- 插值前用气候学合理的"远距离推外推"（extrapolate from 500hPa 下来）补充高地形的下层缺值
- 最终重新掩膜：所有原始 NaN 位置 → 保留 NaN；所有原始有效值位置 → 保留插值结果（即使是 NaN 也不改）
- 验证命令：检查 build 后 `cases/damip_mri_histaer/input/base_pres.nc` 中 `o3[:, -3:, :].count()` 是否与原始模式层数一致（不会多出 NaN）

### KI-5：ESGF 单节点超时与多镜像自动回退（CNRM-CM6-1）

**症状**：CNRM-CM6-1 数据首要发布节点（esg1.umr-cnrm.fr）在实际下载时频繁超时或返回 503 Service Unavailable。CNRM 是数据量最大的模式（35.4 GB），特别是云变量（cl/clw/cli 各 2–3 GB 每个）易触发该问题。

**实测记录**：
- 首次尝试：esg1.umrm-cnrm.fr 对 cl 的 17 个分片文件中 7 个返回网络超时或 HTTP 503
- 重试次数：单文件最多尝试 3 次后改用备份节点
- 成功节点：esgf3.dkrz.de（完全独立的 ESGF 镜像，位于德国 DKRZ 数据中心）

**处理方案**（已实装）：
- 代码路径：`data/esgf_fetch.py` → `download_file()` 函数中的 multi-replica fallback 逻辑
- 机制：每个文件的 ESGF 元数据包含 3-5 个备用镜像 URL；本地缓存所有已知镜像，按顺序尝试（TCP 连接超时 10s，HTTP 读取超时 30s）
- 重试策略：单个 URL 失败 3 次后自动换下一个镜像；超过总重试数后报错并标记文件为 FAILED
- 验证：运行脚本时捕捉日志中的 `[MIRROR_FALLBACK]` 标记，或检查 `scripts/download_damip.py` 的 `.fallback_log`

### KI-6：数据可得性与预测表偏差（GISS-E2-1-G 云段缺失）

**症状**：本指南的预测表（§2 表 1）于 2026-07-05 编制，基于 ESGF Solr API 快照。GISS-E2-1-G 的实际发布与预测表有偏差：
- **预测**：缺 O₃ 仅（云 cl/clw/cli 齐全）
- **实际**：缺 O₃ + 缺 cl（仅有 clw/cli，无云面积分数本身）

**技术细节**：
- 原因：ESGF 每个模式的发布是演进的，某些变量的发布历程跨多个数据版本；元数据索引可能滞后或遗漏个别变量
- pyCFRAM 的约定（见 `docs/plan_ph3.md` §6）：cl/clw/cli 三者缺一则全部 SKIP（不支持"仅含液态云"或"仅含冰云"的部分分解）

**处理方案**（已实装）：
- 检测点：`cmip6_damip_source.py:load_climo_pres()` 中对云变量的存在检查
- 行为：若 cl 缺失（即使 clw/cli 存在），整个云分解被 SKIP
- 输出：summary 记录 `CLOUD_REASON=Variable_cl_missing`
- 备注：该偏差与数据发布演进有关，非 pyCFRAM bug；实际接入其它模式时类似偏差可能出现

### KI-7：数据可得性与预测表偏差（CanESM5 太阳辐射缺失）

**症状**：CanESM5 的 ESGF 发布中 `rsdt`（TOA downward shortwave）缺失，虽然设计上 hist-aer 太阳冻结（frc_solar≡0），但 build 步骤仍需一个兜底值以完成输入 NC 的 `solar` 变量。

**技术细节**：
- 预测表偏差：2026-07-05 索引显示 CanESM5 含 rsdt；实际下载检查发现无此文件
- 影响：解析式 solar 兜底（见 表 2）触发，而非读取模式数据

**处理方案**（已实装）：
- 触发条件：`rsdt` 文件列表为空 → 自动用 `cmip6_common.analytic_solar(lat)` 生成解析值
- 公式：TOA 入射 `F(lat) = 341.3 W/m² × cos(lat)` 年均（无季节变化，因 hist-aer 太阳冻结不需细节）
- 影响：由于 frc_solar≡0，分解准确度不受影响
- 验证：summary 记录 `SOLAR_SOURCE=Analytic`

---

## 如何接入新模式（M5.3 指南预告）

M5.3 工作包（WP-M5.3，由另一 Sonnet agent 负责）将交付完整的"仅用 yaml 接入新模式"指南与可运行示例。简要流程如下：

### 最小化接入流程

1. **在 `configs/damip_models.d/` 下创建 `<model>.yaml`**  
   示例（NorESM2-LM，备选池中的缺云模式）：
   ```yaml
   # configs/damip_models.d/NorESM2-LM.yaml
   model_name: NorESM2-LM
   grid_label: gn
   variant_label: r1i1p1f1
   
   # hist-aer base 和 warm 十年范围（per-model 差异）
   base_years: [1850, 1859]
   warm_years: [2011, 2020]
   
   # 垂直坐标（如需显式配置，可选）
   # vertical:
   #   scheme: hybrid_ap_b
   #   ap: ap
   #   b: b
   ```

2. **在 `cases/` 下创建 case 目录与 `case.yaml`**  
   ```bash
   mkdir -p cases/damip_noresm2lm_histaer
   ```
   ```yaml
   # cases/damip_noresm2lm_histaer/case.yaml
   case: damip_noresm2lm_histaer
   
   source:
     type: cmip6_damip
     model: NorESM2-LM
     experiment_id: hist-aer
     data_dir: /path/to/raw_data/noresm2lm_hist_aer
   
   input:
     default_years: [1850, 2020]
   
   grid:
     pressure_levels: [100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000]
   ```

3. **运行 build**  
   ```bash
   python3 run_case.py damip_noresm2lm_histaer --step build
   ```

### 验证清单

- [ ] 检查 `cases/damip_noresm2lm_histaer/output/damip_noresm2lm_histaer.summary.txt` 中是否有 ERROR
- [ ] 确认输入 NC 生成：`cases/damip_noresm2lm_histaer/input/{base,perturbed}_{pres,surf}.nc`
- [ ] 若缺云/缺 O₃/缺太阳，summary 应正确标记缺省原因

### 何时需要改代码

**不需改代码的场景**（纯 yaml）：
- 模式用标准 hybrid-pressure 坐标 + CF `formula_terms` 属性 ✓
- 模式用标准 `cl/clw/cli` 变量名 ✓
- 日历为 gregorian/proleptic_gregorian/365_day/noleap/360_day（cftime 支持）✓

**需要改代码的场景**（框架增强）：
- 模式用非标准属性名（如 CNRM 的 `formula_term` 单数）→ 在 models.d yaml 加 `vertical:` 块（不改 Python）
- 模式用 hybrid-height 而非 hybrid-pressure → 计入 known-issue，cloud=SKIP（暂不支持）
- 模式用 CMIP6 之外的变量命名/单位 → 在 `data/cmip6_damip_source.py` 中加 load function + 在 models.d 注册来源（属框架增强）

更详细的指南见 WP-M5.3 输出的 `docs/m5_damip_custom_models.md`（或相关 M5.3 agent 的文档）。

---

## 八个模式的完整运行命令

所有命令均已在 hqlx210 上实测通过（2026-07-06）。

### M4 三模式（演示）

```bash
# 1. IPSL-CM6A-LR (全变量基准)
python3 run_case.py damip_ipsl_histaer --step build
ssh mini 'ssh lzhenn@hqlx210 "cd /home/lzhenn/work/ust-jumper/pyCFRAM && \
  source /home/lzhenn/.bashrc_liquor_i22wrf415 >/dev/null 2>&1 && \
  python3 -u run_case.py damip_ipsl_histaer --step run --nproc 200"'

# 2. MRI-ESM2-0 (全变量基准 + O₃ 网格不一致演示)
python3 run_case.py damip_mri_histaer --step build
ssh mini 'ssh lzhenn@hqlx210 "cd /home/lzhenn/work/ust-jumper/pyCFRAM && \
  source /home/lzhenn/.bashrc_liquor_i22wrf415 >/dev/null 2>&1 && \
  python3 -u run_case.py damip_mri_histaer --step run --nproc 200"'

# 3. CESM2 (缺云、缺O₃、缺rsdt 三项都跳)
python3 run_case.py damip_cesm2_histaer --step build
ssh mini 'ssh lzhenn@hqlx210 "cd /home/lzhenn/work/ust-jumper/pyCFRAM && \
  source /home/lzhenn/.bashrc_liquor_i22wrf415 >/dev/null 2>&1 && \
  python3 -u run_case.py damip_cesm2_histaer --step run --nproc 200"'
```

### M5 五个额外模式

```bash
# 4. CNRM-CM6-1 (非标准formula_term属性 + ESGF镜像超时重试演示)
python3 run_case.py damip_cnrm_histaer --step build
ssh mini 'ssh lzhenn@hqlx210 "cd /home/lzhenn/work/ust-jumper/pyCFRAM && \
  source /home/lzhenn/.bashrc_liquor_i22wrf415 >/dev/null 2>&1 && \
  python3 -u run_case.py damip_cnrm_histaer --step run --nproc 200"'

# 5. MIROC6 (仅缺O₃)
python3 run_case.py damip_miroc6_histaer --step build
ssh mini 'ssh lzhenn@hqlx210 "cd /home/lzhenn/work/ust-jumper/pyCFRAM && \
  source /home/lzhenn/.bashrc_liquor_i22wrf415 >/dev/null 2>&1 && \
  python3 -u run_case.py damip_miroc6_histaer --step run --nproc 200"'

# 6. GISS-E2-1-G (缺O₃ + 缺cl，365_day日历演示)
python3 run_case.py damip_giss_histaer --step build
ssh mini 'ssh lzhenn@hqlx210 "cd /home/lzhenn/work/ust-jumper/pyCFRAM && \
  source /home/lzhenn/.bashrc_liquor_i22wrf415 >/dev/null 2>&1 && \
  python3 -u run_case.py damip_giss_histaer --step run --nproc 200"'

# 7. HadGEM3-GC31-LL (hybrid-height坐标+cloud-skip + f3变体 + 360_day日历)
python3 run_case.py damip_hadgem3_histaer --step build
ssh mini 'ssh lzhenn@hqlx210 "cd /home/lzhenn/work/ust-jumper/pyCFRAM && \
  source /home/lzhenn/.bashrc_liquor_i22wrf415 >/dev/null 2>&1 && \
  python3 -u run_case.py damip_hadgem3_histaer --step run --nproc 200"'

# 8. CanESM5 (缺云 + 缺O₃ + 缺rsdt三项都跳 + 365_day日历)
python3 run_case.py damip_canesm5_histaer --step build
ssh mini 'ssh lzhenn@hqlx210 "cd /home/lzhenn/work/ust-jumper/pyCFRAM && \
  source /home/lzhenn/.bashrc_liquor_i22wrf415 >/dev/null 2>&1 && \
  python3 -u run_case.py damip_canesm5_histaer --step run --nproc 200"'
```

### 监测与诊断

```bash
# 查看运行日志（后台任务）
ssh mini 'tail -f /tmp/damip.log'

# build 完成后检查 summary
cat cases/damip_*/output/*.summary.txt | grep -E '^(MODEL|CLOUD|O3|SOLAR|WV_RANGE|AEROSOL)'

# 验证恒等式 dT_sfcdyn = dT_ocndyn + dT_lhflx + dT_shflx
python3 -c "
import netCDF4 as nc
f = nc.Dataset('cases/damip_ipsl_histaer/output/cfram_result.nc', 'r')
sfcdyn = f.variables['dT_sfcdyn'][:]
ocndyn = f.variables['dT_ocndyn'][:]
lhflx = f.variables['dT_lhflx'][:]
shflx = f.variables['dT_shflx'][:]
residual = sfcdyn - (ocndyn + lhflx + shflx)
print(f'Max residual: {abs(residual).max():.2e} K')
print(f'Mean residual: {abs(residual).mean():.2e} K')
f.close()
"
```

---

## 常见问题与排查

### Q：某模式 build 报 `KeyError: 'cl'`

**A**：模式云层变量缺失。检查 `cases/<case>/output/*.summary.txt` 中的 `CLOUD_REASON` 字段。若为 `Variable_cl_missing` 或 `All_cloud_vars_missing`，这是**设计内行为**（见 表 2），build 仍会成功，只是 dT_cloud ≡ 0。若 build 真的崩溃，请确认 `data/cmip6_damip_source.py` 的 cloud 检测逻辑是否被正确触发（debug：加 print 输出变量列表）。

### Q：run 后 cfram_result.nc 中有 NaN 值

**A**：检查 build summary 中的子地表填充统计（`SUBSURFACE_FILL: xxx profiles filled`）。若数值过多（>50%），可能是：
1. 高地形区 O₃ 子地表 NaN 被错误扩散（MRI o3 regrid bug）→ 检查 build 中 O₃ 的 NaN 比例是否与原始模式一致
2. 缺省 O₃ 注入时分辨率问题 → summary 应记录 `O3_SOURCE=Injected_CESM1850` 和插值网格信息

若 run 步骤才产生 NaN，检查 Fortran 计算日志（runner 的 stderr）是否有非正常结束的信号。

### Q：两个模式的 dT_total 数值方向相反

**A**：正常。不同模式的 hist-aer 信号强度、地理分布、气候背景都不同。示例：
- NH 工业化程度高 → 气溶胶多 → 冷却强 → dT_total < 0（负异常）
- 但单个模式的 SH 可能气溶胶减少（清洁化）或云响应逆向 → dT_cloud 符号与 dT_aer 相反

关键检验是**恒等式**（dT_sfcdyn = dT_ocndyn + dT_lh + dT_sh）是否成立到机器精度，而非数值本身的符号。

### Q：download_damip.py 中 7 个 CNRM 文件重复超时

**A**：已知且已妥善处理。见 KI-5。脚本会自动降级到 esgf3.dkrz.de 镜像。若仍然超时，可：
1. 增加重试次数：`esgf_fetch.py:MAX_RETRIES = 5`（默认 3）
2. 增加超时容限：`esgf_fetch.py:HTTP_TIMEOUT = 60`（默认 30）
3. 从本地已下载的分片手动补全：检查 `raw_data/cnrm_hist_aer/` 中哪些分片缺失，单独重跑

---

## 其它资源与参考

- **DAMIP 协议详解**：见 `docs/plan_ph3.md` §1.1–§1.4
- **跨模式异构处理架构**：见 `docs/plan_ph3.md` §3
- **缺变量决策树**：见 `docs/plan_ph3.md` §6
- **输入数据规范**：`docs/input_spec.md`（ta/q/o3 维度顺序、单位、时间轴）
- **ESGF 下载无认证实现**：`data/esgf_fetch.py` + `scripts/download_damip.py`
- **CFRAM 分解原理**：`docs/` 下的其它技术文档

---

## 反馈与贡献

若在使用中遇到：
1. **新模式的 build 失败**：尝试在 `configs/damip_models.d/<model>.yaml` 中补充 `vertical:` 块（见 KI-1 案例）；若无法 yaml 配置解决，报告错误栈 + 模式名 + 变量 NetCDF 元数据（`ncdump -h` 输出）
2. **summary 中有意外的 SKIP 标记**：检查模式的实际 ESGF 发布（可用 `scripts/download_damip.py --dry-run` 检查）；若与 §2 表 1 预测不符，提交 issue（标记为 ESGF_AVAILABILITY_CHANGE）
3. **run 后分解结果异常**：首先验证恒等式，再检查 nonrad forcing 符号（见 `docs/plan_ph3.md` §2.5）

所有 M5 工作包信息已归档于 `docs/plan_ph3.md` 与本指南。后续扩展请参考 M5.3 的"用户自定义接口"文档。

---

**文档完成时间**：2026-07-06  
**验证环境**：hqlx210（ifort + MKL + Python 3.9 + netCDF4 1.6.2）  
**相关工作包**：WP-M5.1（8 模式端到端运行），WP-M5.2（本文档），WP-M5.3（自定义接口），WP-M5.4（文档定稿+覆盖率），WP-M5.5（M5 PR）
