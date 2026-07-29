# Phase 3（M4/M5）完成汇报

**日期**：2026-07-06
**对应合同**：`contract/contract_10.pdf` Phase 3 — M4（DAMIP 单强迫 CMIP6 预处理模块，≥3 模式）+ M5（多模式扩展，≥8 模式接入、≥6 通过）
**执行依据**：`docs/plan_ph3.md`（v2）
**分支**：`feat/m4-m5-damip`，基于 Phase 2 已合并分支 `a87611d` 起共 30 次提交

---

## 一、目标与范围

DAMIP（Detection & Attribution MIP）hist-aer 实验：CO2/O3/太阳辐射冻结在 1850 年水平，只有气溶胶强迫历史变化，用于把 CFRAM 分解结果归因到"纯气溶胶效应"。M4 要求打通≥3个 CMIP6 模式的下载→预处理→CFRAM 运行全链路；M5 要求扩到≥8模式且≥6个真正跑通，并证明新增模式不需要碰 `core/`/`fortran/` 核心代码。

## 二、做了什么

- **`data/esgf_fetch.py`**：纯 stdlib 的 ESGF 检索+下载客户端（Solr 搜索 + THREDDS 下载），免账号认证，多镜像自动重试。
- **`data/cmip6_common.py`**：从既有 CESM2 代码里抽出的模式无关公共层（时间解码、垂直坐标探测、混合坐标转标准气压层、水平重网格、O3 气候态等）。
- **`data/cmip6_damip_source.py`**：新的可插拔数据源插件 `cmip6_damip`，与既有 `era5`/`cesm2_cmip6` 源用同一套 `run_case.py`/`build_case_input.py` 调度框架。
- **`configs/damip_models.d/*.yaml`**：每模式一个"怪癖"配置文件（变量缺失表、混合坐标 scheme、日历类型等），新增模式=加一个 yaml。
- **9 个真实数据 case**（`cases/damip_*_histaer/`）：IPSL-CM6A-LR、MRI-ESM2-0、CESM2、CNRM-CM6-1、MIROC6、GISS-E2-1-G、HadGEM3-GC31-LL、CanESM5（M4/M5 八模式）+ NorESM2-LM（M5.3 自定义接入示例）。
- 配套 `scripts/write_run_summary.py`（自动生成可加性核验+复现命令的运行摘要）、`docs/m4_damip_module.md`、`docs/m5_multimodel_userguide.md`（含7条真实踩坑 KI-1~KI-7）。

## 三、效果如何

用真实下载的 ESGF 数据逐模式跑通后，发现并修复了 9 类真实数据结构问题（而非合成测试能暴露的），包括：IPSL 缺 `formula_terms` 属性、界面层↔层中值坐标不匹配、混合坐标存储方向（地表→TOA 与既有假设相反）；MRI 的 O3 在与其余变量不同的水平网格上、且含地下层 NaN；CNRM 用非标准的单数 `formula_term` 属性名、且官方节点下载超时需多镜像兜底；HadGEM3 是混合高度坐标（非混合气压），通过优雅降级把云处理跳过而非整体失败；`http.client.RemoteDisconnected` 未被原有异常捕获逻辑覆盖。

**结果**：8/8 模式全部构建+运行成功（超过合同"≥6/8"门槛），全部通过物理合理性检验（地表净冷却、北半球比南半球更冷，气溶胶集中在北半球的预期信号）与恒等式检验（`dT_sfcdyn = dT_ocndyn+dT_lhflx+dT_shflx` 残差达机器精度）。NorESM2-LM 作为"新模式接入"示例验证：只加了1个 yaml 模式配置+1个 case.yaml，`core/`、`fortran/` 目录零改动（有专门回归测试 `test_phase3_diff_never_touches_core_or_fortran` 断言）。

## 四、状态

- 测试：194 passed, 1 skipped（跳过项为需真实网络的可选检验）
- 新增模块行覆盖率：`esgf_fetch.py` 98%，`cmip6_common.py` 98%，`cmip6_damip_source.py` 93%，`lr_kernel/lr_attribution` 100%，`kernels.py` 94%，均超合同60%门槛
- 代码、文档、9个 case 已全部提交并推送到 `feat/m4-m5-damip`
- 本次同步创建 PR，待走 code review 合并流程

## 五、复审与修复轮（2026-07-29）

对 Phase 2/3 全部产出做了一轮独立复审——不复用开发期结论，直接从 `cfram_result.nc`
重算验收指标、对活的 ESGF 端点验证下载逻辑、逐模块读代码。结论：**主结论全部成立**
（8/8 模式跑通、物理检验与恒等式检验确实通过、覆盖率数字属实、9 个输出 md5 互不相同），
另发现并修复 3 个真实缺陷 + 2 处交付缺口：

| # | 问题 | 性质 | 处理 |
|---|---|---|---|
| 1 | `esgf_fetch.list_files()` 的 HTTPServer 识别对 ESGF 真实格式 `url\|mime\|SERVICE` **永远不匹配**，静默退化为"取第一个 http 开头的串"——当 OPENDAP 排在前面时会取回 `.nc.html` 浏览页当 NetCDF 下载 | 真实 bug（潜伏）；原测试用的是臆造格式，给了假信心 | 重写为按服务令牌匹配 + 排除 OPENDAP 端点；补 4 个用真实 ESGF 报文格式的测试 |
| 2 | `fetch()` 断点续传：服务器忽略 Range 头返回 200 全量时仍按 `'ab'` 追加 → 生成 (残片+全量) 的损坏文件；有 checksum 时报"莫名不匹配"，无 checksum 时静默损坏 | 真实 bug | 仅在确认 206 Partial Content 时才追加，否则整体覆写；补回归测试 |
| 3 | `cmip6_damip_source.build_states()` 对 base/warm 每个字段都做了 `normalize_grid`，唯独漏了 `nonrad`（lhflx/shflx）；而 writer 是按**归一化后**的 lat/lon 轴写这两个场的 | 真实 bug（潜伏）——当前 9 个模式原生均为 lat S→N、lon 0-360 升序，故已交付结果无一受影响；换一个 N→S 发布的模式则地表通量强迫会在纬向镜像且不报错 | 补上同一置换；测试用 lat 翻转的 fixture 副本验证（去掉修复即失败） |
| 4 | M4 验收门（净冷却、NH<SH、`dT_sfcdyn=dT_ocndyn+dT_lhflx+dT_shflx` 残差 <1e-10）只在开发期手工核过一次，**没有落进任何产物**——审阅者无从复核 | 交付缺口 | `write_run_summary.py` 新增 "Acceptance gates" 段，每次 run 从 `cfram_result.nc` 重算；9 个 case 的 summary 已全部重生成，全部 PASS（恒等式残差 3.6e-15 ~ 1.6e-14 K，地表 NaN 占比 0） |
| 5 | `docs/technical_notes_{en,zh}.md` 有 Phase 3 章节但**完全没有 Phase 2 辐射核模块**（grep `kernel\|lapse` 命中 0），与 M5 交付物 3"文档覆盖 Phase 2 + Phase 3 全部新增功能"不符 | 文档缺口 | 两个语种各补一章（算法要点、免 xesmf 决策、与 ClimKern 的 0.9997/0.9999 一致性、M3 归因、masked-array 陷阱、入口脚本表） |

复审中特别核查、确认**不是** bug 的两点，记录以免后人重复怀疑：

- `_climo_pair_for_variable()` 里的 `np.asarray(nc.variables[...][...])` 确实会把
  masked 元素还原成原始 fill value（~1e20）而非 NaN（已在真实 CMIP6 文件上实测确认），
  但 `annual_climo_from_monthly()` 在入口用 `|value| > 1e15 → NaN` 拦截了，
  `fill_subsurface()` 也有同样的护栏——是有意设计，不是 Phase 2 那个 masked-array 事故的重演。
- 输出 NC 里各 `dT_` 场在 1000/925/850 hPa 有 3.6–3.9% 的 NaN，与 Phase 2 基线
  `cesm2_4xco2_official` 的分布逐层一致（43.4%/15.5%/10.4%…），是地形以下层的正常掩膜，
  不是泄漏；`dT_observed` 地表行 NaN 为 0。

修复后：194 passed / 1 skipped，覆盖率无回退，9 个 case 的 summary 已在 hqlx210 重生成并全部通过验收门。
