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

- 测试：181 passed, 1 skipped（跳过项为需真实网络的可选检验）
- 新增模块行覆盖率：`esgf_fetch.py` 98-99%，`cmip6_common.py` 98%，`cmip6_damip_source.py` 92-93%，均超合同60%门槛
- 代码、文档、9个 case 已全部提交并推送到 `feat/m4-m5-damip`
- 本次同步创建 PR，待走 code review 合并流程
