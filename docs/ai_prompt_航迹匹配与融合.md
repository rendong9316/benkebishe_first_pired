# 航迹匹配与融合实现任务

## 项目背景

双基地外辐射源雷达仿真，两部雷达（R1、R2）独立对3架民航机目标进行UKF跟踪，
已在 `run_simulation.m` 的 Phase 5 完成。现在需要在 Phase 5 之后增加 Phase 6（航迹匹配）和 Phase 7（航迹融合）。

## 开题报告中的匹配算法（双门限法）

来自开题报告 4.3.2 节：

1. **时间对齐预处理**：将R1和R2航迹的状态点通过CV模型外推/内插到统一时间网格（30s整秒），使两源时间戳一致。

2. **双门限关联判决**：
   - 第一门限（距离门限 T₁）：计算时间对齐后两条航迹对应点的欧氏距离或马氏距离，< T₁ 进入候选。
   - 第二门限（次数门限 T₂）：对候选航迹对滑窗统计，连续 N 次判断中满足 T₁ 的次数 > T₂ 时确认关联。
   - 关联质量修正：对短航迹（life < 阈值）适当放宽门限。

3. **全局配对**：多目标场景下，每条R1航迹只能匹配一条R2航迹，需要用代价矩阵+贪心分配（类似JNN思路）解决冲突。

4. **评估指标**：正确匹配概率（目标 > 80%）、漏关联概率、虚关联概率。

## 开题报告中的融合算法（四种）

来自开题报告 4.5 节，需全部实现并对比：

1. **SCC（简单凸组合）**：假设两源估计误差独立，P⁻¹x = P₁⁻¹x₁ + P₂⁻¹x₂, P⁻¹ = P₁⁻¹ + P₂⁻¹
2. **Bar-Shalom-Campo**：考虑互协方差 Pᵢⱼ，需在UKF内部递推 Pᵢⱼ(k|k) = [I-K₁H₁]·[F·Pᵢⱼ·F'+Q]·[I-K₂H₂]'
3. **CI（协方差交叉）**：无需互协方差，保守融合 P⁻¹ = ωP₁⁻¹ + (1-ω)P₂⁻¹，ω 通过最小化 det(P) 或 trace(P) 求取
4. **FCI（快速协方差交叉）**：ω = tr(P₁)⁻¹ / (tr(P₁)⁻¹ + tr(P₂)⁻¹)，无需迭代优化

注意：BC需要互协方差递推，当前UKF代码（`ukf/ukf_filter_update.m`）没有维护此变量，需要修改。

## 现有代码上下文（关键！必须使用这些精确变量名和路径）

### 运行完 Phase 5 后已有的变量：
```matlab
% 时间网格
t1_grid  % R1采样时刻: 0, 30, 60, ...  (1×n_frames)
t2_grid  % R2采样时刻: 13, 43, 73, ...  (1×n_frames)
n_frames % 仿真帧数

% R1/R2最终航迹列表（cell array）
trackList_R1  % 每元素为航迹struct，R1通常有3条RELIABLE航迹
trackList_R2  % R2通常有3-4条航迹（含可能鬼影）

% 航迹struct字段：
%   trk.id      - 航迹编号
%   trk.type    - 1=RELIABLE, 2=MAINTAIN, 6=TEMPORARY, 7=HISTORY
%   trk.quality - 航迹质量分
%   trk.life    - 存在帧数
%   trk.ukf.x   - UKF状态 [lon; lon_dot; lat; lat_dot] (4×1)
%   trk.ukf.P   - UKF协方差 (4×4)
%   trk.lat, trk.lon - 当前位置

% 逐帧快照
trackSnapshots_R1{k}  % struct带字段 frameID, trackList
trackSnapshots_R2{k}

% 真实航迹（用于验证匹配正确性）
truthTrajs  % cell array of struct, 每个含 label/time_sec/lat/lon/lon_rate/lat_rate
```

### 已有但未接入的工具函数：
```
fusion/time_align_tracks.m  — 用CV模型将R2航迹状态回退13s到R1时间网格
                             函数签名: aligned_R2 = time_align_tracks(trackSnapshots_R2, params)
fusion/regularize_cov.m     — 协方差正则化
params.gate_sigma = 3.0     — 关联波门σ
params.dt_sec = 30          — 采样间隔
```

### 航迹struct中的trk.type含义：
- 1 = RELIABLE（稳定跟踪）
- 7 = HISTORY（已终止，应跳过）

### 关于时间对齐：
`time_align_tracks.m` 已实现：对 `trackSnapshots_R2{k}.trackList{t}` 中的每条非HISTORY航迹，
用 `dt = -13s` 构造 CV 状态转移矩阵 F，做 `x_aligned = F * x` 和 `P_aligned = F*P*F' + Q`。

## 验收标准

### 匹配阶段：
1. 程序无报错跑通
2. 打印 R1-R2 航迹配对表，验证3架飞机全部正确匹配（R1#1↔R2#?, R1#2↔R2#?, R1#3↔R2#?）
3. R2的鬼影航迹（HISTORY类型）被正确排除
4. 正确匹配概率 > 80%（目标值）

### 融合阶段：
1. 四种融合算法全部实现
2. 每种算法跑一遍，输出：融合后 vs 单站的RMSE、中位误差对比表
3. 融合航迹误差 < min(R1误差, R2误差)（即融合有正向收益）
4. 对比四种算法的优劣（精度 vs 计算量 vs 假设条件）

### 可视化：
- 在现有地图上叠加图层，用不同颜色显示：R1航迹、R2航迹、匹配后融合航迹
- 每条航迹标注来源（如"R1-航迹#1→飞机A"、"融合航迹-飞机A"）
- 另画误差收敛曲线（横轴=帧号，纵轴=位置误差km），R1/R2/融合三条线对比

## 约束

- 不要修改 `ukf/` 下已有的核心UKF函数（`ukf_filter_update.m`, `ukf_predict_step.m` 等），
  如果BC算法需要互协方差，在 `fusion/` 下新建函数或在航迹struct中附加字段，不要改动现有UKF。
- 新代码全部放在 `run_simulation.m` 的 Phase 6/7（或新建 `fusion/` 下文件 + `tracker/track_matcher.m`）
- 保持与 Phase 1-5 的变量名兼容
- 不要改动 `config/simulation_params.m`
