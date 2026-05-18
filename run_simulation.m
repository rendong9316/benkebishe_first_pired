% =========================================================================
% run_simulation.m
% 双基地外辐射源雷达逐帧仿真主程序
% =========================================================================
% 流程:
%   Phase 0: 场景初始化（参数、航迹）
%   Phase 1: 系统偏差离线标定（标校点 -> LS估计）
%   Phase 2: 逐帧主循环（点迹生成 -> 偏差校正 -> 跟踪处理）
%   Phase 3: 可视化 + 数据输出
% =========================================================================

clear; close all; clc;
addpath('config', 'utils', 'simulation', 'registration', 'ukf', 'fusion', 'visualization', 'io', 'tracker', 'evaluation');

%% ==================== Phase 0: 场景初始化 ====================
fprintf('========== Phase 0: 场景初始化 ==========\n');

params = simulation_params();
rng(params.random_seed);

% 多目标航迹配置
aircraft_labels = {'A', 'B', 'C'};
aircraft_wps = {params.aircraft_A_waypoints, params.aircraft_B_waypoints, ...
                params.aircraft_C_waypoints};
aircraft_spds = [params.aircraft_A_speed_ms, params.aircraft_B_speed_ms, ...
                 params.aircraft_C_speed_ms];

% 生成多目标真实航迹
trajs = cell(params.num_aircraft, 1);
true_tracks = cell(params.num_aircraft, 1);
for a = 1:params.num_aircraft
    trajs{a} = aircraft_trajectory_create(aircraft_wps{a}, aircraft_spds(a), params.dt_sec);
    true_tracks{a} = aircraft_trajectory_generate(trajs{a});
    fprintf('飞机%s航迹: %d 点, 总时长 %.0f s, 速度 %.0f m/s\n', ...
        aircraft_labels{a}, size(true_tracks{a},1), trajs{a}.duration_sec, aircraft_spds(a));
end

% 检查各航迹在双方威力范围内的覆盖
for a = 1:params.num_aircraft
    n_in_r1 = 0; n_in_r2 = 0;
    tt = true_tracks{a};
    for i = 1:size(tt, 1)
        [in1, ~, ~] = radar_coverage_check(params.radar1_lon, params.radar1_lat, ...
            tt(i,1), tt(i,2), params.radar1_beam_center_deg, params);
        [in2, ~, ~] = radar_coverage_check(params.radar2_lon, params.radar2_lat, ...
            tt(i,1), tt(i,2), params.radar2_beam_center_deg, params);
        if in1, n_in_r1 = n_in_r1 + 1; end
        if in2, n_in_r2 = n_in_r2 + 1; end
    end
    fprintf('  飞机%s: 在R1威力内 %d 点, 在R2威力内 %d 点 (共%d点)\n', ...
        aircraft_labels{a}, n_in_r1, n_in_r2, size(tt,1));
end

% 构建采样时间网格（以飞机A的航迹时长为基准，dt=30s对齐）
t1_grid = params.time_offset_radar1_sec : params.dt_sec : trajs{1}.duration_sec;
t2_grid = params.time_offset_radar2_sec : params.dt_sec : trajs{1}.duration_sec;
n_frames = min(length(t1_grid), length(t2_grid));
fprintf('仿真帧数: %d (dt=%.0fs)\n', n_frames, params.dt_sec);

%% ==================== Phase 1: 系统偏差离线标定 ====================
fprintf('\n========== Phase 1: 系统偏差标定 ==========\n');

rng(params.random_seed);  % 标定使用独立随机流

% 使用飞机A的真值航迹进行标定
true_track_A = true_tracks{1};

n_cal = min(30, size(true_track_A, 1));
cal_step = max(1, floor(size(true_track_A, 1) / n_cal));
cal_idxs = 1:cal_step:size(true_track_A, 1);
cal_idxs = cal_idxs(1:min(n_cal, length(cal_idxs)));

dr1_list = []; da1_list = [];
dr2_list = []; da2_list = [];

for idx = cal_idxs
    t_lon = true_track_A(idx,1);  t_lat = true_track_A(idx,2);

    [in1, ~, ~] = radar_coverage_check(params.radar1_lon, params.radar1_lat, ...
        t_lon, t_lat, params.radar1_beam_center_deg, params);
    if in1
        r0 = sphere_utils_haversine_distance(params.radar1_tx_lon, params.radar1_tx_lat, t_lon, t_lat);
        r1 = sphere_utils_haversine_distance(params.radar1_lon, params.radar1_lat, t_lon, t_lat);
        Rg_true = r0 + r1;
        az_true = sphere_utils_azimuth(params.radar1_lon, params.radar1_lat, t_lon, t_lat);
        Rg_meas = Rg_true + params.radar1_range_bias_m + randn() * params.range_noise_std_m;
        az_meas = az_true + params.radar1_azimuth_bias_deg + randn() * params.azimuth_noise_std_deg;
        dr1_list(end+1) = Rg_meas - Rg_true;
        daz = az_meas - az_true;
        if daz > 180, daz = daz - 360; elseif daz < -180, daz = daz + 360; end
        da1_list(end+1) = daz;
    end

    [in2, ~, ~] = radar_coverage_check(params.radar2_lon, params.radar2_lat, ...
        t_lon, t_lat, params.radar2_beam_center_deg, params);
    if in2
        r0 = sphere_utils_haversine_distance(params.radar2_tx_lon, params.radar2_tx_lat, t_lon, t_lat);
        r1 = sphere_utils_haversine_distance(params.radar2_lon, params.radar2_lat, t_lon, t_lat);
        Rg_true = r0 + r1;
        az_true = sphere_utils_azimuth(params.radar2_lon, params.radar2_lat, t_lon, t_lat);
        Rg_meas = Rg_true + params.radar2_range_bias_m + randn() * params.range_noise_std_m;
        az_meas = az_true + params.radar2_azimuth_bias_deg + randn() * params.azimuth_noise_std_deg;
        dr2_list(end+1) = Rg_meas - Rg_true;
        daz = az_meas - az_true;
        if daz > 180, daz = daz - 360; elseif daz < -180, daz = daz + 360; end
        da2_list(end+1) = daz;
    end
end

dr1_est = mean(dr1_list);  da1_est = mean(da1_list);
dr2_est = mean(dr2_list);  da2_est = mean(da2_list);
fprintf('标校点数: R1=%d, R2=%d\n', length(dr1_list), length(dr2_list));
fprintf('R1: dr_est=%.1f (true=%.0f) m, da_est=%.4f (true=%.1f) deg\n', ...
    dr1_est, params.radar1_range_bias_m, da1_est, params.radar1_azimuth_bias_deg);
fprintf('R2: dr_est=%.1f (true=%.0f) m, da_est=%.4f (true=%.1f) deg\n', ...
    dr2_est, params.radar2_range_bias_m, da2_est, params.radar2_azimuth_bias_deg);

%% ==================== Phase 2: 逐帧主循环 ====================
fprintf('\n========== Phase 2: 逐帧主循环 (多目标点迹生成 + 多航迹跟踪) ==========\n');

% 预分配点迹列表
detList_R1 = cell(n_frames, 1);
detList_R2 = cell(n_frames, 1);

% 多目标航迹管理器状态
trackList_R1 = {};  tempPool_R1 = {};
trackList_R2 = {};  tempPool_R2 = {};
trackSnapshots_R1 = cell(n_frames, 1);
trackSnapshots_R2 = cell(n_frames, 1);

% UKF模板（用于新航迹初始化）
ukf1_tpl = ukf_filter(params, params.radar1_lon, params.radar1_lat, ...
    params.radar1_tx_lon, params.radar1_tx_lat, params.dt_sec);
ukf2_tpl = ukf_filter(params, params.radar2_lon, params.radar2_lat, ...
    params.radar2_tx_lon, params.radar2_tx_lat, params.dt_sec);

% 统计每个飞机的检出情况
ac_det_counts_r1 = zeros(params.num_aircraft, 1);
ac_det_counts_r2 = zeros(params.num_aircraft, 1);

for k = 1:n_frames
    % ---- R1 多目标点迹生成 ----
    all_dets_r1 = [];
    for a = 1:params.num_aircraft
        [pos, vel] = aircraft_trajectory_interpolate(trajs{a}, t1_grid(k));
        add_clut = (a == 1);
        rng(params.random_seed + a*1000 + k);

        dets_a = generate_frame_detections(params.radar1_lon, params.radar1_lat, ...
            params.radar1_tx_lon, params.radar1_tx_lat, ...
            pos(1), pos(2), vel(1), vel(2), k, t1_grid(k), ...
            params.radar1_range_bias_m, params.radar1_azimuth_bias_deg, ...
            params.radar1_beam_center_deg, params, add_clut);

        for d = 1:length(dets_a)
            Rgc = dets_a(d).prange - dr1_est;
            azc = dets_a(d).paz - da1_est;
            dets_a(d).drange = Rgc;
            dets_a(d).daz = azc;
            dets_a(d).range_meas = Rgc;
            dets_a(d).azimuth_meas = azc;
            dets_a(d).aircraft_id = a;
            if ~(isfield(dets_a(d), 'lat') && ~isnan(dets_a(d).lat))
                [~, lat_e, lon_e] = bistatic_inverse_solver(Rgc, azc, ...
                    params.radar1_tx_lon, params.radar1_tx_lat, ...
                    params.radar1_lon, params.radar1_lat);
                dets_a(d).lat = lat_e;
                dets_a(d).lon = lon_e;
            end
            [~, raw_lat, raw_lon] = bistatic_inverse_solver(dets_a(d).prange, dets_a(d).paz, ...
                params.radar1_tx_lon, params.radar1_tx_lat, ...
                params.radar1_lon, params.radar1_lat);
            dets_a(d).raw_lat = raw_lat;
            dets_a(d).raw_lon = raw_lon;
            if ~dets_a(d).is_clutter
                ac_det_counts_r1(a) = ac_det_counts_r1(a) + 1;
            end
        end
        all_dets_r1 = [all_dets_r1, dets_a];
    end
    detList_R1{k} = all_dets_r1;

    % ---- R2 多目标点迹生成 ----
    all_dets_r2 = [];
    for a = 1:params.num_aircraft
        [pos, vel] = aircraft_trajectory_interpolate(trajs{a}, t2_grid(k));
        add_clut = (a == 1);
        rng(params.random_seed + 10000 + a*1000 + k);

        dets_a = generate_frame_detections(params.radar2_lon, params.radar2_lat, ...
            params.radar2_tx_lon, params.radar2_tx_lat, ...
            pos(1), pos(2), vel(1), vel(2), k, t2_grid(k), ...
            params.radar2_range_bias_m, params.radar2_azimuth_bias_deg, ...
            params.radar2_beam_center_deg, params, add_clut);

        for d = 1:length(dets_a)
            Rgc = dets_a(d).prange - dr2_est;
            azc = dets_a(d).paz - da2_est;
            dets_a(d).drange = Rgc;
            dets_a(d).daz = azc;
            dets_a(d).range_meas = Rgc;
            dets_a(d).azimuth_meas = azc;
            dets_a(d).aircraft_id = a;
            if ~(isfield(dets_a(d), 'lat') && ~isnan(dets_a(d).lat))
                [~, lat_e, lon_e] = bistatic_inverse_solver(Rgc, azc, ...
                    params.radar2_tx_lon, params.radar2_tx_lat, ...
                    params.radar2_lon, params.radar2_lat);
                dets_a(d).lat = lat_e;
                dets_a(d).lon = lon_e;
            end
            [~, raw_lat, raw_lon] = bistatic_inverse_solver(dets_a(d).prange, dets_a(d).paz, ...
                params.radar2_tx_lon, params.radar2_tx_lat, ...
                params.radar2_lon, params.radar2_lat);
            dets_a(d).raw_lat = raw_lat;
            dets_a(d).raw_lon = raw_lon;
            if ~dets_a(d).is_clutter
                ac_det_counts_r2(a) = ac_det_counts_r2(a) + 1;
            end
        end
        all_dets_r2 = [all_dets_r2, dets_a];
    end
    detList_R2{k} = all_dets_r2;

    % ---- 多目标航迹管理 ----
    [trackList_R1, tempPool_R1, trackSnapshots_R1{k}] = multi_track_manager(...
        trackList_R1, tempPool_R1, all_dets_r1, ukf1_tpl, params, k);
    [trackList_R2, tempPool_R2, trackSnapshots_R2{k}] = multi_track_manager(...
        trackList_R2, tempPool_R2, all_dets_r2, ukf2_tpl, params, k);
end

fprintf('仿真完成: %d 帧处理完毕\n', n_frames);
for a = 1:params.num_aircraft
    fprintf('  飞机%s: R1目标检出=%d, R2目标检出=%d\n', ...
        aircraft_labels{a}, ac_det_counts_r1(a), ac_det_counts_r2(a));
end

% 航迹统计
fprintf('\n--- 航迹统计 ---\n');
fprintf('R1: 共产生 %d 条航迹\n', length(trackList_R1));
for t = 1:length(trackList_R1)
    trk = trackList_R1{t};
    type_str = get_type_str(trk.type);
    fprintf('  R1航迹#%d: type=%s quality=%d life=%d\n', trk.id, type_str, trk.quality, trk.life);
end
fprintf('R2: 共产生 %d 条航迹\n', length(trackList_R2));
for t = 1:length(trackList_R2)
    trk = trackList_R2{t};
    type_str = get_type_str(trk.type);
    fprintf('  R2航迹#%d: type=%s quality=%d life=%d\n', trk.id, type_str, trk.quality, trk.life);
end

%% ==================== Phase 3: 统计 ====================
fprintf('\n========== Phase 3: 统计 ==========\n');

for a = 1:params.num_aircraft
    s_r1 = compute_stats_for_aircraft(detList_R1, a);
    s_r2 = compute_stats_for_aircraft(detList_R2, a);
    fprintf('飞机%s  R1: 目标检出=%d   R2: 目标检出=%d\n', ...
        aircraft_labels{a}, s_r1.target, s_r2.target);
end

% 总体统计
total_r1 = struct('detections', 0, 'target', 0, 'clutter', 0);
total_r2 = struct('detections', 0, 'target', 0, 'clutter', 0);
for k = 1:n_frames
    for d = 1:length(detList_R1{k})
        total_r1.detections = total_r1.detections + 1;
        if detList_R1{k}(d).is_clutter
            total_r1.clutter = total_r1.clutter + 1;
        else
            total_r1.target = total_r1.target + 1;
        end
    end
    for d = 1:length(detList_R2{k})
        total_r2.detections = total_r2.detections + 1;
        if detList_R2{k}(d).is_clutter
            total_r2.clutter = total_r2.clutter + 1;
        else
            total_r2.target = total_r2.target + 1;
        end
    end
end
fprintf('R1总计: 点迹=%d, 目标检出=%d, 杂波=%d\n', ...
    total_r1.detections, total_r1.target, total_r1.clutter);
fprintf('R2总计: 点迹=%d, 目标检出=%d, 杂波=%d\n', ...
    total_r2.detections, total_r2.target, total_r2.clutter);

%% ==================== Phase 3.5: 定量误差评估 ====================
fprintf('\n========== 定量误差评估 ==========\n');

% 构建真值结构体（供误差计算用）
truthTrajs = cell(params.num_aircraft, 1);
for a = 1:params.num_aircraft
    tt = true_tracks{a};
    truthTrajs{a} = struct('label', aircraft_labels{a}, ...
        'speed_ms', aircraft_spds(a), ...
        'time_sec', tt(:,5), 'lat', tt(:,2), 'lon', tt(:,1), ...
        'lon_rate', tt(:,3), 'lat_rate', tt(:,4));
end

errorStats_R1 = compute_tracking_errors(trackSnapshots_R1, detList_R1, ...
    truthTrajs, n_frames, params.dt_sec, 'R1');
errorStats_R2 = compute_tracking_errors(trackSnapshots_R2, detList_R2, ...
    truthTrajs, n_frames, params.dt_sec, 'R2');

% 控制台输出
for es = {errorStats_R1, errorStats_R2}
    e = es{1};
    fprintf('\n--- %s UKF滤波误差 ---\n', e.radar);
    fprintf('%-6s %6s %8s %8s %8s %8s %8s\n', ...
        '飞机', '点数', '中位(km)', '均值(km)', 'RMSE(km)', '95%(km)', 'vs检测');
    for a = 1:length(e.summary)
        s_u = e.summary(a).ukf;
        s_d = e.summary(a).det_calibrated;
        fprintf('飞机%d   %6d %8.1f %8.1f %8.1f %8.1f %7.0f%%\n', ...
            a, s_u.n, s_u.median, s_u.mean, s_u.rms, s_u.pct95, ...
            e.summary(a).ukf_vs_det_pct);
    end
    fprintf('%-6s %6d %8.1f %8.1f %8.1f %8.1f\n', ...
        '总计', e.overall.ukf.n, e.overall.ukf.median, ...
        e.overall.ukf.mean, e.overall.ukf.rms, e.overall.ukf.pct95);

    fprintf('\n--- %s 检测点迹误差(校准后) ---\n', e.radar);
    for a = 1:length(e.summary)
        s = e.summary(a).det_calibrated;
        fprintf('飞机%d   %6d %8.1f %8.1f %8.1f %8.1f\n', ...
            a, s.n, s.median, s.mean, s.rms, s.pct95);
    end
end

%% ==================== Phase 4: 可视化 ====================
fprintf('\n========== Phase 4: 可视化 ==========\n');
if ~exist('results', 'dir'), mkdir('results'); end

% 场景总览（多航迹）
plot_scene_overview_multi(true_tracks, aircraft_labels, params, 'results');
% 三维点迹图
plot_point_cloud_3d(detList_R1, 'R1', 'results/fig2a_R1_point_cloud.png');
plot_point_cloud_3d(detList_R2, 'R2', 'results/fig2b_R2_point_cloud.png');
% 多目标点迹分布
plot_multi_target_detections(true_tracks, aircraft_labels, detList_R1, detList_R2, ...
    params, 'results');
% 多目标跟踪综合图（真值 + 原始/校准点迹 + UKF滤波 + 图层复选框）
plot_multi_track_result(true_tracks, aircraft_labels, detList_R1, detList_R2, ...
    trackSnapshots_R1, trackSnapshots_R2, params, 'results');

% 定量误差分析图（RMSE时间线 + CDF + 汇总表）
plot_error_analysis(errorStats_R1, errorStats_R2, 'results');

%% ==================== Phase 5: 数据保存 ====================
fprintf('\n========== Phase 5: 数据保存 ==========\n');

% 系统参数（对标参考项目 sysPara + 仿真独有参数）
sysPara = struct(...
    'dt_sec', params.dt_sec, ...
    'n_frames', n_frames, ...
    'R1_lon', params.radar1_lon, 'R1_lat', params.radar1_lat, ...
    'R1_tx_lon', params.radar1_tx_lon, 'R1_tx_lat', params.radar1_tx_lat, ...
    'R1_beam_center_deg', params.radar1_beam_center_deg, ...
    'R1_range_bias_m', params.radar1_range_bias_m, ...
    'R1_azimuth_bias_deg', params.radar1_azimuth_bias_deg, ...
    'R2_lon', params.radar2_lon, 'R2_lat', params.radar2_lat, ...
    'R2_tx_lon', params.radar2_tx_lon, 'R2_tx_lat', params.radar2_tx_lat, ...
    'R2_beam_center_deg', params.radar2_beam_center_deg, ...
    'R2_range_bias_m', params.radar2_range_bias_m, ...
    'R2_azimuth_bias_deg', params.radar2_azimuth_bias_deg, ...
    'beam_width_deg', params.beam_width_deg, ...
    'range_km', [params.range_min_km, params.range_max_km], ...
    'detection_probability', params.detection_probability, ...
    'false_alarm_rate', params.false_alarm_rate, ...
    'range_noise_std_m', params.range_noise_std_m, ...
    'azimuth_noise_std_deg', params.azimuth_noise_std_deg, ...
    'radial_vel_noise_std_ms', params.radial_vel_noise_std_ms, ...
    'random_seed', params.random_seed);

% 偏差标定结果（对标参考项目标校模块输出）
calibResult = struct(...
    'dr1_est', dr1_est, 'da1_est', da1_est, ...
    'dr2_est', dr2_est, 'da2_est', da2_est, ...
    'dr1_true', params.radar1_range_bias_m, 'da1_true', params.radar1_azimuth_bias_deg, ...
    'dr2_true', params.radar2_range_bias_m, 'da2_true', params.radar2_azimuth_bias_deg);

% 各雷达数据（对标参考项目逐帧文件）
% R1
R1 = struct();
R1.detList = detList_R1;           % 点迹列表 (cell, 每帧struct array)
R1.trackSnapshots = trackSnapshots_R1;  % 逐帧航迹快照
R1.finalTrackList = trackList_R1;  % 最终航迹列表 (cell, 含全部航迹)
R1.tempTrackList = tempPool_R1;    % 临时航迹池 (对标参考项目tempTrackList)
R1.targetDetCounts = ac_det_counts_r1;  % 各飞机目标检出数
R1.totalDetections = total_r1.detections;
R1.totalClutter = total_r1.clutter;

% R2
R2 = struct();
R2.detList = detList_R2;
R2.trackSnapshots = trackSnapshots_R2;
R2.finalTrackList = trackList_R2;
R2.tempTrackList = tempPool_R2;
R2.targetDetCounts = ac_det_counts_r2;
R2.totalDetections = total_r2.detections;
R2.totalClutter = total_r2.clutter;

outf = fullfile('results', sprintf('simulation_multi_%s.mat', datestr(now, 'yyyymmdd_HHMMSS')));
save(outf, 'sysPara', 'calibResult', 'truthTrajs', 'R1', 'R2', 'params', ...
    'errorStats_R1', 'errorStats_R2');
fprintf('数据已保存: %s\n', outf);
fprintf('\nDone.\n');

% =========================================================================
% 内部函数
% =========================================================================

function s = compute_stats_for_aircraft(detList, aircraft_id)
    s.target = 0;
    for k = 1:length(detList)
        dets = detList{k};
        for d = 1:length(dets)
            if ~dets(d).is_clutter && isfield(dets(d), 'aircraft_id') ...
                    && dets(d).aircraft_id == aircraft_id
                s.target = s.target + 1;
            end
        end
    end
end

function s = get_type_str(t)
    switch t
        case 1, s = 'RELIABLE';
        case 2, s = 'MAINTAIN';
        case 6, s = 'TEMPORARY';
        case 7, s = 'HISTORY';
        otherwise, s = sprintf('UNKNOWN(%d)', t);
    end
end

function t = init_track_state()
    t = struct('status', 'UNINITIATED', 'x', zeros(4,1), 'P', eye(4), ...
        'ukf', [], 'init_window', {{}}, 'init_dets', {{}}, ...
        'first_det', [], 'life', 0, 'missed', 0, 'quality', 0, ...
        'nis_history', []);
end

function [t, snap] = process_one_frame(t, detList, ukf_tpl, params, k, beam_center)
    snap = struct('frameID', k, 'status', t.status, 'lat', NaN, 'lon', NaN, ...
        'det_lat', NaN, 'det_lon', NaN, 'det_raw_lat', NaN, 'det_raw_lon', NaN, ...
        'associated', false, 'assc_is_clutter', false);

    switch t.status
        case {'UNINITIATED', 'INITIATING'}
            % M/N起始：收集每帧所有点迹，在触发时做多假设配对
            if ~isempty(detList)
                t.init_window{end+1} = 1;
                t.init_dets{end+1} = detList;
            else
                t.init_window{end+1} = 0;
                t.init_dets{end+1} = [];
            end
            if length(t.init_window) > params.tracker_N
                t.init_window(1) = [];
                t.init_dets(1) = [];
            end
            n_det = sum(cell2mat(t.init_window));
            if n_det >= params.tracker_M && ~isempty(detList)
                % 多假设共识配对：遍历窗内所有帧-帧检测对
                det_now = detList(1);  % 当前帧首检测

                best_prev = [];
                best_support = -1;

                % 遍历之前各帧的全部点迹，做配对（不检查速度——位置噪声太大时速度不可靠）
                for i = 1:(length(t.init_dets)-1)
                    prev_dets = t.init_dets{i};
                    if isempty(prev_dets), continue; end
                    for p = 1:length(prev_dets)
                        dp = prev_dets(p);
                        if ~isfield(dp, 'lat') || isnan(dp.lat), continue; end
                        if ~isfield(det_now, 'lat') || isnan(det_now.lat), continue; end

                        % 共识评分：统计窗内其他帧有多少点迹同时靠近dp和det_now
                        support = 0;
                        for jj = 1:(length(t.init_dets)-1)
                            if jj == i, continue; end
                            other2 = t.init_dets{jj};
                            if isempty(other2), continue; end
                            for oo = 1:length(other2)
                                do = other2(oo);
                                if ~isfield(do, 'lat') || isnan(do.lat), continue; end
                                d1 = sphere_utils_haversine_distance(dp.lon, dp.lat, do.lon, do.lat);
                                d2 = sphere_utils_haversine_distance(det_now.lon, det_now.lat, do.lon, do.lat);
                                if d1 < 80000 && d2 < 80000
                                    support = support + 1;
                                end
                            end
                        end
                        if support > best_support
                            best_support = support;
                            best_prev = dp;
                        end
                    end
                end

                % 仅当存在共识配对时才起始（至少1个其他帧的点迹同时靠近配对两点）
                if best_support >= 1
                    t.status = 'TRACKING';
                    t.life = 0;  t.missed = 0;  t.quality = 0;
                    t.nis_history = [];
                    t.ukf = ukf_filter_init(ukf_tpl, best_prev, det_now);
                    t.x = t.ukf.x;  t.P = t.ukf.P;
                    snap.lat = t.x(3);  snap.lon = t.x(1);
                    snap.det_lat = det_now.lat;  snap.det_lon = det_now.lon;
                    if isfield(det_now, 'raw_lat')
                        snap.det_raw_lat = det_now.raw_lat;
                        snap.det_raw_lon = det_now.raw_lon;
                    end
                    snap.associated = true;
                    snap.assc_is_clutter = det_now.is_clutter;
                end
                % 无共识配对：暂不起始，继续收集下一帧
            end
            snap.status = t.status;

        case 'TRACKING'
            t.ukf.dt = params.dt_sec;
            [x_pred, P_pred, X_pred, t.ukf] = ukf_predict_step(t.ukf);

            % 计算预测量测及完整新息协方差 P_zz（含状态不确定性投影）
            z_pred = ukf_measurement_model(t.ukf, x_pred);
            Z_pred = zeros(t.ukf.m, 2*t.ukf.n + 1);
            for i = 1:(2*t.ukf.n + 1)
                Z_pred(:,i) = ukf_measurement_model(t.ukf, X_pred(:,i));
            end
            P_zz = t.ukf.R;
            for i = 1:(2*t.ukf.n + 1)
                dz = Z_pred(:,i) - z_pred;
                P_zz = P_zz + t.ukf.Wc(i) * (dz * dz');
            end
            P_zz_2d = P_zz(1:2, 1:2);

            % 数值稳定性守卫
            if any(isnan(P_zz_2d(:))) || any(isnan(z_pred))
                dets_in_gate = {};
            else
                % 航迹确认期：地理距离预筛选
                use_geo_gate = (t.life <= 15);
                gate_threshold = params.gate_sigma^2 * 2;
                dets_in_gate = {};
                for d = 1:length(detList)
                    dp = detList(d);
                    if ~isfield(dp, 'lat') || isnan(dp.lat), continue; end

                    if use_geo_gate
                        geo_err = sphere_utils_haversine_distance(...
                            x_pred(1), x_pred(3), dp.lon, dp.lat);
                        if geo_err > 80000, continue; end
                    end

                    z_meas_2d = [dp.drange; dp.daz];
                    innov_2d = z_meas_2d - z_pred(1:2);
                    if innov_2d(2) > 180, innov_2d(2) = innov_2d(2) - 360;
                    elseif innov_2d(2) < -180, innov_2d(2) = innov_2d(2) + 360; end
                    mahal = innov_2d' * (P_zz_2d \ innov_2d);
                    if mahal < gate_threshold
                        dets_in_gate{end+1} = dp;
                    end
                end
            end

            if ~isempty(dets_in_gate)
                % ---- 栅格加权PDA 或 最近邻UKF更新 ----
                if params.use_pda_weighting && length(dets_in_gate) >= 2
                    % 多量测PDA加权更新
                    [~, ~, t.ukf, best, nis_used] = ukf_pda_update(...
                        t.ukf, dets_in_gate, z_pred, Z_pred, ...
                        x_pred, P_pred, P_zz, params);
                    t.nis_history(end+1) = nis_used;
                else
                    % 单量测标准UKF更新
                    [~, ~, t.ukf] = ukf_filter_update(t.ukf, dets_in_gate{1});
                    best = dets_in_gate{1};
                    % 计算该量测的NIS
                    innov_best = [best.drange; best.daz] - z_pred(1:2);
                    if innov_best(2) > 180, innov_best(2) = innov_best(2) - 360;
                    elseif innov_best(2) < -180, innov_best(2) = innov_best(2) + 360; end
                    nis_val = innov_best' * (P_zz_2d \ innov_best);
                    t.nis_history(end+1) = nis_val;
                end

                t.x = t.ukf.x;  t.P = t.ukf.P;
                t.missed = 0;  t.life = t.life + 1;
                t.quality = min(t.quality + 1, 5);
                snap.lat = t.x(3);  snap.lon = t.x(1);
                snap.det_lat = best.lat;  snap.det_lon = best.lon;
                if isfield(best, 'raw_lat')
                    snap.det_raw_lat = best.raw_lat;
                    snap.det_raw_lon = best.raw_lon;
                end
                snap.associated = true;
                snap.assc_is_clutter = best.is_clutter;

                % ---- 模糊自适应Q调节 ----
                if params.use_fuzzy_adaptive
                    if length(t.nis_history) > params.fuzzy_window_size
                        t.nis_history(1) = [];
                    end
                    t.ukf = ukf_fuzzy_adapt(t.ukf, ...
                        t.nis_history, t.life, params);
                end
            else
                % 波门内无点迹：纯预测
                t.ukf.x = x_pred;  t.ukf.P = P_pred;
                t.missed = t.missed + 1;  t.life = t.life + 1;
                t.quality = max(t.quality - 1, 0);
                snap.lat = x_pred(3);  snap.lon = x_pred(1);
                snap.associated = false;
                if t.missed >= params.tracker_K_loss
                    t.status = 'LOST';
                end
            end
            % 航迹质量确认
            if t.life <= 15 && t.quality <= 0 && t.missed >= 2
                t.status = 'LOST';
            end
            snap.status = t.status;

        case 'LOST'
            % 航迹终止后重新尝试起始
            t.status = 'UNINITIATED';
            t.init_window = {};
            t.init_dets = {};
            t.first_det = [];
            t.quality = 0;
            snap.status = 'LOST';
    end
end

function n = count_valid_frames(stateList)
    n = 0;
    for k = 1:length(stateList)
        s = stateList{k};
        if ~isempty(s) && isfield(s, 'lat') && ~isnan(s.lat)
            n = n + 1;
        end
    end
end

function s = compute_stats(detList, true_track, t_grid, params, tag)
    s.total = 0;  s.target = 0;  s.false = 0;
    for k = 1:length(detList)
        dets = detList{k};
        for d = 1:length(dets)
            s.total = s.total + 1;
            if dets(d).is_clutter
                s.false = s.false + 1;
            else
                s.target = s.target + 1;
            end
        end
    end
end
