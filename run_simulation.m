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
addpath('config', 'utils', 'simulation', 'registration', 'ukf', 'fusion', 'visualization', 'io');

%% ==================== Phase 0: 场景初始化 ====================
fprintf('========== Phase 0: 场景初始化 ==========\n');

params = simulation_params();
rng(params.random_seed);

% 生成真实航迹
traj = aircraft_trajectory_create(params.aircraft_waypoints, ...
    params.aircraft_speed_ms, params.dt_sec);
true_track = aircraft_trajectory_generate(traj);
fprintf('真实航迹: %d 点, 总时长 %.0f s\n', size(true_track, 1), traj.duration_sec);

% 检查航迹在双方威力范围内的覆盖
n_in_r1 = 0; n_in_r2 = 0;
for i = 1:size(true_track, 1)
    [in1, ~, ~] = radar_coverage_check(params.radar1_lon, params.radar1_lat, ...
        true_track(i,1), true_track(i,2), params.radar1_beam_center_deg, params);
    [in2, ~, ~] = radar_coverage_check(params.radar2_lon, params.radar2_lat, ...
        true_track(i,1), true_track(i,2), params.radar2_beam_center_deg, params);
    if in1, n_in_r1 = n_in_r1 + 1; end
    if in2, n_in_r2 = n_in_r2 + 1; end
end
fprintf('  在R1威力内: %d 点, 在R2威力内: %d 点 (共%d点)\n', ...
    n_in_r1, n_in_r2, size(true_track,1));

% 构建采样时间网格
t1_grid = params.time_offset_radar1_sec : params.dt_sec : traj.duration_sec;
t2_grid = params.time_offset_radar2_sec : params.dt_sec : traj.duration_sec;
n_frames = min(length(t1_grid), length(t2_grid));
fprintf('仿真帧数: %d (dt=%.0fs)\n', n_frames, params.dt_sec);

%% ==================== Phase 1: 系统偏差离线标定 ====================
fprintf('\n========== Phase 1: 系统偏差标定 ==========\n');

rng(params.random_seed);  % 标定使用独立随机流

n_cal = min(30, size(true_track, 1));
cal_step = max(1, floor(size(true_track, 1) / n_cal));
cal_idxs = 1:cal_step:size(true_track, 1);
cal_idxs = cal_idxs(1:min(n_cal, length(cal_idxs)));

dr1_list = []; da1_list = [];
dr2_list = []; da2_list = [];

for idx = cal_idxs
    t_lon = true_track(idx,1);  t_lat = true_track(idx,2);

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
fprintf('\n========== Phase 2: 逐帧主循环 ==========\n');

% 预分配
detList_R1 = cell(n_frames, 1);
detList_R2 = cell(n_frames, 1);
trackState_R1 = cell(n_frames, 1);
trackState_R2 = cell(n_frames, 1);

% UKF 模板
ukf1_tpl = ukf_filter(params, params.radar1_lon, params.radar1_lat, ...
    params.radar1_tx_lon, params.radar1_tx_lat, params.dt_sec);
ukf2_tpl = ukf_filter(params, params.radar2_lon, params.radar2_lat, ...
    params.radar2_tx_lon, params.radar2_tx_lat, params.dt_sec);

% 航迹状态初始化
trk1 = init_track_state();
trk2 = init_track_state();

for k = 1:n_frames
    % R1随机流
    rng(params.random_seed + k);

    % 目标真实位置
    [pos, vel] = aircraft_trajectory_interpolate(traj, t1_grid(k));

    % ---- R1 点迹生成 + 偏差校正 ----
    det_r1 = generate_frame_detections(params.radar1_lon, params.radar1_lat, ...
        params.radar1_tx_lon, params.radar1_tx_lat, ...
        pos(1), pos(2), vel(1), vel(2), k, t1_grid(k), ...
        params.radar1_range_bias_m, params.radar1_azimuth_bias_deg, ...
        params.radar1_beam_center_deg, params);
    for d = 1:length(det_r1)
        Rgc = det_r1(d).prange - dr1_est;
        azc = det_r1(d).paz - da1_est;
        det_r1(d).drange = Rgc;
        det_r1(d).daz = azc;
        det_r1(d).range_meas = Rgc;
        det_r1(d).azimuth_meas = azc;
        % 杂波已在(r1,az)空间预设计算好地理坐标，跳过反解
        if isfield(det_r1(d), 'lat') && ~isnan(det_r1(d).lat)
            % lat/lon already set
        else
            [~, lat_e, lon_e] = bistatic_inverse_solver(Rgc, azc, ...
                params.radar1_tx_lon, params.radar1_tx_lat, ...
                params.radar1_lon, params.radar1_lat);
            det_r1(d).lat = lat_e;
            det_r1(d).lon = lon_e;
        end
        % 计算原始（校准前）地理位置，用于对比显示偏差配准效果
        [~, raw_lat, raw_lon] = bistatic_inverse_solver(det_r1(d).prange, det_r1(d).paz, ...
            params.radar1_tx_lon, params.radar1_tx_lat, ...
            params.radar1_lon, params.radar1_lat);
        det_r1(d).raw_lat = raw_lat;
        det_r1(d).raw_lon = raw_lon;
    end
    detList_R1{k} = det_r1;

    % ---- R1 跟踪 ----
    [trk1, trackState_R1{k}] = process_one_frame(trk1, det_r1, ukf1_tpl, params, k, ...
        params.radar1_beam_center_deg);

    % ---- R2 点迹生成 + 偏差校正 ----
    rng(params.random_seed + 10000 + k);
    [pos2, vel2] = aircraft_trajectory_interpolate(traj, t2_grid(k));

    det_r2 = generate_frame_detections(params.radar2_lon, params.radar2_lat, ...
        params.radar2_tx_lon, params.radar2_tx_lat, ...
        pos2(1), pos2(2), vel2(1), vel2(2), k, t2_grid(k), ...
        params.radar2_range_bias_m, params.radar2_azimuth_bias_deg, ...
        params.radar2_beam_center_deg, params);
    for d = 1:length(det_r2)
        Rgc = det_r2(d).prange - dr2_est;
        azc = det_r2(d).paz - da2_est;
        det_r2(d).drange = Rgc;
        det_r2(d).daz = azc;
        det_r2(d).range_meas = Rgc;
        det_r2(d).azimuth_meas = azc;
        % 杂波已在(r1,az)空间预设计算好地理坐标，跳过反解
        if isfield(det_r2(d), 'lat') && ~isnan(det_r2(d).lat)
            % lat/lon already set
        else
            [~, lat_e, lon_e] = bistatic_inverse_solver(Rgc, azc, ...
                params.radar2_tx_lon, params.radar2_tx_lat, ...
                params.radar2_lon, params.radar2_lat);
            det_r2(d).lat = lat_e;
            det_r2(d).lon = lon_e;
        end
        % 计算原始（校准前）地理位置，用于对比显示偏差配准效果
        [~, raw_lat, raw_lon] = bistatic_inverse_solver(det_r2(d).prange, det_r2(d).paz, ...
            params.radar2_tx_lon, params.radar2_tx_lat, ...
            params.radar2_lon, params.radar2_lat);
        det_r2(d).raw_lat = raw_lat;
        det_r2(d).raw_lon = raw_lon;
    end
    detList_R2{k} = det_r2;

    % ---- R2 跟踪 ----
    [trk2, trackState_R2{k}] = process_one_frame(trk2, det_r2, ukf2_tpl, params, k, ...
        params.radar2_beam_center_deg);
end

fprintf('仿真完成: %d 帧处理完毕\n', n_frames);

%% ==================== Phase 3: 统计 ====================
fprintf('\n========== Phase 3: 统计 ==========\n');
% 统计
s1 = compute_stats(detList_R1, true_track, t1_grid, params, 'R1');
s2 = compute_stats(detList_R2, true_track, t2_grid, params, 'R2');

% 滤波航迹有效帧数
n_filt1 = count_valid_frames(trackState_R1);
n_filt2 = count_valid_frames(trackState_R2);
fprintf('R1: 总点迹=%d, 目标检出=%d, 虚警=%d, 滤波有效帧=%d\n', ...
    s1.total, s1.target, s1.false, n_filt1);
fprintf('R2: 总点迹=%d, 目标检出=%d, 虚警=%d, 滤波有效帧=%d\n', ...
    s2.total, s2.target, s2.false, n_filt2);

%% ==================== Phase 4: 可视化 ====================
fprintf('\n========== Phase 4: 可视化 ==========\n');
if ~exist('results', 'dir'), mkdir('results'); end

plot_scene_overview(true_track, params, 'results');
plot_point_cloud_3d(detList_R1, 'R1', 'results/fig2a_R1_point_cloud.png');
plot_point_cloud_3d(detList_R2, 'R2', 'results/fig2b_R2_point_cloud.png');
plot_combined_tracks(true_track, detList_R1, detList_R2, ...
    trackState_R1, trackState_R2, params, 'results');
plot_error_timeline(trackState_R1, trackState_R2, detList_R1, detList_R2, ...
    true_track, t1_grid, t2_grid, params, 'results');

%% ==================== Phase 5: 数据保存 ====================
fprintf('\n========== Phase 5: 数据保存 ==========\n');

calib = struct('dr1_est', dr1_est, 'da1_est', da1_est, ...
    'dr2_est', dr2_est, 'da2_est', da2_est, ...
    'dr1_true', params.radar1_range_bias_m, 'da1_true', params.radar1_azimuth_bias_deg, ...
    'dr2_true', params.radar2_range_bias_m, 'da2_true', params.radar2_azimuth_bias_deg);
scene = struct('R1_lon', params.radar1_lon, 'R1_lat', params.radar1_lat, ...
    'R2_lon', params.radar2_lon, 'R2_lat', params.radar2_lat, ...
    'Tx1_lon', params.radar1_tx_lon, 'Tx1_lat', params.radar1_tx_lat, ...
    'Tx2_lon', params.radar2_tx_lon, 'Tx2_lat', params.radar2_tx_lat);
truth = struct('time_sec', true_track(:,5), 'lat', true_track(:,2), ...
    'lon', true_track(:,1), 'lon_rate', true_track(:,3), 'lat_rate', true_track(:,4));

outf = fullfile('results', sprintf('simulation_%s.mat', datestr(now, 'yyyymmdd_HHMMSS')));
save(outf, 'detList_R1', 'detList_R2', 'trackState_R1', 'trackState_R2', ...
    'calib', 'scene', 'truth', 'params', 's1', 's2');
fprintf('数据已保存: %s\n', outf);
fprintf('\nDone.\n');

% =========================================================================
% 内部函数
% =========================================================================

function t = init_track_state()
    t = struct('status', 'UNINITIATED', 'x', zeros(4,1), 'P', eye(4), ...
        'ukf', [], 'init_window', {{}}, 'init_dets', {{}}, ...
        'first_det', [], 'life', 0, 'missed', 0, 'quality', 0);
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
            P_zz_2d = P_zz(1:2, 1:2);  % 只用距离+方位做波门

            % 数值稳定性守卫：若 P_zz 含 NaN（Sigma点量测异常），退化为纯预测
            if any(isnan(P_zz_2d(:))) || any(isnan(z_pred))
                best = [];
            else
                % 航迹确认期：用地理距离预筛选（UKF不确定性大时P_zz门过宽）
                use_geo_gate = (t.life <= 15);
                best = [];  best_dist = inf;
                for d = 1:length(detList)
                    dp = detList(d);
                    if ~isfield(dp, 'lat') || isnan(dp.lat), continue; end

                    % 确认期地理距离预筛选：目标点迹必须在预测位置80 km内
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
                    if mahal < params.gate_sigma^2 * 2 && mahal < best_dist
                        best_dist = mahal;  best = dp;
                    end
                end
            end

            if ~isempty(best)
                [~, ~, t.ukf] = ukf_filter_update(t.ukf, best);
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
            % 航迹质量确认：起始后15帧内若质量归零且连续2帧无关联，判定起始失败
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
