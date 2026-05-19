% =========================================================================
% run_simulation.m
% 双基地外辐射源雷达仿真主程序
% =========================================================================
% Phase 0: 场景初始化（航迹生成 + 覆盖检查）
% Phase 1: ADS-B系统偏差标定
% Phase 2: 两部雷达各自产生原始点迹（含偏差，不做校正）
% Phase 3: 时间对齐（R2点迹→R1时间网格）
% Phase 4: 偏差校正（几何反解）
% Phase 5: 多目标航迹起始 + UKF滤波跟踪
% Phase 6: 定量误差评估
% Phase 7: 可视化 + 数据保存
% =========================================================================

clear; close all; clc;
addpath(genpath('.'));

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

trajs = cell(params.num_aircraft, 1);
true_tracks = cell(params.num_aircraft, 1);
for a = 1:params.num_aircraft
    trajs{a} = aircraft_trajectory_create(aircraft_wps{a}, aircraft_spds(a), params.dt_sec);
    true_tracks{a} = aircraft_trajectory_generate(trajs{a});
    fprintf('飞机%s航迹: %d 点, 总时长 %.0f s, 速度 %.0f m/s\n', ...
        aircraft_labels{a}, size(true_tracks{a},1), trajs{a}.duration_sec, aircraft_spds(a));
end

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

t1_grid = params.time_offset_radar1_sec : params.dt_sec : trajs{1}.duration_sec;
t2_grid = params.time_offset_radar2_sec : params.dt_sec : trajs{1}.duration_sec;
n_frames = min(length(t1_grid), length(t2_grid));
fprintf('仿真帧数: %d (dt=%.0fs)\n', n_frames, params.dt_sec);

%% ==================== Phase 1: ADS-B系统偏差标定 ====================
fprintf('\n========== Phase 1: ADS-B系统偏差标定 ==========\n');

rng(params.random_seed);

fprintf('加载ADS-B合作目标: %s\n', params.adsb_csv_path);
T_adsb = readtable(params.adsb_csv_path, 'ReadVariableNames', false);
adsb_lat = T_adsb.Var2;
adsb_lon = T_adsb.Var3;

dr1_list = []; da1_list = [];
dr2_list = []; da2_list = [];

n_check = min(5000, height(T_adsb));
cal_step = max(1, floor(height(T_adsb) / n_check));

for idx = 1:cal_step:height(T_adsb)
    t_lon = adsb_lon(idx);  t_lat = adsb_lat(idx);
    if isnan(t_lon) || isnan(t_lat), continue; end

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
fprintf('ADS-B标校点数: R1=%d, R2=%d\n', length(dr1_list), length(dr2_list));
fprintf('R1: dr_est=%.1f (true=%.0f) m, da_est=%.4f (true=%.1f) deg\n', ...
    dr1_est, params.radar1_range_bias_m, da1_est, params.radar1_azimuth_bias_deg);
fprintf('R2: dr_est=%.1f (true=%.0f) m, da_est=%.4f (true=%.1f) deg\n', ...
    dr2_est, params.radar2_range_bias_m, da2_est, params.radar2_azimuth_bias_deg);

%% ==================== Phase 2: 原始点迹生成（含偏差，不校正） ====================
fprintf('\n========== Phase 2: 原始点迹生成 ==========\n');

detRaw_R1 = cell(n_frames, 1);
detRaw_R2 = cell(n_frames, 1);

for k = 1:n_frames
    % R1
    all_raw_r1 = [];
    for a = 1:params.num_aircraft
        [pos, vel] = aircraft_trajectory_interpolate(trajs{a}, t1_grid(k));
        add_clut = (a == 1);
        rng(params.random_seed + a*1000 + k);

        dets = generate_frame_detections(params.radar1_lon, params.radar1_lat, ...
            params.radar1_tx_lon, params.radar1_tx_lat, ...
            pos(1), pos(2), vel(1), vel(2), k, t1_grid(k), ...
            params.radar1_range_bias_m, params.radar1_azimuth_bias_deg, ...
            params.radar1_beam_center_deg, params, add_clut);
        for d = 1:length(dets)
            dets(d).aircraft_id = a;
        end
        all_raw_r1 = [all_raw_r1, dets];
    end
    detRaw_R1{k} = all_raw_r1;

    % R2
    all_raw_r2 = [];
    for a = 1:params.num_aircraft
        [pos, vel] = aircraft_trajectory_interpolate(trajs{a}, t2_grid(k));
        add_clut = (a == 1);
        rng(params.random_seed + 10000 + a*1000 + k);

        dets = generate_frame_detections(params.radar2_lon, params.radar2_lat, ...
            params.radar2_tx_lon, params.radar2_tx_lat, ...
            pos(1), pos(2), vel(1), vel(2), k, t2_grid(k), ...
            params.radar2_range_bias_m, params.radar2_azimuth_bias_deg, ...
            params.radar2_beam_center_deg, params, add_clut);
        for d = 1:length(dets)
            dets(d).aircraft_id = a;
        end
        all_raw_r2 = [all_raw_r2, dets];
    end
    detRaw_R2{k} = all_raw_r2;
end

fprintf('原始点迹生成完成: R1共%d帧, R2共%d帧\n', n_frames, n_frames);

%% ==================== Phase 3: 时间对齐（一阶外推） ====================
fprintf('\n========== Phase 3: 时间对齐（一阶外推） ==========\n');
fprintf('R1采样: 0s/30s/60s/...  R2采样: 13s/43s/73s/...  偏移=%ds\n', ...
    params.time_offset_radar2_sec);
fprintf('算法: 利用双基地径向速度 vd 外推群距离 Rg(t)=Rg(t+dt)-vd*dt\n');

dt_align = params.time_offset_radar2_sec;  % 13s

detAligned_R2 = cell(n_frames, 1);
for k = 1:n_frames
    dets = detRaw_R2{k};
    for d = 1:length(dets)
        % 一阶外推：利用径向速度修正群距离
        % Rg(t) ≈ Rg(t+dt) - pvr * dt
        dets(d).prange = dets(d).prange - dets(d).pvr * dt_align;
        dets(d).time_sec = dets(d).time_sec - dt_align;
    end
    detAligned_R2{k} = dets;
end
fprintf('R2点迹时间对齐+距离外推完成 (%d帧)\n', n_frames);

%% ==================== Phase 4: 偏差校正 + 几何反解 ====================
fprintf('\n========== Phase 4: 偏差校正 ==========\n');

detList_R1 = cell(n_frames, 1);
detList_R2 = cell(n_frames, 1);

for k = 1:n_frames
    % R1: 偏差校正
    dets_r1 = detRaw_R1{k};
    for d = 1:length(dets_r1)
        Rgc = dets_r1(d).prange - dr1_est;
        azc = dets_r1(d).paz - da1_est;
        dets_r1(d).drange = Rgc;
        dets_r1(d).daz = azc;
        dets_r1(d).range_meas = Rgc;
        dets_r1(d).azimuth_meas = azc;

        if ~(isfield(dets_r1(d), 'lat') && ~isnan(dets_r1(d).lat))
            [~, lat_e, lon_e] = bistatic_inverse_solver(Rgc, azc, ...
                params.radar1_tx_lon, params.radar1_tx_lat, ...
                params.radar1_lon, params.radar1_lat);
            dets_r1(d).lat = lat_e;
            dets_r1(d).lon = lon_e;
        end
        [~, raw_lat, raw_lon] = bistatic_inverse_solver(dets_r1(d).prange, dets_r1(d).paz, ...
            params.radar1_tx_lon, params.radar1_tx_lat, ...
            params.radar1_lon, params.radar1_lat);
        dets_r1(d).raw_lat = raw_lat;
        dets_r1(d).raw_lon = raw_lon;
    end
    detList_R1{k} = dets_r1;

    % R2: 偏差校正（使用时间对齐后的点迹）
    dets_r2 = detAligned_R2{k};
    for d = 1:length(dets_r2)
        Rgc = dets_r2(d).prange - dr2_est;
        azc = dets_r2(d).paz - da2_est;
        dets_r2(d).drange = Rgc;
        dets_r2(d).daz = azc;
        dets_r2(d).range_meas = Rgc;
        dets_r2(d).azimuth_meas = azc;

        if ~(isfield(dets_r2(d), 'lat') && ~isnan(dets_r2(d).lat))
            [~, lat_e, lon_e] = bistatic_inverse_solver(Rgc, azc, ...
                params.radar2_tx_lon, params.radar2_tx_lat, ...
                params.radar2_lon, params.radar2_lat);
            dets_r2(d).lat = lat_e;
            dets_r2(d).lon = lon_e;
        end
        [~, raw_lat, raw_lon] = bistatic_inverse_solver(dets_r2(d).prange, dets_r2(d).paz, ...
            params.radar2_tx_lon, params.radar2_tx_lat, ...
            params.radar2_lon, params.radar2_lat);
        dets_r2(d).raw_lat = raw_lat;
        dets_r2(d).raw_lon = raw_lon;
    end
    detList_R2{k} = dets_r2;
end

fprintf('偏差校正完成: R1=%d帧, R2=%d帧\n', n_frames, n_frames);

%% ==================== Phase 5: 多目标航迹跟踪 ====================
fprintf('\n========== Phase 5: 多目标航迹跟踪 ==========\n');

ukf1_tpl = ukf_filter(params, params.radar1_lon, params.radar1_lat, ...
    params.radar1_tx_lon, params.radar1_tx_lat, params.dt_sec);
ukf2_tpl = ukf_filter(params, params.radar2_lon, params.radar2_lat, ...
    params.radar2_tx_lon, params.radar2_tx_lat, params.dt_sec);

trackList_R1 = {};  tempPool_R1 = {};
trackList_R2 = {};  tempPool_R2 = {};
trackSnapshots_R1 = cell(n_frames, 1);
trackSnapshots_R2 = cell(n_frames, 1);

ac_det_counts_r1 = zeros(params.num_aircraft, 1);
ac_det_counts_r2 = zeros(params.num_aircraft, 1);

for k = 1:n_frames
    % 检出统计
    for d = 1:length(detList_R1{k})
        dp = detList_R1{k}(d);
        if ~dp.is_clutter && isfield(dp,'aircraft_id')
            ac_det_counts_r1(dp.aircraft_id) = ac_det_counts_r1(dp.aircraft_id) + 1;
        end
    end
    for d = 1:length(detList_R2{k})
        dp = detList_R2{k}(d);
        if ~dp.is_clutter && isfield(dp,'aircraft_id')
            ac_det_counts_r2(dp.aircraft_id) = ac_det_counts_r2(dp.aircraft_id) + 1;
        end
    end

    % 航迹管理
    [trackList_R1, tempPool_R1, trackSnapshots_R1{k}] = multi_track_manager(...
        trackList_R1, tempPool_R1, detList_R1{k}, ukf1_tpl, params, k);
    [trackList_R2, tempPool_R2, trackSnapshots_R2{k}] = multi_track_manager(...
        trackList_R2, tempPool_R2, detList_R2{k}, ukf2_tpl, params, k);
end

fprintf('跟踪完成: %d 帧\n', n_frames);
for a = 1:params.num_aircraft
    fprintf('  飞机%s: R1目标检出=%d, R2目标检出=%d\n', ...
        aircraft_labels{a}, ac_det_counts_r1(a), ac_det_counts_r2(a));
end

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

%% ==================== Phase 6: 定量误差评估 ====================
fprintf('\n========== Phase 6: 定量误差评估 ==========\n');

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

for es = {errorStats_R1, errorStats_R2}
    e = es{1};
    fprintf('\n--- %s UKF滤波误差 ---\n', e.radar);
    fprintf('%-6s %6s %8s %8s %8s %8s %8s\n', ...
        '飞机', '点数', '中位(km)', '均值(km)', 'RMSE(km)', '95%(km)', 'vs检测');
    for a = 1:length(e.summary)
        s_u = e.summary(a).ukf;
        fprintf('飞机%d   %6d %8.1f %8.1f %8.1f %8.1f %7.0f%%\n', ...
            a, s_u.n, s_u.median, s_u.mean, s_u.rms, s_u.pct95, ...
            e.summary(a).ukf_vs_det_pct);
    end
    fprintf('%-6s %6d %8.1f %8.1f %8.1f %8.1f\n', ...
        '总计', e.overall.ukf.n, e.overall.ukf.median, ...
        e.overall.ukf.mean, e.overall.ukf.rms, e.overall.ukf.pct95);
end

%% ==================== Phase 7: 可视化 + 数据保存 ====================
fprintf('\n========== Phase 7: 可视化 ==========\n');
if ~exist('results', 'dir'), mkdir('results'); end

plot_scene_overview_multi(true_tracks, aircraft_labels, params, 'results');
plot_point_cloud_3d(detList_R1, 'R1', 'results/fig2a_R1_point_cloud.png');
plot_point_cloud_3d(detList_R2, 'R2', 'results/fig2b_R2_point_cloud.png');
plot_multi_target_detections(true_tracks, aircraft_labels, detList_R1, detList_R2, ...
    params, 'results');
plot_multi_track_result(true_tracks, aircraft_labels, detList_R1, detList_R2, ...
    trackSnapshots_R1, trackSnapshots_R2, params, 'results');
plot_error_analysis(errorStats_R1, errorStats_R2, 'results');

fprintf('\n========== Phase 7: 数据保存 ==========\n');
sysPara = struct(...
    'dt_sec', params.dt_sec, 'n_frames', n_frames, ...
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

calibResult = struct(...
    'dr1_est', dr1_est, 'da1_est', da1_est, ...
    'dr2_est', dr2_est, 'da2_est', da2_est, ...
    'dr1_true', params.radar1_range_bias_m, 'da1_true', params.radar1_azimuth_bias_deg, ...
    'dr2_true', params.radar2_range_bias_m, 'da2_true', params.radar2_azimuth_bias_deg, ...
    'n_cal_R1', length(dr1_list), 'n_cal_R2', length(dr2_list));

total_r1 = count_det_stats(detList_R1);
total_r2 = count_det_stats(detList_R2);

R1 = struct('detRaw', {detRaw_R1}, 'detList', {detList_R1}, ...
    'trackSnapshots', {trackSnapshots_R1}, 'finalTrackList', {trackList_R1}, ...
    'tempTrackList', {tempPool_R1}, 'targetDetCounts', ac_det_counts_r1, ...
    'totalDetections', total_r1.detections, 'totalClutter', total_r1.clutter);

R2 = struct('detRaw', {detRaw_R2}, 'detList', {detList_R2}, ...
    'trackSnapshots', {trackSnapshots_R2}, 'finalTrackList', {trackList_R2}, ...
    'tempTrackList', {tempPool_R2}, 'targetDetCounts', ac_det_counts_r2, ...
    'totalDetections', total_r2.detections, 'totalClutter', total_r2.clutter);

outf = fullfile('results', sprintf('simulation_%s.mat', datestr(now, 'yyyymmdd_HHMMSS')));
save(outf, 'sysPara', 'calibResult', 'truthTrajs', 'R1', 'R2', 'params', ...
    'errorStats_R1', 'errorStats_R2');
fprintf('数据已保存: %s\n', outf);
fprintf('\nDone.\n');

% =========================================================================
% 内部函数
% =========================================================================

function s = get_type_str(t)
    switch t
        case 1, s = 'RELIABLE';
        case 2, s = 'MAINTAIN';
        case 6, s = 'TEMPORARY';
        case 7, s = 'HISTORY';
        otherwise, s = 'UNKNOWN';
    end
end

function s = count_det_stats(detList)
    s.detections = 0; s.clutter = 0;
    for k = 1:length(detList)
        dets = detList{k};
        s.detections = s.detections + length(dets);
        for d = 1:length(dets)
            if dets(d).is_clutter, s.clutter = s.clutter + 1; end
        end
    end
end
