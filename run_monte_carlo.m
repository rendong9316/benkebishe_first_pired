% =========================================================================
% run_monte_carlo.m
% 蒙特卡洛仿真: 1000次独立运行, 统计匹配/融合/滤波指标
% =========================================================================
% 优化: Phase0(航迹生成)+Phase1(ADS-B标定) 预计算, parfor并行Phase2-9
% =========================================================================

function run_monte_carlo(N_mc)
    if nargin < 1, N_mc = 1000; end

    addpath(genpath('.'));
    params = simulation_params();

    % ====== 预计算 Phase 0: 场景初始化 (确定性) ======
    fprintf('预计算 Phase 0 (航迹生成)...\n');
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
    end

    t1_grid = params.time_offset_radar1_sec : params.dt_sec : trajs{1}.duration_sec;
    t2_grid = params.time_offset_radar2_sec : params.dt_sec : trajs{1}.duration_sec;
    n_frames = min(length(t1_grid), length(t2_grid));
    fprintf('  帧数: %d\n', n_frames);

    % 真值航迹结构体 (用于误差评估)
    truthTrajs = cell(params.num_aircraft, 1);
    for a = 1:params.num_aircraft
        tt = true_tracks{a};
        truthTrajs{a} = struct('label', aircraft_labels{a}, 'speed_ms', aircraft_spds(a), ...
            'time_sec', tt(:,5), 'lat', tt(:,2), 'lon', tt(:,1), ...
            'lon_rate', tt(:,3), 'lat_rate', tt(:,4));
    end

    % ====== 预计算 Phase 1: ADS-B偏差标定 (大样本均值稳定) ======
    fprintf('预计算 Phase 1 (ADS-B标定)...\n');
    rng(42);  % 固定种子, 标定参数稳定
    T_adsb = readtable(params.adsb_csv_path, 'ReadVariableNames', false);
    adsb_lat = T_adsb.Var2;
    adsb_lon = T_adsb.Var3;
    dr1_list = []; da1_list = [];
    dr2_list = []; da2_list = [];
    n_check = min(5000, height(T_adsb));
    cal_step = max(1, floor(height(T_adsb) / n_check));
    for idx = 1:cal_step:height(T_adsb)
        t_lon = adsb_lon(idx); t_lat = adsb_lat(idx);
        if isnan(t_lon) || isnan(t_lat), continue; end
        [in1, ~, ~] = radar_coverage_check(params.radar1_lon, params.radar1_lat, t_lon, t_lat, params.radar1_beam_center_deg, params);
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
        [in2, ~, ~] = radar_coverage_check(params.radar2_lon, params.radar2_lat, t_lon, t_lat, params.radar2_beam_center_deg, params);
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
    dr1_est = mean(dr1_list); da1_est = mean(da1_list);
    dr2_est = mean(dr2_list); da2_est = mean(da2_list);
    fprintf('  R1: dr=%.1fm da=%.4fdeg | R2: dr=%.1fm da=%.4fdeg\n', dr1_est, da1_est, dr2_est, da2_est);

    % 打包预计算数据为struct (parfor要求单一变量传递)
    pre = struct();
    pre.params = params;
    pre.trajs = {trajs};
    pre.n_frames = n_frames;
    pre.t1_grid = t1_grid;
    pre.t2_grid = t2_grid;
    pre.dr1_est = dr1_est; pre.da1_est = da1_est;
    pre.dr2_est = dr2_est; pre.da2_est = da2_est;
    pre.truthTrajs = {truthTrajs};
    pre.aircraft_labels = {aircraft_labels};
    pre.aircraft_spds = aircraft_spds;

    % ====== Monte Carlo 主循环 (parfor) ======
    method_names = {'SCC', 'BC', 'CI', 'FCI'};
    n_methods = length(method_names);
    n_ac = params.num_aircraft;

    % 预分配结果数组
    match_rate     = zeros(N_mc, 1);
    n_matched_arr  = zeros(N_mc, 1);
    fusion_rms     = zeros(N_mc, 1);
    fusion_med     = zeros(N_mc, 1);
    r1_rms         = zeros(N_mc, 1);
    r1_med         = zeros(N_mc, 1);
    r2_rms         = zeros(N_mc, 1);
    r2_med         = zeros(N_mc, 1);
    improvement    = zeros(N_mc, 1);
    ac_fusion_rms  = zeros(N_mc, n_ac);
    ac_fusion_med  = zeros(N_mc, n_ac);
    ac_r1_rms      = zeros(N_mc, n_ac);
    ac_r1_med      = zeros(N_mc, n_ac);
    ac_r2_rms      = zeros(N_mc, n_ac);
    ac_r2_med      = zeros(N_mc, n_ac);
    total_tracks_r1 = zeros(N_mc, 1);
    total_tracks_r2 = zeros(N_mc, 1);
    best_m_arr      = zeros(N_mc, 1);

    t_start = tic;
    fprintf('\n========== Monte Carlo N=%d (parfor 6 workers) ==========\n', N_mc);

    parfor mc = 1:N_mc
        seed = mc * 100 + 1;  % 确保种子分散
        result = simulate_once_parallel(pre, seed);

        match_rate(mc)      = result.match_correct_rate;
        n_matched_arr(mc)   = result.n_matched;
        fusion_rms(mc)      = result.fusion_rmse;
        fusion_med(mc)      = result.fusion_median;
        r1_rms(mc)          = result.r1_rmse;
        r1_med(mc)          = result.r1_median;
        r2_rms(mc)          = result.r2_rmse;
        r2_med(mc)          = result.r2_median;
        improvement(mc)     = result.improvement_pct;
        ac_fusion_rms(mc,:) = result.ac_fusion_rmse;
        ac_fusion_med(mc,:) = result.ac_fusion_med;
        ac_r1_rms(mc,:)     = result.ac_r1_rmse;
        ac_r1_med(mc,:)     = result.ac_r1_med;
        ac_r2_rms(mc,:)     = result.ac_r2_rmse;
        ac_r2_med(mc,:)     = result.ac_r2_med;
        total_tracks_r1(mc) = result.n_tracks_r1;
        total_tracks_r2(mc) = result.n_tracks_r2;
        best_m_arr(mc)      = result.best_fusion_method;
    end

    elapsed = toc(t_start);
    fprintf('\nMonte Carlo 完成: %d 次 / %.0f 秒 (%.1f 秒/次)\n', N_mc, elapsed, elapsed/N_mc);

    % ====== 统计分析 ======
    fprintf('\n========== Monte Carlo 统计汇总 (N=%d) ==========\n', N_mc);

    function p(label, vals, unit)
        if nargin < 3, unit = ''; end
        v = vals(~isnan(vals) & ~isinf(vals));
        if isempty(v)
            fprintf('  %-25s (无有效数据)\n', label);
        else
            fprintf('  %-25s  %7.2f ± %6.2f %s  [%5.1f, %5.1f]\n', ...
                label, mean(v), std(v), unit, prctile(v,5), prctile(v,95));
        end
    end

    fprintf('\n--- 匹配性能 ---\n');
    p('正确匹配率', match_rate, '');
    n_complete = sum(total_tracks_r1 >= n_ac & total_tracks_r2 >= n_ac);
    fprintf('  完整匹配运行数: %d/%d (%.0f%%)\n', n_complete, N_mc, n_complete/N_mc*100);
    p('匹配对数', n_matched_arr, '对');
    p('R1活跃航迹数', total_tracks_r1, '条');
    p('R2活跃航迹数', total_tracks_r2, '条');

    fprintf('\n--- 融合性能 (总体, 最佳融合算法) ---\n');
    p('融合 RMSE', fusion_rms, 'km');
    p('融合 中位误差', fusion_med, 'km');
    p('R1_only RMSE', r1_rms, 'km');
    p('R2_only RMSE', r2_rms, 'km');
    p('融合 vs 最佳单站改善', improvement, '%');

    % 最佳算法分布
    fprintf('\n--- 最佳融合算法分布 ---\n');
    for m = 1:n_methods
        cnt = sum(best_m_arr == m);
        fprintf('  %s: %d 次 (%.1f%%)\n', method_names{m}, cnt, cnt/N_mc*100);
    end

    fprintf('\n--- 分飞机 RMSE (最佳融合) ---\n');
    fprintf('  %-8s', '算法');
    for a = 1:n_ac
        fprintf('  飞机%-4s  ', pre.aircraft_labels{1}{a});
    end
    fprintf('  总体\n');
    fprintf('  Fusion ');
    for a = 1:n_ac
        v = ac_fusion_rms(:,a); v = v(~isnan(v) & ~isinf(v));
        fprintf('%7.1f±%.1f', mean(v), std(v));
    end
    v = fusion_rms(~isnan(fusion_rms) & ~isinf(fusion_rms));
    fprintf('%7.1f±%.1f\n', mean(v), std(v));
    fprintf('  R1     ');
    for a = 1:n_ac
        v = ac_r1_rms(:,a); v = v(~isnan(v) & ~isinf(v));
        fprintf('%7.1f±%.1f', mean(v), std(v));
    end
    v = r1_rms(~isnan(r1_rms) & ~isinf(r1_rms));
    fprintf('%7.1f±%.1f\n', mean(v), std(v));
    fprintf('  R2     ');
    for a = 1:n_ac
        v = ac_r2_rms(:,a); v = v(~isnan(v) & ~isinf(v));
        fprintf('%7.1f±%.1f', mean(v), std(v));
    end
    v = r2_rms(~isnan(r2_rms) & ~isinf(r2_rms));
    fprintf('%7.1f±%.1f\n', mean(v), std(v));

    fprintf('\n--- 分飞机中位误差 (最佳融合) ---\n');
    fprintf('  Fusion ');
    for a = 1:n_ac
        v = ac_fusion_med(:,a); v = v(~isnan(v) & ~isinf(v));
        fprintf('%7.1f±%.1f', mean(v), std(v));
    end
    v = fusion_med(~isnan(fusion_med) & ~isinf(fusion_med));
    fprintf('%7.1f±%.1f\n', mean(v), std(v));
    fprintf('  R1     ');
    for a = 1:n_ac
        v = ac_r1_med(:,a); v = v(~isnan(v) & ~isinf(v));
        fprintf('%7.1f±%.1f', mean(v), std(v));
    end
    v = r1_med(~isnan(r1_med) & ~isinf(r1_med));
    fprintf('%7.1f±%.1f\n', mean(v), std(v));
    fprintf('  R2     ');
    for a = 1:n_ac
        v = ac_r2_med(:,a); v = v(~isnan(v) & ~isinf(v));
        fprintf('%7.1f±%.1f', mean(v), std(v));
    end
    v = r2_med(~isnan(r2_med) & ~isinf(r2_med));
    fprintf('%7.1f±%.1f\n', mean(v), std(v));

    % ====== 保存 ======
    if ~exist('results', 'dir'), mkdir('results'); end
    outf = fullfile('results', sprintf('monte_carlo_N%d_%s.mat', N_mc, datestr(now, 'yyyymmdd_HHMMSS')));
    save(outf, 'N_mc', 'match_rate', 'n_matched_arr', 'fusion_rms', 'fusion_med', ...
        'r1_rms', 'r1_med', 'r2_rms', 'r2_med', 'improvement', ...
        'ac_fusion_rms', 'ac_fusion_med', 'ac_r1_rms', 'ac_r1_med', 'ac_r2_rms', 'ac_r2_med', ...
        'total_tracks_r1', 'total_tracks_r2', 'best_m_arr', 'method_names', ...
        'aircraft_labels', 'elapsed');
    fprintf('\n数据已保存: %s\n', outf);
    fprintf('Done.\n');
end

% =========================================================================
% simulate_once_parallel: parfor兼容的单次仿真
% =========================================================================
function result = simulate_once_parallel(pre, seed)
    params = pre.params;
    trajs = pre.trajs{1};
    n_frames = pre.n_frames;
    t1_grid = pre.t1_grid;
    t2_grid = pre.t2_grid;
    dr1_est = pre.dr1_est; da1_est = pre.da1_est;
    dr2_est = pre.dr2_est; da2_est = pre.da2_est;
    truthTrajs = pre.truthTrajs{1};
    aircraft_labels = pre.aircraft_labels{1};
    n_ac = params.num_aircraft;

    rng(seed);  % 每次独立随机种子

    % ---- Phase 2: 原始点迹生成 ----
    detRaw_R1 = cell(n_frames, 1);
    detRaw_R2 = cell(n_frames, 1);
    for k = 1:n_frames
        all_raw_r1 = [];
        for a = 1:n_ac
            [pos, vel] = aircraft_trajectory_interpolate(trajs{a}, t1_grid(k));
            add_clut = (a == 1);
            rng(seed + a*1000 + k);
            dets = generate_frame_detections(params.radar1_lon, params.radar1_lat, ...
                params.radar1_tx_lon, params.radar1_tx_lat, ...
                pos(1), pos(2), vel(1), vel(2), k, t1_grid(k), ...
                params.radar1_range_bias_m, params.radar1_azimuth_bias_deg, ...
                params.radar1_beam_center_deg, params, add_clut);
            for d = 1:length(dets), dets(d).aircraft_id = a; end
            all_raw_r1 = [all_raw_r1, dets];
        end
        detRaw_R1{k} = all_raw_r1;

        all_raw_r2 = [];
        for a = 1:n_ac
            [pos, vel] = aircraft_trajectory_interpolate(trajs{a}, t2_grid(k));
            add_clut = (a == 1);
            rng(seed + 10000 + a*1000 + k);
            dets = generate_frame_detections(params.radar2_lon, params.radar2_lat, ...
                params.radar2_tx_lon, params.radar2_tx_lat, ...
                pos(1), pos(2), vel(1), vel(2), k, t2_grid(k), ...
                params.radar2_range_bias_m, params.radar2_azimuth_bias_deg, ...
                params.radar2_beam_center_deg, params, add_clut);
            for d = 1:length(dets), dets(d).aircraft_id = a; end
            all_raw_r2 = [all_raw_r2, dets];
        end
        detRaw_R2{k} = all_raw_r2;
    end

    % ---- Phase 4: 偏差校正 ----
    detList_R1 = cell(n_frames, 1);
    detList_R2 = cell(n_frames, 1);
    for k = 1:n_frames
        dets_r1 = detRaw_R1{k};
        for d = 1:length(dets_r1)
            Rgc = dets_r1(d).prange - dr1_est; azc = dets_r1(d).paz - da1_est;
            dets_r1(d).drange = Rgc; dets_r1(d).daz = azc;
            dets_r1(d).range_meas = Rgc; dets_r1(d).azimuth_meas = azc;
            if ~(isfield(dets_r1(d), 'lat') && ~isnan(dets_r1(d).lat))
                [~, lat_e, lon_e] = bistatic_inverse_solver(Rgc, azc, ...
                    params.radar1_tx_lon, params.radar1_tx_lat, params.radar1_lon, params.radar1_lat);
                dets_r1(d).lat = lat_e; dets_r1(d).lon = lon_e;
            end
            [~, raw_lat, raw_lon] = bistatic_inverse_solver(dets_r1(d).prange, dets_r1(d).paz, ...
                params.radar1_tx_lon, params.radar1_tx_lat, params.radar1_lon, params.radar1_lat);
            dets_r1(d).raw_lat = raw_lat; dets_r1(d).raw_lon = raw_lon;
        end
        detList_R1{k} = dets_r1;

        dets_r2 = detRaw_R2{k};
        for d = 1:length(dets_r2)
            Rgc = dets_r2(d).prange - dr2_est; azc = dets_r2(d).paz - da2_est;
            dets_r2(d).drange = Rgc; dets_r2(d).daz = azc;
            dets_r2(d).range_meas = Rgc; dets_r2(d).azimuth_meas = azc;
            if ~(isfield(dets_r2(d), 'lat') && ~isnan(dets_r2(d).lat))
                [~, lat_e, lon_e] = bistatic_inverse_solver(Rgc, azc, ...
                    params.radar2_tx_lon, params.radar2_tx_lat, params.radar2_lon, params.radar2_lat);
                dets_r2(d).lat = lat_e; dets_r2(d).lon = lon_e;
            end
            [~, raw_lat, raw_lon] = bistatic_inverse_solver(dets_r2(d).prange, dets_r2(d).paz, ...
                params.radar2_tx_lon, params.radar2_tx_lat, params.radar2_lon, params.radar2_lat);
            dets_r2(d).raw_lat = raw_lat; dets_r2(d).raw_lon = raw_lon;
        end
        detList_R2{k} = dets_r2;
    end

    % ---- Phase 5: 多目标航迹跟踪 ----
    ukf1_tpl = ukf_filter(params, params.radar1_lon, params.radar1_lat, ...
        params.radar1_tx_lon, params.radar1_tx_lat, params.dt_sec);
    ukf2_tpl = ukf_filter(params, params.radar2_lon, params.radar2_lat, ...
        params.radar2_tx_lon, params.radar2_tx_lat, params.dt_sec);

    trackList_R1 = {}; tempPool_R1 = {};
    trackList_R2 = {}; tempPool_R2 = {};
    trackSnapshots_R1 = cell(n_frames, 1);
    trackSnapshots_R2 = cell(n_frames, 1);

    for k = 1:n_frames
        [trackList_R1, tempPool_R1, trackSnapshots_R1{k}] = multi_track_manager(...
            trackList_R1, tempPool_R1, detList_R1{k}, ukf1_tpl, params, k);
        [trackList_R2, tempPool_R2, trackSnapshots_R2{k}] = multi_track_manager(...
            trackList_R2, tempPool_R2, detList_R2{k}, ukf2_tpl, params, k);
    end

    % ---- Phase 6: 航迹匹配 ----
    matcher = track_matcher(trackSnapshots_R1, trackSnapshots_R2, params);
    n_matched = length(matcher.matched_pairs);

    % 匹配正确率验证
    matched_ac_r1 = zeros(n_matched, 1);
    matched_ac_r2 = zeros(n_matched, 1);
    for p = 1:n_matched
        mp = matcher.matched_pairs(p);
        r1_idx = find(matcher.r1_ids == mp.R1_track_id, 1);
        r2_idx = find(matcher.r2_ids == mp.R2_track_id, 1);

        best_d = inf; best_a = 1;
        for a = 1:n_ac
            tt = truthTrajs{a};
            t_lat = interp1(tt.time_sec, tt.lat, t1_grid, 'linear', 'extrap');
            t_lon = interp1(tt.time_sec, tt.lon, t1_grid, 'linear', 'extrap');
            if ~isempty(r1_idx)
                r1_lons = squeeze(matcher.r1_pos(r1_idx,:,1))';
                r1_lats = squeeze(matcher.r1_pos(r1_idx,:,2))';
                valid = ~isnan(r1_lons);
                if any(valid)
                    d = mean(haversine_km_vec(r1_lons(valid), r1_lats(valid), t_lon(valid), t_lat(valid)));
                    if ~isnan(d) && d < best_d, best_d = d; best_a = a; end
                end
            end
        end
        matched_ac_r1(p) = best_a;

        best_d = inf; best_a = 1;
        for a = 1:n_ac
            tt = truthTrajs{a};
            t_lat = interp1(tt.time_sec, tt.lat, t1_grid, 'linear', 'extrap');
            t_lon = interp1(tt.time_sec, tt.lon, t1_grid, 'linear', 'extrap');
            if ~isempty(r2_idx)
                r2_lons = squeeze(matcher.r2_pos(r2_idx,:,1))';
                r2_lats = squeeze(matcher.r2_pos(r2_idx,:,2))';
                valid = ~isnan(r2_lons);
                if any(valid)
                    d = mean(haversine_km_vec(r2_lons(valid), r2_lats(valid), t_lon(valid), t_lat(valid)));
                    if ~isnan(d) && d < best_d, best_d = d; best_a = a; end
                end
            end
        end
        matched_ac_r2(p) = best_a;
    end

    n_correct = sum(matched_ac_r1 == matched_ac_r2);
    match_correct_rate = n_correct / max(n_matched, 1);

    % ---- Phase 7: 航迹融合 ----
    method_names = {'SCC', 'BC', 'CI', 'FCI'};
    all_fused_snapshots = cell(length(method_names), 1);
    for m = 1:length(method_names)
        all_fused_snapshots{m} = run_track_fusion(matcher.matched_pairs, ...
            trackSnapshots_R1, matcher.aligned_R2, params, method_names{m});
    end

    % ---- Phase 8: 融合误差评估 ----
    fusion_eval = evaluate_fusion(all_fused_snapshots, method_names, ...
        matcher.matched_pairs, trackSnapshots_R1, trackSnapshots_R2, ...
        truthTrajs, n_frames, params.dt_sec, matcher);

    % ---- 提取指标 (最佳融合算法) ----
    n_methods = length(method_names);
    best_fusion_rmse = inf;
    best_fusion_med = NaN;
    best_method = 1;
    for m = 1:n_methods
        if fusion_eval.overall(m).s.rms < best_fusion_rmse
            best_fusion_rmse = fusion_eval.overall(m).s.rms;
            best_fusion_med  = fusion_eval.overall(m).s.median;
            best_method = m;
        end
    end

    ac_fusion_rmse = zeros(1, n_ac);
    ac_fusion_med = zeros(1, n_ac);
    ac_r1_rmse = zeros(1, n_ac);
    ac_r1_med = zeros(1, n_ac);
    ac_r2_rmse = zeros(1, n_ac);
    ac_r2_med = zeros(1, n_ac);

    for a = 1:n_ac
        idx = (a-1)*(n_methods+2) + best_method;
        s_f = fusion_eval.summary(idx).s;
        ac_fusion_rmse(a) = s_f.rms;
        ac_fusion_med(a) = s_f.median;

        s_r1 = fusion_eval.summary((a-1)*(n_methods+2) + n_methods+1).s;
        s_r2 = fusion_eval.summary((a-1)*(n_methods+2) + n_methods+2).s;
        ac_r1_rmse(a) = s_r1.rms;
        ac_r1_med(a) = s_r1.median;
        ac_r2_rmse(a) = s_r2.rms;
        ac_r2_med(a) = s_r2.median;
    end

    r1_rmse = fusion_eval.overall(n_methods+1).s.rms;
    r1_med  = fusion_eval.overall(n_methods+1).s.median;
    r2_rmse = fusion_eval.overall(n_methods+2).s.rms;
    r2_med  = fusion_eval.overall(n_methods+2).s.median;

    improvement_pct = (1 - best_fusion_rmse/min(r1_rmse, r2_rmse)) * 100;

    % 航迹统计
    n_tracks_r1 = 0; n_tracks_r2 = 0;
    for t = 1:length(trackList_R1)
        if trackList_R1{t}.type ~= 7, n_tracks_r1 = n_tracks_r1 + 1; end
    end
    for t = 1:length(trackList_R2)
        if trackList_R2{t}.type ~= 7, n_tracks_r2 = n_tracks_r2 + 1; end
    end

    result = struct(...
        'seed', seed, ...
        'n_matched', n_matched, ...
        'match_correct_rate', match_correct_rate, ...
        'best_fusion_method', best_method, ...
        'fusion_rmse', best_fusion_rmse, ...
        'fusion_median', best_fusion_med, ...
        'r1_rmse', r1_rmse, 'r1_median', r1_med, ...
        'r2_rmse', r2_rmse, 'r2_median', r2_med, ...
        'improvement_pct', improvement_pct, ...
        'ac_fusion_rmse', ac_fusion_rmse, ...
        'ac_fusion_med', ac_fusion_med, ...
        'ac_r1_rmse', ac_r1_rmse, ...
        'ac_r1_med', ac_r1_med, ...
        'ac_r2_rmse', ac_r2_rmse, ...
        'ac_r2_med', ac_r2_med, ...
        'n_tracks_r1', n_tracks_r1, ...
        'n_tracks_r2', n_tracks_r2);
end

function d_vec = haversine_km_vec(lon1, lat1, lon2, lat2)
    R = 6371;
    d_vec = zeros(size(lon1));
    for i = 1:length(lon1)
        if isnan(lon1(i)) || isnan(lat1(i)) || isnan(lon2(i)) || isnan(lat2(i))
            d_vec(i) = NaN; continue;
        end
        dlat = deg2rad(lat2(i) - lat1(i));
        dlon = deg2rad(lon2(i) - lon1(i));
        a = sin(dlat/2)^2 + cos(deg2rad(lat1(i)))*cos(deg2rad(lat2(i)))*sin(dlon/2)^2;
        a = max(0, min(1, a));
        d_vec(i) = R * 2 * atan2(sqrt(a), sqrt(1 - a));
    end
end
