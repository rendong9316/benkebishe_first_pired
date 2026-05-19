% =========================================================================
% single_track_runner.m
% 单目标简化的逐帧航迹管理 (替代 multi_track_manager)
% =========================================================================
% 不使用 M/N 起始逻辑, 直接从首帧可用点迹初始化 UKF
% 后续帧: UKF 预测 → 最近邻关联 (加地理门限) → PDA 更新或纯预测
% =========================================================================

function [trackSnapshots, finalTrack] = single_track_runner(detList, ukf_tpl, params, n_frames)
    trackSnapshots = cell(n_frames, 1);
    ukf = [];
    track_state = 'UNINIT';
    life = 0;
    missed = 0;
    quality = 0;

    for k = 1:n_frames
        snap = struct('frameID', k, 'trackList', {{}});
        dets = detList{k};

        switch track_state
            case 'UNINIT'
                % 找距离波束中心最近的非杂波点迹作为初始检测
                if ~isempty(dets)
                    best_det = []; best_idx = 1;
                    for d = 1:length(dets)
                        if dets(d).is_clutter, continue; end
                        if isempty(best_det)
                            best_det = dets(d); best_idx = d;
                        end
                    end
                    if ~isempty(best_det)
                        % 单点初始化UKF (零速假设)
                        rng_meas = best_det.range_meas;
                        az_meas = best_det.azimuth_meas;
                        [lon, lat] = ukf_meas_to_latlon(ukf_tpl, rng_meas, az_meas);
                        ukf = ukf_tpl;
                        ukf.x = [lon; 0; lat; 0];
                        ukf.P = diag([params.ukf_P_pos_std^2, params.ukf_P_vel_std^2, ...
                                      params.ukf_P_pos_std^2, params.ukf_P_vel_std^2]);
                        ukf.initialized = true;
                        ukf.Q_base = ukf.Q;
                        ukf.Q_ema = 1.0;
                        track_state = 'TRACKING';
                        life = 1; missed = 0; quality = 5;

                        snap.trackList{1} = make_track_snap(1, 1, ukf.x(3), ukf.x(1), ...
                            ukf, life, quality, 0, best_det);
                    end
                end

            case 'TRACKING'
                ukf.dt = params.dt_sec;
                [x_pred, P_pred, X_pred, ukf] = ukf_predict_step(ukf);

                z_pred = ukf_measurement_model(ukf, x_pred);
                Z_pred = zeros(ukf.m, 2*ukf.n + 1);
                for s = 1:(2*ukf.n + 1)
                    Z_pred(:, s) = ukf_measurement_model(ukf, X_pred(:, s));
                end
                P_zz = ukf.R;
                for s = 1:(2*ukf.n + 1)
                    dz = Z_pred(:, s) - z_pred;
                    P_zz = P_zz + ukf.Wc(s) * (dz * dz');
                end
                if any(isnan(P_zz(:))), P_zz = ukf.R; end

                % 最近邻关联 (地理预筛选 + 马氏距离)
                best_det = [];
                best_mahal = inf;
                geo_gate_m = 120000;  % 初始阶段120km
                if life > 15, geo_gate_m = 60000; end  % 收敛后60km

                for d = 1:length(dets)
                    dp = dets(d);
                    if ~isfield(dp, 'lat') || isnan(dp.lat), continue; end
                    % 地理距离预筛选
                    geo_dist = sphere_utils_haversine_distance(...
                        x_pred(1), x_pred(3), dp.lon, dp.lat);
                    if geo_dist > geo_gate_m, continue; end

                    z_m = [dp.drange; dp.daz];
                    innov = z_m - z_pred(1:2);
                    if innov(2) > 180, innov(2) = innov(2) - 360;
                    elseif innov(2) < -180, innov(2) = innov(2) + 360; end
                    mahal = innov' * (P_zz(1:2,1:2) \ innov);
                    if mahal < params.gate_sigma^2 * 2 && mahal < best_mahal
                        best_mahal = mahal;
                        best_det = dp;
                    end
                end

                if ~isempty(best_det)
                    % PDA 更新
                    dets_in_gate = {best_det};
                    % 收集门内其他点迹 (用于PDA加权)
                    gate_threshold = params.gate_sigma^2 * 2;
                    for d = 1:length(dets)
                        dp = dets(d);
                        if isequal(dp, best_det), continue; end
                        if ~isfield(dp, 'drange') || isnan(dp.drange), continue; end
                        z_m = [dp.drange; dp.daz];
                        innov = z_m - z_pred(1:2);
                        if innov(2) > 180, innov(2) = innov(2) - 360;
                        elseif innov(2) < -180, innov(2) = innov(2) + 360; end
                        if innov' * (P_zz(1:2,1:2) \ innov) < gate_threshold
                            dets_in_gate{end+1} = dp;
                        end
                    end

                    [~, ~, ukf, ~, nis_val] = ukf_pda_update(ukf, dets_in_gate, ...
                        z_pred, Z_pred, X_pred, x_pred, P_pred, P_zz, params);

                    missed = 0;
                    life = life + 1;
                    quality = min(quality + 1, 15);
                    if ~isfield(ukf, 'nis_history'), ukf.nis_history = []; end
                    ukf.nis_history(end+1) = nis_val;
                    if length(ukf.nis_history) > params.fuzzy_window_size
                        ukf.nis_history(1) = [];
                    end
                    if params.use_fuzzy_adaptive && life > 12
                        ukf = ukf_fuzzy_adapt(ukf, ukf.nis_history, life, params);
                    end
                else
                    % 纯预测
                    ukf.x = x_pred;
                    ukf.P = P_pred;
                    missed = missed + 1;
                    life = life + 1;
                    quality = max(quality - 1, 0);
                    best_det = [];
                end

                % 终止检查
                if missed >= params.tracker_K_loss
                    track_state = 'LOST';
                end

                snap.trackList{1} = make_track_snap(1, 1, ukf.x(3), ukf.x(1), ...
                    ukf, life, quality, missed, best_det);

            case 'LOST'
                snap.trackList{1} = make_track_snap(1, 7, NaN, NaN, ukf, life, quality, missed, []);
        end

        trackSnapshots{k} = snap;
    end

    finalTrack = struct('id', 1, 'type', iif(strcmp(track_state,'TRACKING'),1,7), ...
        'quality', quality, 'life', life);
end

function trk = make_track_snap(id, type, lat, lon, ukf, life, quality, missed, det)
    trk.id = id;
    trk.type = type;
    trk.lat = lat;
    trk.lon = lon;
    trk.ukf = ukf;
    trk.life = life;
    trk.quality = quality;
    trk.missed = missed;
    trk.assoc_det = det;
    if ~isempty(det)
        trk.x_pred = [];
        trk.P_pred = [];
    end
end

function v = iif(cond, t, f)
    if cond, v = t; else, v = f; end
end
