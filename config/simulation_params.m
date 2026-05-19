% =========================================================================
% simulation_params.m
% 双基地外辐射源雷达单目标仿真 — 全局参数配置（优化版）
% =========================================================================

function params = simulation_params()
    % ==================== 1. 时间参数 ====================
    params.dt_sec = 30.0;
    params.duration_sec = 3600.0;
    params.ref_start_time = datetime(2026, 4, 27, 9, 30, 0);
    params.time_offset_radar1_sec = 0.0;
    params.time_offset_radar2_sec = 13.0;

    % ==================== 2. 站点几何布局 ====================
    params.radar1_lon = 112.0;   params.radar1_lat = 33.5;
    params.radar2_lon = 114.0;   params.radar2_lat = 31.5;
    params.radar1_tx_lon = 108.0;  params.radar1_tx_lat = 33.5;
    params.radar2_tx_lon = 110.0;  params.radar2_tx_lat = 31.5;

    % ==================== 3. 雷达探测范围约束 ====================
    params.radar1_beam_center_deg = 91.5;
    params.radar2_beam_center_deg = 83.5;
    params.beam_width_deg = 15.0;
    params.range_min_km = 1000.0;
    params.range_max_km = 2000.0;
    params.range_min_m = params.range_min_km * 1000;
    params.range_max_m = params.range_max_km * 1000;

    % ==================== 4. 目标航迹 ====================
    params.aircraft_waypoints = [127.5, 31.0, 0.0; 130.5, 33.0, 0.0];
    params.aircraft_speed_ms = 230.0;
    params.trajectory_mode = "straight";

    % ==================== 5. 量测噪声 ====================
    params.range_noise_std_m = 4000.0;
    params.azimuth_noise_std_deg = 0.4;
    params.radial_vel_noise_std_ms = 0.5;

    % ==================== 6. 系统偏差 ====================
    params.radar1_range_bias_m = 20000.0;
    params.radar1_azimuth_bias_deg = -3.0;
    params.radar2_range_bias_m = -15000.0;
    params.radar2_azimuth_bias_deg = 3.5;
    params.ukf_range_std_m = params.range_noise_std_m;
    params.ukf_azimuth_std_deg = params.azimuth_noise_std_deg;
    params.ukf_rv_std_ms = params.radial_vel_noise_std_ms;

    % ==================== 7. UKF滤波参数 (单目标极简调优) ====================
    params.ukf_alpha = 1e-3;
    params.ukf_beta = 2.0;
    params.ukf_kappa = 0.0;
    params.ukf_Q_scale = 5e4;
    params.ukf_P_pos_std = 0.2;
    params.ukf_P_vel_std = 0.003;
    params.ukf_mode = "standard";

    % ==================== 8. 航迹管理参数 ====================
    params.tracker_M = 4;
    params.tracker_N = 8;
    params.tracker_K_loss = 12;
    params.gate_sigma = 2.0;

    % ==================== 9. 检测概率 & 虚警率 ====================
    params.detection_probability = 0.6;
    params.false_alarm_rate = 0.001;
    params.range_resolution_km = 10.0;
    params.azimuth_resolution_deg = 1.0;
    params.n_resolution_cells = ...
        ((params.range_max_km - params.range_min_km) / params.range_resolution_km) * ...
        (params.beam_width_deg / params.azimuth_resolution_deg);

    % ==================== 10. PDA栅格加权参数 ====================
    params.use_pda_weighting = true;
    params.pda_pd_gate = 0.8647;
    params.pda_clutter_intensity = 1.5 / (2000e3 * 15);

    % ==================== 11. 模糊自适应参数 ====================
    params.use_fuzzy_adaptive = true;
    params.fuzzy_window_size = 8;
    params.fuzzy_Q_min_factor = 0.6;
    params.fuzzy_Q_max_factor = 1.5;

    % ADS-B标定数据
    params.adsb_csv_path = '2026-04-27 09-30-00.csv';

    % ==================== 12. 随机种子 ====================
    params.random_seed = 42;
end
