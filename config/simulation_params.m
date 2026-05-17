% =========================================================================
% simulation_params.m
% 双基地外辐射源雷达仿真 —— 全局参数配置
% =========================================================================
% 参数分类：
%   1. 时间参数
%   2. 站点几何布局（Tx发射站 + Rx接收站）
%   3. 雷达探测范围约束
%   4. 目标航迹
%   5. 量测噪声
%   6. 系统偏差（真值 + UKF用的噪声）
%   7. UKF滤波参数
%   8. 航迹管理参数（M/N起始 + K_loss终止 + 关联波门）
%   9. 检测概率 & 虚警率
%  10. 随机种子
% =========================================================================

function params = simulation_params()
    % ==================== 1. 时间参数 ====================
    params.dt_sec = 30.0;                          % 数据更新间隔 (s)
    params.duration_sec = 3600.0;                  % 仿真总时长 (s), 1小时

    % 参考起始时间（与ADS-B数据时间对齐）
    params.ref_start_time = datetime(2026, 4, 27, 9, 30, 0);

    % 两部雷达异步采样偏移
    params.time_offset_radar1_sec = 0.0;           % R1准点采样
    params.time_offset_radar2_sec = 13.0;          % R2偏移13s

    % ==================== 2. 站点几何布局 ====================
    % 接收站
    params.radar1_lon = 112.0;   params.radar1_lat = 33.5;   % R1
    params.radar2_lon = 114.0;   params.radar2_lat = 31.5;   % R2

    % 照射站（发射站，位于接收站西侧）
    params.radar1_tx_lon = 108.0;  params.radar1_tx_lat = 33.5;  % Tx1
    params.radar2_tx_lon = 110.0;  params.radar2_tx_lat = 31.5;  % Tx2

    % ==================== 3. 雷达探测范围约束 ====================
    % 波束中心分别指向目标区域（~129°E, 32°N），从各站计算：
    %   R1(112,33.5)→目标区域: az ≈ 91.5°
    %   R2(114,31.5)→目标区域: az ≈ 83.5°
    params.radar1_beam_center_deg = 91.5;
    params.radar2_beam_center_deg = 83.5;
    params.beam_width_deg = 15.0;
    params.range_min_km = 1000.0;                  % 探测近界 (km)
    params.range_max_km = 2000.0;                  % 探测远界 (km)
    % 派生：距离范围 (m)
    params.range_min_m = params.range_min_km * 1000;
    params.range_max_m = params.range_max_km * 1000;

    % ==================== 4. 目标航迹 ====================
    % 航迹需在双方雷达探测范围内
    % 从R1(112,33.5): 正东~126-130°E, 32-34°N → azimuth ~85-95°
    % 从R2(114,31.5): 正东~126-130°E, 32-34°N → azimuth ~82-97°
    % 航迹: 东海自南向北飞行
    params.aircraft_waypoints = [ ...
        127.5, 31.0, 0.0; ...   % 起点: 东海南部
        130.5, 33.0, 0.0 ...    % 终点: 济州岛东南
    ];
    params.aircraft_speed_ms = 230.0;              % 巡航速度 ~828 km/h
    params.trajectory_mode = "straight";            % 大圆直飞

    % ==================== 5. 量测噪声 ====================
    params.range_noise_std_m = 10000.0;            % 距离噪声 σ (m)
    params.azimuth_noise_std_deg = 1.0;            % 方位噪声 σ (deg)
    params.radial_vel_noise_std_ms = 0.5;          % 径向速度噪声 σ (m/s)

    % ==================== 6. 系统偏差（真值 + UKF量测噪声） ====================
    params.radar1_range_bias_m = 20000.0;          % R1距离偏置真值 (m)
    params.radar1_azimuth_bias_deg = -3.0;         % R1方位偏置真值 (deg)
    params.radar2_range_bias_m = -15000.0;         % R2距离偏置真值 (m)
    params.radar2_azimuth_bias_deg = 3.5;          % R2方位偏置真值 (deg)

    params.ukf_range_std_m = params.range_noise_std_m;
    params.ukf_azimuth_std_deg = params.azimuth_noise_std_deg;
    params.ukf_rv_std_ms = params.radial_vel_noise_std_ms;

    % ==================== 7. UKF滤波参数 ====================
    params.ukf_alpha = 1e-3;
    params.ukf_beta = 2.0;
    params.ukf_kappa = 0.0;
    params.ukf_Q_scale = 8e5;
    params.ukf_P_pos_std = 0.5;
    params.ukf_P_vel_std = 0.01;
    params.ukf_mode = "standard";

    % ==================== 8. 航迹管理参数 ====================
    params.tracker_M = 3;                          % M/N起始: N帧中至少M帧
    params.tracker_N = 5;
    params.tracker_K_loss = 4;                     % 连续K帧无关联→终止
    params.gate_sigma = 3.0;                       % 关联波门: 3σ椭圆门

    % ==================== 9. 检测概率 & 虚警率 ====================
    params.detection_probability = 0.6;            % Pd = 60%
    params.false_alarm_rate = 0.001;               % Pfa = 10^{-3}
    % 分辨率单元数 = (2000-1000)/10 × 15/1 = 100 × 15 = 1500
    params.range_resolution_km = 10.0;
    params.azimuth_resolution_deg = 1.0;
    params.n_resolution_cells = ...
        ((params.range_max_km - params.range_min_km) / params.range_resolution_km) * ...
        (params.beam_width_deg / params.azimuth_resolution_deg);

    % ==================== 10. 随机种子 ====================
    params.random_seed = 42;
end
