% =========================================================================
% plot_multi_track_result.m
% 多目标跟踪综合图 — 真值 + 原始/校准点迹 + UKF滤波航迹，带图层复选框
% =========================================================================

function plot_multi_track_result(true_tracks, labels, detList_R1, detList_R2, ...
        trackSnapshots_R1, trackSnapshots_R2, params, out_dir)

    ac_colors = {'g', 'm', 'c'};  % A=绿, B=洋红, C=青
    n_ac = length(true_tracks);

    fig = figure('Position', [50, 50, 1400, 750]);

    % ---- 左侧地理图 ----
    ax = geoaxes('Units', 'normalized', 'Position', [0.04, 0.10, 0.70, 0.88]);
    ax.Basemap = 'landcover';
    hold(ax, 'on');

    % ---- Layer 0: 真实航迹 (黑色虚线) ----
    h_truth = gobjects(1, n_ac);
    for a = 1:n_ac
        tt = true_tracks{a};
        col = ac_colors{a};
        h_truth(a) = geoplot(ax, tt(:,2), tt(:,1), '--', 'Color', col, ...
            'LineWidth', 2, 'DisplayName', sprintf('%s 真值', labels{a}));
    end

    % ---- Layer 1: R1 原始点迹（校准前） ----
    [r1_raw_lat, r1_raw_lon] = extract_dets_by_aircraft(detList_R1, 'raw');
    h_r1_raw = gobjects(1, n_ac);
    for a = 1:n_ac
        if ~isempty(r1_raw_lat{a})
            h_r1_raw(a) = geoplot(ax, r1_raw_lat{a}, r1_raw_lon{a}, '--', ...
                'Color', [0.4 0.6 1.0], 'LineWidth', 1.0, 'Marker', 'o', ...
                'MarkerSize', 4, 'MarkerFaceColor', [0.4 0.6 1.0], ...
                'DisplayName', sprintf('R1 %s 原始', labels{a}));
        end
    end

    % ---- Layer 2: R1 校准后点迹 ----
    [r1_cal_lat, r1_cal_lon] = extract_dets_by_aircraft(detList_R1, 'cal');
    h_r1_cal = gobjects(1, n_ac);
    for a = 1:n_ac
        if ~isempty(r1_cal_lat{a})
            h_r1_cal(a) = geoplot(ax, r1_cal_lat{a}, r1_cal_lon{a}, '-', ...
                'Color', [0.0 0.4 1.0], 'LineWidth', 1.2, 'Marker', 'o', ...
                'MarkerSize', 5, 'MarkerFaceColor', 'b', ...
                'DisplayName', sprintf('R1 %s 校准', labels{a}));
        end
    end

    % ---- Layer 3: R1 UKF滤波航迹 ----
    r1_tracks = collect_active_tracks(trackSnapshots_R1);
    h_r1_ukf = gobjects(1, length(r1_tracks));
    for t = 1:length(r1_tracks)
        trk = r1_tracks{t};
        if length(trk.lat_history) > 2
            h_r1_ukf(t) = geoplot(ax, trk.lat_history, trk.lon_history, 'b-', ...
                'LineWidth', 2.5, 'DisplayName', sprintf('R1 UKF#%d', trk.id));
        end
    end

    % ---- Layer 4: R2 原始点迹（校准前） ----
    [r2_raw_lat, r2_raw_lon] = extract_dets_by_aircraft(detList_R2, 'raw');
    h_r2_raw = gobjects(1, n_ac);
    for a = 1:n_ac
        if ~isempty(r2_raw_lat{a})
            h_r2_raw(a) = geoplot(ax, r2_raw_lat{a}, r2_raw_lon{a}, '--', ...
                'Color', [1.0 0.6 0.6], 'LineWidth', 1.0, 'Marker', 'o', ...
                'MarkerSize', 4, 'MarkerFaceColor', [1.0 0.6 0.6], ...
                'DisplayName', sprintf('R2 %s 原始', labels{a}));
        end
    end

    % ---- Layer 5: R2 校准后点迹 ----
    [r2_cal_lat, r2_cal_lon] = extract_dets_by_aircraft(detList_R2, 'cal');
    h_r2_cal = gobjects(1, n_ac);
    for a = 1:n_ac
        if ~isempty(r2_cal_lat{a})
            h_r2_cal(a) = geoplot(ax, r2_cal_lat{a}, r2_cal_lon{a}, '-', ...
                'Color', [1.0 0.2 0.2], 'LineWidth', 1.2, 'Marker', 'o', ...
                'MarkerSize', 5, 'MarkerFaceColor', 'r', ...
                'DisplayName', sprintf('R2 %s 校准', labels{a}));
        end
    end

    % ---- Layer 6: R2 UKF滤波航迹 ----
    r2_tracks = collect_active_tracks(trackSnapshots_R2);
    h_r2_ukf = gobjects(1, length(r2_tracks));
    for t = 1:length(r2_tracks)
        trk = r2_tracks{t};
        if length(trk.lat_history) > 2
            h_r2_ukf(t) = geoplot(ax, trk.lat_history, trk.lon_history, 'r-', ...
                'LineWidth', 2.5, 'DisplayName', sprintf('R2 UKF#%d', trk.id));
        end
    end

    % ---- 站点标记（始终可见） ----
    geoplot(ax, params.radar1_lat, params.radar1_lon, 'bs', ...
        'MarkerSize', 14, 'MarkerFaceColor', 'b', 'DisplayName', 'R1');
    geoplot(ax, params.radar2_lat, params.radar2_lon, 'rs', ...
        'MarkerSize', 14, 'MarkerFaceColor', 'r', 'DisplayName', 'R2');
    geoplot(ax, params.radar1_tx_lat, params.radar1_tx_lon, 'b^', ...
        'MarkerSize', 10, 'DisplayName', 'Tx1');
    geoplot(ax, params.radar2_tx_lat, params.radar2_tx_lon, 'r^', ...
        'MarkerSize', 10, 'DisplayName', 'Tx2');

    title(ax, '多目标双基地雷达航迹综合对比');
    legend(ax, 'Location', 'northeastoutside');

    % ---- 右侧控制面板 ----
    panel = uipanel('Units', 'normalized', 'Position', [0.76, 0.04, 0.22, 0.94], ...
        'Title', '图层显隐控制', 'FontSize', 11);

    % 收集所有图层句柄和名称（配对添加，保证对齐）
    h_all = [];
    layer_names = {};

    for a = 1:n_ac
        h_all(end+1) = h_truth(a);
        layer_names{end+1} = sprintf('%s 真值', labels{a});
    end
    for a = 1:n_ac
        if h_r1_raw(a) ~= 0
            h_all(end+1) = h_r1_raw(a);
            layer_names{end+1} = sprintf('R1 %s 原始(校准前)', labels{a});
        end
    end
    for a = 1:n_ac
        if h_r1_cal(a) ~= 0
            h_all(end+1) = h_r1_cal(a);
            layer_names{end+1} = sprintf('R1 %s 校准后', labels{a});
        end
    end
    for t = 1:length(r1_tracks)
        if h_r1_ukf(t) ~= 0
            h_all(end+1) = h_r1_ukf(t);
            layer_names{end+1} = sprintf('R1 UKF滤波#%d', r1_tracks{t}.id);
        end
    end
    for a = 1:n_ac
        if h_r2_raw(a) ~= 0
            h_all(end+1) = h_r2_raw(a);
            layer_names{end+1} = sprintf('R2 %s 原始(校准前)', labels{a});
        end
    end
    for a = 1:n_ac
        if h_r2_cal(a) ~= 0
            h_all(end+1) = h_r2_cal(a);
            layer_names{end+1} = sprintf('R2 %s 校准后', labels{a});
        end
    end
    for t = 1:length(r2_tracks)
        if h_r2_ukf(t) ~= 0
            h_all(end+1) = h_r2_ukf(t);
            layer_names{end+1} = sprintf('R2 UKF滤波#%d', r2_tracks{t}.id);
        end
    end

    % 创建复选框
    n_layers = length(layer_names);
    cb_handles = gobjects(1, n_layers);
    for i = 1:n_layers
        ypos = 0.92 - (i-1) * 0.035;
        cb_handles(i) = uicontrol('Parent', panel, 'Style', 'checkbox', ...
            'String', layer_names{i}, 'Value', 1, ...
            'Units', 'normalized', 'Position', [0.05, ypos, 0.9, 0.032], ...
            'FontSize', 7, ...
            'Callback', @(src, ~) set(h_all(i), 'Visible', onoff(src.Value)));
    end

    % 一键隐藏/显示按钮
    btn_y = 0.92 - n_layers * 0.035 - 0.02;
    uicontrol('Parent', panel, 'Style', 'pushbutton', ...
        'String', '全部隐藏', ...
        'Units', 'normalized', 'Position', [0.1, btn_y, 0.8, 0.04], ...
        'FontSize', 8, ...
        'Callback', @(src, ~) toggle_all(src, cb_handles, h_all));

    % ---- 底部统计信息 ----
    n_r1_filt = length(r1_tracks);
    n_r2_filt = length(r2_tracks);
    n_r1_active = sum_active(trackSnapshots_R1);
    n_r2_active = sum_active(trackSnapshots_R2);
    uicontrol('Parent', panel, 'Style', 'text', ...
        'Units', 'normalized', 'Position', [0.03, 0.005, 0.94, 0.04], ...
        'String', sprintf('R1:%d航迹(%d帧)  R2:%d航迹(%d帧) | Pd=%.0f%% Pfa=%.3f', ...
        n_r1_filt, n_r1_active, n_r2_filt, n_r2_active, ...
        params.detection_probability*100, params.false_alarm_rate), ...
        'FontSize', 7, 'HorizontalAlignment', 'center');

    saveas(fig, fullfile(out_dir, 'fig4_multi_track_result.png'));
    fprintf('  多目标跟踪综合图已保存: fig4_multi_track_result.png\n');
end

function v = onoff(val)
    if val, v = 'on'; else, v = 'off'; end
end

function toggle_all(btn, cb_handles, h_all)
    % 检查当前状态：如果大部分可见，则全部隐藏；反之全部显示
    if strcmp(btn.String, '全部隐藏')
        new_val = 0;
        btn.String = '全部显示';
    else
        new_val = 1;
        btn.String = '全部隐藏';
    end
    for i = 1:length(cb_handles)
        set(cb_handles(i), 'Value', new_val);
        set(h_all(i), 'Visible', onoff(new_val));
    end
end

% ---- 从点迹列表中按飞机提取原始/校准后地理位置 ----
function [lats_cell, lons_cell] = extract_dets_by_aircraft(detList, mode)
    n_ac = 3;  % 硬编码3架飞机
    lats_cell = cell(1, n_ac);
    lons_cell = cell(1, n_ac);
    for k = 1:length(detList)
        dets = detList{k};
        for d = 1:length(dets)
            dp = dets(d);
            if dp.is_clutter, continue; end
            if ~isfield(dp, 'aircraft_id'), continue; end
            aid = dp.aircraft_id;
            if aid < 1 || aid > n_ac, continue; end

            if strcmp(mode, 'raw') && isfield(dp, 'raw_lat') && ~isnan(dp.raw_lat)
                lats_cell{aid}(end+1) = dp.raw_lat;
                lons_cell{aid}(end+1) = dp.raw_lon;
            elseif strcmp(mode, 'cal') && isfield(dp, 'lat') && ~isnan(dp.lat)
                lats_cell{aid}(end+1) = dp.lat;
                lons_cell{aid}(end+1) = dp.lon;
            end
        end
    end
end

% ---- 收集活跃航迹（非HISTORY）的历史位置 ----
function tracks = collect_active_tracks(snapshots)
    track_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    for k = 1:length(snapshots)
        snap = snapshots{k};
        if isempty(snap.trackList), continue; end
        for t = 1:length(snap.trackList)
            trk = snap.trackList{t};
            if trk.type == 7, continue; end  % skip HISTORY
            tid = trk.id;
            if ~track_map.isKey(tid)
                track_map(tid) = struct('id', tid, 'lat_history', [], ...
                    'lon_history', [], 'life', 0, 'quality', 0, 'n_points', 0);
            end
            rec = track_map(tid);
            rec.lat_history(end+1) = trk.lat;
            rec.lon_history(end+1) = trk.lon;
            rec.life = trk.life;
            rec.quality = trk.quality;
            rec.n_points = rec.n_points + 1;
            track_map(tid) = rec;
        end
    end
    tracks = values(track_map);
end

function n = sum_active(snapshots)
    n = 0;
    for k = 1:length(snapshots)
        snap = snapshots{k};
        if ~isempty(snap.trackList)
            for t = 1:length(snap.trackList)
                if snap.trackList{t}.type ~= 7, n = n + 1; end
            end
        end
    end
end
