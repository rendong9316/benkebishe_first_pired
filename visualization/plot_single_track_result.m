% =========================================================================
% plot_single_track_result.m
% 单目标跟踪综合图 — 真值 + 原始/校准点迹 + UKF滤波航迹，带图层复选框
% =========================================================================

function plot_single_track_result(true_track, detList_R1, detList_R2, ...
        trackSnapshots_R1, trackSnapshots_R2, params, out_dir)

    fig = figure('Position', [50, 50, 1400, 750]);
    ax = geoaxes('Units', 'normalized', 'Position', [0.04, 0.10, 0.70, 0.88]);
    ax.Basemap = 'landcover';
    hold(ax, 'on');

    h_all = [];
    layer_names = {};

    % ---- 真值 (绿色虚线+方块) ----
    h_truth = geoplot(ax, true_track(:,2), true_track(:,1), '--s', ...
        'Color', 'g', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'g', ...
        'DisplayName', '真值');
    h_all(end+1) = h_truth;
    layer_names{end+1} = '真值航迹';

    % ---- R1 原始点迹 ----
    [r1_raw_lat, r1_raw_lon] = extract_dets(detList_R1, 'raw');
    if ~isempty(r1_raw_lat)
        h = geoplot(ax, r1_raw_lat, r1_raw_lon, '--o', ...
            'Color', [0.4 0.6 1.0], 'LineWidth', 1.0, 'MarkerSize', 4, ...
            'MarkerFaceColor', [0.4 0.6 1.0], 'DisplayName', 'R1 原始点迹');
        h_all(end+1) = h;
        layer_names{end+1} = 'R1 原始点迹(校准前)';
    end

    % ---- R1 校准点迹 ----
    [r1_cal_lat, r1_cal_lon] = extract_dets(detList_R1, 'cal');
    if ~isempty(r1_cal_lat)
        h = geoplot(ax, r1_cal_lat, r1_cal_lon, '-o', ...
            'Color', [0.0 0.4 1.0], 'LineWidth', 1.2, 'MarkerSize', 5, ...
            'MarkerFaceColor', 'b', 'DisplayName', 'R1 校准点迹');
        h_all(end+1) = h;
        layer_names{end+1} = 'R1 校准后点迹';
    end

    % ---- R1 UKF航迹 ----
    r1_tracks = collect_active_tracks(trackSnapshots_R1);
    for t = 1:length(r1_tracks)
        trk = r1_tracks{t};
        if length(trk.lat_history) > 2
            h = geoplot(ax, trk.lat_history, trk.lon_history, 'b-o', ...
                'LineWidth', 2.0, 'MarkerSize', 5, 'MarkerFaceColor', 'b', ...
                'DisplayName', sprintf('R1 UKF#%d', trk.id));
            h_all(end+1) = h;
            layer_names{end+1} = sprintf('R1 UKF航迹#%d', trk.id);
        end
    end

    % ---- R2 原始点迹 ----
    [r2_raw_lat, r2_raw_lon] = extract_dets(detList_R2, 'raw');
    if ~isempty(r2_raw_lat)
        h = geoplot(ax, r2_raw_lat, r2_raw_lon, '--o', ...
            'Color', [1.0 0.6 0.6], 'LineWidth', 1.0, 'MarkerSize', 4, ...
            'MarkerFaceColor', [1.0 0.6 0.6], 'DisplayName', 'R2 原始点迹');
        h_all(end+1) = h;
        layer_names{end+1} = 'R2 原始点迹(校准前)';
    end

    % ---- R2 校准点迹 ----
    [r2_cal_lat, r2_cal_lon] = extract_dets(detList_R2, 'cal');
    if ~isempty(r2_cal_lat)
        h = geoplot(ax, r2_cal_lat, r2_cal_lon, '-o', ...
            'Color', [1.0 0.2 0.2], 'LineWidth', 1.2, 'MarkerSize', 5, ...
            'MarkerFaceColor', 'r', 'DisplayName', 'R2 校准点迹');
        h_all(end+1) = h;
        layer_names{end+1} = 'R2 校准后点迹';
    end

    % ---- R2 UKF航迹 ----
    r2_tracks = collect_active_tracks(trackSnapshots_R2);
    for t = 1:length(r2_tracks)
        trk = r2_tracks{t};
        if length(trk.lat_history) > 2
            h = geoplot(ax, trk.lat_history, trk.lon_history, 'r-^', ...
                'LineWidth', 2.0, 'MarkerSize', 5, 'MarkerFaceColor', 'r', ...
                'DisplayName', sprintf('R2 UKF#%d', trk.id));
            h_all(end+1) = h;
            layer_names{end+1} = sprintf('R2 UKF航迹#%d', trk.id);
        end
    end

    % ---- 站点标记 ----
    geoplot(ax, params.radar1_lat, params.radar1_lon, 'bs', ...
        'MarkerSize', 14, 'MarkerFaceColor', 'b', 'DisplayName', 'R1');
    geoplot(ax, params.radar2_lat, params.radar2_lon, 'rs', ...
        'MarkerSize', 14, 'MarkerFaceColor', 'r', 'DisplayName', 'R2');
    geoplot(ax, params.radar1_tx_lat, params.radar1_tx_lon, 'b^', ...
        'MarkerSize', 10, 'DisplayName', 'Tx1');
    geoplot(ax, params.radar2_tx_lat, params.radar2_tx_lon, 'r^', ...
        'MarkerSize', 10, 'DisplayName', 'Tx2');

    title(ax, '单目标双基地雷达航迹综合对比');
    legend(ax, 'Location', 'northeastoutside');

    % ---- 右侧控制面板 ----
    panel = uipanel('Units', 'normalized', 'Position', [0.76, 0.04, 0.22, 0.94], ...
        'Title', '图层显隐控制', 'FontSize', 11);

    n_layers = length(layer_names);
    cb_handles = gobjects(1, n_layers);
    row_height = 0.035;
    for i = 1:n_layers
        ypos = 0.92 - (i-1) * row_height;
        if ypos < 0.05, break; end
        cb_handles(i) = uicontrol('Parent', panel, 'Style', 'checkbox', ...
            'String', layer_names{i}, 'Value', 1, ...
            'Units', 'normalized', 'Position', [0.05, ypos, 0.9, row_height-0.002], ...
            'FontSize', 7, ...
            'Callback', @(src, ~) set_layer_visible(h_all(i), src.Value));
    end

    btn_y = 0.92 - n_layers * row_height - 0.02;
    if btn_y > 0.02
        uicontrol('Parent', panel, 'Style', 'pushbutton', ...
            'String', '全部隐藏', ...
            'Units', 'normalized', 'Position', [0.1, btn_y, 0.8, 0.04], ...
            'FontSize', 8, ...
            'Callback', @(src, ~) toggle_all_layers(src, cb_handles, h_all));
    end

    % 底部统计
    uicontrol('Parent', panel, 'Style', 'text', ...
        'Units', 'normalized', 'Position', [0.03, 0.005, 0.94, 0.04], ...
        'String', sprintf('R1:%d航迹  R2:%d航迹 | Pd=%.0f%% Pfa=%.3f', ...
        length(r1_tracks), length(r2_tracks), ...
        params.detection_probability*100, params.false_alarm_rate), ...
        'FontSize', 7, 'HorizontalAlignment', 'center');

    saveas(fig, fullfile(out_dir, 'fig4_single_track_result.png'));
    fprintf('  单目标跟踪综合图已保存: fig4_single_track_result.png\n');
end

function [lats, lons] = extract_dets(detList, mode)
    lats = []; lons = [];
    for k = 1:length(detList)
        dets = detList{k};
        for d = 1:length(dets)
            dp = dets(d);
            if dp.is_clutter, continue; end
            if strcmp(mode, 'raw') && isfield(dp, 'raw_lat') && ~isnan(dp.raw_lat)
                lats(end+1) = dp.raw_lat;
                lons(end+1) = dp.raw_lon;
            elseif strcmp(mode, 'cal') && isfield(dp, 'lat') && ~isnan(dp.lat)
                lats(end+1) = dp.lat;
                lons(end+1) = dp.lon;
            end
        end
    end
end

function tracks = collect_active_tracks(snapshots)
    track_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    for k = 1:length(snapshots)
        snap = snapshots{k};
        if isempty(snap.trackList), continue; end
        for t = 1:length(snap.trackList)
            trk = snap.trackList{t};
            if trk.type == 7, continue; end
            tid = trk.id;
            if ~track_map.isKey(tid)
                track_map(tid) = struct('id', tid, 'lat_history', [], 'lon_history', []);
            end
            rec = track_map(tid);
            rec.lat_history(end+1) = trk.lat;
            rec.lon_history(end+1) = trk.lon;
            track_map(tid) = rec;
        end
    end
    tracks = values(track_map);
end

function v = onoff(val)
    if val, v = 'on'; else, v = 'off'; end
end

function set_layer_visible(h, val)
    try
        set(h, 'Visible', onoff(val));
    catch
    end
end

function toggle_all_layers(btn, cb_handles, h_all)
    if strcmp(btn.String, '全部隐藏')
        new_val = 0; btn.String = '全部显示';
    else
        new_val = 1; btn.String = '全部隐藏';
    end
    for i = 1:length(cb_handles)
        if cb_handles(i) ~= 0
            set(cb_handles(i), 'Value', new_val);
        end
        try
            set(h_all(i), 'Visible', onoff(new_val));
        catch
        end
    end
end
