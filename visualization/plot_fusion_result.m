% =========================================================================
% plot_fusion_result.m
% 融合航迹可视化: 地图叠加 + 误差收敛曲线 + 算法对比表
% 所有图层支持复选框独立显隐控制
% =========================================================================

function plot_fusion_result(true_tracks, labels, trackSnapshots_R1, trackSnapshots_R2, ...
        all_fused_snapshots, method_names, matched_pairs, fusion_eval, truthTrajs, params, out_dir)

    ac_colors = {'g', 'm', 'c'};
    fusion_colors_ac = {[0 0.5 0], [0.5 0 0.5], [0 0.5 0.5]};
    n_ac = length(true_tracks);
    n_methods = length(method_names);

    % 找出最佳融合方法 (RMSE最小)
    best_idx = 1;
    best_rmse = inf;
    for m = 1:n_methods
        if fusion_eval.overall(m).s.rms < best_rmse
            best_rmse = fusion_eval.overall(m).s.rms;
            best_idx = m;
        end
    end

    % 收集所有航迹数据（用于绘图）
    r1_tracks = collect_track_positions(trackSnapshots_R1);
    r2_tracks = collect_track_positions(trackSnapshots_R2);
    best_snaps = all_fused_snapshots{best_idx};
    fused_tracks = collect_fused_positions(best_snaps, matched_pairs, fusion_eval.pair_to_aircraft);

    %% ===== Figure 1: 地图叠加 (带图层复选框) =====
    fig1 = figure('Position', [50, 50, 1400, 750]);
    ax = geoaxes('Units', 'normalized', 'Position', [0.04, 0.10, 0.68, 0.88]);
    ax.Basemap = 'landcover';
    hold(ax, 'on');

    h_all = [];      % 图层句柄数组
    layer_names = {}; % 图层名称

    % ---- Layer: 真值 (虚线) ----
    h_truth = gobjects(1, n_ac);
    for a = 1:n_ac
        tt = true_tracks{a};
        h_truth(a) = geoplot(ax, tt(:,2), tt(:,1), '--s', 'Color', ac_colors{a}, ...
            'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', ac_colors{a}, ...
            'DisplayName', sprintf('%s 真值', labels{a}));
        h_all(end+1) = h_truth(a);
        layer_names{end+1} = sprintf('%s 真值', labels{a});
    end

    % ---- Layer: R1 UKF航迹 (蓝色) ----
    h_r1_ukf = gobjects(1, length(r1_tracks));
    for t = 1:length(r1_tracks)
        trk = r1_tracks{t};
        if length(trk.lat_history) > 2
            h_r1_ukf(t) = geoplot(ax, trk.lat_history, trk.lon_history, 'b-o', ...
                'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'b', ...
                'DisplayName', sprintf('R1 UKF#%d', trk.id));
            h_all(end+1) = h_r1_ukf(t);
            layer_names{end+1} = sprintf('R1 UKF航迹#%d', trk.id);
        end
    end

    % ---- Layer: R2 UKF航迹 (红色) ----
    h_r2_ukf = gobjects(1, length(r2_tracks));
    for t = 1:length(r2_tracks)
        trk = r2_tracks{t};
        if length(trk.lat_history) > 2
            h_r2_ukf(t) = geoplot(ax, trk.lat_history, trk.lon_history, 'r-^', ...
                'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'r', ...
                'DisplayName', sprintf('R2 UKF#%d', trk.id));
            h_all(end+1) = h_r2_ukf(t);
            layer_names{end+1} = sprintf('R2 UKF航迹#%d', trk.id);
        end
    end

    % ---- Layer: 四种融合航迹 (粗实线) ----
    h_fused_methods = cell(n_methods, 1);
    method_colors_map = {[0 0.5 0], [0.8 0.4 0], [0 0 0.8], [0.6 0 0.6]};  % SCC/BC/CI/FCI
    for m = 1:n_methods
        h_fused_methods{m} = gobjects(1, n_ac);
        snaps_m = all_fused_snapshots{m};
        fused_m = collect_fused_positions(snaps_m, matched_pairs, fusion_eval.pair_to_aircraft);
        for t = 1:length(fused_m)
            ft = fused_m{t};
            if length(ft.lat_history) > 2 && ft.aircraft > 0 && ft.aircraft <= n_ac
                lw = 3.0;
                if m == best_idx, lw = 3.5; end
                ac_label = labels{ft.aircraft};
                h_fused_methods{m}(ft.aircraft) = geoplot(ax, ft.lat_history, ft.lon_history, '-d', ...
                    'Color', method_colors_map{m}, 'LineWidth', lw, ...
                    'MarkerSize', 5, 'MarkerFaceColor', method_colors_map{m}, ...
                    'DisplayName', sprintf('%s-%s', method_names{m}, ac_label));
                h_all(end+1) = h_fused_methods{m}(ft.aircraft);
                layer_names{end+1} = sprintf('%s 融合-%s', method_names{m}, ac_label);
            end
        end
    end

    % ---- 站点标记 (始终可见, 不加入图层控制) ----
    geoplot(ax, params.radar1_lat, params.radar1_lon, 'bs', ...
        'MarkerSize', 14, 'MarkerFaceColor', 'b', 'DisplayName', 'R1');
    geoplot(ax, params.radar2_lat, params.radar2_lon, 'rs', ...
        'MarkerSize', 14, 'MarkerFaceColor', 'r', 'DisplayName', 'R2');
    geoplot(ax, params.radar1_tx_lat, params.radar1_tx_lon, 'b^', ...
        'MarkerSize', 10, 'DisplayName', 'Tx1');
    geoplot(ax, params.radar2_tx_lat, params.radar2_tx_lon, 'r^', ...
        'MarkerSize', 10, 'DisplayName', 'Tx2');

    % ---- 右侧控制面板 ----
    panel = uipanel('Units', 'normalized', 'Position', [0.74, 0.04, 0.24, 0.94], ...
        'Title', '图层显隐控制', 'FontSize', 11);

    n_layers = length(layer_names);
    cb_handles = gobjects(1, n_layers);
    row_height = 0.028;
    max_rows = 30;
    for i = 1:n_layers
        ypos = 0.92 - (i-1) * row_height;
        if ypos < 0.05, break; end
        cb_handles(i) = uicontrol('Parent', panel, 'Style', 'checkbox', ...
            'String', layer_names{i}, 'Value', 1, ...
            'Units', 'normalized', 'Position', [0.05, ypos, 0.9, row_height-0.002], ...
            'FontSize', 6.5, ...
            'Callback', @(src, ~) set_layer_visible(h_all(i), src.Value));
    end

    % 一键隐藏/显示按钮
    btn_y = 0.92 - n_layers * row_height - 0.02;
    if btn_y > 0.02
        uicontrol('Parent', panel, 'Style', 'pushbutton', ...
            'String', '全部隐藏', ...
            'Units', 'normalized', 'Position', [0.1, btn_y, 0.8, 0.04], ...
            'FontSize', 8, ...
            'Callback', @(src, ~) toggle_all_layers(src, cb_handles, h_all));
    end

    % 底部统计信息
    uicontrol('Parent', panel, 'Style', 'text', ...
        'Units', 'normalized', 'Position', [0.03, 0.005, 0.94, 0.035], ...
        'String', sprintf('R1:%d | R2:%d | 匹配:%d对 | %s RMSE=%.1fkm', ...
        length(r1_tracks), length(r2_tracks), length(matched_pairs), ...
        method_names{best_idx}, best_rmse), ...
        'FontSize', 7, 'HorizontalAlignment', 'center');

    % 左上角匹配标注
    dim = [0.05, 0.75, 0.25, 0.12];
    match_str = sprintf('匹配: R1↔R2 (对=%d)', length(matched_pairs));
    for p = 1:length(matched_pairs)
        ac = fusion_eval.pair_to_aircraft(p);
        match_str = sprintf('%s\n  R1#%d↔R2#%d→%s', match_str, ...
            matched_pairs(p).R1_track_id, matched_pairs(p).R2_track_id, labels{ac});
    end
    annotation('textbox', dim, 'String', match_str, ...
        'FontSize', 7, 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
        'VerticalAlignment', 'top');

    title(ax, sprintf('双基地雷达航迹融合结果 (%s最优)', method_names{best_idx}));
    legend(ax, 'Location', 'northeastoutside');

    saveas(fig1, fullfile(out_dir, 'fig6_fusion_map.png'));
    fprintf('  融合地图已保存: fig6_fusion_map.png\n');

    %% ===== Figure 2: 误差收敛曲线 (带图层复选框) =====
    fig2 = figure('Position', [50, 50, 1400, 750]);

    frame_times = (0:length(trackSnapshots_R1)-1) * params.dt_sec;
    method_ls = {'-', '--', '-.', ':'};
    method_clr = {[0 0 0], [1 0 0], [0 0 1], [0 0.7 0]};
    win = 10;

    % 收集所有轴和线句柄
    ax_err = gobjects(1, 4);
    all_h_lines = gobjects(0);  % 所有线条句柄
    all_line_names = {};        % 对应名称

    for a = 1:n_ac
        ax_err(a) = subplot(2, 2, a);
        hold on;
        grid on;

        for m = 1:n_methods
            frame_errs = build_frame_errors(all_fused_snapshots{m}, ...
                fusion_eval.pair_to_aircraft, a, truthTrajs{a}, frame_times);
            if length(frame_errs) >= win
                smoothed = movmean(frame_errs, win, 'omitnan');
                h = plot(frame_times, smoothed, 'LineStyle', method_ls{m}, ...
                    'Color', method_clr{m}, 'LineWidth', 2);
                all_h_lines(end+1) = h;
                all_line_names{end+1} = sprintf('飞机%s-%s', labels{a}, method_names{m});
            end
        end

        % R1误差
        r1_fe = build_frame_errors_single(trackSnapshots_R1, ...
            fusion_eval.pair_to_aircraft, a, truthTrajs{a}, frame_times, ...
            matched_pairs, 'R1');
        if length(r1_fe) >= win
            h = plot(frame_times, movmean(r1_fe, win, 'omitnan'), ...
                'LineStyle', ':', 'Color', [0 0 0.7], 'LineWidth', 1.5);
            all_h_lines(end+1) = h;
            all_line_names{end+1} = sprintf('飞机%s-R1 UKF', labels{a});
        end

        % R2误差
        r2_fe = build_frame_errors_single(trackSnapshots_R2, ...
            fusion_eval.pair_to_aircraft, a, truthTrajs{a}, frame_times, ...
            matched_pairs, 'R2');
        if length(r2_fe) >= win
            h = plot(frame_times, movmean(r2_fe, win, 'omitnan'), ...
                'LineStyle', ':', 'Color', [0.7 0 0], 'LineWidth', 1.5);
            all_h_lines(end+1) = h;
            all_line_names{end+1} = sprintf('飞机%s-R2 UKF', labels{a});
        end

        xlabel('时间 (s)'); ylabel('位置误差 (km)');
        title(sprintf('飞机%s 误差收敛曲线', labels{a}));
        legend('FontSize', 6, 'Location', 'best');
    end

    % 总体CDF对比
    ax_err(4) = subplot(2, 2, 4);
    hold on;
    grid on;
    for m = 1:n_methods
        combined = [];
        for a = 1:n_ac
            combined = [combined, fusion_eval.fusion_errors{m, a}];
        end
        if ~isempty(combined)
            [f, x] = ecdf(combined);
            h = plot(x, f*100, 'LineStyle', method_ls{m}, ...
                'Color', method_clr{m}, 'LineWidth', 2);
            all_h_lines(end+1) = h;
            all_line_names{end+1} = sprintf('总体-%s', method_names{m});
        end
    end
    combined_r1 = []; combined_r2 = [];
    for a = 1:n_ac
        combined_r1 = [combined_r1, fusion_eval.r1_errors{a}];
        combined_r2 = [combined_r2, fusion_eval.r2_errors{a}];
    end
    if ~isempty(combined_r1)
        [f, x] = ecdf(combined_r1);
        h = plot(x, f*100, ':', 'Color', [0 0 0.7], 'LineWidth', 1.5);
        all_h_lines(end+1) = h;
        all_line_names{end+1} = '总体-R1 UKF';
    end
    if ~isempty(combined_r2)
        [f, x] = ecdf(combined_r2);
        h = plot(x, f*100, ':', 'Color', [0.7 0 0], 'LineWidth', 1.5);
        all_h_lines(end+1) = h;
        all_line_names{end+1} = '总体-R2 UKF';
    end
    xlabel('位置误差 (km)'); ylabel('累积概率 (%)');
    title('误差CDF对比');
    legend('FontSize', 6, 'Location', 'southeast');

    % ---- Figure 2 右侧复选框面板 ----
    panel2 = uipanel('Units', 'normalized', 'Position', [0.76, 0.04, 0.22, 0.94], ...
        'Title', '曲线显隐控制', 'FontSize', 11);

    n_lines2 = length(all_line_names);
    cb2 = gobjects(1, n_lines2);
    row_h2 = 0.028;
    for i = 1:n_lines2
        ypos = 0.92 - (i-1) * row_h2;
        if ypos < 0.05, break; end
        cb2(i) = uicontrol('Parent', panel2, 'Style', 'checkbox', ...
            'String', all_line_names{i}, 'Value', 1, ...
            'Units', 'normalized', 'Position', [0.05, ypos, 0.9, row_h2-0.002], ...
            'FontSize', 6.5, ...
            'Callback', @(src, ~) set(all_h_lines(i), 'Visible', onoff(src.Value)));
    end

    btn_y2 = 0.92 - n_lines2 * row_h2 - 0.02;
    if btn_y2 > 0.02
        uicontrol('Parent', panel2, 'Style', 'pushbutton', ...
            'String', '全部隐藏', ...
            'Units', 'normalized', 'Position', [0.1, btn_y2, 0.8, 0.04], ...
            'FontSize', 8, ...
            'Callback', @(src, ~) toggle_all_lines(src, cb2, all_h_lines));
    end

    saveas(fig2, fullfile(out_dir, 'fig7_fusion_error.png'));
    fprintf('  误差收敛曲线已保存: fig7_fusion_error.png\n');

    %% ===== Figure 3: 算法性能对比柱状图 (不含复选框，结构简单) =====
    fig3 = figure('Position', [50, 50, 1200, 500]);

    all_method_labels = [method_names, {'R1_only', 'R2_only'}];
    n_rows = length(all_method_labels);
    rmse_vals = zeros(n_rows, n_ac + 1);
    median_vals = zeros(n_rows, n_ac + 1);
    for m = 1:n_rows
        for a = 1:n_ac
            idx = (a-1)*n_rows + m;
            if idx <= size(fusion_eval.summary, 1)
                s = fusion_eval.summary(idx).s;
                rmse_vals(m, a) = s.rms;
                median_vals(m, a) = s.median;
            end
        end
        rmse_vals(m, n_ac+1) = fusion_eval.overall(m).s.rms;
        median_vals(m, n_ac+1) = fusion_eval.overall(m).s.median;
    end

    ac_legends = cell(1, n_ac + 1);
    for a = 1:n_ac
        ac_legends{a} = sprintf('飞机%s', labels{a});
    end
    ac_legends{n_ac+1} = 'Overall';

    subplot(1, 2, 1);
    b = bar(rmse_vals);
    set(gca, 'XTickLabel', all_method_labels);
    ylabel('RMSE (km)');
    title('各算法RMSE对比');
    legend(b, ac_legends, 'FontSize', 7, 'Location', 'best');
    grid on;
    xtickangle(45);

    subplot(1, 2, 2);
    b2 = bar(median_vals);
    set(gca, 'XTickLabel', all_method_labels);
    ylabel('中位误差 (km)');
    title('各算法中位误差对比');
    legend(b2, ac_legends, 'FontSize', 7, 'Location', 'best');
    grid on;
    xtickangle(45);

    saveas(fig3, fullfile(out_dir, 'fig8_fusion_comparison.png'));
    fprintf('  融合算法对比图已保存: fig8_fusion_comparison.png\n');
end

% =========================================================================
% 图层显隐工具函数
% =========================================================================

function v = onoff(val)
    if val, v = 'on'; else, v = 'off'; end
end

function set_layer_visible(h, val)
    try
        set(h, 'Visible', onoff(val));
    catch
        % GraphicsPlaceholder or deleted handle
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

function toggle_all_lines(btn, cb_handles, h_all)
    if strcmp(btn.String, '全部隐藏')
        new_val = 0; btn.String = '全部显示';
    else
        new_val = 1; btn.String = '全部隐藏';
    end
    for i = 1:length(cb_handles)
        if cb_handles(i) ~= 0
            set(cb_handles(i), 'Value', new_val);
            set(h_all(i), 'Visible', onoff(new_val));
        end
    end
end

% =========================================================================
% 数据提取辅助函数
% =========================================================================

function tracks = collect_track_positions(snapshots)
    track_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    for k = 1:length(snapshots)
        snap = snapshots{k};
        if isempty(snap) || ~isfield(snap, 'trackList'), continue; end
        for t = 1:length(snap.trackList)
            trk = snap.trackList{t};
            if isfield(trk, 'type') && trk.type == 7, continue; end
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

function tracks = collect_fused_positions(snapshots, matched_pairs, pair_to_ac)
    track_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    for k = 1:length(snapshots)
        snap = snapshots{k};
        if isempty(snap) || ~isfield(snap, 'trackList'), continue; end
        for t = 1:length(snap.trackList)
            ft = snap.trackList{t};
            pid = ft.id;
            if ~track_map.isKey(pid)
                ac = 0;
                if pid <= length(pair_to_ac), ac = pair_to_ac(pid); end
                track_map(pid) = struct('id', pid, 'lat_history', [], ...
                    'lon_history', [], 'aircraft', ac);
            end
            rec = track_map(pid);
            rec.lat_history(end+1) = ft.lat;
            rec.lon_history(end+1) = ft.lon;
            track_map(pid) = rec;
        end
    end
    tracks = values(track_map);
end

function frame_errs = build_frame_errors(fused_snaps, pair_to_ac, ac, truth, frame_times)
    n_frames = length(fused_snaps);
    frame_errs = nan(1, n_frames);
    for k = 1:n_frames
        t_lat = interp1(truth.time_sec, truth.lat, frame_times(k), 'linear', 'extrap');
        t_lon = interp1(truth.time_sec, truth.lon, frame_times(k), 'linear', 'extrap');
        if isnan(t_lat), continue; end
        snap = fused_snaps{k};
        if isempty(snap) || ~isfield(snap, 'trackList'), continue; end
        best_d = inf;
        for t = 1:length(snap.trackList)
            ft = snap.trackList{t};
            if ft.id > length(pair_to_ac) || pair_to_ac(ft.id) ~= ac, continue; end
            d = haversine_km(ft.lon, ft.lat, t_lon, t_lat);
            if d < best_d, best_d = d; end
        end
        if best_d < inf, frame_errs(k) = best_d; end
    end
end

function frame_errs = build_frame_errors_single(snapshots, pair_to_ac, ac, truth, ...
        frame_times, matched_pairs, which)
    n_frames = length(snapshots);
    frame_errs = nan(1, n_frames);
    target_id = 0;
    for p = 1:length(matched_pairs)
        if p <= length(pair_to_ac) && pair_to_ac(p) == ac
            if strcmp(which, 'R1')
                target_id = matched_pairs(p).R1_track_id;
            else
                target_id = matched_pairs(p).R2_track_id;
            end
            break;
        end
    end
    if target_id == 0, return; end
    for k = 1:n_frames
        t_lat = interp1(truth.time_sec, truth.lat, frame_times(k), 'linear', 'extrap');
        t_lon = interp1(truth.time_sec, truth.lon, frame_times(k), 'linear', 'extrap');
        if isnan(t_lat), continue; end
        snap = snapshots{k};
        if isempty(snap) || ~isfield(snap, 'trackList'), continue; end
        for t = 1:length(snap.trackList)
            trk = snap.trackList{t};
            if trk.id == target_id
                frame_errs(k) = haversine_km(trk.lon, trk.lat, t_lon, t_lat);
                break;
            end
        end
    end
end

function d = haversine_km(lon1, lat1, lon2, lat2)
    R = 6371;
    dlat = deg2rad(lat2 - lat1);
    dlon = deg2rad(lon2 - lon1);
    a = sin(dlat/2)^2 + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dlon/2)^2;
    a = max(0, min(1, a));
    d = R * 2 * atan2(sqrt(a), sqrt(1 - a));
end
