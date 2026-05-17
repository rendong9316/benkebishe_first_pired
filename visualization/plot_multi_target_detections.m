% =========================================================================
% plot_multi_target_detections.m
% 多目标点迹图 — 真实航迹 + 各飞机检测点迹（按aircraft_id着色）
% =========================================================================

function plot_multi_target_detections(true_tracks, labels, detList_R1, detList_R2, ...
        params, out_dir)

    ac_colors = {'g', 'm', 'c'};  % A=绿, B=洋红, C=青

    fig = figure('Position', [50, 50, 1400, 750]);

    % ---- 左侧地理图 ----
    ax = geoaxes('Units', 'normalized', 'Position', [0.05, 0.12, 0.70, 0.85]);
    ax.Basemap = 'landcover';
    hold(ax, 'on');

    % 各飞机真值航迹 (虚线)
    for a = 1:length(true_tracks)
        tt = true_tracks{a};
        col = ac_colors{a};
        geoplot(ax, tt(:,2), tt(:,1), '--', 'Color', col, ...
            'LineWidth', 1.5, 'DisplayName', sprintf('飞机%s 真值', labels{a}));
    end

    % R1 校准后关联点迹 (按aircraft_id着色)
    r1_h = [];
    for a = 1:length(true_tracks)
        [lats, lons] = extract_aircraft_dets(detList_R1, a);
        if ~isempty(lats)
            h = geoplot(ax, lats, lons, 'o', 'Color', ac_colors{a}, ...
                'MarkerSize', 4, 'MarkerFaceColor', ac_colors{a}, ...
                'DisplayName', sprintf('R1 飞机%s 检出', labels{a}));
            r1_h = [r1_h, h];
        end
    end

    % R2 校准后关联点迹 (三角标记，与R1区分)
    r2_h = [];
    for a = 1:length(true_tracks)
        [lats, lons] = extract_aircraft_dets(detList_R2, a);
        if ~isempty(lats)
            h = geoplot(ax, lats, lons, '^', 'Color', ac_colors{a}, ...
                'MarkerSize', 4, 'MarkerFaceColor', ac_colors{a}, ...
                'DisplayName', sprintf('R2 飞机%s 检出', labels{a}));
            r2_h = [r2_h, h];
        end
    end

    % 站点标记
    geoplot(ax, params.radar1_lat, params.radar1_lon, 'bs', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'b', 'DisplayName', 'R1');
    geoplot(ax, params.radar2_lat, params.radar2_lon, 'rs', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', 'R2');

    title(ax, '多目标双基地雷达点迹分布');
    legend(ax, 'Location', 'northeast');

    % ---- 右侧统计面板 ----
    panel = uipanel('Units', 'normalized', 'Position', [0.77, 0.08, 0.21, 0.88], ...
        'Title', '检出统计', 'FontSize', 11);

    ypos = 0.85;
    for a = 1:length(true_tracks)
        n_r1 = count_aircraft_dets(detList_R1, a);
        n_r2 = count_aircraft_dets(detList_R2, a);
        uicontrol('Parent', panel, 'Style', 'text', ...
            'Units', 'normalized', 'Position', [0.05, ypos, 0.9, 0.12], ...
            'String', sprintf('飞机%s\nR1检出:%d  R2检出:%d', ...
            labels{a}, n_r1, n_r2), ...
            'FontSize', 9, 'HorizontalAlignment', 'left');
        ypos = ypos - 0.14;
    end

    % 总杂波统计
    n_clut_r1 = count_clutter(detList_R1);
    n_clut_r2 = count_clutter(detList_R2);
    uicontrol('Parent', panel, 'Style', 'text', ...
        'Units', 'normalized', 'Position', [0.05, ypos, 0.9, 0.10], ...
        'String', sprintf('杂波 R1:%d  R2:%d', n_clut_r1, n_clut_r2), ...
        'FontSize', 8, 'HorizontalAlignment', 'left');

    saveas(fig, fullfile(out_dir, 'fig3_multi_target_detections.png'));
    fprintf('  多目标点迹图已保存: fig3_multi_target_detections.png\n');
end

function [lats, lons] = extract_aircraft_dets(detList, aircraft_id)
    lats = []; lons = [];
    for k = 1:length(detList)
        dets = detList{k};
        for d = 1:length(dets)
            dp = dets(d);
            if ~dp.is_clutter && isfield(dp, 'aircraft_id') ...
                    && dp.aircraft_id == aircraft_id ...
                    && isfield(dp, 'lat') && ~isnan(dp.lat)
                lats(end+1) = dp.lat;
                lons(end+1) = dp.lon;
            end
        end
    end
end

function n = count_aircraft_dets(detList, aircraft_id)
    n = 0;
    for k = 1:length(detList)
        dets = detList{k};
        for d = 1:length(dets)
            dp = dets(d);
            if ~dp.is_clutter && isfield(dp, 'aircraft_id') ...
                    && dp.aircraft_id == aircraft_id
                n = n + 1;
            end
        end
    end
end

function n = count_clutter(detList)
    n = 0;
    for k = 1:length(detList)
        dets = detList{k};
        for d = 1:length(dets)
            if dets(d).is_clutter, n = n + 1; end
        end
    end
end
