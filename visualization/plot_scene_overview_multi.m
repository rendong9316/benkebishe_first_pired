% =========================================================================
% plot_scene_overview_multi.m
% 多目标场景总览 — 显示站点、波束扇区、多飞机真实航迹
% =========================================================================

function plot_scene_overview_multi(true_tracks, labels, params, out_dir)
    colors = {'g', 'm', 'c'};  % A=绿, B=洋红, C=青
    figure('Position', [50, 50, 1400, 750]);
    geoaxes('Basemap', 'landcover');
    hold on;

    % 接收站
    geoplot(params.radar1_lat, params.radar1_lon, 'bs', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'b', 'DisplayName', 'R1');
    geoplot(params.radar2_lat, params.radar2_lon, 'rs', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', 'R2');

    % 发射站
    geoplot(params.radar1_tx_lat, params.radar1_tx_lon, 'b^', ...
        'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Tx1');
    geoplot(params.radar2_tx_lat, params.radar2_tx_lon, 'r^', ...
        'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Tx2');

    % 波束扇区
    draw_beam_sector(params.radar1_lat, params.radar1_lon, ...
        params.radar1_beam_center_deg, params.beam_width_deg, ...
        params.range_min_m, params.range_max_m, [0 0 1]);
    draw_beam_sector(params.radar2_lat, params.radar2_lon, ...
        params.radar2_beam_center_deg, params.beam_width_deg, ...
        params.range_min_m, params.range_max_m, [1 0 0]);

    % 各飞机真实航迹
    for a = 1:length(true_tracks)
        tt = true_tracks{a};
        col = colors{a};
        geoplot(tt(:,2), tt(:,1), '-s', 'Color', col, 'LineWidth', 1.5, ...
            'MarkerSize', 4, 'MarkerFaceColor', col, ...
            'DisplayName', sprintf('飞机%s 真值', labels{a}));
        geoplot(tt(1,2), tt(1,1), 'o', 'Color', col, ...
            'MarkerSize', 8, 'MarkerFaceColor', col, ...
            'DisplayName', sprintf('飞机%s 起点', labels{a}));
        geoplot(tt(end,2), tt(end,1), 'x', 'Color', col, ...
            'MarkerSize', 10, 'LineWidth', 2, ...
            'DisplayName', sprintf('飞机%s 终点', labels{a}));
    end

    title(sprintf('多目标双基地雷达仿真场景 (%d架飞机)', length(true_tracks)));
    subtitle(sprintf('Pd=%.0f%%, Pfa=%.3f, dt=%.0fs, 波束15°, %d-%d km', ...
        params.detection_probability*100, params.false_alarm_rate, ...
        params.dt_sec, params.range_min_km, params.range_max_km));
    legend('Location', 'best');

    saveas(gcf, fullfile(out_dir, 'fig1_scene_overview.png'));
    fprintf('  图1 已保存: fig1_scene_overview.png\n');
end

function draw_beam_sector(rx_lat, rx_lon, center_az, width, r_min, r_max, color)
    az_edges = linspace(center_az - width/2, center_az + width/2, 20);
    lats_inner = zeros(1, length(az_edges));
    lons_inner = zeros(1, length(az_edges));
    lats_outer = zeros(1, length(az_edges));
    lons_outer = zeros(1, length(az_edges));
    for i = 1:length(az_edges)
        [lons_inner(i), lats_inner(i)] = sphere_utils_destination_point(...
            rx_lon, rx_lat, r_min, az_edges(i));
        [lons_outer(i), lats_outer(i)] = sphere_utils_destination_point(...
            rx_lon, rx_lat, r_max, az_edges(i));
    end
    geoplot(lats_inner, lons_inner, '--', 'Color', [color 0.5], 'LineWidth', 1);
    geoplot(lats_outer, lons_outer, '--', 'Color', [color 0.5], 'LineWidth', 1);
    for az_edge = [center_az - width/2, center_az + width/2]
        [lon1, lat1] = sphere_utils_destination_point(rx_lon, rx_lat, r_min, az_edge);
        [lon2, lat2] = sphere_utils_destination_point(rx_lon, rx_lat, r_max, az_edge);
        geoplot([lat1 lat2], [lon1 lon2], '-', 'Color', [color 0.3], 'LineWidth', 1);
    end
end
