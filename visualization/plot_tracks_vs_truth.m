% =========================================================================
% plot_tracks_vs_truth.m
% 图3: UKF滤波航迹 vs 真实航迹 (geoplot)
% =========================================================================

function plot_tracks_vs_truth(trackState_R1, trackState_R2, true_track, params, out_dir)
    figure('Position', [50, 50, 1400, 750]);

    % ---- R1 ----
    subplot(1, 2, 1);
    geoaxes('Basemap', 'landcover');
    hold on;
    title('R1 UKF滤波航迹');

    plot_track_on_map(trackState_R1, true_track, params.radar1_lat, params.radar1_lon);

    % ---- R2 ----
    subplot(1, 2, 2);
    geoaxes('Basemap', 'landcover');
    hold on;
    title('R2 UKF滤波航迹');

    plot_track_on_map(trackState_R2, true_track, params.radar2_lat, params.radar2_lon);

    sgtitle('UKF滤波航迹 vs 真实航迹');
    saveas(gcf, fullfile(out_dir, 'fig3_tracks_vs_truth.png'));
    fprintf('  图3 已保存: fig3_tracks_vs_truth.png\n');
end

function plot_track_on_map(stateList, true_track, rx_lat, rx_lon)
    % 真实航迹
    geoplot(true_track(:,2), true_track(:,1), 'k--', 'LineWidth', 1.5, ...
        'DisplayName', '真实航迹');

    % 滤波航迹
    lats = []; lons = [];
    for k = 1:length(stateList)
        s = stateList{k};
        if isempty(s) || ~isfield(s, 'lat') || isnan(s.lat), continue; end
        lats(end+1) = s.lat;
        lons(end+1) = s.lon;
    end
    if ~isempty(lats)
        geoplot(lats, lons, 'b-', 'LineWidth', 1.5, 'DisplayName', 'UKF滤波');
    end

    % 接收站
    geoplot(rx_lat, rx_lon, 'rs', 'MarkerSize', 10, ...
        'MarkerFaceColor', 'r', 'DisplayName', '接收站');

    legend('Location', 'best');
end
