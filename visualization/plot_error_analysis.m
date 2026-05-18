% =========================================================================
% plot_error_analysis.m
% 定量误差分析图 — RMSE时间线 + 误差CDF分布 + 逐飞机汇总表
% =========================================================================

function plot_error_analysis(es_r1, es_r2, out_dir)
    n_ac = length(es_r1.summary);
    ac_colors = {[0 0.6 0], [0.8 0 0.8], [0 0.7 0.7]};  % A=绿, B=紫, C=青
    ac_labels = {'A', 'B', 'C'};

    fig = figure('Position', [50, 50, 1400, 750]);

    % ---- 子图1: R1误差时间线 ----
    ax1 = subplot(2, 3, 1);
    hold on;
    for a = 1:n_ac
        n_pts = length(es_r1.ukf_errors_km{a});
        if n_pts > 0
            plot(1:n_pts, movmean(es_r1.ukf_errors_km{a}, 5), '-', ...
                'Color', ac_colors{a}, 'LineWidth', 1.5);
        end
    end
    yline(es_r1.overall.ukf.median, 'b--', 'LineWidth', 1);
    title(sprintf('R1 UKF位置误差 (中位=%.1f km)', es_r1.overall.ukf.median));
    xlabel('帧序号'); ylabel('误差 (km)'); grid on;

    % ---- 子图2: R2误差时间线 ----
    ax2 = subplot(2, 3, 4);
    hold on;
    for a = 1:n_ac
        n_pts = length(es_r2.ukf_errors_km{a});
        if n_pts > 0
            plot(1:n_pts, movmean(es_r2.ukf_errors_km{a}, 5), '-', ...
                'Color', ac_colors{a}, 'LineWidth', 1.5);
        end
    end
    yline(es_r2.overall.ukf.median, 'r--', 'LineWidth', 1);
    title(sprintf('R2 UKF位置误差 (中位=%.1f km)', es_r2.overall.ukf.median));
    xlabel('帧序号'); ylabel('误差 (km)'); grid on;

    % ---- 子图3: R1误差CDF ----
    ax3 = subplot(2, 3, 2);
    hold on;
    x_max = 0;
    for a = 1:n_ac
        errs = es_r1.ukf_errors_km{a};
        if length(errs) < 3, continue; end
        [f, x] = ecdf(errs);
        plot(x, f*100, '-', 'Color', ac_colors{a}, 'LineWidth', 1.5);
        x_max = max(x_max, max(x));
    end
    % 叠加检测误差CDF（灰色虚线）
    all_det = []; for a = 1:n_ac, all_det = [all_det, es_r1.det_errors_km{a}]; end
    [f_d, x_d] = ecdf(all_det);
    plot(x_d, f_d*100, 'k--', 'LineWidth', 1);
    xlim([0, max(x_max, 50)]);
    ylim([0, 100]);
    title('R1 误差累积分布 (CDF)');
    xlabel('误差 (km)'); ylabel('累积概率 (%)'); grid on;
    legend([ac_labels, {'检测'}], 'Location', 'southeast');

    % ---- 子图4: R2误差CDF ----
    ax4 = subplot(2, 3, 5);
    hold on;
    for a = 1:n_ac
        errs = es_r2.ukf_errors_km{a};
        if length(errs) < 3, continue; end
        [f, x] = ecdf(errs);
        plot(x, f*100, '-', 'Color', ac_colors{a}, 'LineWidth', 1.5);
    end
    all_det2 = []; for a = 1:n_ac, all_det2 = [all_det2, es_r2.det_errors_km{a}]; end
    [f_d2, x_d2] = ecdf(all_det2);
    plot(x_d2, f_d2*100, 'k--', 'LineWidth', 1);
    xlim([0, max(x_d2(end), 50)]);
    ylim([0, 100]);
    title('R2 误差累积分布 (CDF)');
    xlabel('误差 (km)'); ylabel('累积概率 (%)'); grid on;

    % ---- 子图3+6: 汇总表 ----
    ax_table = subplot(2, 3, [3, 6]);
    axis off;
    hold on;

    % 构建表格数据
    tbl_data = {};
    row = 1;
    for radar_idx = 1:2
        if radar_idx == 1, es = es_r1; tag = 'R1'; else, es = es_r2; tag = 'R2'; end
        tbl_data{row, 1} = sprintf('=== %s ===', tag);
        row = row + 1;
        tbl_data{row, 1} = '飞机'; tbl_data{row, 2} = '类型';
        tbl_data{row, 3} = '点数'; tbl_data{row, 4} = '中位(km)';
        tbl_data{row, 5} = '均值(km)'; tbl_data{row, 6} = 'RMSE(km)';
        tbl_data{row, 7} = '95%(km)'; tbl_data{row, 8} = 'vs检测';
        row = row + 1;

        for a = 1:n_ac
            s_ukf = es.summary(a).ukf;
            s_det = es.summary(a).det_calibrated;
            imp = es.summary(a).ukf_vs_det_pct;
            tbl_data{row, 1} = ac_labels{a};
            tbl_data{row, 2} = 'UKF';
            tbl_data{row, 3} = sprintf('%d', s_ukf.n);
            tbl_data{row, 4} = sprintf('%.1f', s_ukf.median);
            tbl_data{row, 5} = sprintf('%.1f', s_ukf.mean);
            tbl_data{row, 6} = sprintf('%.1f', s_ukf.rms);
            tbl_data{row, 7} = sprintf('%.1f', s_ukf.pct95);
            tbl_data{row, 8} = sprintf('%.0f%%', imp);
            row = row + 1;

            tbl_data{row, 2} = '检测';
            tbl_data{row, 3} = sprintf('%d', s_det.n);
            tbl_data{row, 4} = sprintf('%.1f', s_det.median);
            tbl_data{row, 5} = sprintf('%.1f', s_det.mean);
            tbl_data{row, 6} = sprintf('%.1f', s_det.rms);
            tbl_data{row, 7} = sprintf('%.1f', s_det.pct95);
            tbl_data{row, 8} = '-';
            row = row + 1;
        end
        % 总体行
        s = es.overall.ukf;
        tbl_data{row, 1} = '总计'; tbl_data{row, 2} = 'UKF';
        tbl_data{row, 3} = sprintf('%d', s.n);
        tbl_data{row, 4} = sprintf('%.1f', s.median);
        tbl_data{row, 5} = sprintf('%.1f', s.mean);
        tbl_data{row, 6} = sprintf('%.1f', s.rms);
        tbl_data{row, 7} = sprintf('%.1f', s.pct95);
        tbl_data{row, 8} = '-';
        row = row + 2;
    end

    t = uitable('Data', tbl_data, 'Units', 'normalized', ...
        'Position', [0.71, 0.05, 0.28, 0.90], ...
        'ColumnWidth', {40, 40, 45, 60, 60, 60, 55, 55}, ...
        'FontSize', 8);

    sgtitle(sprintf('UKF滤波误差分析 | Pd=60%%, Pfa=10^{-3}, \\sigma_{range}=4km, \\sigma_{az}=0.4deg'));
    saveas(fig, fullfile(out_dir, 'fig5_error_analysis.png'));
    fprintf('  误差分析图已保存: fig5_error_analysis.png\n');
end
