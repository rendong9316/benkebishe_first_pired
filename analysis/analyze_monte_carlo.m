% =========================================================================
% analyze_monte_carlo.m
% 蒙特卡洛实验数据分析与报告生成
% =========================================================================

function analyze_monte_carlo(datafile)
    if nargin < 1
        % 自动找最新的 MC 结果文件
        files = dir('results/monte_carlo_N*.mat');
        if isempty(files)
            error('未找到 MC 结果文件');
        end
        [~, idx] = max([files.datenum]);
        datafile = fullfile('results', files(idx).name);
    end

    fprintf('加载: %s\n', datafile);
    S = load(datafile);

    N_mc = S.N_mc;
    match_rate = S.match_rate;
    n_matched_arr = S.n_matched_arr;
    fusion_rms = S.fusion_rms;
    fusion_med = S.fusion_med;
    r1_rms = S.r1_rms;
    r1_med = S.r1_med;
    r2_rms = S.r2_rms;
    r2_med = S.r2_med;
    improvement = S.improvement;
    ac_fusion_rms = S.ac_fusion_rms;
    ac_fusion_med = S.ac_fusion_med;
    ac_r1_rms = S.ac_r1_rms;
    ac_r1_med = S.ac_r1_med;
    ac_r2_rms = S.ac_r2_rms;
    ac_r2_med = S.ac_r2_med;
    total_tracks_r1 = S.total_tracks_r1;
    total_tracks_r2 = S.total_tracks_r2;
    best_m_arr = S.best_m_arr;
    method_names = S.method_names;
    elapsed = S.elapsed;

    if ~exist('out_dir', 'dir'), out_dir = 'results'; end

    % 过滤有效数据
    valid = ~isnan(fusion_rms) & ~isinf(fusion_rms);
    n_valid = sum(valid);
    fprintf('有效运行: %d / %d\n', n_valid, N_mc);

    %% ===== Figure 1: 融合RMSE分布直方图 + CDF =====
    figure('Position', [50, 50, 1400, 600]);

    subplot(1, 3, 1);
    histogram(fusion_rms(valid), 50, 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none');
    hold on;
    xline(mean(fusion_rms(valid)), 'r--', 'LineWidth', 2);
    xline(median(fusion_rms(valid)), 'g--', 'LineWidth', 2);
    xlabel('RMSE (km)'); ylabel('频次');
    title(sprintf('融合RMSE分布 (N=%d)\n均值=%.1f 中位=%.1f km', ...
        n_valid, mean(fusion_rms(valid)), median(fusion_rms(valid))));
    legend({'RMSE', '均值', '中位'}, 'FontSize', 7);
    grid on;

    subplot(1, 3, 2);
    [f, x] = ecdf(fusion_rms(valid));
    plot(x, f*100, 'b-', 'LineWidth', 2); hold on;
    [f1, x1] = ecdf(r1_rms(valid));
    plot(x1, f1*100, 'Color', [0 0 0.7], 'LineWidth', 1.5);
    [f2, x2] = ecdf(r2_rms(valid));
    plot(x2, f2*100, 'Color', [0.7 0 0], 'LineWidth', 1.5);
    xlabel('RMSE (km)'); ylabel('累积概率 (%)');
    title('RMSE CDF 对比');
    legend({'融合', 'R1', 'R2'}, 'FontSize', 7, 'Location', 'southeast');
    grid on;

    subplot(1, 3, 3);
    histogram(improvement(valid & ~isinf(improvement)), 50, ...
        'FaceColor', [0.8 0.4 0.2], 'EdgeColor', 'none');
    hold on;
    xline(0, 'k-', 'LineWidth', 1.5);
    xline(mean(improvement(valid & ~isinf(improvement))), 'r--', 'LineWidth', 2);
    xlabel('融合改善 (%)'); ylabel('频次');
    title(sprintf('融合改善分布\n均值=%.1f%%', ...
        mean(improvement(valid & ~isinf(improvement)))));
    grid on;

    saveas(gcf, fullfile(out_dir, 'mc_fig1_rmse_distribution.png'));
    fprintf('  图1 已保存: mc_fig1_rmse_distribution.png\n');

    %% ===== Figure 2: 分飞机 RMSE vs 运行次数 =====
    figure('Position', [50, 50, 1400, 600]);

    aircraft_labels_ = {'A', 'B', 'C'};
    n_ac = size(ac_fusion_rms, 2);
    colors = {[0.2 0.4 0.8], [0.8 0.4 0.2], [0.2 0.8 0.4]};

    for a = 1:n_ac
        subplot(1, 3, a);
        v_f = ac_fusion_rms(:,a); v_f = v_f(~isnan(v_f) & ~isinf(v_f));
        v1 = ac_r1_rms(:,a); v1 = v1(~isnan(v1) & ~isinf(v1));
        v2 = ac_r2_rms(:,a); v2 = v2(~isnan(v2) & ~isinf(v2));

        [f_f, x_f] = ecdf(v_f);
        [f1, x1] = ecdf(v1);
        [f2, x2] = ecdf(v2);

        plot(x_f, f_f*100, '-', 'Color', colors{a}, 'LineWidth', 2); hold on;
        plot(x1, f1*100, '--', 'Color', colors{a}*0.6, 'LineWidth', 1.5);
        plot(x2, f2*100, ':', 'Color', colors{a}*0.3, 'LineWidth', 1.5);

        xlabel('RMSE (km)'); ylabel('累积概率 (%)');
        title(sprintf('飞机%s RMSE CDF\n融合: %.1f±%.1f | R1: %.1f±%.1f | R2: %.1f±%.1f', ...
            aircraft_labels_{a}, mean(v_f), std(v_f), mean(v1), std(v1), mean(v2), std(v2)));
        legend({'融合', 'R1', 'R2'}, 'FontSize', 7, 'Location', 'southeast');
        grid on;
    end

    saveas(gcf, fullfile(out_dir, 'mc_fig2_per_aircraft_cdf.png'));
    fprintf('  图2 已保存: mc_fig2_per_aircraft_cdf.png\n');

    %% ===== Figure 3: 算法选择分布 + 完整匹配率 =====
    figure('Position', [50, 50, 1200, 500]);

    subplot(1, 3, 1);
    method_counts = zeros(1, length(method_names));
    for m = 1:length(method_names)
        method_counts(m) = sum(best_m_arr(valid) == m);
    end
    bar(method_counts);
    set(gca, 'XTickLabel', method_names);
    ylabel('次数');
    title(sprintf('最佳融合算法分布 (N=%d)', n_valid));
    for m = 1:length(method_names)
        text(m, method_counts(m) + 10, sprintf('%.1f%%', method_counts(m)/n_valid*100), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
    end
    grid on;

    subplot(1, 3, 2);
    match_bins = [0, 1, 2, 3, 4];
    h = histogram(n_matched_arr(valid), match_bins, 'FaceColor', [0.3 0.6 0.3]);
    xlabel('匹配对数'); ylabel('频次');
    title(sprintf('匹配对数分布\n均值=%.2f ± %.2f', ...
        mean(n_matched_arr(valid)), std(n_matched_arr(valid))));
    grid on;

    subplot(1, 3, 3);
    p_complete = sum(total_tracks_r1 >= n_ac & total_tracks_r2 >= n_ac) / n_valid;
    bins = 0:6;
    h1 = histogram(total_tracks_r1(valid), bins, 'FaceColor', [0 0 0.7], 'FaceAlpha', 0.5);
    hold on;
    h2 = histogram(total_tracks_r2(valid), bins, 'FaceColor', [0.7 0 0], 'FaceAlpha', 0.5);
    xlabel('活跃航迹数'); ylabel('频次');
    title(sprintf('航迹数分布 (完整匹配=%.0f%%)', p_complete*100));
    legend({'R1', 'R2'}, 'FontSize', 7);
    grid on;

    saveas(gcf, fullfile(out_dir, 'mc_fig3_method_distribution.png'));
    fprintf('  图3 已保存: mc_fig3_method_distribution.png\n');

    %% ===== Figure 4: RMSE vs 改善散点图 (诊断用) =====
    figure('Position', [50, 50, 1200, 500]);

    % 先提取有效数据
    r1_v = r1_rms(valid); r2_v = r2_rms(valid);
    fus_v = fusion_rms(valid); imp_v = improvement(valid);

    subplot(1, 2, 1);
    scatter(r1_v, fus_v, 15, imp_v, 'filled');
    hold on;
    max_val = max([r1_v; fus_v]);
    plot([0 max_val], [0 max_val], 'k--');
    xlabel('R1 RMSE (km)'); ylabel('融合 RMSE (km)');
    title('融合 vs R1 (颜色=改善率)');
    colorbar; grid on;
    axis equal;

    subplot(1, 2, 2);
    idx_good = imp_v > 0;
    n_good = sum(idx_good); n_bad = sum(~idx_good);
    scatter(r1_v(idx_good), fus_v(idx_good), 15, 'g', 'filled', ...
        'DisplayName', '融合改善');
    hold on;
    scatter(r1_v(~idx_good), fus_v(~idx_good), 15, 'r', 'filled', ...
        'DisplayName', '融合劣化');
    plot([0 max_val], [0 max_val], 'k--');
    xlabel('R1 RMSE (km)'); ylabel('融合 RMSE (km)');
    title(sprintf('融合收益条件\n改善=%d次(%.0f%%) 劣化=%d次(%.0f%%)', ...
        n_good, n_good/n_valid*100, n_bad, n_bad/n_valid*100));
    legend('FontSize', 7, 'Location', 'best');
    grid on;
    axis equal;

    saveas(gcf, fullfile(out_dir, 'mc_fig4_scatter_diagnosis.png'));
    fprintf('  图4 已保存: mc_fig4_scatter_diagnosis.png\n');

    %% ===== 统计报告输出 =====
    fprintf('\n');
    fprintf('================================================================================\n');
    fprintf('  双基地外辐射源雷达多目标跟踪融合 — 蒙特卡洛实验统计报告\n');
    fprintf('================================================================================\n');
    fprintf('  实验日期: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf('  数据文件: %s\n', datafile);
    fprintf('  仿真次数: N = %d  (有效 %d)\n', N_mc, n_valid);
    fprintf('  总耗时: %.0f 秒 (%.2f 秒/次)\n\n', elapsed, elapsed/N_mc);

    fprintf('  ┌─────────────────────────────────────────────────────────────┐\n');
    fprintf('  │  第一部分: 匹配性能                                         │\n');
    fprintf('  ├─────────────────────────────────────────────────────────────┤\n');
    fprintf('  │  正确匹配率 (均值±std)        │  %5.1f%% ± %4.1f%%              │\n', ...
        mean(match_rate(valid))*100, std(match_rate(valid))*100);
    fprintf('  │  完整3对匹配率                │  %d/%d (%.0f%%)                 │\n', ...
        sum(total_tracks_r1>=n_ac & total_tracks_r2>=n_ac), n_valid, p_complete*100);
    fprintf('  │  平均匹配对数                  │  %.2f ± %.2f                  │\n', ...
        mean(n_matched_arr(valid)), std(n_matched_arr(valid)));
    fprintf('  │  R1 平均活跃航迹数            │  %.2f ± %.2f                  │\n', ...
        mean(total_tracks_r1(valid)), std(total_tracks_r1(valid)));
    fprintf('  │  R2 平均活跃航迹数            │  %.2f ± %.2f                  │\n', ...
        mean(total_tracks_r2(valid)), std(total_tracks_r2(valid)));
    fprintf('  └─────────────────────────────────────────────────────────────┘\n\n');

    fprintf('  ┌─────────────────────────────────────────────────────────────┐\n');
    fprintf('  │  第二部分: 跟踪与融合性能 (总体)                             │\n');
    fprintf('  ├─────────────────────────────────────────────────────────────┤\n');

    function print_row(label, vals, unit)
        v = vals(valid & ~isnan(vals) & ~isinf(vals));
        fprintf('  │  %-30s │  %6.1f ± %5.1f %-3s [%5.1f, %5.1f] │\n', ...
            label, mean(v), std(v), unit, prctile(v,5), prctile(v,95));
    end

    print_row('融合 RMSE', fusion_rms, 'km');
    print_row('融合 中位误差', fusion_med, 'km');
    print_row('R1 RMSE', r1_rms, 'km');
    print_row('R1 中位误差', r1_med, 'km');
    print_row('R2 RMSE', r2_rms, 'km');
    print_row('R2 中位误差', r2_med, 'km');
    fprintf('  │  %-30s │  %6.1f ± %5.1f %-3s [%5.1f, %5.1f] │\n', ...
        '融合改善率', mean(improvement(valid)), std(improvement(valid)), '%', ...
        prctile(improvement(valid),5), prctile(improvement(valid),95));
    fprintf('  └─────────────────────────────────────────────────────────────┘\n\n');

    fprintf('  ┌─────────────────────────────────────────────────────────────┐\n');
    fprintf('  │  第三部分: 分飞机性能 (最佳融合算法, omitnan)                 │\n');
    fprintf('  ├─────────────────────────────────────────────────────────────┤\n');
    for a = 1:n_ac
        v_f = ac_fusion_rms(:,a); v_f = v_f(~isnan(v_f) & ~isinf(v_f));
        v1 = ac_r1_rms(:,a); v1 = v1(~isnan(v1) & ~isinf(v1));
        v2 = ac_r2_rms(:,a); v2 = v2(~isnan(v2) & ~isinf(v2));
        fprintf('  │  飞机%s                                        有效样本: %d  │\n', ...
            aircraft_labels_{a}, length(v_f));
        fprintf('  │    融合 RMSE:  %6.1f ± %5.1f km  [%5.1f, %5.1f]            │\n', ...
            mean(v_f), std(v_f), prctile(v_f,5), prctile(v_f,95));
        fprintf('  │    R1   RMSE:  %6.1f ± %5.1f km  [%5.1f, %5.1f]            │\n', ...
            mean(v1), std(v1), prctile(v1,5), prctile(v1,95));
        fprintf('  │    R2   RMSE:  %6.1f ± %5.1f km  [%5.1f, %5.1f]            │\n', ...
            mean(v2), std(v2), prctile(v2,5), prctile(v2,95));
    end
    fprintf('  └─────────────────────────────────────────────────────────────┘\n\n');

    fprintf('  ┌─────────────────────────────────────────────────────────────┐\n');
    fprintf('  │  第四部分: 融合算法选择分布                                  │\n');
    fprintf('  ├─────────────────────────────────────────────────────────────┤\n');
    for m = 1:length(method_names)
        fprintf('  │  %-4s  │  %4d 次 (%.1f%%)                                      │\n', ...
            method_names{m}, method_counts(m), method_counts(m)/n_valid*100);
    end
    fprintf('  └─────────────────────────────────────────────────────────────┘\n\n');

    % 融合正向收益的条件分析
    imp_v = improvement(valid); r1_v = r1_rms(valid); r2_v = r2_rms(valid);
    idx_good = imp_v > 5;
    idx_bad = imp_v < -5;
    idx_neutral = abs(imp_v) <= 5;
    n_v = length(imp_v);

    fprintf('  ┌─────────────────────────────────────────────────────────────┐\n');
    fprintf('  │  第五部分: 融合收益条件分析                                  │\n');
    fprintf('  ├─────────────────────────────────────────────────────────────┤\n');
    fprintf('  │  融合显著改善 (>5%%):  %d 次 (%.0f%%)                          │\n', ...
        sum(idx_good), sum(idx_good)/n_v*100);
    fprintf('  │  融合显著劣化 (<-5%%): %d 次 (%.0f%%)                          │\n', ...
        sum(idx_bad), sum(idx_bad)/n_v*100);
    fprintf('  │  融合无明显变化:      %d 次 (%.0f%%)                          │\n', ...
        sum(idx_neutral), sum(idx_neutral)/n_v*100);

    if sum(idx_good) > 0
        fprintf('  │  改善时 R1-R2 精度差: %.1f km (R2差于R1)                     │\n', ...
            mean(r2_v(idx_good) - r1_v(idx_good)));
    end
    if sum(idx_bad) > 0
        fprintf('  │  劣化时 R1-R2 精度差: %.1f km                                │\n', ...
            mean(r2_v(idx_bad) - r1_v(idx_bad)));
    end
    fprintf('  └─────────────────────────────────────────────────────────────┘\n\n');

    % 保存统计结果
    stats = struct();
    stats.N_mc = N_mc;
    stats.n_valid = n_valid;
    stats.complete_match_rate = p_complete;
    stats.match_rate_mean = mean(match_rate(valid));
    stats.match_rate_std = std(match_rate(valid));
    stats.fusion_rmse_mean = mean(fusion_rms(valid));
    stats.fusion_rmse_std = std(fusion_rms(valid));
    stats.fusion_med_mean = mean(fusion_med(valid));
    stats.fusion_med_std = std(fusion_med(valid));
    stats.r1_rmse_mean = mean(r1_rms(valid));
    stats.r1_rmse_std = std(r1_rms(valid));
    stats.r2_rmse_mean = mean(r2_rms(valid));
    stats.r2_rmse_std = std(r2_rms(valid));
    stats.improvement_mean = mean(improvement(valid));
    stats.improvement_std = std(improvement(valid));
    stats.improvement_p5 = prctile(improvement(valid), 5);
    stats.improvement_p95 = prctile(improvement(valid), 95);
    stats.method_distribution = method_counts;
    stats.p_good = sum(idx_good) / n_v;
    stats.p_bad = sum(idx_bad) / n_v;
    stats.p_neutral = sum(idx_neutral) / n_v;
    stats.elapsed = elapsed;

    statfile = fullfile(out_dir, sprintf('mc_analysis_%s.mat', datestr(now, 'yyyymmdd_HHMMSS')));
    save(statfile, 'stats');
    fprintf('  统计数据已保存: %s\n', statfile);
    fprintf('Done.\n');
end
