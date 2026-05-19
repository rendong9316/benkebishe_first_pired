% =========================================================================
% track_matcher.m
% 双门限航迹匹配: 将R1和R2航迹进行全局配对
% =========================================================================
% 算法 (来自开题报告 4.3.2 节):
%   1. 时间对齐: 将R2航迹快照回退到R1时间网格
%   2. 第一门限 T1 (距离): 逐帧计算对齐后位置距离, < T1 计入匹配计数
%   3. 第二门限 T2 (次数): 滑窗内满足T1的帧数 >= T2 确认关联
%   4. 全局配对: 代价矩阵 + 贪心分配 (一R1对一R2)
%   5. 短航迹修正: life < 阈值时放宽T2比例
% =========================================================================

function matcher = track_matcher(trackSnapshots_R1, trackSnapshots_R2, params)
    n_frames = length(trackSnapshots_R1);

    % 航迹级时间对齐: R2航迹由CV模型全状态回退13s到R1时间网格
    % 使用UKF滤波后的2D速度矢量, 优于点迹级的单分量多普勒外推
    aligned_R2 = time_align_tracks(trackSnapshots_R2, params);

    % ---- 提取每条航迹的逐帧位置历史 ----
    [r1_pos, r1_ids, r1_lives] = extract_track_positions(trackSnapshots_R1, n_frames);
    [r2_pos, r2_ids, r2_lives] = extract_track_positions(aligned_R2, n_frames);

    n_r1 = length(r1_ids);
    n_r2 = length(r2_ids);

    fprintf('\n--- 航迹匹配: R1共%d条, R2共%d条活跃航迹 ---\n', n_r1, n_r2);
    for i = 1:n_r1
        fprintf('  R1#%d: life=%d\n', r1_ids(i), r1_lives(i));
    end
    for j = 1:n_r2
        fprintf('  R2#%d: life=%d\n', r2_ids(j), r2_lives(j));
    end

    % ---- 参数 ----
    T1_km = 80;           % 第一门限: 位置距离 (km)
    T2_ratio = 0.40;      % 第二门限: 匹配帧数/共存帧数 最低比例
    T2_abs_min = 6;       % 第二门限: 绝对最低匹配帧数
    short_life_thresh = 20;  % 短航迹阈值
    T2_ratio_short = 0.30;   % 短航迹放宽比例

    % ---- 逐帧计算距离并累积匹配计数 ----
    match_count = zeros(n_r1, n_r2);     % 绝对匹配帧数
    coexist_count = zeros(n_r1, n_r2);   % 共存帧数
    dist_history = cell(n_r1, n_r2);     % 逐帧距离记录
    for i = 1:n_r1
        for j = 1:n_r2
            dist_history{i,j} = [];
        end
    end

    for k = 1:n_frames
        for i = 1:n_r1
            if isnan(r1_pos(i,k,1)), continue; end
            for j = 1:n_r2
                if isnan(r2_pos(j,k,1)), continue; end
                coexist_count(i,j) = coexist_count(i,j) + 1;
                d = haversine_km(r1_pos(i,k,1), r1_pos(i,k,2), ...
                                 r2_pos(j,k,1), r2_pos(j,k,2));
                dist_history{i,j}(end+1) = d;
                if d < T1_km
                    match_count(i,j) = match_count(i,j) + 1;
                end
            end
        end
    end

    % ---- 打印匹配矩阵 ----
    fprintf('\n匹配计数矩阵 (R1行 × R2列):\n');
    fprintf('        ');
    for j = 1:n_r2, fprintf('R2#%-4d ', r2_ids(j)); end
    fprintf('  | coexist\n');
    for i = 1:n_r1
        fprintf('R1#%-3d: ', r1_ids(i));
        for j = 1:n_r2
            fprintf('%3d/%-3d ', match_count(i,j), coexist_count(i,j));
        end
        fprintf('\n');
    end

    % ---- 候选匹配判定 (双门限) ----
    is_candidate = false(n_r1, n_r2);
    match_quality = zeros(n_r1, n_r2);  % 匹配质量分

    for i = 1:n_r1
        life_i = r1_lives(i);
        for j = 1:n_r2
            if coexist_count(i,j) == 0, continue; end

            % 根据航迹长度选择T2比例
            if life_i < short_life_thresh || r2_lives(j) < short_life_thresh
                t2_ratio = T2_ratio_short;
            else
                t2_ratio = T2_ratio;
            end

            match_ratio = match_count(i,j) / coexist_count(i,j);
            meets_ratio = match_ratio >= t2_ratio;
            meets_abs = match_count(i,j) >= T2_abs_min;

            % 额外检查: 连续匹配帧数 (防止零星匹配)
            dhist = dist_history{i,j};
            consec_count = max_consecutive_below(dhist, T1_km);

            is_candidate(i,j) = meets_ratio && meets_abs && (consec_count >= 4);
            if is_candidate(i,j)
                match_quality(i,j) = match_ratio * 100 + consec_count;
            end
        end
    end

    % ---- 全局贪心分配 ----
    assigned_r1 = false(n_r1, 1);
    assigned_r2 = false(n_r2, 1);
    pairs = [];  % [r1_idx, r2_idx, quality]

    remaining = find(any(is_candidate, 2) & ~assigned_r1);
    while ~isempty(remaining)
        % 找当前最佳匹配
        best_q = -inf;
        best_i = 0; best_j = 0;
        for ii = 1:length(remaining)
            i = remaining(ii);
            for j = 1:n_r2
                if assigned_r2(j) || ~is_candidate(i,j), continue; end
                if match_quality(i,j) > best_q
                    best_q = match_quality(i,j);
                    best_i = i; best_j = j;
                end
            end
        end

        if best_i == 0, break; end

        pairs(end+1, :) = [best_i, best_j, best_q];
        assigned_r1(best_i) = true;
        assigned_r2(best_j) = true;

        remaining = find(any(is_candidate, 2) & ~assigned_r1);
    end

    % ---- 构建输出 ----
    matched_pairs = struct('R1_track_id', {}, 'R2_track_id', {}, ...
        'match_count', {}, 'coexist_count', {}, 'match_ratio', {}, ...
        'mean_dist_km', {}, 'quality', {});
    for p = 1:size(pairs, 1)
        i = pairs(p,1); j = pairs(p,2);
        dhist = dist_history{i,j};
        matched_pairs(p).R1_track_id = r1_ids(i);
        matched_pairs(p).R2_track_id = r2_ids(j);
        matched_pairs(p).match_count = match_count(i,j);
        matched_pairs(p).coexist_count = coexist_count(i,j);
        matched_pairs(p).match_ratio = match_count(i,j) / coexist_count(i,j);
        matched_pairs(p).mean_dist_km = mean(dhist);
        matched_pairs(p).quality = pairs(p,3);
    end

    % 未匹配航迹
    unmatched_r1 = r1_ids(~assigned_r1);
    unmatched_r2 = r2_ids(~assigned_r2);

    matcher = struct(...
        'matched_pairs', matched_pairs, ...
        'unmatched_R1', unmatched_r1, ...
        'unmatched_R2', unmatched_r2, ...
        'match_count', match_count, ...
        'coexist_count', coexist_count, ...
        'match_quality', match_quality, ...
        'r1_ids', r1_ids, ...
        'r2_ids', r2_ids, ...
        'r1_pos', r1_pos, ...
        'r2_pos', r2_pos, ...
        'aligned_R2', {aligned_R2});

    % ---- 打印匹配结果 ----
    fprintf('\n========== 匹配结果 ==========\n');
    for p = 1:length(matched_pairs)
        mp = matched_pairs(p);
        fprintf('  R1#%d <-> R2#%d: match=%d/%d (%.0f%%), mean_dist=%.1fkm, quality=%.1f\n', ...
            mp.R1_track_id, mp.R2_track_id, mp.match_count, mp.coexist_count, ...
            mp.match_ratio*100, mp.mean_dist_km, mp.quality);
    end
    if ~isempty(unmatched_r1)
        fprintf('  未匹配R1航迹: %s\n', mat2str(unmatched_r1));
    end
    if ~isempty(unmatched_r2)
        fprintf('  未匹配R2航迹: %s\n', mat2str(unmatched_r2));
    end
end

% =========================================================================
% 辅助函数
% =========================================================================

function [pos, ids, lives] = extract_track_positions(snapshots, n_frames)
    % 从快照中提取每条航迹的逐帧位置历史
    % pos: n_tracks × n_frames × 2 (lon, lat), NaN表示不存在

    % 第一遍: 收集所有航迹ID
    all_ids = [];
    for k = 1:n_frames
        snap = snapshots{k};
        if isempty(snap) || ~isfield(snap, 'trackList'), continue; end
        for t = 1:length(snap.trackList)
            trk = snap.trackList{t};
            if trk.type == 7, continue; end
            all_ids = [all_ids, trk.id];
        end
    end
    ids = unique(all_ids);

    n_tracks = length(ids);
    pos = nan(n_tracks, n_frames, 2);
    lives = zeros(n_tracks, 1);

    for k = 1:n_frames
        snap = snapshots{k};
        if isempty(snap) || ~isfield(snap, 'trackList'), continue; end
        for t = 1:length(snap.trackList)
            trk = snap.trackList{t};
            if trk.type == 7, continue; end
            idx = find(ids == trk.id, 1);
            if ~isempty(idx)
                pos(idx, k, 1) = trk.lon;
                pos(idx, k, 2) = trk.lat;
                lives(idx) = lives(idx) + 1;
            end
        end
    end
end

function n = max_consecutive_below(values, threshold)
    % 计算序列中连续低于阈值的最大长度
    n = 0;
    current = 0;
    for i = 1:length(values)
        if values(i) < threshold
            current = current + 1;
            n = max(n, current);
        else
            current = 0;
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
