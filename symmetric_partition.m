function result = symmetric_partition(n, num_parts)
    % n: 要素の数 (整数)
    % num_parts: パートの数（整数）
    % result: 要素が1と0で構成された配列

    if num_parts > n
        error('パート数は要素数以下である必要があります');
    end

    % 初期化
    result = zeros(1, n);

    % 中心点の計算
    center = ceil(n / 2);
    half_parts = floor(num_parts / 2);

    % 対称的に1を配置
    for i = 1:half_parts
        offset = round((i - 0.5) * n / num_parts); % 偏差を計算

        % 左右のインデックスを計算
        left_index = center - offset;
        right_index = center + offset;

        % 配置（境界チェック込み）
        if left_index > 0
            result(left_index) = 1;
        end
        if right_index <= n
            result(right_index) = 1;
        end
    end

    % 奇数の場合は中央に1を追加
    if mod(num_parts, 2) == 1
        result(center) = 1;
    end
end
