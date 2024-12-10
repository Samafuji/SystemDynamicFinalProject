function result = symmetric_partition(n, num_parts)
    % n: 要素の数 (整数)
    % num_parts: パートの数（整数）
    % result: 要素が1と0で構成された配列

    if num_parts > n
        error('num_parts need to be smaller than n');
    end

    result = zeros(1, n);

    center = ceil(n / 2);
    half_parts = floor(num_parts / 2);

    for i = 1:half_parts
        offset = round((i - 0.5) * n / num_parts);

        left_index = center - offset;
        right_index = center + offset;

        if left_index > 0
            result(left_index) = 1;
        end
        if right_index <= n
            result(right_index) = 1;
        end
    end

    if mod(num_parts, 2) == 1
        result(center) = 1;
    end
end
