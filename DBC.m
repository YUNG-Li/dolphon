function CENTRIOD = DBC(images)
    
    % 读取当前图像
    I = imread(images);
  
    % 获取当前图像的尺寸
    [height, width, ~] = size(I);
    
    % 对图像进行裁剪（保持整个图像，没有实际裁剪）
    I = imcrop(I, [1 1 width height]);

    

    % 如果图像有多个通道，将其转换为灰度图
    [height, width, chennal] = size(I);
    if chennal > 1
        I = rgb2gray(I);  % 将 RGB 图像转换为灰度图
    end

    % 使用标准差滤波器进行边缘检测
    Zs = stdfilt(I);  % 计算标准差
    Z = double(I);  % 将图像转换为 double 类型以便后续处理
    Zs = Zs .* Z;  % 元素逐次相乘，增强强度变化
    
    %% 使用梯度阈值粗略地选取候选点
    rate = 0.5;  
    max_rate = 0.9;  
    mean_std = mean2(Zs); 
    max_std = max(Zs(:));  

    % 找到强度大于阈值的候选点
    [col_candi, row_candi] = find(Zs > rate * max_std); 
    points_candi = [row_candi, col_candi];  % 候选像素点的位置

    % DBSCAN 聚类参数设置
    epsilon = 10;  % 邻域半径
    MinPts = 1;  % 形成簇的最小点数
    IDX_candi = DBSCAN(points_candi, epsilon, MinPts);  % 对候选点进行聚类
    k_candi = max(IDX_candi);  % 聚类数量

    % 初始化每个聚类的平均位置
    mean_points_candi = zeros(k_candi, 2);  
    retina_pixel_FOV = 0.1;  
    small_rect = (1.2 / retina_pixel_FOV / 2);  
    large_rect = round(3 / retina_pixel_FOV / 2); 

    cumulative_points = [];  % 用于进一步分析的累积点
    min_width_large = zeros(k_candi, 1);
    max_width_large = zeros(k_candi, 1);  
    min_heigh_large = zeros(k_candi, 1);  
    max_heigh_large = zeros(k_candi, 1);  

  
    for j = 1:k_candi
        % 计算聚类的平均位置
        if size(points_candi(IDX_candi == j, :), 1) == 1
            mean_points_candi(j, :) = points_candi(IDX_candi == j, :);
        else
            mean_points_candi(j, :) = mean(points_candi(IDX_candi == j, :));
        end

        
        min_width_large(j) = round(min(mean_points_candi(j, 1))) - large_rect; % 行
        max_width_large(j) = round(max(mean_points_candi(j, 1))) + large_rect;
        min_heigh_large(j) = round(min(mean_points_candi(j, 2))) - large_rect; % 列
        max_heigh_large(j) = round(max(mean_points_candi(j, 2))) + large_rect;

        
        min_width_large(j) = max(min_width_large(j), 1);
        max_width_large(j) = min(max_width_large(j), width);
        min_heigh_large(j) = max(min_heigh_large(j), 1);
        max_heigh_large(j) = min(max_heigh_large(j), height);
        
        
        mask_candi = zeros(height, width);
        mask_candi(min_heigh_large(j):max_heigh_large(j), min_width_large(j):max_width_large(j)) = 1;
        
        
        meanz = mean2(mask_candi .* Z);
        maxz = max(max(mask_candi .* Z));
        
        
        [col{j}, row{j}] = find(mask_candi .* Z > meanz + max_rate * (maxz - meanz));

        points{j} = [row{j} col{j}];  % 存储聚类的像素点
        cumulative_points = [cumulative_points; points{j}];  % 累积点用于进一步分析
    end

    % 初始化一些变量
    count = 0;
    cum_IDX = DBSCAN(cumulative_points, epsilon / 2, MinPts * 2);  % 对累积点进行更严格的 DBSCAN 聚类
    k_obj = max(cum_IDX);  
    new_points = [];  

    
    if k_obj == 0
        new_points = cumulative_points;
        k_obj = size(new_points, 1);
    end

    
    for j = 1:k_obj
       
        if cum_IDX(j) ~= 0 || max(cum_IDX) ~= 0
            new_points(j, :) = mean(cumulative_points(cum_IDX == j, :));  % 计算平均位置
            new_points_min(j, 1) = min(cumulative_points(cum_IDX == j, 1));  % 获取最小X坐标
            new_points_min(j, 2) = min(cumulative_points(cum_IDX == j, 2));  % 获取最小Y坐标
            new_points_max(j, 1) = max(cumulative_points(cum_IDX == j, 1));  % 获取最大X坐标
            new_points_max(j, 2) = max(cumulative_points(cum_IDX == j, 2));  % 获取最大Y坐标
        else
          
            new_points_min(j, 1) = min(new_points(:, 1));
            new_points_min(j, 2) = min(new_points(:, 2));
            new_points_max(j, 1) = max(new_points(:, 1));
            new_points_max(j, 2) = max(new_points(:, 2));
        end


% 初始化目标点的数量和位置
n = 0;  % 目标点数量初始化
target_points = zeros(k_obj, 2);  % 目标点坐标初始化
for j = 1:k_obj
    
    if new_points_max(j, 1) - new_points_min(j, 1) <= small_rect * 3 && new_points_max(j, 2) - new_points_min(j, 2) <= small_rect * 3
        n = n + 1;
        target_points(n, :) = new_points(j, :);  % 保存目标点位置
    end
    
   
    if j == k_obj && n == 0
        n = k_obj;  % 目标点数量设为总对象数
        target_points(n, :) = new_points(n, :);  % 保存最后一个聚类作为目标
    end
end

% 初始化目标点的矩形框边界
min_width = zeros(n, 1);
max_width = zeros(n, 1);
min_heigh = zeros(n, 1);
max_heigh = zeros(n, 1);


patch = zeros(round(2 * small_rect), round(2 * small_rect));
number = 0;  % 目标的计数
len = 4;  
nr = 3;  
nc = 3; 
leny = len * nr;  
lenx = len * nc;  
op = zeros(leny, lenx, nr * nc); 
% 生成操作矩阵（分割为9个小矩形）
for ii = 1:nr * nc
    temp = zeros(len * nr, len * nc);  % 创建临时矩阵
    [r, c] = ind2sub([nr, nc], ii);  % 获取行列索引
    temp((r - 1) * len + 1:r * len, (c - 1) * len + 1:c * len) = 1;  % 设置小矩形区域
    temp = temp';  % 转置矩阵
    op(:, :, ii) = temp;  % 保存到操作矩阵
end


small_rect = 6;  % 统一的小矩形大小
min_size = 12;   % 确保裁剪区域的最小尺寸

% 遍历每个目标点
for k = 1:n
    % 定义小矩形边界
    min_width = round(target_points(k, 1)) - small_rect;  
    max_width = round(target_points(k, 1)) + small_rect - 1;  
    min_heigh = round(target_points(k, 2)) - small_rect;  
    max_heigh = round(target_points(k, 2)) + small_rect - 1;  

    % 确保裁剪区域在图像范围内
    min_width = max(min_width, 1);
    max_width = min(max_width, width - 1);
    min_heigh = max(min_heigh, 1);
    max_heigh = min(max_heigh, height - 1);

    % 裁剪并确保大小一致
    patch = imcrop(I, [min_width min_heigh 2 * small_rect - 1 2 * small_rect - 1]);
    if size(patch, 1) ~= min_size || size(patch, 2) ~= min_size
        % warning('Patch size is inconsistent for image %d', i);
        continue;  % 跳过此帧的处理
    end

    % 计算每个小矩形的强度
    for ii = 1:nr * nc
        gimg(:, :, ii) = double(patch) .* op(:, :, ii);  % 将补丁与操作矩阵相乘
        intensity(k, ii) = max(max(gimg(:, :, ii)));  % 找到每个小矩形中的最大强度
    end

    % 判断目标点是否符合强度要求
    if intensity(k, 5) - max(intensity(k, [1:4, 6:end])) >= 5
        number = number + 1;
        target_position{i}(number, :) = round(target_points(k, :));  % 保存目标点的位置
    end
    
    % 如果没有符合条件的目标点，选取最后一个点作为目标
    if number == 0 && n == 1
        number = 1;
        target_position{i}(number, :) = round(target_points(n, :));
    elseif number == 0
        target_position{i} = [];  % 无符合条件的目标点
    end
end
CENTRIOD = target_position;




