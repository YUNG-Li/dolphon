clear; clc; close all;  % 清除工作空间，命令窗口，关闭所有图形窗口

% 设置输入图像文件夹路径
input_folder = 'D:\dolphin_dataset\handled\oringin\detect-train-1-3\images';  
ext = {'*.jpg', '*.png', '*.bmp'};  % 支持的图像文件扩展名
images = [];  % 初始化图像列表

% 创建变量以存储上一个图像的尺寸
prev_height = -1;
prev_width = -1;

% 初始化视频段计数器
videoSegmentCounter = 0;

% 遍历每种扩展名，找到输入文件夹中的所有图像
for i = 1:length(ext)
    images = [images; dir(fullfile(input_folder, ext{i}))];  % 将找到的图像文件添加到列表中
end
% 遍历所有图像  
tic
for i = 1:numel(images)
    tic  % 开始计时
    % 提取文件名（不包括扩展名）
    [~, name, exte] = fileparts(images(i).name);
    
    % 读取当前图像
    I = imread(fullfile(input_folder, images(i).name));
    elapsedTime1 = toc;  % 结束计时并计算读取时间
    fprintf('读图经过的时间: %.4f 秒\n', elapsedTime1);  % 输出读取时间

    tic  % 重新开始计时
    % 获取当前图像的尺寸
    [height, width, ~] = size(I);
    
    % 对图像进行裁剪（保持整个图像，没有实际裁剪）
    I = imcrop(I, [1 1 width height]);
%%
    % 检查当前图像的尺寸是否与前一个图像不同
%     if height ~= prev_height || width ~= prev_width
%         % 如果是新的尺寸，关闭之前的视频文件（如果已打开）
%         if i > 1  % 确保不是第一张图像
%             close(outputVideo);
%             fprintf('开始新的视频段：%d x %d\n', height, width);  % 输出新的视频段信息
%         end
%         
%         % 更新视频段计数器
%         videoSegmentCounter = videoSegmentCounter + 1;
%         
%         % 创建新的视频对象，使用计数器确保文件名唯一
%         outputVideo = VideoWriter(sprintf('output_video_%dx%d_segment%d.mp4', width, height, videoSegmentCounter), 'MPEG-4'); 
%         open(outputVideo);  % 打开视频写入
%         
%         % 更新之前图像的尺寸
%         prev_height = height;
%         prev_width = width;
%     end
%%
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
    rate = 0.5;  % 梯度阈值
    max_rate = 0.9;  % 最大阈值
    mean_std = mean2(Zs);  % 过滤后的图像的均值
    max_std = max(Zs(:));  % 过滤后图像的最大强度

    % 找到强度大于阈值的候选点
    [col_candi, row_candi] = find(Zs > rate * max_std); 
    points_candi = [row_candi, col_candi];  % 候选像素点的位置

    % DBSCAN 聚类参数设置
    epsilon = 10;  % 邻域半径
    MinPts = 1;  % 形成簇的最小点数
    tic
    IDX_candi = DBSCAN(points_candi, epsilon, MinPts);  % 对候选点进行聚类
    elapsedTime3 = toc;  % 记录处理当前帧所用的时间
    fprintf('DBSCAN处理时间为: %.4f 秒\n', elapsedTime3)
    k_candi = max(IDX_candi);  % 聚类数量

    % 初始化每个聚类的平均位置
    mean_points_candi = zeros(k_candi, 2);  
    retina_pixel_FOV = 0.1;  % 每像素的视场角（单位：毫弧度）
    small_rect = (1.2 / retina_pixel_FOV / 2);  % 小物体的矩形框大小
    large_rect = round(3 / retina_pixel_FOV / 2);  % 大物体的矩形框大小

    cumulative_points = [];  % 用于进一步分析的累积点
    min_width_large = zeros(k_candi, 1);  % 初始化最小行边界
    max_width_large = zeros(k_candi, 1);  % 初始化最大行边界
    min_heigh_large = zeros(k_candi, 1);  % 初始化最小列边界
    max_heigh_large = zeros(k_candi, 1);  % 初始化最大列边界

    % 遍历每个聚类
    for j = 1:k_candi
        % 计算聚类的平均位置
        if size(points_candi(IDX_candi == j, :), 1) == 1
            mean_points_candi(j, :) = points_candi(IDX_candi == j, :);
        else
            mean_points_candi(j, :) = mean(points_candi(IDX_candi == j, :));
        end

        % 定义聚类的矩形边界
        min_width_large(j) = round(min(mean_points_candi(j, 1))) - large_rect; % 行
        max_width_large(j) = round(max(mean_points_candi(j, 1))) + large_rect;
        min_heigh_large(j) = round(min(mean_points_candi(j, 2))) - large_rect; % 列
        max_heigh_large(j) = round(max(mean_points_candi(j, 2))) + large_rect;

        % 确保边界在图像范围内
        min_width_large(j) = max(min_width_large(j), 1);
        max_width_large(j) = min(max_width_large(j), width);
        min_heigh_large(j) = max(min_heigh_large(j), 1);
        max_heigh_large(j) = min(max_heigh_large(j), height);
        
        % 创建候选区域的掩膜
        mask_candi = zeros(height, width);
        mask_candi(min_heigh_large(j):max_heigh_large(j), min_width_large(j):max_width_large(j)) = 1;
        
        % 计算候选区域的平均和最大强度
        meanz = mean2(mask_candi .* Z);
        maxz = max(max(mask_candi .* Z));
        
        % 找到强度高于平均值的像素
        [col{j}, row{j}] = find(mask_candi .* Z > meanz + max_rate * (maxz - meanz));

        points{j} = [row{j} col{j}];  % 存储聚类的像素点
        cumulative_points = [cumulative_points; points{j}];  % 累积点用于进一步分析
    end

    % 初始化一些变量
    count = 0;
    cum_IDX = DBSCAN(cumulative_points, epsilon / 2, MinPts * 2);  % 对累积点进行更严格的 DBSCAN 聚类
    k_obj = max(cum_IDX);  % 找到对象的数量
    new_points = [];  % 初始化新聚类中心

    % 如果没有找到聚类，则将所有累积点视为一个新点
    if k_obj == 0
        new_points = cumulative_points;
        k_obj = size(new_points, 1);
    end

    % 遍历每个聚类
    for j = 1:k_obj
        % 如果对象存在，则计算聚类的平均位置
        if cum_IDX(j) ~= 0 || max(cum_IDX) ~= 0
            new_points(j, :) = mean(cumulative_points(cum_IDX == j, :));  % 计算平均位置
            new_points_min(j, 1) = min(cumulative_points(cum_IDX == j, 1));  % 获取最小X坐标
            new_points_min(j, 2) = min(cumulative_points(cum_IDX == j, 2));  % 获取最小Y坐标
            new_points_max(j, 1) = max(cumulative_points(cum_IDX == j, 1));  % 获取最大X坐标
            new_points_max(j, 2) = max(cumulative_points(cum_IDX == j, 2));  % 获取最大Y坐标
        else
            % 如果没有找到聚类，取累积点的最小和最大0
            % 如果没有找到聚类，取累积点的最小和最大值
            new_points_min(j, 1) = min(new_points(:, 1));
            new_points_min(j, 2) = min(new_points(:, 2));
            new_points_max(j, 1) = max(new_points(:, 1));
            new_points_max(j, 2) = max(new_points(:, 2));
        end
    end 

% 初始化目标点的数量和位置
n = 0;  % 目标点数量初始化
target_points = zeros(k_obj, 2);  % 目标点坐标初始化
for j = 1:k_obj
    % 如果对象满足小矩形的大小约束，则认为它是一个目标点
    if new_points_max(j, 1) - new_points_min(j, 1) <= small_rect * 3 && new_points_max(j, 2) - new_points_min(j, 2) <= small_rect * 3
        n = n + 1;
        target_points(n, :) = new_points(j, :);  % 保存目标点位置
    end
    
    % 如果遍历完所有对象没有目标点，则选择最后一个聚类作为目标
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

% 初始化用于裁剪的小矩形补丁
patch = zeros(round(2 * small_rect), round(2 * small_rect));
number = 0;  % 目标的计数
len = 4;  % 矩形分割参数
nr = 3;  % 行数
nc = 3;  % 列数
leny = len * nr;  % 矩形高度
lenx = len * nc;  % 矩形宽度
op = zeros(leny, lenx, nr * nc);  % 初始化操作矩阵

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
        warning('Patch size is inconsistent for image %d', i);
        patch = imresize(patch, [min_size, min_size]);  % 调整为 min_size
%         continue;  % 跳过此帧的处理
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

%%
% % 绘制最终确认的目标点
% figure('Visible', 'off');  % 不显示窗口以提高处理速度
% imshow(I); hold on;
% for check = 1:n
%     % 绘制每个目标点的矩形边界
%     rectangle('Position', [target_points(check, 1) - small_rect, target_points(check, 2) - small_rect, 2 * small_rect, 2 * small_rect], 'EdgeColor', 'r', 'LineWidth', 2);
% end
% hold off;
elapsedTime2 = toc;  % 记录处理当前帧所用的时间
fprintf('总处理经过的时间: %.4f 秒\n', elapsedTime2);
% % 写入当前帧到视频
% frame = getframe(gcf);  % 获取当前帧
% writeVideo(outputVideo, frame);  % 将帧写入视频
% close;  % 关闭图形窗口

%%



end
toc
elapsedTime4 = toc;  % 记录处理当前帧所用的时间
fprintf('总时间为: %.4f 秒\n', elapsedTime4)
% % 关闭视频文件（在最后一帧处理完后）
% if exist('outputVideo', 'var') && isvalid(outputVideo)
%     close(outputVideo);  % 关闭视频文件
% end
% disp('处理完成，视频已保存。');  % 输出处理完成的提示

