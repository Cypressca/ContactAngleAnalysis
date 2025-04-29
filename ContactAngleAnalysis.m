function ContactAngleAnalysis()
    % 选择文件夹
    folder = uigetdir('选择包含TIFF图像的文件夹');
    if folder == 0, return; end
    
    % 获取所有TIFF文件
    files = dir(fullfile(folder, '*.jpg'));
    if isempty(files)
        errordlg('未找到TIFF文件');
        return;
    end
    
    % 创建结果保存路径
    result_dir = fullfile(folder, 'results');
    if ~exist(result_dir, 'dir')
        mkdir(result_dir);
    end
    
    % 处理每个图像
    for i = 1:length(files)
        processImage(fullfile(folder, files(i).name), result_dir);
    end
end

function processImage(imgPath, result_dir)
    % 读取图像
    try
        I = imread(imgPath);
    catch
        return;
    end
    
    % 转换为灰度图像
    if size(I,3) == 3
        Igray = rgb2gray(I);
    else
        Igray = I;
    end
    
    % 交互式阈值设置
    h = figure('Name','设置阈值');
    imshow(Igray);
    title('灰度图像 - 设置阈值后关闭窗口');
    
    % 创建阈值调节界面
    uicontrol('Style','slider',...
        'Position',[20 20 300 20],...
        'Min',0,'Max',1,'Value',0.5,...
        'Callback',@updateThreshold);
    
    threshold = 0.6;
    function updateThreshold(src,~)
        threshold = get(src,'Value');
        BW = imbinarize(Igray, threshold);
        subplot(1,1,1);
        imshowpair(Igray, BW, 'montage');
        title(['当前阈值: ' num2str(threshold,'%.2f')]);
    end
    
    waitfor(h);
    
    % 应用最终阈值
    BW = imbinarize(Igray, threshold);
    BW = imfill(BW, 'holes');
    
    % 形态学操作
    BW = imopen(BW, strel('disk',3));
    
    % 找到最大连通区域
    [BW_clean, boundary] = findLargestRegion(BW);
    
    % 检测角点
    corners = detectCorners(boundary);
    
    % 筛选四个纵坐标最低的点，并选择中间的两个
    ellipsePoints = [];
    if size(corners, 1) >= 4
        [~, idx] = sort(corners(:, 2), 'descend');
        sorted_corners = corners(idx, :);
        lowest_four = sorted_corners(1:4, :);
        
        x_coords = lowest_four(:, 1);
        median_x = median(x_coords);
        
        distances = abs(x_coords - median_x);
        [~, idx_sorted] = sort(distances);
        selected_indices = idx_sorted(1:2);
        central_two = lowest_four(selected_indices, :);
        
        [~, x_order] = sort(central_two(:, 1));
        central_two = central_two(x_order, :);
        
        % 计算直线方程
        x1 = central_two(1, 1); y1 = central_two(1, 2);
        x2 = central_two(2, 1); y2 = central_two(2, 2);
        
        if x2 ~= x1
            m = (y2 - y1)/(x2 - x1);
            b = y1 - m*x1;
            line_eq = [m, b];
        else
            m = Inf;
            b = x1;
            line_eq = [m, b];
        end
            
        % 椭圆拟合
        try
            % 查找边界点索引
            pt1 = [x1, y1]; pt2 = [x2, y2];
            boundary_xy = fliplr(boundary); % 转换为[x,y]格式
            
            % 寻找最近点索引
            [~, idx1] = min(sum((boundary_xy - pt1).^2, 2));
            [~, idx2] = min(sum((boundary_xy - pt2).^2, 2));
            
            % 提取两点间较短的轮廓段
            N = size(boundary_xy, 1);
            if abs(idx1 - idx2) < N/2
                contour_segment = boundary_xy(min(idx1,idx2):max(idx1,idx2), :);
            else
                contour_segment = boundary_xy([max(idx1,idx2):end, 1:min(idx1,idx2)], :);
            end
            
            % 椭圆拟合
            if size(contour_segment,1) >= 5
                ellipse = fitEllipse(contour_segment);
                
                % 生成椭圆点
                theta = ellipse.theta;
                h = ellipse.center(1); k = ellipse.center(2);
                a = ellipse.major_axis; b = ellipse.minor_axis;
                
                t = linspace(0, 2*pi, 100);
                x_ellipse = h + a*cos(t)*cos(theta) - b*sin(t)*sin(theta);
                y_ellipse = k + a*cos(t)*sin(theta) + b*sin(t)*cos(theta);
                ellipsePoints = [x_ellipse', y_ellipse'];
                
                % 计算R²值
                x_data = contour_segment(:,1);
                y_data = contour_segment(:,2);
                
                % 获取椭圆方程系数
                A = ellipse.A; B = ellipse.B; C = ellipse.C;
                D = ellipse.D; E = ellipse.E; F = ellipse.F;
                
                % 计算残差平方和
                residuals = A*x_data.^2 + B*x_data.*y_data + C*y_data.^2 + D*x_data + E*y_data + F;
                SS_res = sum(residuals.^2);
                
                % 计算总平方和（基于中心点）
                SS_tot = sum((x_data - h).^2 + (y_data - k).^2);
                
                % 计算R²
                R_squared = 1 - (SS_res / SS_tot);
                
                % 保存椭圆参数和R²
                [~, name] = fileparts(imgPath);
                save_path = fullfile(result_dir, [name '_ellipse.txt']);
                fid = fopen(save_path, 'w');
                fprintf(fid, '椭圆标准方程:\n');
                fprintf(fid, '((x-%.4f)cos(%.4f) + (y-%.4f)sin(%.4f))^2/%.4f^2 + ', h, theta, k, theta, a);
                fprintf(fid, '((x-%.4f)sin(%.4f) - (y-%.4f)cos(%.4f))^2/%.4f^2 = 1\n', h, theta, k, theta, b);
                fprintf(fid, 'R²值: %.4f\n', R_squared);
                fclose(fid);
            else
                ellipsePoints = [];
                R_squared = NaN;
            end
        catch
            ellipsePoints = [];
            R_squared = NaN;
        end
    else
        central_two = [];
        line_eq = [];
        R_squared = NaN;
    end

    % 计算接触角（新增部分）
    contact_angles = [];
    if ~isempty(central_two) && ~isempty(ellipsePoints)
        % 提取椭圆参数
        A = ellipse.A; B = ellipse.B; C = ellipse.C;
        D = ellipse.D; E = ellipse.E; 
        % 基线参数
        m_line = line_eq(1); 
        
        % 遍历每个接触点
        for i = 1:size(central_two,1)
            x0 = central_two(i,1);
            y0 = central_two(i,2);
            
            % 计算椭圆在该点的切线斜率
            numerator = -(2*A*x0 + B*y0 + D);
            denominator = B*x0 + 2*C*y0 + E;
            if denominator == 0
                m_tan = Inf; % 垂直切线
            else
                m_tan = numerator / denominator;
            end
            
            % 计算与基线夹角
            if m_line == Inf || m_line == -Inf
                % 基线垂直
                angle_rad = abs(pi/2 - atan(m_tan));
            elseif abs(m_tan) > 1e10
                % 切线垂直
                angle_rad = abs(pi/2 - atan(m_line));
            else
                angle_rad = atan(abs((m_tan - m_line)/(1 + m_tan*m_line)));
            end
            contact_angles = [contact_angles; rad2deg(angle_rad)];
        end
    end

    visualizeResults(I, Igray, BW_clean, boundary, imgPath, result_dir,...
        corners, central_two, line_eq, ellipsePoints, R_squared, contact_angles);
end


function visualizeResults(I, Igray, BW, boundary, imgPath, save_dir,...
        corners, selectedPoints, lineParams, ellipsePoints, R_squared, contact_angles)
    fig = figure('Visible', 'off', 'Position', [100 100 1400 800]);
    
    % 获取图像尺寸
    img_height = size(BW, 1);
    img_width = size(BW, 2);
    
    % 子图1：原始图像
    subplot(1, 3, 1);
    imshow(I);
    title('原始图像');
    
    % 子图2：灰度图像
    subplot(1, 3, 2);
    imshow(Igray);
    title('灰度图像');
    
    % 子图3：绘制拟合椭圆和接触角
    subplot(1, 3, 3);
    imshow(BW);
    hold on;

    % 绘制角点
    if ~isempty(corners)
        plot(corners(:, 1), corners(:, 2), 'ro', 'MarkerSize', 10);
    end
    
    % 绘制拟合椭圆和R²值
    if ~isempty(ellipsePoints)
        plot(ellipsePoints(:,1), ellipsePoints(:,2), 'm-', 'LineWidth', 2);
        % 在左上角显示R²值
        text(10, 20, sprintf('R² = %.3f', R_squared),...
            'Color', 'white', 'BackgroundColor', 'black',...
            'VerticalAlignment', 'top', 'FontSize', 10);
    end

    % 绘制选中点和基线
    if ~isempty(selectedPoints)
        plot(selectedPoints(:, 1), selectedPoints(:, 2), 'g*', 'MarkerSize', 12);
        x_lim = xlim; y_lim = ylim;
        if lineParams(1) ~= Inf
            x_fit = linspace(x_lim(1), x_lim(2), 100);
            y_fit = lineParams(1)*x_fit + lineParams(2);
        else
            y_fit = linspace(y_lim(1), y_lim(2), 100);
            x_fit = lineParams(2)*ones(size(y_fit));
        end
        plot(x_fit, y_fit, 'c--', 'LineWidth', 2);
    end
    
    % 绘制拟合椭圆
    if ~isempty(ellipsePoints)
        plot(ellipsePoints(:,1), ellipsePoints(:,2), 'm-', 'LineWidth', 2);
    end
    title('二值化、角点和拟合结果');

    if ~isempty(contact_angles)
        % 主标注：底部中央
        text_str = sprintf('接触角: %.1f° | %.1f°',...
            contact_angles(1), contact_angles(2));
        text(img_width/2, img_height-30, text_str,...
            'Color', [1 1 0],...       % 黄色文字
            'FontSize', 14,...
            'FontWeight', 'bold',...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'bottom',...
            'BackgroundColor', [1 0 0],... % 红色背景
            'Margin', 3);
        
        % 切线标注：接触点旁
        for i = 1:length(contact_angles)
            % 箭头起点（接触点）
            x_point = selectedPoints(i, 1);
            y_point = selectedPoints(i, 2);
            
            % 绘制指示箭头
            annotation('arrow',...
                [x_point/img_width, (x_point+50)/img_width],... % X方向偏移
                [y_point/img_height, (y_point-80)/img_height],... % Y方向偏移
                'Color', [0 1 0], 'LineWidth', 1.5);
            
            % 角度值标注
            text(x_point+20, y_point-60,...
                sprintf('%.1f°', contact_angles(i)),...
                'Color', [0 1 1],...    % 青色文字
                'FontSize', 12,...
                'FontWeight', 'bold');
        end
    end
    
    % 保存结果（保持原有代码不变）
    [~, name] = fileparts(imgPath);
    save_path = fullfile(save_dir, [name '_result.png']);
    print(fig, save_path, '-dpng', '-r300');
    close(fig);
end
function [BW_clean, boundary] = findLargestRegion(BW)
    CC = bwconncomp(BW);
    numPixels = cellfun(@numel, CC.PixelIdxList);
    [~, idx] = max(numPixels);
    BW_clean = false(size(BW));
    BW_clean(CC.PixelIdxList{idx}) = true;
    boundary = bwboundaries(BW_clean);
    boundary = boundary{1};
end

function corners = detectCorners(boundary)
    x = boundary(:,2); y = boundary(:,1);
    points = [x, y];
    
    min_pt = min(points); max_pt = max(points);
    bbox_diag = norm(max_pt - min_pt);
    
    dx = diff(x); dy = diff(y);
    distances = sqrt(dx.^2 + dy.^2);
    perimeter = sum(distances);
    
    tol_ratio = 0.01;
    tol_pixels = perimeter * tol_ratio;
    tol = max(min(tol_pixels / bbox_diag, 1), 0.001);
    
    simplified_points = reducepoly(points, tol);
    [~, idx] = unique(simplified_points, 'rows', 'stable');
    corners = simplified_points(sort(idx), :);
end

function ellipse = fitEllipse(points)
    % 椭圆拟合函数（基于直接最小二乘法）
    x = points(:,1); y = points(:,2);
    
    % 构建设计矩阵
    D = [x.^2, x.*y, y.^2, x, y, ones(size(x))];
    S = D'*D;
    
    % 构建约束矩阵
    C = zeros(6);
    C(1,3) = 2; C(2,2) = -1; C(3,1) = 2;
    
    % 求解广义特征值
    [V, ~] = eig(S, C);
    
    % 选择正确特征向量
    cond = 4*V(1,:).*V(3,:) - V(2,:).^2;
    idx = find(cond > 0 & ~isinf(cond));
    if isempty(idx), error('无法拟合椭圆'); end
    params = V(:,idx(1));
    
    % 解析椭圆参数
    A = params(1); B = params(2); C = params(3);
    D = params(4); E = params(5); F = params(6);
    
    % 计算中心坐标
    denominator = B^2 - 4*A*C;
    h = (2*C*D - B*E)/denominator;
    k = (2*A*E - B*D)/denominator;
    
    % 计算轴长和角度
    term = sqrt((A - C)^2 + B^2);
    a = sqrt(-(F - (A*h^2 + B*h*k + C*k^2))/((A + C + term)/2));
    b = sqrt(-(F - (A*h^2 + B*h*k + C*k^2))/((A + C - term)/2));
    theta = 0.5*atan2(B, A - C);
    
    % 返回结构体（新增方程系数存储）
    ellipse = struct('center', [h, k], 'major_axis', a,...
        'minor_axis', b, 'theta', theta,...
        'A', A, 'B', B, 'C', C, 'D', D, 'E', E, 'F', F);
end