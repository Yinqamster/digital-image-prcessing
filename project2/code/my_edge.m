function output = my_edge(input_image)
%in this function, you should finish the edge detection utility.
%the input parameter is a matrix of a gray image
%the output parameter is a matrix contains the edge index using 0 and 1
%the entries with 1 in the matrix shows that point is on the edge of the
%image
%you can use different methods to complete the edge detection function
%the better the final result and the more methods you have used, you will get higher scores  

%%
% roberts
[m,n]=size(input_image);
outputRoberts=zeros(m,n);
robertsNum=0; %经roberts算子计算得到的每个像素的值
for j=1:m-1 
    for k=1:n-1
        robertsNum = abs(input_image(j,k)-input_image(j+1,k+1)) + abs(input_image(j+1,k)-input_image(j,k+1));
        roberts(j,k)=robertsNum;
    end
end
roberts = mat2gray(roberts);  % 将梯度矩阵转换为灰度图像
robertThreshold=graythresh(roberts);%?计算灰度阈值
for j=1:m-1 %进行边界提取
    for k=1:n-1
        if(roberts(j,k) >= robertThreshold)
            outputRoberts(j,k)=255;
        else
            outputRoberts(j,k)=0;
        end
    end
end
figure;clf;
imshow(outputRoberts);title('roberts');

%%
%Prewitt
[m,n]=size(input_image);
outputPrewitt=zeros(m,n);
PrewittNum=0; %经Prewitt算子计算得到的每个像素的值
for j=2:m-1 %进行边界提取
    for k=2:n-1
        PrewittNum=abs(input_image(j-1,k+1)-input_image(j+1,k+1)+input_image(j-1,k)-input_image(j+1,k)+input_image(j-1,k-1)-input_image(j+1,k-1))+abs(input_image(j-1,k+1)+input_image(j,k+1)+input_image(j+1,k+1)-input_image(j-1,k-1)-input_image(j,k-1)-input_image(j+1,k-1));
        Prewitt(j,k)=PrewittNum;
    end
end
Prewitt = mat2gray(Prewitt);  % 将梯度矩阵转换为灰度图像
PrewittThreshold=graythresh(Prewitt);%?计算灰度阈值
for j=1:m-1 %进行边界提取
    for k=1:n-1
        if(Prewitt(j,k) > PrewittThreshold)
            outputPrewitt(j,k)=255;
        else
            outputPrewitt(j,k)=0;
        end
    end
end
figure;clf;
imshow(outputPrewitt);title('Prewitt');

%%
%sobel
[m,n]=size(input_image);
outputSobel=zeros(m,n);
sobelNum=0; %经sobel算子计算得到的每个像素的值
for j=2:m-1 %进行边界提取
    for k=2:n-1
        sobelNum=abs(input_image(j-1,k+1)+2*input_image(j,k+1)+input_image(j+1,k+1)-input_image(j-1,k-1)-2*input_image(j,k-1)-input_image(j+1,k-1))+abs(input_image(j-1,k-1)+2*input_image(j-1,k)+input_image(j-1,k+1)-input_image(j+1,k-1)-2*input_image(j+1,k)-input_image(j+1,k+1));
        sobel(j,k)=sobelNum;
    end
end

sobel = mat2gray(sobel);  % 将梯度矩阵转换为灰度图像
sobelThreshold=graythresh(sobel);%?计算灰度阈值

for j=1:m-1 %进行边界提取
    for k=1:n-1
        if(sobel(j,k) >= sobelThreshold)
            outputSobel(j,k)=255;
        else
            outputSobel(j,k)=0;
        end
    end
end
figure;clf;
imshow(outputSobel);title('sobel');




%%
%Laplacian
[m,n]=size(input_image);
outputLaplician=zeros(m,n);
LaplacianNum=0; %经Laplician算子计算得到的每个像素的值
for j=2:m-1 %进行边界提取
    for k=2:n-1
        LaplacianNum=abs(4*input_image(j,k)-input_image(j-1,k)-input_image(j+1,k)-input_image(j,k+1)-input_image(j,k-1));
        Laplacian(j,k)=LaplacianNum;
    end
end

Laplacian = mat2gray(Laplacian);  % 将梯度矩阵转换为灰度图像
LaplacianThreshold=graythresh(Laplacian);%?计算灰度阈值

for j=1:m-1 %进行边界提取
    for k=1:n-1
        if(Laplacian(j,k) >= LaplacianThreshold)
            outputLaplacian(j,k)=255;
        else
            outputLaplacian(j,k)=0;
        end
    end
end
figure;clf;
imshow(outputLaplacian);title('Laplacian');


%%
% Marr
sigma=1.0;
[m,n]=size(input_image);
e=repmat(logical(uint8(0)),m,n);
rr=2:m-1;
cc=2:n-1;
fsize=ceil(sigma*7)*6+1;
op=fspecial('log',fsize,sigma);
op=op-sum(op(:))/numel(op);
b=filter2(op,input_image);
thresh=1.75*mean2(abs(b));   %thresh=.75*mean2(abs(b));
% thresh=graythresh(abs(b));
[rx,cx]=find(b(rr,cc)<0&b(rr,cc+1)>0&abs(b(rr,cc)-b(rr,cc+1))>thresh);
e((rx+1)+cx*m)=1;
[rx,cx]=find(b(rr,cc-1)>0&b(rr,cc)<0&abs(b(rr,cc-1)-b(rr,cc))>thresh);
e((rx+1)+cx*m)=1;
[rx,cx]=find(b(rr,cc)<0&b(rr+1,cc)>0&abs(b(rr,cc)-b(rr+1,cc))>thresh);
e((rx+1)+cx*m)=1;
[rx,cx]=find(b(rr-1,cc)>0&b(rr,cc)<0&abs(b(rr-1,cc)-b(rr,cc))>thresh);
e((rx+1)+cx*m)=1;
[rz,cz]=find(b(rr,cc)==0);
if isempty(rz)
    zero=(rz+1)+cz*m;
    zz=find(b(zero-1)<0&b(zero+1)>0&abs(b(zero-1)-b(zero+1))>2*thresh);
    e(zero(zz))=1;
    zz=find(b(zero-1)>0&b(zero+1)<0&abs(b(zero-1)-b(zero+1))>2*thresh);
    e(zero(zz))=1;
    zz=find(b(zero-m)<0&b(zero+m)>0&abs(b(zero-m)-b(zero+m))>2*thresh);
    e(zero(zz))=1;
    zz=find(b(zero-m)>0&b(zero+m)<0&abs(b(zero-m)-b(zero+m))>2*thresh);
    e(zero(zz))=1;
end
figure;clf;
imshow(e);title('Marr');

%%
%Marr
% sigma = 1.0;
% [m,n]=size(input_image);
% outputMarr=zeros(m,n);
% fsize=ceil(sigma*7)*6+1;
% gausFilter = fspecial('log', fsize, sigma);
% gausFilter=gausFilter-sum(gausFilter(:))/numel(gausFilter);
% Marr=filter2(gausFilter,input_image);
[m,n]=size(input_image);
sigma = 1.0;
fsize=ceil(sigma*7)*6+1;
gausFilter = fspecial('log', fsize, sigma);
Marr= imfilter(input_image, gausFilter, 'replicate');
MarrThreshold=graythresh(abs(Marr));%?计算灰度阈值
outputMarr=zeros(m,n);

for j=2:m-1 %进行边界提取
    for k=2:n-1
        if(Marr(j,k)<0 & Marr(j,k+1)>0 & abs(Marr(j,k)-Marr(j,k+1))>MarrThreshold)
            outputMarr(j,k)=255;
        end
        if(Marr(j,k)<0 & Marr(j,k-1)>0 &  abs(Marr(j,k)-Marr(j,k-1))>MarrThreshold)
            outputMarr(j,k)=255;
        end
        if(Marr(j,k)<0 & Marr(j+1,k)>0 & abs(Marr(j,k)-Marr(j+1,k))>MarrThreshold)
            outputMarr(j,k)=255;
        end
        if(Marr(j,k)<0 & Marr(j-1,k)>0 &  abs(Marr(j,k)-Marr(j-1,k))>MarrThreshold)
            outputMarr(j,k)=255;
        end 
    end
end
figure;clf;
imshow(outputMarr,[]);
title('Marr')

%%
%Canny
%高斯滤波
[m,n] = size(input_image);
canny=input_image;
outputCanny=zeros(m,n);

conv = zeros(5,5);%高斯卷积核
sigma = 1;%方差
sigma_2 = sigma * sigma;%临时变量
sum_ = 0;
for i = 1:5
    for j = 1:5
        conv(i,j) = exp((-(i - 3) * (i - 3) - (j - 3) * (j - 3)) / (2 * sigma_2)) / (2 * 3.14 * sigma_2);%高斯公式
        sum_ = sum_ + conv(i,j);
    end
end
conv = conv./sum_;%标准化

%对图像实施高斯滤波
for i = 1:m
    for j = 1:n
        sum_ = 0;%临时变量
        for k = 1:5
            for t = 1:5
                if (i - 3 + k) > 0 && (i - 3 + k) <= m && (j - 3 + t) > 0 && (j - 3 + t) < n
                    sum_ = sum_ + conv(k,t) * input_image(i - 3 + k,j - 3 + t);
                end
            end
        end
        canny(i,j) = sum_;
    end
end
% figure;clf;
% imshow(canny,[]);
% title('高斯滤波后的结果')

%求梯度
dx = zeros(m,n);%x方向梯度
dy = zeros(m,n);%y方向梯度
d = zeros(m,n);
for i = 1:m - 1
    for j = 1:n - 1
        dx(i,j) = canny(i,j + 1) - canny(i,j);
        dy(i,j) = canny(i + 1,j) - canny(i,j);
        d(i,j) = sqrt(dx(i,j) * dx(i,j) + dy(i,j) * dy(i,j));
    end
end
% figure;clf;
% imshow(d,[]);
% title('求梯度后的结果')

%局部非极大值抑制
K = d;%记录进行非极大值抑制后的梯度
for j = 1:n
    K(1,j) = 0;
end
for j = 1:n
    K(m,j) = 0;
end
for i = 2:m - 1
    K(i,1) = 0;
end
for i = 2:m - 1
    K(i,n) = 0;
end

for i = 2:m - 1
    for j = 2:n - 1
        %当前像素点的梯度值为0，则一定不是边缘点
        if d(i,j) == 0
            K(i,j) = 0;
        else
            gradX = dx(i,j);%当前点x方向导数
            gradY = dy(i,j);%当前点y方向导数
            gradTemp = d(i,j);%当前点梯度
            %如果Y方向幅度值较大
            if abs(gradY) > abs(gradX)
                weight = abs(gradX) / abs(gradY);%权重
                grad2 = d(i - 1,j);
                grad4 = d(i + 1,j);
                %如果x、y方向导数符号相同
                %像素点位置关系
                %g1 g2
                %   C
                %   g4 g3
                if gradX * gradY > 0
                    grad1 = d(i - 1,j - 1);
                    grad3 = d(i + 1,j + 1);
                else
                    %如果x、y方向导数符号反
                    %像素点位置关系
                    %   g2 g1
                    %   C
                    %g3 g4
                    grad1 = d(i - 1,j + 1);
                    grad3 = d(i + 1,j - 1);
                end
            %如果X方向幅度值较大
            else
                weight = abs(gradY) / abs(gradX);%权重
                grad2 = d(i,j - 1);
                grad4 = d(i,j + 1);
                %如果x、y方向导数符号相同
                %像素点位置关系
                %g3
                %g4 C g2
                %     g1
                if gradX * gradY > 0
                    grad1 = d(i + 1,j + 1);
                    grad3 = d(i - 1,j - 1);
                else
                    %如果x、y方向导数符号反
                    %像素点位置关系
                    %     g1
                    %g4 C g2
                    %g3
                    grad1 = d(i - 1,j + 1);
                    grad3 = d(i + 1,j - 1);
                end
            end
            %利用grad1-grad4对梯度进行插值
            gradTemp1 = weight * grad1 + (1 - weight) * grad2;
            gradTemp2 = weight * grad3 + (1 - weight) * grad4;
            %当前像素的梯度是局部的最大值，可能是边缘点
            if gradTemp >= gradTemp1 && gradTemp >= gradTemp2
                K(i,j) = gradTemp;
            else
                %不可能是边缘点
                K(i,j) = 0;
            end
        end
    end
end
% figure;clf;
% imshow(K,[]);
% title('非极大值抑制后的结果')

% 定义双阈值：EP_MIN、EP_MAX，且EP_MAX = 2 * EP_MIN
% EP_MIN = 0.02;
% EP_MAX = EP_MIN * 2.5;
% for i = 1:m
%     for j = 1:n
%         if K(i,j) >= EP_MAX%小于小阈值，不可能为边缘点
%             outputCanny(i,j) = 255;
%         end
%     end
% end

%定义双阈值：EP_MIN、EP_MAX，且EP_MAX = 2 * EP_MIN
EP_MIN = 0.02;
EP_MAX = EP_MIN * 2;
EdgeLarge = zeros(m,n);%记录真边缘
EdgeBetween = zeros(m,n);%记录可能的边缘点
for i = 1:m
    for j = 1:n
        if K(i,j) >= EP_MAX%小于小阈值，不可能为边缘点
            EdgeLarge(i,j) = K(i,j);
        else if K(i,j) >= EP_MIN
                EdgeBetween(i,j) = K(i,j);
            end
        end
    end
end

%把EdgeLarge的边缘连成连续的轮廓
MAXSIZE = 999999;
Queue = zeros(MAXSIZE,2);%用数组模拟队列
front = 1;%队头
rear = 1;%队尾
for i = 1:m
    for j = 1:n
        if EdgeLarge(i,j) > 0
            %强点入队
            Queue(rear,1) = i;
            Queue(rear,2) = j;
            rear = rear + 1;
            outputCanny(i,j) = 255;
            EdgeLarge(i,j) = 0;%避免重复计算
        end
        while front ~= rear%队不空
            %队头出队
            temp_i = Queue(front,1);
            temp_j = Queue(front,2);
            front = front + 1;
            %8-连通域寻找可能的边缘点
            %左上方
            if EdgeBetween(temp_i - 1,temp_j - 1) > 0%把在强点周围的弱点变为强点
                EdgeLarge(temp_i - 1,temp_j - 1) = K(temp_i - 1,temp_j - 1);
                EdgeBetween(temp_i - 1,temp_j - 1) = 0;%避免重复计算
                %入队
                Queue(rear,1) = temp_i - 1;
                Queue(rear,2) = temp_j - 1;
                rear = rear + 1;
            end
            %正上方
            if EdgeBetween(temp_i - 1,temp_j) > 0%把在强点周围的弱点变为强点
                EdgeLarge(temp_i - 1,temp_j) = K(temp_i - 1,temp_j);
                EdgeBetween(temp_i - 1,temp_j) = 0;
                %入队
                Queue(rear,1) = temp_i - 1;
                Queue(rear,2) = temp_j;
                rear = rear + 1;
            end
            %右上方
            if EdgeBetween(temp_i - 1,temp_j + 1) > 0%把在强点周围的弱点变为强点
                EdgeLarge(temp_i - 1,temp_j + 1) = K(temp_i - 1,temp_j + 1);
                EdgeBetween(temp_i - 1,temp_j + 1) = 0;
                %入队
                Queue(rear,1) = temp_i - 1;
                Queue(rear,2) = temp_j + 1;
                rear = rear + 1;
            end
            %正左方
            if EdgeBetween(temp_i,temp_j - 1) > 0%把在强点周围的弱点变为强点
                EdgeLarge(temp_i,temp_j - 1) = K(temp_i,temp_j - 1);
                EdgeBetween(temp_i,temp_j - 1) = 0;
                %入队
                Queue(rear,1) = temp_i;
                Queue(rear,2) = temp_j - 1;
                rear = rear + 1;
            end
            %正右方
            if EdgeBetween(temp_i,temp_j + 1) > 0%把在强点周围的弱点变为强点
                EdgeLarge(temp_i,temp_j + 1) = K(temp_i,temp_j + 1);
                EdgeBetween(temp_i,temp_j + 1) = 0;
                %入队
                Queue(rear,1) = temp_i;
                Queue(rear,2) = temp_j + 1;
                rear = rear + 1;
            end
            %左下方
            if EdgeBetween(temp_i + 1,temp_j - 1) > 0%把在强点周围的弱点变为强点
                EdgeLarge(temp_i + 1,temp_j - 1) = K(temp_i + 1,temp_j - 1);
                EdgeBetween(temp_i + 1,temp_j - 1) = 0;
                %入队
                Queue(rear,1) = temp_i + 1;
                Queue(rear,2) = temp_j - 1;
                rear = rear + 1;
            end
            %正下方
            if EdgeBetween(temp_i + 1,temp_j) > 0%把在强点周围的弱点变为强点
                EdgeLarge(temp_i + 1,temp_j) = K(temp_i + 1,temp_j);
                EdgeBetween(temp_i + 1,temp_j) = 0;
                %入队
                Queue(rear,1) = temp_i + 1;
                Queue(rear,2) = temp_j;
                rear = rear + 1;
            end
            %右下方
            if EdgeBetween(temp_i + 1,temp_j + 1) > 0%把在强点周围的弱点变为强点
                EdgeLarge(temp_i + 1,temp_j + 1) = K(temp_i + 1,temp_j + 1);
                EdgeBetween(temp_i + 1,temp_j + 1) = 0;
                %入队
                Queue(rear,1) = temp_i + 1;
                Queue(rear,2) = temp_j + 1;
                rear = rear + 1;
            end
        end
    end
end

figure;clf;
imshow(outputCanny,[]);
title('Canny')


%%
%对比图
figure;clf;
subplot(2,4,1), imshow(input_image);
title('原图')
subplot(2,4,2), imshow(outputRoberts);
title('roberts')
subplot(2,4,3), imshow(outputPrewitt);
title('prewitt');
subplot(2,4,4), imshow(outputSobel);
title('sobel');
subplot(2,4,5), imshow(outputLaplacian);
title('laplacian');
subplot(2,4,6), imshow(e);
title('marr_1');
subplot(2,4,7), imshow(outputMarr);
title('marr');
subplot(2,4,8), imshow(outputCanny);
title('canny');

%%
output=outputRoberts;
% output=outputPrewitt;
% output=outputSobel;
% output=outputLaplacian;
% output=e;
% output=outputMarr;
% output=outputCanny;
