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
robertsNum=0; %��roberts���Ӽ���õ���ÿ�����ص�ֵ
for j=1:m-1 
    for k=1:n-1
        robertsNum = abs(input_image(j,k)-input_image(j+1,k+1)) + abs(input_image(j+1,k)-input_image(j,k+1));
        roberts(j,k)=robertsNum;
    end
end
roberts = mat2gray(roberts);  % ���ݶȾ���ת��Ϊ�Ҷ�ͼ��
robertThreshold=graythresh(roberts);%?����Ҷ���ֵ
for j=1:m-1 %���б߽���ȡ
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
PrewittNum=0; %��Prewitt���Ӽ���õ���ÿ�����ص�ֵ
for j=2:m-1 %���б߽���ȡ
    for k=2:n-1
        PrewittNum=abs(input_image(j-1,k+1)-input_image(j+1,k+1)+input_image(j-1,k)-input_image(j+1,k)+input_image(j-1,k-1)-input_image(j+1,k-1))+abs(input_image(j-1,k+1)+input_image(j,k+1)+input_image(j+1,k+1)-input_image(j-1,k-1)-input_image(j,k-1)-input_image(j+1,k-1));
        Prewitt(j,k)=PrewittNum;
    end
end
Prewitt = mat2gray(Prewitt);  % ���ݶȾ���ת��Ϊ�Ҷ�ͼ��
PrewittThreshold=graythresh(Prewitt);%?����Ҷ���ֵ
for j=1:m-1 %���б߽���ȡ
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
sobelNum=0; %��sobel���Ӽ���õ���ÿ�����ص�ֵ
for j=2:m-1 %���б߽���ȡ
    for k=2:n-1
        sobelNum=abs(input_image(j-1,k+1)+2*input_image(j,k+1)+input_image(j+1,k+1)-input_image(j-1,k-1)-2*input_image(j,k-1)-input_image(j+1,k-1))+abs(input_image(j-1,k-1)+2*input_image(j-1,k)+input_image(j-1,k+1)-input_image(j+1,k-1)-2*input_image(j+1,k)-input_image(j+1,k+1));
        sobel(j,k)=sobelNum;
    end
end

sobel = mat2gray(sobel);  % ���ݶȾ���ת��Ϊ�Ҷ�ͼ��
sobelThreshold=graythresh(sobel);%?����Ҷ���ֵ

for j=1:m-1 %���б߽���ȡ
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
LaplacianNum=0; %��Laplician���Ӽ���õ���ÿ�����ص�ֵ
for j=2:m-1 %���б߽���ȡ
    for k=2:n-1
        LaplacianNum=abs(4*input_image(j,k)-input_image(j-1,k)-input_image(j+1,k)-input_image(j,k+1)-input_image(j,k-1));
        Laplacian(j,k)=LaplacianNum;
    end
end

Laplacian = mat2gray(Laplacian);  % ���ݶȾ���ת��Ϊ�Ҷ�ͼ��
LaplacianThreshold=graythresh(Laplacian);%?����Ҷ���ֵ

for j=1:m-1 %���б߽���ȡ
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
MarrThreshold=graythresh(abs(Marr));%?����Ҷ���ֵ
outputMarr=zeros(m,n);

for j=2:m-1 %���б߽���ȡ
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
%��˹�˲�
[m,n] = size(input_image);
canny=input_image;
outputCanny=zeros(m,n);

conv = zeros(5,5);%��˹�����
sigma = 1;%����
sigma_2 = sigma * sigma;%��ʱ����
sum_ = 0;
for i = 1:5
    for j = 1:5
        conv(i,j) = exp((-(i - 3) * (i - 3) - (j - 3) * (j - 3)) / (2 * sigma_2)) / (2 * 3.14 * sigma_2);%��˹��ʽ
        sum_ = sum_ + conv(i,j);
    end
end
conv = conv./sum_;%��׼��

%��ͼ��ʵʩ��˹�˲�
for i = 1:m
    for j = 1:n
        sum_ = 0;%��ʱ����
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
% title('��˹�˲���Ľ��')

%���ݶ�
dx = zeros(m,n);%x�����ݶ�
dy = zeros(m,n);%y�����ݶ�
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
% title('���ݶȺ�Ľ��')

%�ֲ��Ǽ���ֵ����
K = d;%��¼���зǼ���ֵ���ƺ���ݶ�
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
        %��ǰ���ص���ݶ�ֵΪ0����һ�����Ǳ�Ե��
        if d(i,j) == 0
            K(i,j) = 0;
        else
            gradX = dx(i,j);%��ǰ��x������
            gradY = dy(i,j);%��ǰ��y������
            gradTemp = d(i,j);%��ǰ���ݶ�
            %���Y�������ֵ�ϴ�
            if abs(gradY) > abs(gradX)
                weight = abs(gradX) / abs(gradY);%Ȩ��
                grad2 = d(i - 1,j);
                grad4 = d(i + 1,j);
                %���x��y������������ͬ
                %���ص�λ�ù�ϵ
                %g1 g2
                %   C
                %   g4 g3
                if gradX * gradY > 0
                    grad1 = d(i - 1,j - 1);
                    grad3 = d(i + 1,j + 1);
                else
                    %���x��y���������ŷ�
                    %���ص�λ�ù�ϵ
                    %   g2 g1
                    %   C
                    %g3 g4
                    grad1 = d(i - 1,j + 1);
                    grad3 = d(i + 1,j - 1);
                end
            %���X�������ֵ�ϴ�
            else
                weight = abs(gradY) / abs(gradX);%Ȩ��
                grad2 = d(i,j - 1);
                grad4 = d(i,j + 1);
                %���x��y������������ͬ
                %���ص�λ�ù�ϵ
                %g3
                %g4 C g2
                %     g1
                if gradX * gradY > 0
                    grad1 = d(i + 1,j + 1);
                    grad3 = d(i - 1,j - 1);
                else
                    %���x��y���������ŷ�
                    %���ص�λ�ù�ϵ
                    %     g1
                    %g4 C g2
                    %g3
                    grad1 = d(i - 1,j + 1);
                    grad3 = d(i + 1,j - 1);
                end
            end
            %����grad1-grad4���ݶȽ��в�ֵ
            gradTemp1 = weight * grad1 + (1 - weight) * grad2;
            gradTemp2 = weight * grad3 + (1 - weight) * grad4;
            %��ǰ���ص��ݶ��Ǿֲ������ֵ�������Ǳ�Ե��
            if gradTemp >= gradTemp1 && gradTemp >= gradTemp2
                K(i,j) = gradTemp;
            else
                %�������Ǳ�Ե��
                K(i,j) = 0;
            end
        end
    end
end
% figure;clf;
% imshow(K,[]);
% title('�Ǽ���ֵ���ƺ�Ľ��')

% ����˫��ֵ��EP_MIN��EP_MAX����EP_MAX = 2 * EP_MIN
% EP_MIN = 0.02;
% EP_MAX = EP_MIN * 2.5;
% for i = 1:m
%     for j = 1:n
%         if K(i,j) >= EP_MAX%С��С��ֵ��������Ϊ��Ե��
%             outputCanny(i,j) = 255;
%         end
%     end
% end

%����˫��ֵ��EP_MIN��EP_MAX����EP_MAX = 2 * EP_MIN
EP_MIN = 0.02;
EP_MAX = EP_MIN * 2;
EdgeLarge = zeros(m,n);%��¼���Ե
EdgeBetween = zeros(m,n);%��¼���ܵı�Ե��
for i = 1:m
    for j = 1:n
        if K(i,j) >= EP_MAX%С��С��ֵ��������Ϊ��Ե��
            EdgeLarge(i,j) = K(i,j);
        else if K(i,j) >= EP_MIN
                EdgeBetween(i,j) = K(i,j);
            end
        end
    end
end

%��EdgeLarge�ı�Ե��������������
MAXSIZE = 999999;
Queue = zeros(MAXSIZE,2);%������ģ�����
front = 1;%��ͷ
rear = 1;%��β
for i = 1:m
    for j = 1:n
        if EdgeLarge(i,j) > 0
            %ǿ�����
            Queue(rear,1) = i;
            Queue(rear,2) = j;
            rear = rear + 1;
            outputCanny(i,j) = 255;
            EdgeLarge(i,j) = 0;%�����ظ�����
        end
        while front ~= rear%�Ӳ���
            %��ͷ����
            temp_i = Queue(front,1);
            temp_j = Queue(front,2);
            front = front + 1;
            %8-��ͨ��Ѱ�ҿ��ܵı�Ե��
            %���Ϸ�
            if EdgeBetween(temp_i - 1,temp_j - 1) > 0%����ǿ����Χ�������Ϊǿ��
                EdgeLarge(temp_i - 1,temp_j - 1) = K(temp_i - 1,temp_j - 1);
                EdgeBetween(temp_i - 1,temp_j - 1) = 0;%�����ظ�����
                %���
                Queue(rear,1) = temp_i - 1;
                Queue(rear,2) = temp_j - 1;
                rear = rear + 1;
            end
            %���Ϸ�
            if EdgeBetween(temp_i - 1,temp_j) > 0%����ǿ����Χ�������Ϊǿ��
                EdgeLarge(temp_i - 1,temp_j) = K(temp_i - 1,temp_j);
                EdgeBetween(temp_i - 1,temp_j) = 0;
                %���
                Queue(rear,1) = temp_i - 1;
                Queue(rear,2) = temp_j;
                rear = rear + 1;
            end
            %���Ϸ�
            if EdgeBetween(temp_i - 1,temp_j + 1) > 0%����ǿ����Χ�������Ϊǿ��
                EdgeLarge(temp_i - 1,temp_j + 1) = K(temp_i - 1,temp_j + 1);
                EdgeBetween(temp_i - 1,temp_j + 1) = 0;
                %���
                Queue(rear,1) = temp_i - 1;
                Queue(rear,2) = temp_j + 1;
                rear = rear + 1;
            end
            %����
            if EdgeBetween(temp_i,temp_j - 1) > 0%����ǿ����Χ�������Ϊǿ��
                EdgeLarge(temp_i,temp_j - 1) = K(temp_i,temp_j - 1);
                EdgeBetween(temp_i,temp_j - 1) = 0;
                %���
                Queue(rear,1) = temp_i;
                Queue(rear,2) = temp_j - 1;
                rear = rear + 1;
            end
            %���ҷ�
            if EdgeBetween(temp_i,temp_j + 1) > 0%����ǿ����Χ�������Ϊǿ��
                EdgeLarge(temp_i,temp_j + 1) = K(temp_i,temp_j + 1);
                EdgeBetween(temp_i,temp_j + 1) = 0;
                %���
                Queue(rear,1) = temp_i;
                Queue(rear,2) = temp_j + 1;
                rear = rear + 1;
            end
            %���·�
            if EdgeBetween(temp_i + 1,temp_j - 1) > 0%����ǿ����Χ�������Ϊǿ��
                EdgeLarge(temp_i + 1,temp_j - 1) = K(temp_i + 1,temp_j - 1);
                EdgeBetween(temp_i + 1,temp_j - 1) = 0;
                %���
                Queue(rear,1) = temp_i + 1;
                Queue(rear,2) = temp_j - 1;
                rear = rear + 1;
            end
            %���·�
            if EdgeBetween(temp_i + 1,temp_j) > 0%����ǿ����Χ�������Ϊǿ��
                EdgeLarge(temp_i + 1,temp_j) = K(temp_i + 1,temp_j);
                EdgeBetween(temp_i + 1,temp_j) = 0;
                %���
                Queue(rear,1) = temp_i + 1;
                Queue(rear,2) = temp_j;
                rear = rear + 1;
            end
            %���·�
            if EdgeBetween(temp_i + 1,temp_j + 1) > 0%����ǿ����Χ�������Ϊǿ��
                EdgeLarge(temp_i + 1,temp_j + 1) = K(temp_i + 1,temp_j + 1);
                EdgeBetween(temp_i + 1,temp_j + 1) = 0;
                %���
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
%�Ա�ͼ
figure;clf;
subplot(2,4,1), imshow(input_image);
title('ԭͼ')
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
