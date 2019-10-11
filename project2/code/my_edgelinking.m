function output = my_edgelinking(binary_image, row, col)
%in this function, you should finish the edge linking utility.
%the input parameters are a matrix of a binary image containing the edge
%information and coordinates of one of the edge points of a obeject
%boundary, you should run this function multiple times to find different
%object boundaries
%the output parameter is a Q-by-2 matrix, where Q is the number of boundary 
%pixels. B holds the row and column coordinates of the boundary pixels.
%you can use different methods to complete the edge linking function
%the better the quality of object boundary and the more the object boundaries, you will get higher scores  

% [m,n]=size(binary_image);
MAXSIZE = 999999;
Queue = zeros(MAXSIZE,2);%用数组模拟队列
front = 1;%队头
rear = 1;%队尾
% EdgeLarge=bwmorph(binary_image,'bridge',Inf);
EdgeLarge=bwmorph(binary_image,'thin',Inf);

% imtool(EdgeLarge);
% size(EdgeLarge)
i=row;
j=col;
pos=1;
count=0;
        if EdgeLarge(i,j) > 0
            Queue(rear,1) = i;
            Queue(rear,2) = j;
            rear = rear + 1;
            output(pos,1) = i;
            output(pos,2) = j;
            pos=pos+1;
%             EdgeLarge(i,j) = 0;%避免重复计算
        end
        while front ~= rear%队不空
            count = count + 1;
            %队头出队
            temp_i = Queue(front,1);
            temp_j = Queue(front,2);
            front = front + 1;
            %8-连通域寻找可能的边缘点
            %左上方
            if EdgeLarge(temp_i - 1,temp_j - 1) > 0
%                 output(temp_i - 1,temp_j - 1) = 255;
                output(pos,1) = temp_i - 1;
                output(pos,2) = temp_j - 1;
                pos=pos+1;
                EdgeLarge(temp_i - 1,temp_j - 1) = 0;%避免重复计算
                %入队
                Queue(rear,1) = temp_i - 1;
                Queue(rear,2) = temp_j - 1;
                rear = rear + 1;
                continue;
            end
            %正上方
            if EdgeLarge(temp_i - 1,temp_j) > 0
%                 output(temp_i - 1,temp_j) = 255;
                output(pos,1) = temp_i - 1;
                output(pos,2) = temp_j;
                pos=pos+1;
                EdgeLarge(temp_i - 1,temp_j) = 0;
                %入队
                Queue(rear,1) = temp_i - 1;
                Queue(rear,2) = temp_j;
                rear = rear + 1;
                continue;
            end
            %右上方
            if EdgeLarge(temp_i - 1,temp_j + 1) > 0
%                 output(temp_i - 1,temp_j + 1) = 255;
                output(pos,1) = temp_i - 1;
                output(pos,2) = temp_j + 1;
                pos=pos+1;
                EdgeLarge(temp_i - 1,temp_j + 1) = 0;
                %入队
                Queue(rear,1) = temp_i - 1;
                Queue(rear,2) = temp_j + 1;
                rear = rear + 1;
                continue;
            end
            %正左方
            if EdgeLarge(temp_i,temp_j - 1) > 0
%                 output(temp_i,temp_j - 1) = 255;
                output(pos,1) = temp_i;
                output(pos,2) = temp_j - 1;
                pos=pos+1;
                EdgeLarge(temp_i,temp_j - 1) = 0;
                %入队
                Queue(rear,1) = temp_i;
                Queue(rear,2) = temp_j - 1;
                rear = rear + 1;
                continue;
            end
            %正右方
            if EdgeLarge(temp_i,temp_j + 1) > 0
%                 output(temp_i,temp_j + 1) = 255;
                output(pos,1) = temp_i;
                output(pos,2) = temp_j + 1;
                pos=pos+1;
                EdgeLarge(temp_i,temp_j + 1) = 0;
                %入队
                Queue(rear,1) = temp_i;
                Queue(rear,2) = temp_j + 1;
                rear = rear + 1;
                continue;
            end
            %左下方
            if EdgeLarge(temp_i + 1,temp_j - 1) > 0
%                 output(temp_i + 1,temp_j - 1) = 255;
                output(pos,1) = temp_i + 1;
                output(pos,2) = temp_j - 1;
                pos=pos+1;
                EdgeLarge(temp_i + 1,temp_j - 1) = 0;
                %入队
                Queue(rear,1) = temp_i + 1;
                Queue(rear,2) = temp_j - 1;
                rear = rear + 1;
                continue;
            end
            %正下方
            if EdgeLarge(temp_i + 1,temp_j) > 0
%                 output(temp_i + 1,temp_j) = 255;
                output(pos,1) = temp_i + 1;
                output(pos,2) = temp_j;
                pos=pos+1;
                EdgeLarge(temp_i + 1,temp_j) = 0;
                %入队
                Queue(rear,1) = temp_i + 1;
                Queue(rear,2) = temp_j;
                rear = rear + 1;
                continue;
            end
            %右下方
            if EdgeLarge(temp_i + 1,temp_j + 1) > 0
%                 output(temp_i + 1,temp_j + 1) = 255;
                output(pos,1) = temp_i + 1;
                output(pos,2) = temp_j + 1;
                pos=pos+1;
                EdgeLarge(temp_i + 1,temp_j + 1) = 0;
                %入队
                Queue(rear,1) = temp_i + 1;
                Queue(rear,2) = temp_j + 1;
                rear = rear + 1;
                continue;
            end
        end
% figure;clf;
% imshow(output);



% Isrc = imread('rubberband_cap.png');
% if ndims(Isrc) == 3
%     I = rgb2gray(Isrc);
% else
%     I = Isrc;
% end
% I = im2bw(I,0.5);
% Isrc=binary_image;
% I = Isrc;
% %外围扩展边界
% [H,W] = size(I);
% I = [zeros(2,W);I;zeros(2,W)];
% I = [zeros(H+4,2),I,zeros(H+4,2)];
% [H,W] = size(I);
% 
% for k = 1:10 %迭代次数
% p = zeros(8,1);
% NzI = zeros(H,W);
% ZI = zeros(H,W);
% P024I = zeros(H,W);
% P246I = zeros(H,W);
% I2 = I;
% for i = 2:H-1
%     for j = 2:W-1
%         if I(i,j) 
%         p(1) = I(i,j+1);
%         p(2) = I(i-1,j+1);
%         p(3) = I(i-1,j);
%         p(4) = I(i-1,j-1);
%         p(5) = I(i,j-1);
%         p(6) = I(i+1,j-1);
%         p(7) = I(i+1,j);
%         p(8) = I(i+1,j+1);
% 
%         NzI(i,j) = sum(p);
%         ZI(i,j) = sum(abs(diff(p)))+abs(p(8)-p(1));
% 
%         P024I(i,j) = p(1) * p(3) * p(5);
%         P246I(i,j) = p(3) * p(5) * p(7);
%         end
%     end
% end
% for i = 3:H-1
%     for j = 3:W-1
%         if I(i,j)
%         if ((NzI(i,j) >= 2) && (NzI(i,j) <= 6)) && ...
%              (ZI(i,j) == 2) && ...
%              ((P024I(i,j)==0)||(ZI(i,j-1)~=2)) && ...
%              ((P246I(i,j)==0)||(ZI(i-1,j)~=2)) 
% 
%          I2(i,j) = 0;  
%         end 
% 
%         end
%     end
% end
% I = I2;
% end
% figure(1);
% subplot(121);
% imshow(Isrc);
% subplot(122);
% imshow(I2);