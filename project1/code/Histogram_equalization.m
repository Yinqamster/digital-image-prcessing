function [output] = Histogram_equalization(input_image)
%first test the image is a RGB or gray image
if numel(size(input_image)) == 3
    %this is a RGB image
    %here is just one method, if you have other ways to do the
    %equalization, you can change the following code
    r=input_image(:,:,1);
    v=input_image(:,:,2);
    b=input_image(:,:,3);
%     r1 = histeq(r);
%     v1 = histeq(v);
%     b1 = histeq(b);

%     r1 = adapthisteq(r);
%     v1 = adapthisteq(v);
%     b1 = adapthisteq(b);

%     r1 = imadjust(r);
%     v1 = imadjust(v);
%     b1 = imadjust(b);
    
    r1 = hist_equal(r);
    v1 = hist_equal(v);
    b1 = hist_equal(b);
    
%     r1 = hist_equal1(r);
%     v1 = hist_equal1(v);
%     b1 = hist_equal1(b);
    
%     r1 = hist_equal2(r);
%     v1 = hist_equal2(v);
%     b1 = hist_equal2(b);

    output = cat(3,r1,v1,b1);    
else
    %this is a gray image
%     [output] = hist_equal(input_image);
%   [output] = hist_equal1(input_image);
  [output] = hist_equal2(input_image);
%   [output] = adapthisteq(input_image);
%   [output] = imadjust(input_image);
    
end

    %ֱ��ͼ���⻯  imhist����
    function [output2] = hist_equal(input_channel)
        %you should complete this sub-function
        input_channel = im2double(input_channel);
        output2 = input_channel;
        [M,N]=size(input_channel);
        [counts,x]=imhist(input_channel);%H�Ƕ�ȡ��ͼ��imhist�Ƕ�ͼ��ֱ��ͼ����ͳ�ƣ�����count����ÿ���Ҷ�ֵ�ø�����x����Ҷ�ֵ
        location=find(counts~=0);%�ҵ��������ظ�����Ϊ0�ĻҶȼ�
        MinCDF=min(counts(location));%�ҵ������������ٵĻҶȼ�
        for  j=1:length(location)
          CDF=sum(counts(location(1:j)));%��������Ҷȼ����ظ����ۼƷֲ�
          P= input_channel==x(location(j));%�ҵ�ͼ���е���ĳ���Ҷȼ��������ص�����λ��
          output2(P)=(CDF-MinCDF)/(M*N-MinCDF);%%���ûҶȻ��㹫ʽ���޸�����λ���ϵ�����ֵ
        end
    end

    %�ֶ�ʵ��ֱ��ͼ���⻯
    function [output3] = hist_equal1(input_channel)
        %you should complete this sub-function
        [R, C, K] = size(input_channel); % ������K��ʾ��ɫͨ����
        % output3 = input_channel;
        % ͳ��ÿ������ֵ���ִ���
        cnt = zeros(K, 256);
        for i = 1 : R
            for j = 1 : C
                for k = 1 : K
                    cnt(k, input_channel(i, j, k) + 1) = cnt(k, input_channel(i, j, k) + 1) + 1;
                end
            end
        end
        f = zeros(3, 256);
        f = double(f); cnt = double(cnt);
        % ͳ��ÿ������ֵ���ֵĸ��ʣ� �õ�����ֱ��ͼ
        for k = 1 : K
            for i = 1 : 256
                f(k, i) = cnt(k, i) / (R * C);
            end
        end
        % ���ۼƸ��ʣ��õ��ۼ�ֱ��ͼ
        for k = 1 : K
            for i = 2 : 256
                f(k, i) = f(k, i - 1) + f(k, i);
            end
        end
        % ��f����ʵ������ֵ[0, 255]��ӳ�䡣 
        for k = 1 : K
            for i = 1 : 256
                f(k, i) = f(k, i) * 255;
            end
        end
        % ���ÿ�����ص��ӳ��
        output3 = input_channel;
        % output3 = double(output3);
        for i = 1 : R
            for j = 1 : C
                for k = 1 : K
                    output3(i, j, k) = f(k, output3(i, j, k) + 1);
                end
            end
        end    
    end

    %������
    function [output4] = hist_equal2(input_channel)
        input_channel=im2double(input_channel);%ת���������ͣ���uint8ͼ��תΪdouble���ͣ���ΧΪ0-1
        [height, width, R]=size(input_channel);%���ؾ���I������
     %   output = input_channel;
        for i=2:height-1
            for j=2:width-1
                R(i,j)=abs(input_channel(i+1,j+1)-input_channel(i,j))+abs(input_channel(i+1,j)-input_channel(i,j+1));
            end
        end
        T=R;
        for i=1:height-1
            for j=1:width-1
                if (R(i,j)<0.25)
                    R(i,j)=1;
                else
                    R(i,j)=0;
                end
            end
        end
        
        [m,n]=size(T);
        output4(1:m,1:n)=input_channel(1:m,1:n)+T(1:m,1:n);
    end
end
