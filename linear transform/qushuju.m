function [m1,xv,yv]=qushuju(m)
%QUSHUJU:人工截取数据的某一部分，先将数据通过函数进行显示，然后通过图片进行人机交互来取出数据的某一部分
%输入：
%m：待截取的完整数据
%输出：
%m1:截取的部分数据
%xv,yv：截取数据对应的横纵坐标
[xv,yv] = getline;                            %用鼠标在当前坐标中选出一些线段，[xv,yv]为这些线段对应的横纵坐标
xv=[xv;xv(1)];yv=[yv;yv(1)];                  %将最后一个点与第一个点连接起来，形成一个封闭的多边形区域
[nt,nx]=size(m);                              %输入矩阵的规模
x=[]; 
y=[];

%产生网格
for i=1:nt                                    
    x=[x;[1:nx]'];                            %网格中每个元素对应的横坐标
end
for j=0:nt-1
    y=[y;j*ones(nx,1)];                       %网格中每个元素对应的纵坐标
end

in = inpolygon(x,y,xv,yv);                    %判断数据中哪些值落在多边形区域之内，落在多边形区域之内的数对应的in的值为1，否则为0
x1=x(in);                                     
y1=y(in);                                     
length_xy=length(x1);
m1=zeros(size(m));
for i=1:length_xy
    m1(y1(i)+1,x1(i))=m(y1(i)+1,x1(i));       %截取的数据矩阵
end