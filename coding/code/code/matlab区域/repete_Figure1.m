%绘制图Figure 1.A，展示模长收敛方向
x=0:0.05:15;
lenx = length(x);
a = zeros(1,lenx);
b = zeros(1,lenx);
c = zeros(1,lenx);
for i = 1:lenx
    a(i)=iteration_q(4.0,0,x(i),1);
    b(i)=iteration_q(2.5,0,x(i),1);
    c(i)=iteration_q(1.3,0,x(i),1);
end
figure(1)
plot(x,a,'r-',x,b,'g-',x,c,'b-',x,x,'k:')
title('sigma_b=0,sigma_w=4.0-red,2.5-green,1.3-blue')
xlabel('迭代中输入的i-1值')
ylabel('迭代中输出的i值')

%绘制Figure 1.B，展示模长的收敛速度
a = zeros(3,7);
a(1,1) = 12;
a(2,1) = 12;
a(3,1) = 12;
b = zeros(3,7);
b(1,1) = 8;
b(2,1) = 8;
b(3,1) = 8;
c = zeros(3,7);
c(1,1) = 4;
c(2,1) = 4;
c(3,1) = 4;

for i = 2:7
    a(1,i) = iteration_q(1.3,0,a(1,i-1),1);
    a(2,i) = iteration_q(2.5,0,a(2,i-1),1);
    a(3,i) = iteration_q(4.0,0,a(3,i-1),1);
    b(1,i) = iteration_q(1.3,0,b(1,i-1),1);
    b(2,i) = iteration_q(2.5,0,b(2,i-1),1);
    b(3,i) = iteration_q(4.0,0,b(3,i-1),1);
    c(1,i) = iteration_q(1.3,0,c(1,i-1),1);
    c(2,i) = iteration_q(2.5,0,c(2,i-1),1);
    c(3,i) = iteration_q(4.0,0,c(3,i-1),1);
end
x=0:6;
figure(2);
plot(x,a(1,:),'b',x,a(2,:),'g',x,a(3,:),'r',x,b(1,:),'b',x,b(2,:),'g',x,b(3,:),'r',x,c(1,:),'b',x,c(2,:),'g',x,c(3,:),'r')
xlabel('迭代的次数')
ylabel('模长的数值')
title('sigma_b=0,sigma_w=4.0-red,2.5-green,1.3-blue')

%复现Figure 1.C和Figure 1.D，更全面展示收敛情况
sigma_b = 0:0.05:4;
sigma_w = 0:0.05:5;
lenb = length(sigma_b);
lenw = length(sigma_w);
q_fixed = ones(lenw,lenb);
error = ones(lenw,lenb);
steps = zeros(lenw,lenb);
for i = 1:lenw
    for j = 1:lenb
        [q_fixed(i,j),error(i,j),steps(i,j)] = iteration_q(sigma_w(i),sigma_b(j));
    end
end
figure(3)
surf(sigma_b,sigma_w,q_fixed)
xlabel('b的方差')
ylabel('w的方差')
figure(4)
surf(sigma_b,sigma_w,steps)%图4在sigma_b=0,sigma_w=1附近有突变，收敛速度急剧变慢。
xlabel('b的方差')
ylabel('w的方差')





