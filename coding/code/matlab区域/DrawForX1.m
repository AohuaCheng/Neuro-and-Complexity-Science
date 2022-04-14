sigma_b = 0:0.05:4;
sigma_w = 0:0.05:5;
lenb = length(sigma_b);
lenw = length(sigma_w);
q_fixed = 0.5.*ones(lenw,lenb);
for i = 1:lenw
    for j = 1:lenb
        q_fixed(i,j) = iteration_q(sigma_w(i),sigma_b(j));
    end
end

ki1 = get_X1(q_fixed,sigma_w);
figure(1)
surf(sigma_b,sigma_w,q_fixed)
figure(2)
surf(sigma_b,sigma_w,ki1)
title('X1的二维分布图')
xlabel('b的方差')
ylabel('w的方差')
figure(3)
contour(sigma_b,sigma_w,ki1,[0,0.5,1,1.5,2])
title('X1的等高线图，绿色为X1=1')
xlabel('b的方差')
ylabel('w的方差')
%hold on
%temp = ones(lenw,lenb);
%surf(sigma_b,sigma_w,temp)
