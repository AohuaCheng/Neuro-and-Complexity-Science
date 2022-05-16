function [X1] = compute_X1(sigma_w,sigma_b)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
q11 = iteration_q(sigma_w,sigma_b);
X1 = get_X1(q11,sigma_w);
end

