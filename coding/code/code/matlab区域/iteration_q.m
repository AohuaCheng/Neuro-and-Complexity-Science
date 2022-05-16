function [qo,erroro,stepso] = iteration_q(weight_sigma,bias_sigma,q0,steps)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
l = length(weight_sigma);
qo = zeros(l,1);
erroro = zeros(l,1);
stepso = zeros(l,1);
if nargin<4
for i = 1:l
    error = 100;
    if nargin ==2
        q0=0.5;
    end
    sp = 0;
    while abs(error)>0.000005 %|| steps<=10
        funx = @(x)(1./sqrt(2*pi)).*exp(-0.5.*x.^2).*( tanh(sqrt(q0).*x) ).^2;%(exp(sqrt(q).*x)-exp(-sqrt(q).*x))./(exp(sqrt(q).*x)+exp(-sqrt(q).*x)) ).^2;
        qn = weight_sigma(i).^2 * integral(funx,-inf,+inf) + bias_sigma.^2;
        error = qn - q0;
        q0 = qn;
       sp = sp + 1;
    end
    qo(i)=q0;
    erroro(i)=error;
    stepso(i)=sp;
end
elseif nargin==4
for i = 1:l
    error = 100;
    sp = 0;
    while sp<steps
        funx = @(x)(1./sqrt(2*pi)).*exp(-0.5.*x.^2).*( tanh(sqrt(q0).*x) ).^2;%(exp(sqrt(q).*x)-exp(-sqrt(q).*x))./(exp(sqrt(q).*x)+exp(-sqrt(q).*x)) ).^2;
        qn = weight_sigma(i).^2 * integral(funx,-inf,+inf) + bias_sigma.^2;
        error = qn - q0;
        q0 = qn;
       sp = sp + 1;
    end
    qo(i)=q0;
    erroro(i)=error;
    stepso(i)=sp;
end
end

end

