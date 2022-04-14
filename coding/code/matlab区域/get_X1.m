function [ki1] = get_X1(q11,sigma_w)
%用来获取判断临界线的参数ki1(ki1<1-->stable,ki1>1-->chaos)。
%对输入的q11和sigma_w做批量处理，输出的ki1中不同行代表不同的sigma_w，不同列代表不同的q11。
syms x
syms q
[lenw,lenq] = size(q11);
if lenw~=length(sigma_w)
    error('输入的q11与sigma_w维数不匹配。正确：q11的行数是w，列数是b');
end

%创建作为操作对象的表达式
funt = tanh(sqrt(q).*x);
funt1 = diff(funt);
funt2 =(1./sqrt(2*pi)).*exp(-0.5.*x.^2).*funt1.^2;
funx = matlabFunction(funt2);

ki1 = zeros(lenw,lenq);
for i = 1:lenw
    for j = 1:lenq
        if q11(i,j)==0
            ki1(i,j)=0;
            continue
        end
        ki1(i,j) = sigma_w(i).^2.*integral(@(x)funx(q11(i,j),x),-inf,+inf);
    end
end
end

