load matr.mat

g = [1 2 4]';
m = [3 5 6]';

Amm = zeros(numel(g));
Amg = zeros(numel(m), numel(g));

for i = 1:numel(m)
    for j = 1:numel(m)
        Amm(i,j) = A(m(i), m(j));
    end
end
for i = 1:numel(m)
    for j = 1:numel(g)
        Amg(i,j) = A(m(i), g(j));
    end
end

%%
for i = 1:10
    p(1,i) = mean(snp(:,i));
end
for i = 1:10
    z(:,i) = snp(:,i) - p(1,i);
end

%%
Mm = -inv(Amm)*Amg*snp
Mm2 = -inv(Amm)*Amg*z

for i = 1:10
    z2(:,i) = Mm(:,i) - p(1,i);
end
