
% gblup - solution of GBLUP
% sstep - solution of SSTEP

sol_sstep = zeros(size(gblup)); % sstep solutions only for genotyped ids

for i = 1:size(gid,1) % extract solutions
    f = find( sstep(:,1) == gid(i,1) );
    sol_sstep(i,:) = sstep(f,:);
end

corr(sol_sstep(:,2), gblup(:,2))