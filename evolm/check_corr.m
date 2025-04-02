clear;

load solutions.mat;

% gblup - solution of GBLUP
% sstep - solution of SSTEP
% snp   - solution to SSTEP-SNP

%sol_sstep = zeros(size(gblup)); % sstep solutions only for genotyped ids
%sol_snp = zeros(size(gblup)); % sstep solutions only for genotyped ids

%gid = gblup(:,1);
gid = sstep(:,1);

not_found = 0;
for i = 1:size(gid,1) % extract solutions
    f1 = find( sstep(:,1) == gid(i,1) );
    f2 = find( snp(:,1) == gid(i,1) );
    %sol_sstep(i,:) = sstep(f1,:);
    if isempty(f2)
        %sol_snp(i,:) = 0;
        disp(gid(i,1));
        not_found = not_found + 1;
    else
        sol_sstep(i,:) = sstep(f1,:);
        sol_snp(i,:) = snp(f2,:);
    end
end

% corr(sol_sstep(:,2), gblup(:,2))
% corr(sol_snp(:,2), gblup(:,2))

corr(sol_sstep(:,2), sol_snp(:,2))