clear;

load solutions.mat;

% gblup - solution of GBLUP
% sstep - solution of SSTEP
% snp   - solution to SSTEP-SNP

%%
sol_snpsstep = zeros(size(sstep)); % sstep solutions only for genotyped ids

gid = sstep(:,1);

not_found = 0;
for i = 1:size(gid,1) % extract solutions
    f1 = find( snpsstep(:,1) == gid(i,1) );
    sol_snpsstep(i,:) = snpsstep(f1,:);
end

corr(sol_snpsstep(:,2), sstep(:,2))

figure;
plot(sol_snpsstep(:,2), sstep(:,2), '*')

%%
sol_snpsstep2 = zeros(size(gblup)); % sstep solutions only for genotyped ids

gid = gblup(:,1);

not_found = 0;
for i = 1:size(gid,1) % extract solutions
    f1 = find( snpsstep(:,1) == gid(i,1) );
    sol_snpsstep2(i,:) = snpsstep(f1,:);
end

corr(sol_snpsstep2(:,2), gblup(:,2))

figure;
plot(sol_snpsstep2(:,2), gblup(:,2), '*')
%%
corr(sol_snpsstep(:,2), sstep(:,2))

figure;
plot(sol_snpsstep(:,2), sstep(:,2), '*')
% corr(sol_snp(:,2), gblup(:,2))

% corr(snpblup, gblup(:,2))
% corr(sol_sstep(:,2), snpblup)
% corr(sol_sstep(:,2), gblup(:,2))
%%
figure;
hold on
%subplot(2,1,1)
plot(sol_sstep(:,1), sol_sstep(:,2), 'o', 'MarkerFaceColor','b')
%subplot(2,1,2)
plot(gblup(:,1), gblup(:,2), 'o', 'MarkerFaceColor','r')

%% get some statistiocs

inum = 0;
for i = 1:59200
    if obs(i,1) ~= -9999
        inum = inum + 1;
    end
end
disp(["non missing obs =", num2str(inum)]);

inum = 0;
for i = 1:59200
    if obs_id(i,1) ~= -9999
        inum = inum + 1;
    end
end
disp(["non missing ids =", num2str(inum)]);

inum = 0;
for i = 1:59200
    if obs_ngid(i,1) ~= -9999
        inum = inum + 1;
    end
end
disp(["non missing ngids =", num2str(inum)]);

%%

inum = 0;
for i = 1:size(gid,1)
    f = find( obs_id == gid(i) );
    if ( obs(f) ~= -9999 )
        inum = inum + 1;
    end
end
disp(["In obs: non missing gids =", num2str(inum)]);

inum = 0;
for i = 1:size(ngid,1)
    f = find( obs_id == ngid(i) );
    if ( obs(f) ~= -9999 )
        inum = inum + 1;
    end
end
disp(["In obs: non missing ngids =", num2str(inum)]);

inum = 0;
for i = 1:size(refid,1)
    f = find( obs_id == refid(i) );
    if ( obs(f) ~= -9999 )
        inum = inum + 1;
    end
end
disp(["In obs: non missing refid =", num2str(inum)]);