% clear
% 
% seed = 1011;
% rng(seed,"twister");
% s = rng;


%
%rng(s);
clear
rng("default");
% GENERATE POPULATION

n_animals = 250;
ref_all_prop = 0.7;
n_snp = 10000;
n_ploidy = 2;

[ haplotypes, genotypes ] = sim_genotypes( n_animals, n_snp, n_ploidy, ref_all_prop );

% -------------------------------------------------------------------------
% correlation matrix
r = [ 1.0 0.4 -0.2 0.3;
      0.4 1.0 -0.6 -0.4;
      -0.2 -0.6 1.0 0.1;
      0.3 -0.4 0.1 1.0 ];

% generate positions and effects
n_qtl = 100;

opt_rounds = 100;
opt_res2 = zeros(opt_rounds,1);

alpha_min = 1;
alpha_max = 1.5;
alpha = alpha_min + ( alpha_max - alpha_min ).*rand( opt_rounds,1 );

best_res = 10;
best_alpha = 0.0;

for opt = 1:opt_rounds % start of optimization

    qtl_wnd = floor( alpha(opt,1)*(n_snp / n_qtl) );
    
    loci = sim_positions( n_qtl, n_snp, r, qtl_wnd );
    a_map = sim_effects_map( r, loci );
    
    % calculate trait
    
    t2 = zeros(n_animals,size(r,1));
    
    for l = 1:n_animals % loop over individuals
        for i = 1:size(r,1) % loop over traits
            for j = 1:n_qtl % loop over QTLs
                iqtl = loci(j,i);
                effects = a_map(iqtl);
                t2(l,i) = t2(l,i) + genotypes( l,iqtl ) * effects( i );
            end
        end
    end
    
    opt_res2(opt,1) = norm( r * inv(corr(t2)) ); % quality measure
    
    if opt_res2(opt,1) < best_res
        best_res = opt_res2(opt,1);
        best_alpha = alpha(opt,1);
        a_best = a_map;
        loci_best = loci;
    end

end % end of optimization

% Recalculate trait with best parameter

t2 = zeros(n_animals,size(r,1));

for l = 1:n_animals % loop over individuals
    for i = 1:size(r,1) % loop over traits
        for j = 1:n_qtl % loop over QTLs
            iqtl = loci_best(j,i);
            effects = a_best(iqtl);
            t2(l,i) = t2(l,i) + genotypes( l,iqtl ) * effects( i );
        end
    end
end

disp("uniformly distributed positions, corr(loci):")
disp( corr(loci) );

% disp("normaly distributed effects, corr(a):")
% disp( corr(a) );

disp("Requested trait correlation, r:")
disp( r );

disp("Trait correlation, corr(t2):")
disp( corr(t2) );

figure(1);
clf;
hold on
plot(alpha, opt_res2, 'o','MarkerFaceColor','b')
hold off

figure(2);
clf;
subplot(1,6,1)
plot(t2(:,1), t2(:,2), 'o')
subplot(1,6,2)
plot(t2(:,1), t2(:,3), 'o')
subplot(1,6,3)
plot(t2(:,1), t2(:,4), 'o')
subplot(1,6,4)
plot(t2(:,2), t2(:,3), 'o')
subplot(1,6,5)
plot(t2(:,2), t2(:,4), 'o')
subplot(1,6,6)
plot(t2(:,3), t2(:,4), 'o')

%
figure(3);
hold on;
plot(loci(:,1),1, 'o-','MarkerFaceColor','k');
plot(loci(:,2),1.5,'d-','MarkerFaceColor','b');
plot(loci(:,3), 2, 's-','MarkerFaceColor','r');
plot(loci(:,3), 2.5, '*-','MarkerFaceColor','g');
hold off;
grid on
ylim([0 4])
% -------------------------------------------------------------------------

function loci = sim_positions( n_qtl, n_snp, correlation, wnd )
    % sample genes positions
    loci_n = randn(n_qtl,size(correlation,1));
    % make positions correlated
    loci_n = loci_n * chol(correlation);
    % convert sampled values from normal to uniform distribution
    loci_u = cdf('Normal',loci_n,0,1);
    % rescale to the range [1 n_snp]
    loci_u = loci_u * (n_snp-1) + 1;
    % move closely located loci (within wnd range) to the same positions    
    if wnd ~= 0
        loci = ( ceil( loci_u./wnd ) ).*wnd;
    else
        loci = floor(loci_u);
    end
    % if any locus index higher than max n_snp, adjust
    for i = 1:size(correlation,1)
        ii = find (loci(:,i) > n_snp);
        loci(ii,i) = n_snp;
    end
end

function [eff_map,eff] = sim_effects_map( correlation, loci )
    u = unique(loci(:,1));
    for i = 2:size(loci,2)
        iu = unique(loci(:,i));
        u = cat(1,u,iu);
    end
    u = unique(u);
    
    % sample qtl from normal distribution
    eff = randn( size(u,1),size(correlation,1) );
    
    % scale qtl positions according to correlation matrix
    eff = eff * chol( correlation );
        
    eff_cell = cell(size(u,1),1);
    for i = 1:size(u,1)
        eff_cell{i,1} = eff(i,:);
    end
   
    eff_map = containers.Map(u,eff_cell);
end

function [ haplotypes, genotypes ] = sim_genotypes( n_animals, n_snp, n_ploidy, ref_all_prop )

    haplotypes = zeros(n_animals*n_ploidy, n_snp);
    genotypes = zeros(n_animals, n_snp);
    
    % Simulate haplotypes
    for i = 1:n_animals*n_ploidy
        for j = 1:n_snp
            rnum = randi([1 100],1,1);
            if rnum <= ref_all_prop*100
                val = 1;
            else
                val = 0;
            end
            haplotypes(i,j) = val;
        end
    end
    % Calculate genotypes based from haplotypes
    i2 = 0;
    for i = 1:n_animals
        for j = 1:n_snp
            for l = 0:n_ploidy-1
                genotypes(i,j) = genotypes(i,j) + haplotypes(i+i2+l,j);
            end
            % in case of scaling:
            %genotypes(i,j) = genotypes(i,j) - 1; % only for diploids
        end
        i2 = i2 + n_ploidy-1;
    end

end
