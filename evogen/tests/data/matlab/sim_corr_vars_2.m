%% Simulate base population

clear;
rng("default");

disp("NEW ROUND:");

% -------------------------------------------------------------------------
% POPULATION PARAMETERS

n_animals = 500;
ref_all_prop = 0.7;
n_snp = 5500;%1000;
n_ploidy = 2;

% -------------------------------------------------------------------------
% GENERATE POPULATION

[ haplotypes, genotypes ] = sim_genotypes( n_animals, n_snp, n_ploidy, ref_all_prop );

% -------------------------------------------------------------------------
% TRAIT MODEL

% Tr = mu + ploidy_effect * sum_qtl( dom_cases * (1 + k) * a * q ) + sum_qtl( e * q ),
%
% a := qtl effect;
% q := genotype value (number of copies of ref allele);
% k := dominance degree;
%      (i) -1 < k < 1 (sampling from normal dist. or uniform dist.);
%      (ii) 0 < k < 1 (sampling from gamma dist. or uniform dist.).
%
% ASSUMPTIONS:
% In case of correlated traits we assume equal number of QTLs for each
% trait.
% There are two models of trait: in the model 1, all correlated traits have
% (high number) the same trait-responsible loci; traits differ by additive
% and environmental effects, dominance effects are the same by the
% assumption that dominance is the local locus property;
% in the model 2, each trait characterised by specific QTLs and additive
% and environmental effects; correlation between traits realised by shared
% QTLs and the respective effects; in this model the following assumption
% is accepted, (i) all effects (additive, environmental and dominance) are the
% local locus properties, (ii) correlated traits share some number of gene
% regulatory (methabolic, biochemical etc.) pathways, therefore, depending
% on the level of correlation share specific number of QTLs with all
% properties associated with them (effects).

% -------------------------------------------------------------------------
% Summarise input of trait parameters:

n_trait = 3; % number of correlated traits

% (1) Determine the number of QTLs responsible for correlated traits:
snp_prop_tr = 0.65;                   % proportion of all available SNPs

n_qtl = floor( n_snp * snp_prop_tr ); % number of QTLs for each trait

% (2) Expected value of traits in base pop
mu = [ 40.0;   % for the trait_1
       5.0;    % for the trait_2
       0.5; ]; % for the trait_3

% (3) Control of an environment (accounting its changes):
% the range: -1 < env < +1
env = [ 0.0;   % for the trait_1
        0.0;   % for the trait_2
        0.0 ]; % for the trait_3

% (4) Parameters for dominance (only heterozygots)
%
% (i) For uniform distribution, model 1
k_range_U = [-1 1]; % range of k parameter
% (ii) For normal distribution, model 2
k_range_N = [0 0.5]; % mean & std
% (iii) For gamma distribution, mdodel 3
k_range_G = [0.05]; % expected value of k in loci; is '1/b' param. in the distribution

% (4) Traits correlaatiions in terms of genetic effects
corr_G = [ 1.0 0.5 0.7;
           0.5 1.0 0.2;
           0.7 0.2 1.0 ];

% (5) Traits correlaatiions in terms of residual (environmental) effects
corr_E = [ 1.0 0.3 0.5;
           0.3 1.0 0.4;
           0.5 0.4 1.0 ];

% (6) Vector of genetic variances
var_G = [ 100.0; % genetic variance trate_1
          10.0;  % genetic variance trate_2
          0.1 ]; % genetic variance trate_3

% (7) Vector of residual variances
var_E = [ 200.0;  % environmental variance trait_1
          20.0;   % environmental variance trait_2
          0.3  ]; % environmental variance trait_3

disp("Required correlation of additive effects: a");
disp( (corr_G) );
disp("Required correlation of residual effects: e");
disp( (corr_E) );

% -------------------------------------------------------------------------
% SAMPLE GENES (QTLs) RESPONSIBLE FOR TRAITS
%correllation = [1 0.5 0.7; 0.5 1 0.2; 0.7 0.2 1];

[genes] = select_qtls(n_snp, n_qtl, 1); % Model 1

% -------------------------------------------------------------------------
% CALCULATE COVARIANCE MATRICES AND THEIR CHOLESKY DECOMPOSITION

% Does not used here. DO NOT DELETE!!!
%covar_G = sqrt( diag(var_G) ) * corr_G * sqrt( diag(var_G) );
%covar_E = sqrt( diag(var_E) ) * corr_E * sqrt( diag(var_E) );

% -------------------------------------------------------------------------
% CALCULATE CHOLESKY DECOMPOSITION OF CORRELATION MATRICES

Lg = chol(corr_G); % gives upper triangular matrix
Le = chol(corr_E); % gives upper triangular matrix

% -------------------------------------------------------------------------
% SAMPLE EFFECTS

% (1) Sample additive QTL effect
a = sample_effects(n_trait, n_qtl);

% (2) Sample GxE QTL effect
e = sample_effects(n_trait, n_qtl);

% (3) Sample dominance effects
% assume the same for all traits because a local propertiy of a locus
%k = sample_dom( n_snp, k_range_G, 3, size(genes,2) );
k = sample_dom( n_qtl, k_range_U, 1, size(genes,2) );

disp("Correlation of sampled effects before scalling: a");
disp( corr(a) );
disp("Correlation of sampled effects before scalling: e");
disp( corr(e) );
if size(k,2) > 1
    disp("Correlation of sampled effects before scalling: k");
    disp( corr(k) );
end

% -------------------------------------------------------------------------
% Adjust effects acording to correlation matrices

a = a*Lg;
e = e*Le;

if size(k,2) > 1
    k = k*Lg;
end

disp("Correlation of sampled effects after scalling: a");
disp( corr(a) );
disp("Correlation of sampled effects after scalling: e");
disp( corr(e) );
if size(k,2) > 1
    disp("Correlation of sampled effects after scalling: k");
    disp( corr(k) );
end
% -------------------------------------------------------------------------
% SIMULATE TRAITS

% (1) Calculate (un-adjusted to required variances) traits in population
[ tr_G, tr_E ] = sim_trait( n_animals, genotypes, genes, k, a, e, env, n_ploidy, haplotypes );
%
%[ ta2, te2 ] = sim_trait( n_animals, genotypes, genes, k1, a1, e1, env, n_ploidy, haplotypes );
%
disp("Correlation between traits components before variance adjustment: G");
disp( corr( tr_G ) );

disp("Correlation between traits components before variance adjustment: E");
disp( corr( tr_E ) );

figure(1); clf; subplot(2,2,1), plot(a(:,1), a(:,2), 'o'), title('a'), subplot(2,2,2), plot(e(:,1), e(:,2), 'o'), title('e'),
                subplot(2,2,3), plot(tr_G(:,1), tr_G(:,2), 'o'), title('tr_G'), subplot(2,2,4), plot(tr_E(:,1), tr_E(:,2), 'o'), title('tr_E'),

% -------------------------------------------------------------------------
% CALCULATE SCALER FOR RE-ADJUSTMENT 

% (1) Calculate scaling diagonal matrices
% square root of diag_matr of required variances * square root of inverse
% of diag_matr of current variances:
scal_a = sqrt( diag(var_G) ) * sqrt( diag( 1./(var(tr_G)) ) );
scal_e = sqrt( diag(var_E) ) * sqrt( diag( 1./(var(tr_E)) ) );

% (2) Apply re-adjustment
a = a * scal_a;
e = e * scal_e;

%disp("Correlation of sampled effects after re-scalling:");
%disp([corr(a(:,1), a(:,2)) corr(e(:,1), e(:,2))]);

% -------------------------------------------------------------------------
% Adjust QTL effects to match required variance
[ tr_G, tr_E ] = sim_trait( n_animals, genotypes, genes, k, a, e, env, n_ploidy, haplotypes );

%disp("Correlation between traits components after variance adjustment:");
%disp([corr(tr_G(:,1), tr_G(:,2)) corr(tr_E(:,1), tr_E(:,2))]);

disp("Variances of traits components after variance adjustment: G & E");
disp([var(tr_G) var(tr_E)]);

disp("Correlation between traits components after variance adjustment: G");
disp( corr( tr_G ) );

disp("Correlation between traits components after variance adjustment: E");
disp( corr( tr_E ) );

figure(2); clf; subplot(2,2,1), plot(a(:,1), a(:,2), 'o'), title('a'), subplot(2,2,2), plot(e(:,1), e(:,2), 'o'), title('e'),
                subplot(2,2,3), plot(tr_G(:,1), tr_G(:,2), 'o'), title('tr_G'), subplot(2,2,4), plot(tr_E(:,1), tr_E(:,2), 'o'), title('tr_E'),
% -------------------------------------------------------------------------
% Calculate final traits
traits = tr_G + tr_E;

traits_vars = var(traits);
traits_means = mean(traits);

% Correction term, constant value defined within the base population;
% to be applied for all further trait calculations (new, derived renotypes)
correct_const = zeros(size(corr_G,1),1);

for i = 1:size(corr_G,1)
    correct_const(i,1) = mu(i) - traits_means(1,i);
end
%correct_const(2,1) = mu(2) - traits_means(1,2);

% -------------------------------------------------------------------------

disp("Variances of integral trait (not adjusted to the required mean):");
disp(traits_vars);

% disp("Means of integral trait (not adjusted to the required mean):");
% disp(traits_means);

traits_adj = traits;
traits_adj(:,1) = traits_adj(:,1) + correct_const(1,1);
traits_adj(:,2) = traits_adj(:,2) + correct_const(2,1);

%disp("Variances of integral trait (adjusted to the required mean):");
%disp(var(traits_adj));

disp("Means of integral traits (adjusted to the required mean):");
disp(mean(traits_adj));

disp("Realised heritability:");
disp( var(tr_G)./var(traits) );

disp("Expected heritability:");
disp([ var_G(1)/(var_E(1)+var_G(1)) var_G(2)/(var_E(2)+var_G(2)) var_G(3)/(var_E(3)+var_G(3)) ]);

figure(3); clf; subplot(2,2,1), histogram(traits(:,1)), title('tr1'), subplot(2,2,2), histogram(traits(:,2)), title('tr2'),
                subplot(2,2,3), histogram(traits_adj(:,1)), title('tr1_{adj}'), subplot(2,2,4), histogram(traits_adj(:,2)), title('tr2_{adj}'),

% -------------------------------------------------------------------------
% FUNCTIONS
% -------------------------------------------------------------------------

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

% ------------------------------------------------------------------------

function [dom_coef, dom_in_locus] = eval_dom_in_locus( ianimal, iqtl, n_ploidy, haplotypes )

dom_coef = 0; % how many dominance cases in a locus; 1 < dom_coef < inf
dom_in_locus = false; % if a locus is heterozygous

alleles = 0;
ihaplotype = (ianimal-1) * n_ploidy + 1; % find a very first haplotype
                                         % for a specific animal

for l = 0:n_ploidy-1
    alleles = alleles + haplotypes(ihaplotype+l, iqtl );
    pair = mod(l+1,2); % evaluate dominance condition on pair of haplotypes
    if pair == 0
        if alleles == 1
            dom_coef = dom_coef + 1;
            dom_in_locus = true;
        end
        alleles = 0;
    end    
end
% adjust the dom_coeff to non-linear ploidy effect,
% for dipoids we assume dom_coef = 1, for 1 < degree < inf;
% for polyploids 2 (degree = 0) > dom_coef > 1 (degree -> inf)
%dom_coef = dom_coef * ploidy_effect( n_ploidy, 2 );

end

% ------------------------------------------------------------------------

function eff = sample_effects(n_traits, n_genes)
    eff = randn( n_genes,n_traits );
end

% ------------------------------------------------------------------------

function ploidy_eff = ploidy_effect( n_ploidy, degree )
% coefficient accountin a non-linear effect of ploidy,
% for diploids it equals 1, for polyploids slowly decreases,
% depending on degree parameter

    ploidy_eff = 1 / ( 1 + ( 1 - 2/n_ploidy )^degree );
end

% ------------------------------------------------------------------------

function k = sample_dom( n_genes, k_range, option_1, option_2 )   
    % n_genes - number of sampled values
    % option_1 - selects the sampling method:
    %            1 := from uniform distribution,
    %                 k_range := [lower_bound upper_bound];
    %            2 := from normal distribution,
    %                 k_range := [mean std];
    %            3 := from gamma distribution,
    %                 k_range := [k_av], which is expected value of
    %                 dominance degree in loci; is '1/b' param in the
    %                 distribution.
    % option_2 - number of columns of different effects; corresponds to
    %            the model 2.

    k = zeros(n_genes,option_2);
    
    range_sz = size(k_range,2);

    if range_sz == 1 && option_1 ~= 3
        disp("sample_dom(): Wrong paarameters!")
        return;
    end
    
    switch option_1
        case 1
            % Case (i) Uniform distribution
            k = k_range(1) + ( k_range(2) - k_range(1) ).*rand( n_genes,option_2 );
        case 2
            % Case (ii) Normal distribution
            k = normrnd(k_range(1),k_range(2),[n_genes,option_2]);
        case 3
            % Case (iii) Gamma distribution
            pdgamma = makedist('Gamma','a',2,'b',k_range(1));
            k = random(pdgamma,n_genes,option_2);
        otherwise
            disp('other value')
    end

end

% ------------------------------------------------------------------------

function [ G, E ] = sim_trait( n_animals, genotype, qtl_id, k, a, e, env, n_ploidy, haplotypes )
% n_animals - number of individuals in the genotype array
% genotype - all genotypes of selected population
% qtl_id - the list of IDs (position in the genome) of genes associated
%          with traits (note, we assuming all traits have the same amount
%          of genes);
% k - vector of sampled dominance effects;
% a - matrix of sampled additive effects;
% e - matrix of sampled residual (environmental) effects;
% env - vector of environmental changes coefficients.

n_trait = size(a,2);
n_qtl = size(qtl_id,1);

diff_qtl = size( qtl_id,2 );

G = zeros( n_animals, n_trait );
E = zeros( n_animals, n_trait );

% coefficient accountin a non-linear effect of ploidy,
% for diploids it equals 1, for polyploids slowly decreases,
% depending on degree parameter, 1 < degree < inf for diploids
% for popyploids it is allowed the degree bellow 1, 0 < degree < inf;
% 0.5 (degree == 0) < p_eff < 1 (degree -> inf)
p_eff = ploidy_effect( n_ploidy, 2 );

for l = 1:n_trait

    for inividual = 1:n_animals
        for j = 1:n_qtl
            
            if diff_qtl > 1        % if model 2, QTLs are different for different traits
                iqtl = qtl_id(j,l);
                ik = k(j,l);
            else
                iqtl = qtl_id(j,1);
                ik = k(j,1);
            end

            igene = genotype(inividual, iqtl );

            [dom_coef, dom_in_locus] = eval_dom_in_locus( inividual, iqtl, n_ploidy, haplotypes );

            if dom_in_locus % Adjust k-values for homozygosity
                G(inividual,l) = G(inividual,l) + p_eff * ( 1 + dom_coef * ik ) * a(j,l) * igene;
            else
                G(inividual,l) = G(inividual,l) + p_eff * a(j,l) * igene;
            end
            E(inividual,l) = E(inividual,l) + p_eff * ( e(j,l) + 0.5 * env(l) * abs(e(j,l)) ) * igene;
        end
    end
    % end of trait
end

% end of function
end

% ------------------------------------------------------------------------

function [genes] = select_qtls(n_snp, n_qtls, model)
% n_snp - max number of indexed SNPs (or associated genes) in a genome
% n_genes - number of genes associated with a trait;
%           note: in this version we assume all traits have the same
%                 number of associated genes;
% correlation - between-traits correlation matrix
% model - which model of trait to apply;
%         model = 1: many QTLs which are the same for all traits,
%                    correlation between traits realized by scalling;
%         model = 2: arbitrary number of QTLs, different for each
%                    trait, correlation is due to shared QTLs between
%                    traits.

if model == 1
    genes = randi( [1 n_snp], n_qtls,1 );
    return;
end

if model == 2

    disp("Does not implemented yet!");
    return;
    % end of method 2
else
    disp("Provided incorrect model parameter!");
    return;
end

% end of function
end

% ------------------------------------------------------------------------