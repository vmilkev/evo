rng(0,'twister');

% no. chromosomes = 4; 
chr = zeros(4,1);
chr(1) = 10000; % 10K bp
chr(2) = 20000; % 20K bp
chr(3) = 10000; % 10K bp
chr(4) = 15000; % 15K bp

% animals = 5; chromosomes = 4; snp distance = 1K bp; file name
%make_haplotypes(5, chr, 1000, true, 'haplotypes_pop1.dat'); % with pedigree and sex
%[hpl, snps] = make_haplotypes(5, chr.*1, 1000, false, 'haplotypes_pop2.dat'); % only haplotypes
[hpl, snps] = make_haplotypes(1000, chr.*5000, 1000, false, 'haplotypes_pop_large2.dat'); % only haplotypes

disp([hpl snps]);

% -------------------------------------------------------------------
function [ haplotypes, snp_variants ] = make_haplotypes(individuals, chromosome_length, snp_density, with_pedigree, fname)
    snp_variants = sum(chromosome_length/snp_density);
    H = zeros(individuals*2,snp_variants);
    haplotypes = individuals*2;
    if with_pedigree
        P = zeros(individuals*2,3); % cols: ind_id, sex, sire/dam 
    end
    for i = 1:2:individuals*2
        r1 = randi([0 1],1,snp_variants);
        r2 = randi([0 1],1,snp_variants);
        H(i,:) = r1; % father
        H(i+1,:) = r2; % mather
        if with_pedigree
            if randi([1 100],1,1) <= 50
                sex = 1;
            else
                sex = 0;
            end
            id1 = randi([100 1000],1,1);
            id2 = randi([100 1000],1,1);
            id3 = randi([100 1000],1,1);
            P(i,1) = id1; % individual
            P(i,3) = id2; % father
            P(i,2) = sex;
            P(i+1,1) = id1; % individual
            P(i+1,3) = id3; % mather
            P(i+1,2) = sex;
        end
    end
    if with_pedigree
        T = table(P,H);
        writetable(T,fname, 'Delimiter', ' ', 'WriteVariableNames',false);
    else
        writematrix(H,fname, 'Delimiter', ' ');
    end
    name1 = 'struct_';
    fname2 = strcat(name1,fname);
    genotype_structure(chromosome_length, snp_density, fname2);
end

function genotype_structure(chromosome_length, snp_density,fname)
    chromosomes = numel(chromosome_length);
    s = zeros(chromosomes,2);
    for i = 1:chromosomes
        s(i,1) = chromosome_length(i);
        s(i,2) = snp_density;
    end
    writematrix(s,fname, 'Delimiter', ' ');
end
