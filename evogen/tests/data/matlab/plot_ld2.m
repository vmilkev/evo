clear;

file_ld1 = ["ld_at_0_chr_0.ld", "ld_at_10_chr_0.ld", "ld_at_20_chr_0.ld", "ld_at_30_chr_0.ld", "ld_at_40_chr_0.ld", ...
           "ld_at_50_chr_0.ld", "ld_at_60_chr_0.ld", "ld_at_70_chr_0.ld", "ld_at_80_chr_0.ld", "ld_at_90_chr_0.ld", ...
           "ld_at_99_chr_0.ld"];

file_ld2 = ["ld_at_0_chr_0.ld", "ld_at_100_chr_0.ld", "ld_at_200_chr_0.ld", "ld_at_300_chr_0.ld", "ld_at_400_chr_0.ld", ...
           "ld_at_500_chr_0.ld", "ld_at_600_chr_0.ld", "ld_at_700_chr_0.ld", "ld_at_800_chr_0.ld", "ld_at_900_chr_0.ld", ...
           "ld_at_999_chr_0.ld"];

figure(1);
clf;
set(gcf,'Color','white');
hold on;

for j = 1:size(file_ld2,2)
    ld = import_ld(file_ld1(j));
    plot(ld(:,1), ld(:,2), '*-');
end

function ld = import_ld(filename)
ld = readtable(filename, "FileType","text");
ld = table2array(ld);
end
