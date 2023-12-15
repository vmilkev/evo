clear;

x = 0:0.01:5;
a = 1; % expected number of events/crossovers; shape parameter
%mu = 2/max(x); % average intensity, events per chromosome
mu = 2;

y = gampdf(x, a, 1/mu);

pd = makedist('Gamma','a',a,'b',1/mu);
pd2 = makedist('Binomial','N',2,'p',0.09);

rng('default')  % For reproducibility
r = random(pd2,10000,1);

clf;
figure(1);
%plot(x./max(x),y, 'o-');
%plot(x,y, 'o-');
histogram(r)
hold on;
%histogram(r,100,'Normalization','pdf');
hold off;

%% simulation coincedence
clear

intervals = 20;
sim_rounds = 10;
chromosome_u = zeros(intervals,sim_rounds);
chromosome_g = zeros(intervals,sim_rounds);
cross_av = 3;

% Sampling
for i = 1:sim_rounds
    i2 = randi([1 intervals],1,1);
    chromosome_g(i2,i) = 1;
    dir_right = true;
    for j = 1:cross_av
        %pdgamma = makedist('Gamma','a',j+1,'b',1/(cross_av/intervals));
        pdgamma = makedist('Gamma','a',2,'b',1/(cross_av/intervals));
        r2 = round( random(pdgamma,1,1) );
        if ( i2 + r2 ) <= intervals && dir_right
            i2 = i2 + r2;
            chromosome_g(i2,i) = 1;
            dir_right = true;
        elseif ( i2 - r2 ) >= 1
            i2 = i2 - r2;
            chromosome_g(i2,i) = 1;
            dir_right = false;
        end
        i1 = randi([1 intervals],1,1); % interval of crossover event
        chromosome_u(i1,i) = 1;
    end
end
%%
figure(2);
plot(chromosome_g(:,10), '*-');
%%
% freq for interval L = 1
for i = 1:intervals

    f1_u = sum( chromosome_u(i,:) );
    f2_u = 0;

    f1_g = sum( chromosome_g(i,:) )
    f2_g = 0;

    for j = 1:sim_rounds
        if i == intervals
            break
        end
        s_u = chromosome_u(i,j) + chromosome_u(i+1,j);
        s_g = chromosome_g(i,j) + chromosome_g(i+1,j);

        if s_u == 2
            f2_u = f2_u + 1;
        end

        if s_g == 2
            f2_g = f2_g + 1;
        end
    end

    disp([i f1_g f2_g]);

    %f2_u = f2_u/sim_rounds;
    f_u(i,1) = (f1_u)/sim_rounds;
    f_u(i,2) = f2_u/sim_rounds;

    %f2_g = f2_g/sim_rounds;
    f_g(i,1) = (f1_g)/sim_rounds;
    f_g(i,2) = f2_g/sim_rounds;

end
for i = 1:intervals-1
    f_u(i,3) = f_u(i,2)/( ( f_u(i,1) + f_u(i,2) )*( f_u(i+1,1) + f_u(i,2) ) );
    f_g(i,3) = f_g(i,2)/( ( f_g(i,1) + f_g(i,2) )*( f_g(i+1,1) + f_g(i,2) ) );
end


% for j = 0:floor(intervals/2)
%     for i = 1:intervals
%         f1 = chromosome(i:i+j,1)./50;
%         c(j+1,1) = chromosome(i,1)
%     end
% end