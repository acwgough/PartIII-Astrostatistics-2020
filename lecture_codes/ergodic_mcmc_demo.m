%simulation illustrating ergodic convergence of MCMC

clear all
close all

fs = 16;

% log of the target probability density

% a standard Gaussian
log_p = @(xx) -0.5*xx.^2 -0.5*log(2*pi);

% a mixture of Gaussians
%log_p = @(xx) log(0.5*exp(-0.5*(xx+6).^2)/sqrt(2*pi) + ...
%    0.35*exp(-0.5*(xx+2).^2)/sqrt(2*pi) + 0.15*exp(-0.15*(xx-2).^2)/sqrt(2*pi)) ;

sigma_prop = 1.5;

n_mc = 11999;
n_chains = 1e4;

mc = zeros(n_mc,n_chains);
acc = zeros(n_chains,1);

xgrid = (-5:0.1:5)' * 2;

% initial values
% uniformly distributed between -5 and 5
% x0s = 10*rand(n_chains,1) - 5;

% all starting at x0 = -3
% x0s = zeros(n_chains,1) - 3;

% starting at -4 or 4
 x0s = 4*(2*binornd(1,0.5,n_chains,1)-1);

% starting at an off-centered Gaussians centered at + 3
% x0s = 7 + randn(n_chains,1)*1;

%histogram of starting values
figure(1)
histogram(x0s,'Normalization','pdf')
title(['Histogram of initial valuess'],'FontSize',fs);
hold on
plot(xgrid,exp(log_p(xgrid)),'LineWidth',3)
hold off
xlabel('x','FontSize',fs);
ylabel('pdf','FontSize',fs);
set(gca,'FontSize',fs);
%set(gca,'Xtick',[-4:2:4])

% keep track of current values
xs = x0s;
log_p_curr = log_p(xs);

%% run a simple 1D Metropolis accept/reject algorithm

disp('Begin MCMCs...');
tic;

parfor c=1:n_chains
    
    %disp(['c_mc = ' num2str(c)])
    
    for t=1:n_mc
        
        %if mod(i,1000) == 0
            %disp(['i_mc = ' num2str(i)])
        %end
        
        x_prop = xs(c) + sigma_prop*randn;
        
        log_p_prop = log_p(x_prop);
        
        logr = log_p_prop - log_p_curr(c);
        
        if log(rand) < logr
            xs(c) = x_prop;
            acc(c) = acc(c) + 1;
            log_p_curr(c) = log_p_prop;
        end
        
        mc(t,c) = xs(c);
        
    end
end

acc = acc/n_mc

runtime = toc

% include intial values in chain
mc = [x0s';mc];
n_mc = size(mc,1);

save('ergodic_mcmc.mat');

disp('MCMC done');

%% load last run
%load('ergodic_mcmc_uniform.mat');
%load('ergodic_mcmc_point3p0.mat');
%load('ergodic_mcmc_bimodal.mat');
%load('ergodic_mcmc_offc.mat');
%load('ergodic_mcmc_3gaussian.mat');

%% Calculate diagnostics

% Gelman-Rubin Ratio
gr = mcmc_calcrhat_1d(mc)

% autocorrelation function for chain #1
figure(2)
autocorr(mc(:,1))

% calculate autocorrelation time
tau = 1 + 2*sum(autocorr(mc(:,1)))

%% look at time distribution of a random chain
figure(3)
c = randi(n_chains);
plot(mc(:,c))
ylabel('x','FontSize',fs);
xlabel('MCMC sample t','FontSize',fs);
set(gca,'FontSize',fs);
title(['Trace Plot of Chain #' num2str(c)],'FontSize',fs);
xlim([0,n_mc])

figure(4)
histogram(mc(:,c),'Normalization','pdf')
hold on
plot(xgrid,exp(log_p(xgrid)),'LineWidth',3)
hold off
xlabel('x','FontSize',fs);
set(gca,'FontSize',fs);
title(['Time-Histogram Chain #' num2str(c)],'FontSize',fs);
%set(gca,'Xtick',[-4:2:4])

%% Look at trace paths for a handful of chains for 100 steps
figure(5)
for i=1:10
    c = randi(n_chains);
    plot(mc(1:50,c))
    hold on
    ylabel('x','FontSize',fs);
    xlabel('MCMC sample t','FontSize',fs);
    set(gca,'FontSize',fs);
    title(['Trace Plot of 10 Chains'],'FontSize',fs);
end
hold off

%% Look at initial and final chain-ensemble distributions
figure(6)
subplot(1,2,1)
histogram(mc(1,:),'Normalization','pdf');
xlabel('x','FontSize',fs);
set(gca,'FontSize',fs);
xlim([-10,10])
title(['Ensemble at step t = ' num2str(1)],'FontSize',fs);
hold on
plot(xgrid,exp(log_p(xgrid)),'LineWidth',3)
hold off
ylim([0,0.45]);
%set(gca,'Xtick',[-4:2:4])

subplot(1,2,2)
histogram(mc(100,:),'Normalization','pdf');
xlabel('x','FontSize',fs);
set(gca,'FontSize',fs);
xlim([-10,10])
title(['Ensemble at step t = ' num2str(100)],'FontSize',fs);
hold on
plot(xgrid,exp(log_p(xgrid)),'LineWidth',3)
hold off
ylim([0,0.45]);
%set(gca,'Xtick',[-4:2:4])


%% Look at chain-ensemble distribution over time until t_max

t_max = 50;

figure(7)
pause(3)
for t=1:t_max
    histogram(mc(t,:),'Normalization','pdf')
    hold on
    plot(xgrid,exp(log_p(xgrid)),'LineWidth',3)
    hold off
    xlabel('x','FontSize',fs);
    set(gca,'FontSize',fs);
    xlim([-10,10])
    %set(gca,'Xtick',[-4:2:4])
    ylim([0,0.45])
    title([num2str(n_chains) '-Chain Ensemble Histogram at step t = ' num2str(t)],'FontSize',fs);
    F(t) = getframe;
    pause(0.15)
end

        