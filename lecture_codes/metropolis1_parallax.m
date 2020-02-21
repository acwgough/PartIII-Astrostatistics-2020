% basic Metropolis algorithm

clear all
close all

fs = 6

omega = 0.01;
ferr = 0.33;
somega = ferr*omega;
Lgal = 1000;

logposterior = @(r) logpdf_parallax(r,omega,somega,Lgal);

rgrid = 1:10000;

figure(1)
hpost = plot(rgrid,exp(logposterior(rgrid)-max(logposterior(rgrid))),'LineWidth',2);
hold on;
xlim([1,10000])
ylabel('Unnormalised Density','FontSize',fs);
xlabel('Distance r (parsec)','FontSize',fs);
title('Parallax Posterior for Distance','FontSize',fs)

text(5000,0.6,['\omega = ' num2str(omega)],'Color','b','FontSize',fs)
text(5000,0.5,['\sigma_\omega = ' num2str(ferr) ' \times \omega'],'Color','b','FontSize',fs)
text(5000,0.4,['L_{gal} = ' num2str(Lgal)],'Color','b','FontSize',fs)
set(gca,'XScale','Linear')
set(gca,'FontSize',fs)

%%

n_mc = 1e5;
tau = 1000;

rs_mc = zeros(n_mc,1);

rs_mc(1) = 2000;

acceptances = 0;

for i=1:n_mc
    
    r_curr = rs_mc(i);
    
    r_prop = r_curr + randn*tau;
    
    logpost_curr = logposterior(r_curr);
    
    logpost_prop = logposterior(r_prop);
    
    log_mhr = logpost_prop - logpost_curr;
    
    % pick a uniform random variate between 0 and 1
    u = rand;
    
    if log(rand) < log_mhr
        rs_mc(i+1) = r_prop;
        acceptances = acceptances + 1;
    else
        rs_mc(i+1) = r_curr;
    end
    
end

acc_ratio = acceptances/n_mc

%%
%rs_mc = rs_mc(n_mc/5 : 2:  end);

figure(1)

subplot(2,1,1)
plot(rs_mc,'LineWidth',0.5)
ylabel('Parameter Value \mu','FontSize',fs)
%xlabel('Chain step','FontSize',20)
title(['Markov Chain Trace Plot: acc ratio = ' num2str(acc_ratio,'%.2f')],'FontSize',fs)
set(gca,'FontSize',fs)
xlim([1,n_mc])

subplot(2,1,2)
plot(logposterior(rs_mc),'LineWidth',2)
ylabel('Log Post Log P(\mu | y)','FontSize',fs)
xlabel('Chain step','FontSize',fs)
set(gca,'FontSize',fs)
xlim([1,n_mc])

%% plot histogram

%rs_mc = rs_mc(1 : end);
%rs_mc = rs_mc(n_mc/2 : 2:  end);

figure(2)

h=histogram(rs_mc,'Normalization','pdf');
hold on
%h2=plot(x, normpdf(x,0,1/sqrt(N)),'LineWidth',2);
h2=plot(rgrid, 3.5e-4*exp(logposterior(rgrid)-max(logposterior(rgrid))),'LineWidth',2)

xlabel('distance r','FontSize',fs)
ylabel('Posterior P(r | \omega)','FontSize',fs)
legend({'Posterior Histogram','Analytic Posterior'})
title(['N_{mc} = ' num2str(n_mc,'%.0f')],'FontSize',fs)
set(gca,'FontSize',fs)


