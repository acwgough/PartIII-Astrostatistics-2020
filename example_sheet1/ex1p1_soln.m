clear all
close all

% read Calibrator sample data
fid = fopen('dhawan_table1_abridged.txt');
table = textscan(fid,'%s %s %f %f %f %f','CommentStyle','#');
fclose(fid);

mhatk = table{3};
sigma_mk = table{4};
muhatCk = table{5};
sigma_Ck = table{6};

K = length(mhatk);

% read Hubble Flow sample data
fid = fopen('dhawan_table2_abridged.txt');
table = textscan(fid,'%s %s %f %f %f','CommentStyle','#');

zi = table{3};
mhati = table{4};
sigma_mi = table{5};

%% evaluate MLE estimators

sigma_int = 0.11;

tauk = sqrt(sigma_int.^2 + sigma_mk.^2 + sigma_Ck.^2);
sigmai = sqrt(sigma_int.^2 + sigma_mi.^2);

Mk_hat = mhatk - muhatCk;
Mi_tilde = mhati - 25 - 5*log10(3e5 * zi/100);

M0_mle = sum( tauk.^-2 .* Mk_hat)/sum(tauk.^-2)
var_M0 = 1/sum(tauk.^-2)
std_M0 = sqrt(var_M0)

theta_mle = M0_mle - sum(sigmai.^-2 .* Mi_tilde)/sum(sigmai.^-2)
var_theta = 1/sum(tauk.^-2) + 1/sum(sigmai.^-2)
std_theta = sqrt(var_theta)

h_mle = 10^(theta_mle/5)

frac_std_h_mle = 0.46 * sqrt(var_theta)

%% optimise likelihood over M0, theta, and sigma_int;

negloglkhd = @(pars) 0.5 * sum( (Mk_hat - pars(1)).^2 ./(pars(3)^2 + sigma_mk.^2 + sigma_Ck.^2) ) ...
    + 0.5 * sum(log(2*pi*(pars(3)^2 + sigma_mk.^2 + sigma_Ck.^2))) ...
    + 0.5 * sum( (Mi_tilde - pars(1) + pars(2) ).^2 ./(pars(3)^2 + sigma_mi.^2)) ...
    + 0.5 * sum(log(2*pi*(pars(3)^2 + sigma_mi.^2)));

start = [M0_mle; theta_mle; 0.20]
    
[out,fval,exitflag,output,lambda,grad,hessian] = fmincon(negloglkhd,start,[],[],[],[],[],[],[])

% observed Fisher matrix is the inverse of the Hessian of the negative log
% likelihood

mle = out

sigmas_mle = sqrt(diag(inv(hessian)))

