function l = loglkhd_timedelay(y1s,y1errs,y2s,y2errs,tobs,dt,dm,c,A,tau)

if (A <= 0 || tau <= 0)
    l = -1e9;
    return
else
    
    
    Ys = [y1s; y2s - dm];
    ts = [tobs; tobs - dt];
    one = Ys*0 + 1;
    
    W = diag([y1errs; y2errs].^2);
    
    T_grid_mat = ts(:,ones(1,length(ts)));
    
    totalcov = A^2 *exp(-abs(T_grid_mat-T_grid_mat') /tau) + W;
    
    %detcov = det(totalcov)
    %invcov = inv(totalcov);
    cholL = chol(totalcov,'lower');
    logdetcov = 2 * sum(log(diag(cholL)));
    
    resid = (Ys - one*c);
    
    % solve cholL*z = r for z using forward substitution
    z = cholL\resid;
    
    %log likelihood
    l = -0.5*logdetcov - 0.5*(z'*z);
    
    % equivalent to:
    %l = -0.5*logdetcov - 0.5*resid' * (cholL' \(cholL\(resid)));
    %l = -0.5*logdetcov - 0.5*resid'*inv(totalcov)*resid;
end


end