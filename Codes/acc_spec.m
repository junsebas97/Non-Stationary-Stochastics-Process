function [T, Sa] = acc_spec(as, Tmin, Tmax, dt)
%{
This functions produce an acceleration spectre, with equivalent damping of
5%, in a period interval for a given accelaration time series
time series 
    
Inputs:
    as:   a column vector containing the accelaration time series
    Tmin: bottom period considered
    Tmax: top period considered
    dt:   delta time
    
Outputs:
    T:  an array containig the periods
    Sa: an array containing the spectre evaluations
    
NOTE: this code is based on desplin of INTRODUCCIÓN A LA DINÁMICA DE
      ESTRUCTURAS - Jorge Eduardo Hurtado
%}
    
    xi = 0.05;    % equivalent damping
    
    % the period and frequency domain are created
    m  = 200;     % number of evaluations 
    T  = linspace(Tmin, Tmax, m);
    W  = 2*pi./T;
    
    Sa = nan([m, 1]);
    
    % for each period the structure's response is calculated, but on ly the
    % maximum acceleration is considered and returned
    for i = 1:m
        w = W(i);
        [~, ~, acc] = response(as, 1, w, xi, dt);
        Sa(i, 1)    = max(abs(as + acc));
    end
    
end