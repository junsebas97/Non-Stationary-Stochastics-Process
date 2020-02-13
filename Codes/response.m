function [d, v, a] = response(as, m, w, xi, dt)
%{
This function calcules the response of a structure to a given seismic record
   
Inputs:
    as: a column vector containing the seismic record
    m:  mass of the structure
    w:  frequency of the strucutre
    xi: equivalent dumping
    dt: time between records
    
Outputs:
    d: structure's displacements
    v: structure's velocities
    a. structure's acceleration
    
NOTE: this code is based on dmaclin of INTRODUCCIÓN A LA DINÁMICA DE
      ESTRUCTURAS - Jorge Eduardo Hurtado
%}
    
    p    = -m*as;        % the force due to the quake
    n    = length(p);    % number of seismic evaluation
    
    % it's assumed that the structure is in repose 
    d  = zeros([n, 1]);    d0 = 0;    
    v  = zeros([n, 1]);    v0 = 0;                
    a  = zeros([n, 1]);    a0 = 0;
    
    k = m*w^2;                         % structure's stiffness
    c = 2*m*w*xi;                      % structure's damping
    kbar = k + 3*c/dt + 6*m/(dt^2);    % equivalent stiffness
    
    % in each instant:
    for i = 1:n
        
        % 1) the effective acting force is assessed
        p1 = p(i, :);
        dp = m*(6*d0/(dt^2) + 6*v0/dt + 2*a0) + c*(3*d0/dt + 2*v0 + dt*a0/2);
        pbar = p1 + dp;
        
        % 2) the produced displacemet, velocity and acceleration are
        %    calculated
        d1 = pbar/kbar;
        v1 = 3*(d1 - d0)/dt - 2*v0 - dt*a0/2;
        a1 = 6*(d1 - d0)/dt^2 - 6*v0/dt - 2*a0;
        
        % 3) those parameters are stored and will be used in the next
        %    iteration
        d(i, 1) = d1;        d0 = d1; 
        v(i, 1) = v1;        v0 = v1;
        a(i, 1) = a1;        a0 = a1;      
    end
end