function [ wMinMin, q_min, p_min ] = minimize_wigner( rho, q_0, p_0 )
%MINIMIZE_WIGNER finds the minimum of the wigner function
%   
%   W_MIN = MINIMIZE_WIGNER( RHO, Q_0, P_0) calculates the minimum of the
%   Wigner function of density matrix RHO.  It takes as its initialization
%   point (Q_0, P_0) in phase space.
%
%   wMinMin = minimize_wigner(rho) tries initialization points at the
%   origin and at the corners of the unit square.
%
%   [W_MIN, Q_MIN, P_MIN] = ... also returns the point Q_MIN and P_MIN
%   where the minimum was found.
%
%   Warning, the minimization will get stuck in local minima, so the point
%   (Q_0, P_0) must be chose carefully.


options = optimset('Display','off','LargeScale','off','HessUpdate','dfp');

minimize_me = @(q_and_p) wigner(rho, q_and_p(1), q_and_p(2));

if exist('q_0','var') && exist('p_0','var')
    [q_and_p_min, wMinMin(1)] = fminsearch(minimize_me,[q_0,p_0],options);
    iBest=1;
else
    [q_and_p_min(1,:), w_min(1)] = fminsearch(minimize_me,[0,0],options);
    [q_and_p_min(2,:), w_min(2)] = fminsearch(minimize_me,[1,1],options);
    [q_and_p_min(3,:), w_min(3)] = fminsearch(minimize_me,[-1,1],options);
    [q_and_p_min(4,:), w_min(4)] = fminsearch(minimize_me,[1,-1],options);
    [q_and_p_min(5,:), w_min(5)] = fminsearch(minimize_me,[-1,-1],options);
    [wMinMin, iBest] = min(w_min);
end

if nargout >= 2
    q_min = q_and_p_min(iBest,1);
    p_min = q_and_p_min(iBest,2);
end

end

