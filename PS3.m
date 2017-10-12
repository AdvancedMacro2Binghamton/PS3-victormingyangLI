% PROGRAM NAME: ps3huggett
clear, clc

% PARAMETERS
beta = .9932; % discount factor 
sigma = 1.5; % coefficient of risk aversion
b = 0.5; % replacement ratio (unemployment benefits)
y_s = [1, b]; % endowment in employment states
PI = [.97 .03; .5 .5]; % transition matrix

% ASSET VECTOR
% Why do we have the lower and upper bounds of asset?
a_lo = -2; % lower bound of grid points
a_hi = 5; % upper bound of grid points
% How to decide the size of asset vector?
num_a = 1000;

a = linspace(a_lo, a_hi, num_a); % asset (row) vector

% INITIAL GUESS FOR q
% Why is this the range of q?
q_min = 0.98;
q_max = 1;
q_guess = (q_min + q_max) / 2;

% ITERATE OVER ASSET PRICES
% What is the variable aggsav? Excess demand for asset? Yes.
aggsav = 1 ;
% Loop counter
% i = 0;
while abs(aggsav) >= 0.01 ;
    q_guess = (q_min + q_max) / 2;
    % CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus, a', q_guess * a); % a' is today's asset while a is tomorrow's
    cons = bsxfun(@plus, cons, permute(y_s, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret(cons < 0) = -Inf;
    
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(2, num_a);
    
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    while v_tol > 0.0001;
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        value_mat_e = ret(:,:,1) + beta * ...
                     (PI(1,1) * repmat(v_guess(1,:), [num_a 1]) + ...
                      PI(1,2) * repmat(v_guess(2,:), [num_a 1]));
        value_mat_u = ret(:,:,2) + beta * ...
                     (PI(2,1) * repmat(v_guess(1,:), [num_a 1]) + ...
                      PI(2,2) * repmat(v_guess(2,:), [num_a 1]));
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        [vfn_e, pol_indx_e] = max(value_mat_e, [], 2);
        [vfn_u, pol_indx_u] = max(value_mat_u, [], 2);
        vfn = [vfn_e, vfn_u];
        vfn = permute(vfn, [2 1 3]);
        
        % what is the distance between current guess and value function
        v_tol = max(max(abs(vfn - v_guess)));
        
        % if distance is larger than tolerance, update current guess and
        % continue, otherwise exit the loop
        v_guess = vfn;
    end;
    
    % KEEP DECSISION RULE
    pol_fn_e = a(pol_indx_e);
    pol_fn_u = a(pol_indx_u);
    pol_fn = [pol_fn_e; pol_fn_u];
    pol_indx = [pol_indx_e pol_indx_u]';
    
    % SET UP INITITAL DISTRIBUTION
    Mu = zeros(2, num_a) + 1 / (2 * num_a);
    dis = 1;
    % ITERATE OVER DISTRIBUTIONS
    while dis > 0.0001
        [emp_ind, a_ind, mass] = find(Mu > 0); % find non-zero indices
        MuNew = zeros(size(Mu));
        for ii = 1:length(emp_ind)
            apr_ind = pol_indx(emp_ind(ii), a_ind(ii)); % which a prime does the policy fn prescribe?
            MuNew(:, apr_ind) = MuNew(:, apr_ind) + ... % which mass of households goes to which exogenous state?
            [PI(emp_ind(ii), 1) * Mu(emp_ind(ii), a_ind(ii)); 
             PI(emp_ind(ii), 2) * Mu(emp_ind(ii), a_ind(ii))];
        end
        dis = max(max(abs(MuNew - Mu)));
        Mu = MuNew;
    end
    % i = i + 1
    aggsav = Mu(1, :) * a' + Mu(2, :) * a';
    if aggsav >= 0;
        q_min = q_guess;
    else
        q_max = q_guess;
    end
end

% Compute the risk-free interest rate
r = 1 / q_guess;

% Plot the value function of the Employed and Unemployed
plot(a, v_guess);
legend('Employed', 'Unemployed', 'location', 'northwest');
title('Value Function of the Employed and Unemployed');

% Plot the Lorenz curves and compute the Gini coefficients
Mu_trans = Mu';
pop = [Mu_trans(:, 1); Mu_trans(:, 2)];
val_earnings = [repmat(y_s(1), [num_a, 1]); repmat(y_s(2), [num_a, 1])];
figure
[gini_earnings, l_earnings] = gini(pop, val_earnings, true);
title(['Lorenz Curve of Earnings with Gini Coefficient = ', num2str(gini_earnings)])
val_wealth = [a' + y_s(1); a' + y_s(2)];
val_wealth(val_wealth < 0) = 0;
figure
[gini_wealth, l_wealth] = gini(pop, val_wealth, true);
title(['Lorenz Curve of Wealth with Gini Coefficient = ', num2str(gini_wealth)])