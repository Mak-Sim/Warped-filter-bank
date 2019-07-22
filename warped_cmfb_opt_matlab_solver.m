function h_opt = warped_cmfb_opt_matlab_solver(h0,alpha,Sk,mode,max_itr,Npt, psi_term, gamma)
%% Optimization function of warped cosine modulated filter bank.  A short algorithm description
% is given in Vashkevich, M. Wan, W. and Petrovsky, A., "Practical design of multi-channel oversampled 
% warped cosine-modulated filter banks" Proc. of IET International Communication Conference on 
% Wireless Mobile and Computing (CCWMC-2011), pp. 44-49, Shanghai, China,
% 14-16 Nov., 2011. 
% See also: http://arxiv.org/abs/1111.0737

% Parameters:
% h0    -- initial filter-prototype
% alpha -- warping coefficient
% Sk    -- subsampling ratios
% mode  -- Should be set to 'new' when the fucntion called first time with
%          particular values of alpha, Sk, N, and Npt. Otherwise 'old'. This allow
%          to save time.
% max_itr -- Maximum number of iteration
% Npt     -- Number of point in frequency grid

%% 1. Filter bank parameters
N = length(h0); % Filter-prototype order
N2 = N/2;       

% 2. Design algorithm parameters 
w  = linspace(0, pi, Npt);      % Frequency grid
h_vec = h0(N/2+1:end)';         % Vector of filter prototype coefficients
B = ones(1,Npt)/sqrt(Npt);      % Weighting function
gamma = gamma*ones(1,Npt);      % Ripple ratio function

R1 = zeros(N2,N2);
I1 = zeros(N2,N2);
U_matr = zeros(length(w),N2,N2);
U_alias = zeros(length(w),N2,N2);
R1_all = zeros(length(w),N2,N2);
R2_all = zeros(length(w),N2,N2);
I1_all = zeros(length(w),N2,N2);
I2_all = zeros(length(w),N2,N2);
betta = zeros(1,length(w));
teta = 1.05;                    % factor that affect the convergence
psi=1;
Err = ones(1,length(w));
E_min = 10e10;

% 2.1 Pre-calculation
switch (mode)
    case 'new'
        disp('Precalculation started...');
        for i=1:length(w)
            U = U_mat(w(i),Sk,N,alpha,'T_alias');            
            U1 = U_mat(w(i),Sk,N,alpha,'T_dist');
            U_matr(i,:,:) = U_mat(w(i),Sk,N,alpha,'T_all');
            U_alias(i,:,:) = U;
            R1 = gamma(i)*real(U) + real(U1);   
            I1 = gamma(i)*imag(U) + imag(U1);
            R2 = R1+R1';    I2 = I1+I1';
        
            R1_all(i,:,:) = R1;     R2_all(i,:,:) = R2;
            I1_all(i,:,:) = I1;     I2_all(i,:,:) = I2;
        end;
        save R1_all_cur R1_all;
        save R2_all_cur R2_all;
        save I1_all_cur I1_all;
        save I2_all_cur I2_all;
        save U_matr_cur U_matr;
        save U_alias_cur U_alias;
        disp('Precalculation done...');
    case 'old'
        load R1_all_cur;
        load R2_all_cur;
        load I1_all_cur;
        load I2_all_cur;
        load U_matr_cur;
        load U_alias_cur;
        disp('Precalculated data is loaded...');
    otherwise
        disp('Unknown mode!');
        return;
end

% 3. Main optimization cycle
options = optimset('GradObj','on',...
                   'Hessian','on',...
                   'Display','iter-detailed',...
                   'MaxIter',100,...
                   'TolX', 10^-10,...
                   'TypicalX',h_vec+0.001*randn(size(h_vec)));
              
itr =1;
while((psi>psi_term) && itr<max_itr)
            
    cost_fun_val = wcmfb_cost_function(h_vec, w, B, N2);    
    fprintf(1,'Cost function before value:  %012.13f\n', cost_fun_val);    
    cost_fun = @(h_vec)wcmfb_cost_function(h_vec, w, B, N2);
    [h_vec,cost_fun_val] = fminunc(cost_fun, h_vec, options);
    
    fprintf(1,'Cost function after value:  %015.16f\n', cost_fun_val);
    
    
    % Calculation of error function
    for i=1:length(w)
        R1(:,:) = R1_all(i,:,:);
        I1(:,:) = I1_all(i,:,:);
        Err(i) = (h_vec'*R1*h_vec).^2 + (h_vec'*I1*h_vec).^2 - 1;
    end
    
    fprintf(1,'***************** END OF ITERATION #%d*****************\n', itr);
    itr = itr + 1;
    h_i = [flipud(h_vec)' h_vec'];    
    E = abs(Err);
    E_cur = max(E);
    if E_cur<E_min,
        h_opt = h_i;
        E_min = E_cur;
    end
    
    % Now we need to evoluate envelope function V(w)
    [V, nV] = getMax(E,w);
    I = 1;
    for i=1:length(w),
        if w(i)> w(nV(I+1)),
            I=I+1;
        end
        betta(i) = ((w(i)-w(nV(I)))*V(I+1)+(w(nV(I+1))-w(i))*V(I)) / (w(nV(I+1)) - w(nV(I)));
    end
    
    %  Optional information
    if (mod(itr,4)==1)
        figure('Units','pixels', 'Position',[200 200 1000 300]);
        axes('Position', [0.1 0.175 0.85 0.725])
        hold on;
        hb = plot(w/(2*pi), betta,'-.');
        hE = plot(w/(2*pi),E); 
        hTitle = title('Error function / Envelope function');
        ylabel('$|E_{\mu,opt}(\omega)|,\,\beta_{\mu}(\omega)$', 'Interpreter', 'Latex','FontSize',14);   % 
        xlabel('$\omega / 2\pi$', 'Interpreter', 'Latex','FontSize',14);
        ...xlim();
        set (gca                    ,...
            'FontName', 'Helvetica');
        
        set(hTitle                  ,...
            'FontSize'  , 14        ,...
            'FontName', 'AvantGarde');
        
        set(hb                      ,...
            'LineWidth', 1.5        ,...
            'LineStyle','-.'        ,...
            'Color', [5 5 5]/255);        
        
        set(hE,'LineWidth',1.5);        
                
        set(gca, ...
            'Box'       ,'off'          ,...
            'TickDir'   ,'in'           ,...
            'TickLength',[0.02,0.02]    ,...
            'XMinorTick','on'           ,...
            'YMinorTick','on'           ,...
            'YGrid'     ,'on'           ,...
            'XGrid'     ,'on'           ,...
            'XColor'    ,[0.1 0.1 0.1]  ,...
            'YColor'    ,[0.1 0.1 0.1]  ,...
            'XTick'     ,0:0.1:0.5      ,...
            ...'YTick'     ,[-0.005:0.005:0.016],...
            ...'YLim'      ,[-0.000, 0.015],...
            'FontSize', 12              ,...
            'LineWidth' ,1.5              );      
    end
    
    % Update weighted function
    min_b = min(betta);
    max_b = max(betta);
    betta_t = betta.^teta;
    B = B.*betta_t;
    B = B./norm(B);
    psi = (max_b-min_b)/(max_b+min_b);
    fprintf(1,'PSI = %03.4f \n', psi);
end

