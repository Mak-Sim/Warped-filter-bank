function U = U_mat(w,Sk,N,alpha,opt)
%% Function returns matrix U(w) for calculation of: 
% 1) T_all  (opt='T_all')
% 2) T_dist (opt='T_dist')
% 3) T_alias(opt='T_alias')
% Parameters:
% w  -- frequency value
% Sk -- subsampling ratios
% N  -- filter-prototye order
% alpha -- warping coefficient

% Filter bank parameters
M = length(Sk);     % Number of channel
U = zeros(N/2,N/2); 
Hk = zeros(1,N/2);  % k-th analysis filter
Fk = zeros(1,N/2);  % k-th synthesis filter
switch opt
    case 'T_all'
        for k=0:M-1
            ak = exp(1j*(pi/4)*(-1)^k);
            bk = exp(-1j*2*pi*((N-1)/2)*(k+0.5)/(2*M));
            gamma_1_0_k = -freq_warp(w,alpha) - 2*pi*(k+0.5)/(2*M);
            gamma_2_0_k = -freq_warp(w,alpha) + 2*pi*(k+0.5)/(2*M);
            Fk = (conj(ak) * bk * exp(-1j*(N-1)*gamma_1_0_k/2) * C_matr(gamma_1_0_k,N)' + ...
                  ak * conj(bk) * exp(-1j*(N-1)*gamma_2_0_k/2) * C_matr(gamma_2_0_k,N)');
            for l=0:Sk(k+1)-1
                gamma_1_l_k = -freq_warp(w + 2*pi*l/Sk(k+1),alpha) - 2*pi*(k+0.5)/(2*M);
                gamma_2_l_k = -freq_warp(w + 2*pi*l/Sk(k+1),alpha) + 2*pi*(k+0.5)/(2*M);
                Hk = Hk + (ak * bk * exp(-1j*(N-1)*gamma_1_l_k/2) * C_matr(gamma_1_l_k,N)' + ...
                           conj(ak) * conj(bk) * exp(-1j*(N-1)*gamma_2_l_k/2) * C_matr(gamma_2_l_k,N)');                
            end
            U = U + Fk.'*Hk;
            Hk = zeros(1,N/2);   
        end
        
    case 'T_dist'
        for k=0:M-1
            ak = exp(1j*(pi/4)*(-1)^k);
            bk = exp(-1j*2*pi*((N-1)/2)*(k+0.5)/(2*M));
            gamma_1_0_k = -freq_warp(w,alpha) - 2*pi*(k+0.5)/(2*M);
            gamma_2_0_k = -freq_warp(w,alpha) + 2*pi*(k+0.5)/(2*M);
            
            Fk = (conj(ak) * bk * exp(-1j*(N-1)*gamma_1_0_k/2) * C_matr(gamma_1_0_k,N)' + ...
                  ak * conj(bk) * exp(-1j*(N-1)*gamma_2_0_k/2) * C_matr(gamma_2_0_k,N)');
              
            U = U + Fk.'*(ak * bk * exp(-1j*(N-1)*gamma_1_0_k/2) * C_matr(gamma_1_0_k,N)' + ...
                           conj(ak) * conj(bk) * exp(-1j*(N-1)*gamma_2_0_k/2) * C_matr(gamma_2_0_k,N)'); 
        end
        
    case 'T_alias'
        for k=0:M-1
            ak = exp(1j*(pi/4)*(-1)^k);
            bk = exp(-1j*2*pi*((N-1)/2)*(k+0.5)/(2*M));
            gamma_1_0_k = -freq_warp(w,alpha) - pi*(2*k+1)/(2*M);
            gamma_2_0_k =- freq_warp(w,alpha) + pi*(2*k+1)/(2*M);
            Fk = (conj(ak) * bk * exp(-1j*(N-1)*gamma_1_0_k/2) * C_matr(gamma_1_0_k,N)' + ...
                  ak * conj(bk) * exp(-1j*(N-1)*gamma_2_0_k/2) * C_matr(gamma_2_0_k,N)');
            for l=1:Sk(k+1)-1
                gamma_1_l_k = -freq_warp(w + 2*pi*l/Sk(k+1),alpha) - pi*(2*k+1)/(2*M);
                gamma_2_l_k = -freq_warp(w + 2*pi*l/Sk(k+1),alpha) + pi*(2*k+1)/(2*M);
                Hk = Hk + (ak * bk * exp(-1j*(N-1)*gamma_1_l_k/2) * C_matr(gamma_1_l_k,N)' + ...
                           conj(ak) * conj(bk) * exp(-1j*(N-1)*gamma_2_l_k/2) * C_matr(gamma_2_l_k,N)');                
            end
            U = U + Fk.'*Hk;
            Hk = zeros(1,N/2);   
        end
        
    otherwise
        disp('Unknown method called (function U_vect).')
end
