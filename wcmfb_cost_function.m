function [ e, grad, Hess ] = wcmfb_cost_function(h_vec, w, B, N2)
%WCMFB_COST_FUNCTION -- function used by matlab solver to optimize filter
% prototype.

load R1_all_cur;
load R2_all_cur;
load I1_all_cur;
load I2_all_cur;

grad= zeros(N2,1);
Hess= zeros(N2,N2);
e   = 0.0;
        
for i=1:length(w)
    R1(:,:) = R1_all(i,:,:);
    R2(:,:) = R2_all(i,:,:);
    I1(:,:) = I1_all(i,:,:);
    I2(:,:) = I2_all(i,:,:);
    
    r_val   = h_vec'*R1*h_vec;
    i_val   = h_vec'*I1*h_vec;
    r_grad  = R2*h_vec;
    i_grad  = I2*h_vec;
    r_hess  = R2;
    i_hess  = I2;
    f   = (r_val).^2 + (i_val).^2 - 1;
    
    
    % Function value
    val = f^2;
    e       = e + B(i)*val;
    
    if nargout>=2
        % Gradient
        gr  = 2*f*(2*r_val*r_grad + 2*i_val*i_grad);
        grad    = grad  + B(i)*gr;
        if nargout==3
            % Hessian
            H   = 2*f*(2*(r_grad*r_grad') + 2*(i_grad*i_grad') + 2*r_hess*r_val + 2*i_hess*i_val)+...
                2*(gr*gr');
            
            Hess    = Hess + 2*B(i)*H;
        end
    end
    
    

    
    
    
end

end