function w_def=freq_warp(w,alpha)
%% Mapping w-> \varphi(w)
w_def = - w +2*atan(alpha*sin(w)./(alpha*cos(w)-1));
end