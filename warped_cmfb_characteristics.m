function [xr] = warped_cmfb_characteristics(h0, h_opt, alpha, M, lang, Npt)
%% Filter bank characterestics
% Parameters:
% h0    -- initial filter prototype
% h_opt -- optimized filter-prototype
% lang  -- language

N = length(h0);             % Filter-prototype order
N2 = N/2;                   
w  = linspace(0, pi, Npt);  % Frequency grid

h_new = h_opt(N/2+1:end)';
h_old = h0(N/2+1:end)';

load U_matr_cur;
load U_alias_cur;
        
T_all_after    = zeros(1,length(w));
T_all_before   = zeros(1,length(w));
T_alias_after  = zeros(1,length(w));
T_alias_before = zeros(1,length(w));

for i=1:length(w)    
    U(:,:) = U_matr(i,:,:);
    T_all_after(i)  = abs(h_new'*U*h_new);
    T_all_before(i) = abs(h_old'*U*h_old) ;
    % Alias components
    U(:,:) = U_alias(i,:,:);
    T_alias_after(i)  = abs(h_new'*U*h_new);
    T_alias_before(i) = abs(h_old'*U*h_old);
end

% Transfer fucntions
figure('Position',[300 300 600 500]);
axes('Position',[0.12, 0.55, 0.85, 0.4]);
hold on;
plot(w/(2*pi),20*log10(abs(T_all_before)),'LineWidth',1.5,'color',[0.04 0.52 0.78]);
plot(w/(2*pi),20*log10(abs(T_all_after)),'LineWidth',1.5,'color',[0.00 0.12 0.38]);
ylabel('$20\log_{10}|\mathrm{T}_{\mathrm{all}}(\omega)|,\, \mathrm{dB}$', 'Interpreter', 'Latex','FontSize',14)
xlabel('$\omega / 2\pi$', 'Interpreter', 'Latex','FontSize',14);
xlim([0 w(end)/(2*pi)]);
grid on;
switch lang
    case 'eng'
        legend({'T_{all}(e^{j\omega}) funciton of initial warped CMFB',...
               'T_{all}(e^{j\omega}) function of optimized warped CMFB'},'FontSize',12);
    otherwise
        legend({'Функция T_{all}(e^{j\omega}) исходного НКМБФ',...
               'Функция T_{all}(e^{j\omega}) оптимизированного НКМБФ'},'FontSize',12);
end


% Aliasing components transfer functions
axes('Position',[0.12, 0.1, 0.85, 0.4]);
hold on;
plot(w/(2*pi),20*log10(abs(T_alias_before)),'LineWidth',1.5,'color',[0.04 0.52 0.78]);
plot(w/(2*pi),20*log10(abs(T_alias_after)),'LineWidth',1.5,'color',[0.00 0.12 0.38]);
switch lang
    case 'eng'
        legend({'T_{alias}(e^{j\omega}) funciton of initial warped CMFB',...
               'T_{alias}(e^{j\omega}) funciton of optimizedwarped CMFB'}, 'FontSize',12);
    otherwise
        legend({'Функция T_{alias}(e^{j\omega}) исходного НКМБФ',...
            'Функция T_{alias}(e^{j\omega}) оптимизированного НКМБФ'},'FontSize',12);
end
ylabel('$20\log_{10}|\mathrm{T}_{\mathrm{alias}}(\omega)|, \, \mathrm{dB}$', 'Interpreter', 'Latex','FontSize',14)
xlabel('$\omega / 2\pi$', 'Interpreter', 'Latex','FontSize',14);
xlim([0 w(end)/(2*pi)]);
grid on;
%% 

figure('Units','pixels', 'Position',[200 200 1000 400]);
f_range=w/(2*pi);
h_i = [flipud(h_opt)' h_opt'];
G = my_zerophase(h_i,f_range,N);
plot(f_range, 20*log10(abs(G)),'LineWidth',2.5,'color',[0.0 0.0 0.0]);
xlim([0 0.5]); ylim([-130 5]); grid on; 
hold on;
G = my_zerophase(h0,f_range,N);
plot(f_range, 20*log10(abs(G)),'LineWidth',2.0,'color',[0.04 0.52 0.78],'LineStyle','--');
switch lang
    case 'eng'
        legend('Optimized filter prototype','Initial filter prototype');
    otherwise
        legend('Оптимизированный фильтр-прототип','Исходный фильтр-прототип');
end
ylabel('$H(e^{j\omega}), \mathrm{dB}$', 'Interpreter', 'Latex','FontSize',14)
xlabel('$\omega / 2\pi$', 'Interpreter', 'Latex','FontSize',14);

% Filter bank magnitude response
xr=cmfb_system(h_opt,M,[1 zeros(1,1535)],alpha);


end

