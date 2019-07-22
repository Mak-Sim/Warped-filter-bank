function [ x ] = warped_cmfb_synthesis(wcmfb_params, X)
%warped_cmfb_synthesis_by_sum -- Выполняет реконструкцию сигнала по его  
%   субполосным компонентам путем суммирования.

% 1. Параметры
eq_filt = wcmfb_params.PhaseCorr;   % Фильтр коррекции ФЧХ
h = wcmfb_params.CoeffOfFP;          % Коэффициенты фильтра-прототипа
M = wcmfb_params.NumberOfChannel;    % Число каналов банка фильтров
N = wcmfb_params.OrderOfFP;          % Порядок фильтра-прототипа (ФП)
m = wcmfb_params.OrderOfPPComponents;% Порядок полифазных компонент ФП
alpha = wcmfb_params.Alpha;

Sig_length = size(X,2);
APC_delay_1 = zeros(N,1);
APC_delay_2 = zeros(N,1);
x = zeros(Sig_length,1);
c = zeros(M,2*M);

% 2. Коэффициенты обобщенной полифазной структуры
h = reshape(h,2*M,length(h)/(2*M));
h(:,2:2:size(h,2))=-h(:,2:2:size(h,2));     %G(z)=G(-z)
h = h(:);

% 3. Подготовка модуляционной матрицы
for k =0:M-1,
    for l=0:2*M-1,
        c(k+1,l+1) = 2*cos( (2*k+1)*(pi/(2*M))*(l - (N-1)/2) + (-1)^(k)*(pi/4) );
    end;
end;

% 4. Реализация банка фильтров синтеза
for n=1:Sig_length
    % 4.1 Косинусная демодуляция
    Curr_X = c'*X(:,n);     

    % 4.2 Полифазная фильтрация
    Xk = repmat(Curr_X, m, 1);
    Xk = Xk.*h;
    
    % 4.3 Цепочка фазовых звеньев
    R1 = Xk(1);
    for i=1:N-1
        R3 = APC_delay_1(i);
        APC_delay_1(i) = R1;
        R2 = APC_delay_2(i);
        R1 = (R2-R1)*alpha + R3;
        APC_delay_2(i) = R1;
        
        R1 = R1 + Xk(i+1);
    end
    
    x(n) = R1;
end
% 5. Выравнивание ФЧХ
x = filter(eq_filt,1,x);           
end

