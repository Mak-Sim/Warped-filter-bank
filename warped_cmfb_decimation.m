function [Y, sb_pt] = warped_cmfb_decimation(wcmfb_params, X)
%warped_cmfb_decimation Выполняет децимацию канальный сигналов банка
%фильтров

Sk = wcmfb_params.SubsamplingFactors;   % Коэффициенты децимации
M = wcmfb_params.NumberOfChannel;       % Число каналов банка фильтров
sb_pt = zeros(1,M);                     % Количество отсчетов в каждой полосе
Y = zeros(size(X));

for i=1:M
    tmp = downsample(X(i,:),Sk(i));
    sb_pt(i) = length(tmp);
    Y(i,1:sb_pt(i)) = tmp;
end

end

