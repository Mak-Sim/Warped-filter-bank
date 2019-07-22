function [Y] = warped_cmfb_interpolation(wcmfb_params, X)
%warped_cmfb_decimation Выполняет интерполяцию канальный сигналов банка
%фильтров

sb_pt = wcmfb_params.SubbandSamples;
Sk = wcmfb_params.SubsamplingFactors;   % Коэффициенты децимации
M = wcmfb_params.NumberOfChannel;       % Число каналов банка фильтров
Y = zeros(size(X));
for i=1:M
    tmp = Sk(i)*upsample(X(i,1:sb_pt(i)),Sk(i));
    len_tmp = length(tmp);
    len_tmp = min(len_tmp, size(Y,2));
    Y(i,1:len_tmp) = tmp(1:len_tmp);
end

end