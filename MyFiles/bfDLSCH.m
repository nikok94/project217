function cwout = bfDLSCH(trblk, modulation, nlayers, rv, outlen, codeRate)
% out представлено в ячейки которая содержит одно или два кодовых слова
% представленных в виде вектора столбца, кодовые слова могут иметь
% различную длину
% nlayers количесво слоев [1 : 8] TS 38.214 п. 5.1.1.1
% indata содержит 1 или 2 транспортных блока, если nlayers > 4, то  trblk
% состоит из 2 транспортных блоков, иначе из 1-го
% modulation ячейка содержащая один или два вида модуляции для формирования
% каждого кодового слова
% rv ячейка из 1:2 версий резервирования выхода, каждая из которых может 
% быть равна  RV [0:3]
if nlayers > 8
    error('bfDLSCH: invalid nlayers');
end

if nlayers > 4
    if length(trblk) < 2 || length(modulation) < 2
        error('bfDLSCH: invalid input data size');
    end
    is2CW = true;
    nl1 = floor(nlayers/2);
    nl2 = ceil(nlayers/2);
else
    is2CW = false;
    nl1 = nlayers;
end

  if is2CW
    cw1 = bfGetDLSCHCodeword(trblk{1}, codeRate{1}, rv{1}, modulation{1}, nl1, outlen{1}); 
    cw2 = bfGetDLSCHCodeword(trblk{2}, codeRate{2}, rv{2}, modulation{2}, nl2, outlen{2});
    cwout = {cw1, cw2};
  else    % одно кодовое слово
    cwout = {bfGetDLSCHCodeword(trblk{1}, codeRate{1}, rv{1}, modulation{1}, nl1, outlen{1})};
  end

end