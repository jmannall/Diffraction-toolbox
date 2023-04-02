function idx = CreateNNIdx(hP, tP, nP)
    nP = rmfield(nP, 'seed');
    % nP = rmfield(nP, 'savePath');
    idx = DataHash({hP, tP, nP});
end