function idx = CreateNNIdx(hP, tP, nP)
    nP = rmfield(nP, 'seed');
    nP = rmfield(nP, 'savePath');
    tP = rmfield(tP, 'lossFunc');
    tP = rmfield(tP, 'dataFunc');
    tP = rmfield(tP, 'testFunc');
    idx = DataHash({hP, tP, nP});
end