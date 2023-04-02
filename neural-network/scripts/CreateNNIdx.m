function idx = CreateNNIdx(hP, tP, nP)
    nP = rmfield(nP, 'seed');
    idx = DataHash({hP, tP, nP});
end