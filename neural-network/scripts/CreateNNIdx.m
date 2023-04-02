function idx = CreateNNIdx(hP, tP, nP)
    nP = rmfield(nP, 'seed');
    nP = rmfield(nP, 'savePath');
    idx = DataHash({hP, tP, nP});
    disp(['hP: ' DataHash({hP})])
    disp(['tP: ' DataHash({tP})])
    disp(['nP: ' DataHash({nP})])
end