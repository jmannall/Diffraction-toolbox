function idx = CreateNNIdx(hP, tP, nP)
    nP = rmfield(nP, 'seed');
    nP = rmfield(nP, 'savePath');
    % idx = DataHash({hP, tP, nP});
    disp(['tP Loss: ' DataHash({tP.lossFunc})])
    disp(['tP Data: ' DataHash({tP.dataFunc})])
    disp(['tP Test: ' DataHash({tP.testFunc})])
    tP = rmfield(tP, 'lossFunc');
    tP = rmfield(tP, 'dataFunc');
    tP = rmfield(tP, 'testFunc');
    idx = DataHash({hP, tP, nP});
end