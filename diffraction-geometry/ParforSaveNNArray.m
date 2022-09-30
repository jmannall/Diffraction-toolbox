%% Save the results for neural network geometry array from a parfor loop

function ParforSaveNNArray(fname, result, geometry, NNinput)
  save(fname, 'result', 'geometry', 'NNinput');
end