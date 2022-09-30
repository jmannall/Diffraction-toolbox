%% Save the results for geometry array from a parfor loop

function ParforSaveArray(fname, result, geometry)
  save(fname, 'result', 'geometry');
end