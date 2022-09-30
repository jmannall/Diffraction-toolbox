function vari = ReshapeForParfor(var, extra, filesPerSave)
    var = [var; zeros(extra, 1)];
    vari = reshape(var,filesPerSave,[]);
end