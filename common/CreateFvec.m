function fvec = CreateFvec(fs, nfft)
    fvec = fs/nfft*(0:nfft/2-1);
end