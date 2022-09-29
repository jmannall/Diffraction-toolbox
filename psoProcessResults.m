function [filterStable, tfiir, fveciir, b, a] = psoProcessResults(BestSol, fs)

    z = [BestSol.Position(1);BestSol.Position(3)];
    p = [BestSol.Position(2);BestSol.Position(4)];
    k = BestSol.Position(5);
    
    [b, a] = zp2tf(z, p, k);
    [tfiir, fveciir] = IIRFilter(z, p, k, fs);
    
    stabilityCheck = abs(p);
    
    if stabilityCheck <= 1
        filterStable = true;
        str = 'is';
    else
        filterStable = false;
        str = 'is not';
    end
    
    disp(['Filter ' str ' stable.']);
end