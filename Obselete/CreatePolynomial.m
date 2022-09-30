function poly = CreatePolynomial(roots)
    numInputs = size(roots,2);

    % Find and initialise correct order polynomial
    n = size(roots,1);
    poly = dlarray(zeros(n+1,numInputs));
    
    % Vector for root indices
    vector = 1:n;

    % Find polynomial coeffieicnts
    poly(1,:) = 1;
    poly(2,:) = -sum(roots,1);
    for i = 2:n
        idx = nchoosek(vector, i);
        evenOdd = 1 - 2 * rem(i,2);
        test = prod(roots(idx,:),1);
        poly(i + 1,:) = evenOdd * sum(prod(roots(idx,:),1), 1);
    end
end
