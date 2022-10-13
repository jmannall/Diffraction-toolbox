clear all
close all

set(0, 'DefaultLineLineWidth', 2);

nfft = 8192;
fs = 48000;

wallHeightStore = [1, 2, 4, 8, 16];
wallThicknessStore = [0.125, 0.25, 0.5, 1, 1.5, 2, 4, 8, 16];
radiiStore = [0.125, 0.25, 0.5, 1, 1.5, 2, 4, 8, 16];
% wallThicknessStore = 2;
% radiiStore = 1;
numRadii = length(radiiStore);
numWalls = length(wallThicknessStore);
numHeights = length(wallHeightStore);
Error = zeros(1, numWalls);
input = zeros(numWalls * numRadii * numHeights, 4);
for m = 1:numHeights
for k = 1:numRadii
    for j = 1:numWalls
        radius = radiiStore(k);
        wallThickness = wallThicknessStore(j);
        wallHeight = wallHeightStore(m);
        controlparameters = struct('nfft', nfft, 'fs', fs, 'difforder', 2);
        radiusS = radius;
        radiusR = radius;
        radiusTotal = radiusS + radiusR;
        [zS, zR] = deal(wallHeight / 2);
            
        minSource = asind((wallThickness / 2) / (radiusR + wallThickness / 2));
        wallIndex =  360 - minSource;
        thetaS = 10;
        %thetaS = 89.99;
        thetaR = 360 - thetaS;
        %thetaR = 270.01;
        % thetaS = 30;
        % thetaR = 200;
        
        % radiusR = radiusR + wallThickness / 2;
        % radiusS = radiusS + wallThickness / 2;
        
        [ir, tfmag, tvec, fvec, tfcomplex] = SingleWall(wallHeight, wallThickness, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, false);
        
        wedgeIndex = 270;
        
        edge1 = [0, 0];
        edge2 = [0, -wallThickness];
        
        source = [radiusS * cosd(thetaS), radiusS * sind(thetaS)];
        receiver = [radiusR * cosd(thetaR), radiusR * sind(thetaR) - wallThickness];
        
        radiusS = norm(source - edge1);
        radiusR1 = wallThickness;
        radiusR2 = wallThickness / 2;
        
        thetavR = wedgeIndex - 0.01;
        
        controlparameters .difforder = 1;
        
        [ir1.one, tfmag1.one, tvec1.one, fvec1.one, tfcomplex1.one] = SingleWedge(wallHeight, wedgeIndex, thetaS, thetavR, radiusS, radiusR1, zS, zR, controlparameters, false);
        [ir1.two, tfmag1.two, tvec1.two, fvec1.two, tfcomplex1.two] = SingleWedge(wallHeight, wedgeIndex, thetaS, thetavR, radiusS, radiusR2, zS, zR, controlparameters, false);
        
        radiusS1 = wallThickness;
        radiusS2 = wallThickness / 2;
        radiusR = norm(receiver - edge2);
        
        thetavS = 0.01;
        thetaR = 270 - thetaS;
        
        [ir2.one, tfmag2.one, tvec2.one, fvec2.one, tfcomplex2.one] = SingleWedge(wallHeight, wedgeIndex, thetavS, thetaR, radiusS1, radiusR, zS, zR, controlparameters, false);
        [ir2.two, tfmag2.two, tvec2.two, fvec2.two, tfcomplex2.two] = SingleWedge(wallHeight, wedgeIndex, thetavS, thetaR, radiusS2, radiusR, zS, zR, controlparameters, false);
        
        diff2 = tfmag.diff2;
        diffPart1 = tfmag1.two.diff1;
        diffPart2 = tfmag2.two.diff1;
        factor1 = (radiusS + radiusR1 + radiusS1 + radiusR) / (radiusS * radiusR1 * radiusS1 * radiusR);
        factor2 = (radiusS + radiusR2 + radiusS2 + radiusR) / (radiusS * radiusR2 * radiusS2 * radiusR);
        diff1_1 = mag2db(abs(tfcomplex1.one.diff1 .* tfcomplex2.one.diff1));
        diff1_2 = mag2db(abs(tfcomplex1.two.diff1 .* tfcomplex2.two.diff1));
        
        n = 100;
        freqs = logspace(log10(20), log10(20e3), n);
        idx = zeros(1, n);
        for i = 1:n
            store = find(fvec > freqs(i));
            idx(i) = store(1);
        end
        
        Error(j, k) = mean(diff1_2(idx) - diff2(idx));
        input(numHeights * (m - 1) + numRadii * (j - 1) + k, :) = [log(wallThickness), log(radiusTotal), log((radiusTotal * wallThickness) / (radiusTotal + wallThickness)), wallHeight];
    end
end
end

figure
semilogx(fvec, [tfmag.diff2])
hold
semilogx(fvec, diff1_2)
semilogx(fvec, [tfmag1.two.diff1])
semilogx(fvec, [tfmag2.two.diff1])
xlim([20 20000])
ylim([-100 0])

%% Plot
close all

X = [ones(numWalls * numRadii * numHeights, 1), input];
Y = reshape(Error, [], 1);
disp('break');
[beta,Sigma,E,CovB,logL] = mvregress(X,Y);

averageError = sum(E.^ 2,1) / length(E);
worstError = max(E,[],1);

B = beta;
xx = ([linspace(0.1, 20, 100); linspace(0.1, 20, 100); linspace(0.1, 20, 100); linspace(0.5, 20, 100)])';
test = [ones(length(xx), 1),xx];
fits = zeros(size(Error));
for i = 1:numRadii
    fits(:,i) = X((i - 1) * numWalls + 1:i * numWalls,:)*B;
end

figure
h = plot(log(wallThicknessStore),Error,'x',log(wallThicknessStore),fits,'-');
for i = 1:numWalls
    set(h(numWalls + i),'color',get(h(i),'color'))
end

scale = ([1, log(wallThickness), log(radiusTotal), log((wallThickness * radiusTotal) / (wallThickness + radiusTotal), wallHeight)])*B;
diff1_2 = diff1_2 - scale;

figure
semilogx(fvec, diff2)
hold on
semilogx(fvec, diff1_1)
semilogx(fvec, diff1_2)
semilogx(fvec, diffPart1,'--')
semilogx(fvec, diffPart2,'--')
legend('2nd Order BTM', '1st Order BTM in series - nontrue path length', '1st Order BTM in series - true path length', '1st Order no. 1', '1st Order no. 2', 'Location','southwest')

PlotSpectogram([tfcomplex.diff2, tfcomplex.diff2], fvec, [1, 2], [-70 0], 'SecondOrderDiff', true, true);
PlotSpectogram([tfcomplex1.one.diff1 .* tfcomplex2.one.diff1, tfcomplex1.one.diff1 .* tfcomplex2.one.diff1], fvec, [1, 2], [-70 0], 'FirstOrderDiff - non true path length', true, true);
PlotSpectogram([tfcomplex1.two.diff1 .* tfcomplex2.two.diff1, tfcomplex1.two.diff1 .* tfcomplex2.two.diff1], fvec, [1, 2], [-70 0], 'FirstOrderDiff - true path length', true, true);