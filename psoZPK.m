%
% Copyright (c) 2016, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YTEA101
% Project Title: Particle Swarm Optimization Video Tutorial
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer and Instructor: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

clc;
clear;
close all;
set(0, 'DefaultLineLineWidth', 2);
set(groot, 'defaultLineMarkerSize', 10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SingleWedgeArray input data

disp('Create training data');

wedgeLength = 20;
radiusS = 1;
radiusR = 1;
zS = 10;
zR = 10;
fs = 96000;
wedge = 190:10:350;
bendingAngle = 190:10:350;
minAngle = 0:10:180;

[result, geometry] = SingleWedgeArray(wedgeLength, radiusS, radiusR, zS, zR, fs, wedge, bendingAngle, minAngle);

result = rmfield(result,'i');

%% Check for previous attempts

numResults = length(result);
count = 0;

%% Problem definition

problem.VarMin =  -0.1;  % Lower Bound of Decision Variables
problem.VarMax =  0.1;   % Upper Bound of Decision Variables

%% Parameters of PSO

params.MaxIt = 20;         % Maximum Number of Iterations
params.nPop = 100;          % Population Size (Swarm Size)
params.w = 0.5;             % Intertia Coefficient
params.wdamp = 0.99;        % Damping Ratio of Inertia Coefficient
params.c1 = 2;              % Personal Acceleration Coefficient
params.c2 = 2;              % Social Acceleration Coefficient
params.ShowIterInfo = true; % Flag for Showing Iteration Information

% Create file info
mFile = mfilename('fullpath');
[inFilePath,fileStem] = fileparts(mFile);

modelNum = 10;

index = DataHash(result);
saveStem = [fileStem, '_model', num2str(modelNum)];
savepath = ['results\', saveStem];
loadpath = ['results\', saveStem, '.mat'];

w = geometry.wedge;
bA = geometry.bendingAngle;
mA = geometry.minAngle;

x = [w, bA, mA, w.*mA, w.*bA, bA.*mA, tanh(1 ./ (2 * pi * (mA + 1))), tanh(180 ./ (2 * pi * w)), tanh(180 ./ (2 * pi * bA)), bA ./ (w - mA - 180)];
y = [    7.47360379643697	5.21859751278169	1.30264563917495	1.40344344687535	-1.15608132279468
-0.0136718499236927	-0.00421070249479635	-0.00103645437641610	-0.000546982701864072	0.000165234228446824
-0.0182970058551008	-0.0115102977216140	-0.00104841952414615	-0.00117980859117657	0.00217525661816178
0.0188605667809560	0.00838797712547043	0.00177129985516914	0.00121963366427006	-0.00183418195958101
-1.61944455294202e-05	-9.44159784517010e-06	-2.12089700195233e-06	-1.45910450890173e-06	2.09351394135796e-06
4.89114958507842e-05	1.87454076366832e-05	3.72831993872321e-06	2.40267617722296e-06	-7.67915729917792e-07
-6.10414975043342e-05	-2.20999599689203e-05	-4.42854793322601e-06	-2.96005607081956e-06	4.58711226956381e-06
-1.06084891499443	-0.294874547015583	-0.0864336048687234	-0.0439328512504777	0.0506537647109190
-4.23830564797583	1.62179964177513	-0.259047902582607	0.231863815023599	0.362521834016905
-14.8497052410362	-15.5186431779468	-0.176231505270719	-1.32868385334778	5.24002459389383
-0.0282545839451643	-0.0162775916345191	-0.00421591977421231	-0.00286122363632939	0.00413399808629926];
x = [w, bA, mA];
y = [-0.383913610540019	0.433687677402994	0.918783004755412	0.943480259257778	0.252639338727523
-0.000619985065659293	-9.81938775405013e-05	2.44892735842296e-05	7.07950477238117e-06	-0.000290577612012910
0.00313936734008067	0.00153598423066938	0.000100687433299577	0.000156240483506246	-0.000506838114513924
-9.83120110689203e-05	-3.98674043972430e-05	-1.87024877200234e-05	-4.78452303549349e-06	8.14515403266362e-06];
x = [sind(w), sind(bA), cosd(bA)];
y = [-0.0866180061247715	0.617922228642913	0.940501680651416	0.967775671149180	0.0851291904082404
-0.0247461048070945	-0.00819944028513898	-0.000645580571043332	-0.000230669623200549	-0.0162988089184426
-0.263122087684439	-0.192052895385786	-0.00752202541860864	-0.0191995378786221	0.0822392922211800
0.0762725942990091	0.0194817656754405	0.00490062494510180	0.00242288059272727	0.00318165841847493];
x = [sind(w), sind(bA), cosd(mA), tanh(1 ./ (2 * pi * mA + 1)), tanh(180 ./ (2 * pi * w)), tanh(180 ./ (2 * pi * bA)), bA ./ (w - mA - 180)];
y = [1.21367066470855	0.936853198146257	1.02184352746259	1.00808593704766	-0.0579097950726242
-0.103723448445735	-0.0235533451946467	-0.00767144127547034	-0.00466017378370826	0.0115221891705918
0.186841299597486	-0.0606780022755254	0.0219655522284815	-0.00181113791158392	0.0357374904039806
0.000980716584579317	-0.0367734753896261	-0.00290104565339755	-0.00512181380294063	0.0114176799179053
-0.200166263068927	-0.0696927227841497	-0.0199947005605236	-0.0102672340745358	0.00428933862943399
2.19430670247603	1.65136544973068	0.140407171661555	0.209683795768729	0.505721752476406
-9.10729878394253	-2.56695520435767	-0.515690739260244	-0.326182069870898	0.521072430507535
-0.0165976928856368	-0.0134103263870067	-0.00116974242158492	-0.00117417806084326	0.00394762599128860];
% var = [90, 0.2, 28.6478, 28.6478];
% tanhFunctions = @(x) TanhFunctions(x, w, bA, mA);
% x = [sind(w), sind(bA), cosd(mA), tanh(1 ./ (var(1) * mA + var(2))), tanh(var(3) ./ w), tanh(var(4) ./ bA), bA ./ (w - mA - 180)];
% startingValues = [1.21112877699547	0.936720343905359	1.02210231863349	1.01116021727733	-0.0531032494830244
% -0.103056624504718	-0.0259022275392574	-0.00977316093211119	-0.00392100975579621	0.0109865004870745
% 0.187720233098162	-0.0610349710368644	0.0217619989783887	-0.00134641349449773	0.0447064930319670
% 0.00985848461414291	-0.0242323243710640	-0.00259284144200428	-0.00428356493291683	0.0106990597016850
% -0.151095117878778	-0.0518214618805442	-0.0144442546721419	-0.00802025463124257	0.00893746740818370
% 2.19412372046271	1.63786902942382	0.139356932980586	0.209998191933889	0.509315255157357
% -9.10539565330663	-2.56369115310443	-0.517627818719466	-0.328561639534056	0.520155211138544
% -0.0152328592703426	-0.0117657748854685	-0.00274611249888769	-0.00203689505906371	0.00354185465096816
% var(1) * ones(1,5)
% var(2) * ones(1,5)
% var(3) * ones(1,5)
% var(4) * ones(1,5)];
x = [tanh(180 ./ (2 * pi * bA)), tanh(180 ./ (2 * pi * w)), sind(bA), sind(bA) .* cosd(mA), sind(w) .* sind(bA), cosd(mA), bA, sind(w), mA, bA ./ (w - mA - 180), w, bA.^2, bA .* mA, w .* mA, w .* mA .* bA];
y = [-19.1787376604614	-8.76921870927808	0.707374965801628	1.62285686387075	0.242931350864402
0	-12.3317197063857	2.26038889103474	-1.27815582738716	6.32998249847313
0	1.55981938825958	-2.08325040939969	-1.61173418690732	5.75634832306518
2.41543657418864	1.40623598114196	-0.0393610958375375	-0.0161240020319533	-0.340234588990170
-0.550322649409118	-0.0268529871734618	0.0129806627875203	0.0127554510309902	-0.0185225281600983
-0.410348740826129	-0.0404042039018342	-0.00903038329740990	-0.00191269497093013	-0.0225606835373325
0.115086778185577	0.162540186674393	0.0756719574990103	0.0422122763520750	-0.0581829503348712
0.159490859777717	0.0931847981979880	0.00206198282471877	-0.000622985814704817	-0.0201005224934598
-0.163690983391891	-0.00951171064822728	0.0163548872358231	0.0187665018987107	-0.0571261736077179
-0.0138383814837911	0.0303526783975991	0.0101142519882353	0.00746146828577909	-0.0101811348104293
-0.0273748852085714	-0.0134478415425385	-0.00392350662638207	-0.00262223733542155	0.00299665963755582
-0.00211449771666270	-0.000258812372455974	-0.00103356268265163	-0.000801242499991702	0.00211944558308655
-0.000288854161768315	-0.000179580844505504	-1.56440941441797e-06	5.52312318502170e-07	4.14256021239727e-05
0.000127574189082770	-0.000123056525320265	-4.31874631963212e-05	-3.22974362582385e-05	4.60168682017226e-05
5.58000499242310e-05	-7.23230900448472e-05	-2.23062272362070e-05	-1.74102553717663e-05	2.39481014995961e-05
-3.67508109470548e-07	3.20212114238252e-07	1.06453404896073e-07	8.16039590077753e-08	-1.17761186608160e-07];

input = [ones(numResults,1), tanh(180 ./ (2 * pi * bA)), tanh(180 ./ (2 * pi * w)), sind(bA), sind(bA) .* cosd(mA), sind(w) .* sind(bA), cosd(mA), bA / 10, sind(w), mA, bA ./ (w - mA - 180), w / 10, (bA.^2) / 1e4, (bA .* mA) / 1e3, (w .* mA) / 1e3, (w .* mA .* bA) / 1e5];
startingValues = [-19.1787376604614	-8.76921870927808	0.707374965801628	1.62285686387075	0.242931350864402
0	-12.3317197063857	2.26038889103474	-1.27815582738716	6.32998249847313
0	1.55981938825958	-2.08325040939969	-1.61173418690732	5.75634832306518
2.41543657418864	1.40623598114196	-0.0393610958375375	-0.0161240020319533	-0.340234588990170
-0.550322649409118	-0.0268529871734618	0.0129806627875203	0.0127554510309902	-0.0185225281600983
-0.410348740826129	-0.0404042039018342	-0.00903038329740990	-0.00191269497093013	-0.0225606835373325
0.115086778185577	0.162540186674393	0.0756719574990103	0.0422122763520750	-0.0581829503348712
1.59490859777717	0.931847981979880	0.0206198282471877	-0.00622985814704817	-0.201005224934598
-0.163690983391891	-0.00951171064822728	0.0163548872358231	0.0187665018987107	-0.0571261736077179
-0.0138383814837911	0.0303526783975991	0.0101142519882353	0.00746146828577909	-0.0101811348104293
-0.0273748852085714	-0.0134478415425385	-0.00392350662638207	-0.00262223733542155	0.00299665963755582
-0.0211449771666270	-0.00258812372455974	-0.0103356268265163	-0.00801242499991702	0.0211944558308655
-2.88854161768315	-1.79580844505504	-1.56440941441797e-02	5.52312318502170e-03	4.14256021239727e-01
0.127574189082770	-0.123056525320265	-4.31874631963212e-02	-3.22974362582385e-02	4.60168682017226e-02
5.58000499242310e-02	-7.23230900448472e-02	-2.23062272362070e-02	-1.74102553717663e-02	2.39481014995961e-02
-3.67508109470548e-02	3.20212114238252e-02	1.06453404896073e-02	8.16039590077753e-03	-1.17761186608160e-02];

% startIndex = abs(y) >= 1e-2;
% input = {[x(:,startIndex(:,1))], [x(:,startIndex(:,2))], [x(:,startIndex(:,3))], [x(:,startIndex(:,4))], [x(:,startIndex(:,5))]};
% startingValues = {[y(startIndex(:,1),1)], [y(startIndex(:,2),2)], [y(startIndex(:,3),3)], [y(startIndex(:,4),4)], [y(startIndex(:,5),5)]};

problem.nVar = [size(startingValues)];       % Number of Unknown (Decision) Variables

% z = [(input{1} * startingValues{1})'; (input{3} * startingValues{3})'];
% p = [(input{2} * startingValues{2})'; (input{4} * startingValues{4})'];
% k = (input{5} * startingValues{5})';

test = exist(loadpath, "file");

if test == 2
    load(loadpath, "BestSol", "BestCosts");
    count = numResults;
    disp('Result loaded from save');
else
    disp('Begin PSO');
    problem.CostFunction = @(x) ZPKError(x, input, [result.tfmag], fs);  % Cost Function
    problem.ConstraintFunction = @(x) ZPKConstraints(x, input); % Constraint Function
    problem.InitialiseFunction = @(x) ZPKInitialise(x, startingValues); % Initialisation Function
    
    %% Calling PSO
    
    out = PSOWithConstraints(problem, params);
    
    BestSol = out.BestSol;
    BestCosts = out.BestCosts;

    save(savepath, "BestSol", "BestCosts", "input")
end

f1 = figure(1);
movegui(f1,'southwest');
% plot(BestCosts, 'LineWidth', 2);
semilogy(BestCosts);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

prediction = (input*BestSol.Position)';

z = [prediction(1,:); prediction(3,:)];
p = [prediction(2,:); prediction(4,:)];
k = [prediction(5,:)];

for i = 1:numResults
    [tfmagModel, fvecModel] = IIRFilter(z(:,i), p(:,i), k(i), fs);

    f2 = figure(2);
    movegui(f2,'south');
    semilogx(fvecModel, tfmagModel)
    hold on
    semilogx(result(i).fvec, result(i).tfmag)
    hold off
    xlim([20, 20000])
    ylim([-40 0])
        
    [b, a] = zp2tf(z(:,i), p(:,i), k(i));
    ts = 1 / fs;
    sys = tf(b, a, ts);

    f3 = figure(3);
    movegui(f3,'southeast')
    h = pzplot(sys);
    grid on
end