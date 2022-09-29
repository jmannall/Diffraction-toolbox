%% Radius example (+ change with fs)
close all

% S and R reversible, however not simply S + R;
% Example of how result changes with fs.
SensitivityRadius(48000)
SensitivityRadius(96000)
SensitivityRadius(128000)
SensitivityRadius(256000)

%% Wedge length example (+ change with fs)
close all

% Identical from 10m (plot 10 - 100), very similar from 5m (plot 1 - 100),
% consider logarithmically spaced (0.1 - 100)
SensitivityWedgeLength(48000)
SensitivityWedgeLength(96000)
SensitivityWedgeLength(128000)
SensitivityWedgeLength(256000)

%% zR example
close all

% Symmetrical in the z axis from the midpoint of the wedge
% Removed directsound line 238 EDfindconvexGApaths due to bug in EDpoinpla
% when planes with different numbers of edges, 256000 requires 8192
% frequencies in the fft when zR < 0
% Decrease in level as delta z increases. Change in shape once zR < 0
% (doesn't mean apex no longer on edge)
SensitivityzR(48000)
SensitivityzR(96000)
SensitivityzR(128000)
SensitivityzR(256000)

%% zRzSsame example
close all

% No change until apex no longer on the edge when same. When a constant
% delta z, change in shape and drop in LF when zR < 0. Constant shape,
% decreasing fc and decrease in level once apex point no longer on the
% wedge.
SensitivityzRzSsame(48000)
SensitivityzRzSsame(96000)
SensitivityzRzSsame(128000)
SensitivityzRzSsame(256000)

%% zRzSopp example
close all

% Increasing the number of frequencies in the BTM fft increases the number
% of frequency responses shown at high fs. As before, decrease in level as 
% delta z increases. Change in shape once delta z > wedgeLength (zR and zS
% outside z bounds of the wedge)
SensitivityzRzSopp(48000)
SensitivityzRzSopp(96000)
SensitivityzRzSopp(128000)
SensitivityzRzSopp(256000)

%% Suggested parameters

%{

Wedge length
Wedge index
Minimum angle
Bending angle
delta z
z midpoint (or apex point?) - run some freq tests to find out
radiusR and radiusS

Find the min / max values of each to be considered (either where no change
happens or level negilible) - look up a dB value where magnitude can be
perceptually considered neglible.

}%
