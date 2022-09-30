%% Create interactive IIR filter plot - IS it possible to generalise these functions

function CreateIIR(z, p, k, update, fs)
    
    S.fs = fs;

    % Format figure
    S.fh = figure('units','pixels', 'position',[10 10 1900 1020],...
                  'menubar','none',...
                  'name','slider_plot',...
                  'numbertitle','off',...
                  'resize','off');    
    S.ax = axes('units', 'pixels', 'position',[100 250 1000 700]);
    
    % Inputs
    S.fcn = @(z, p, k) IIRFilter(z, p, k, S.fs);
    S.z = z;
    S.p = p;
    S.k   = k;

    % Initialise first response
    [tfmag, fvec] = S.fcn(S.z, S.p, S.k);

    % Initialise figure
    figure(S.fh);
    S.x   = fvec;
    S.LN  = semilogx(fvec,tfmag);
    hold on
    S.fc = 1000;
    S.highGain = -80;
    S.FCy = semilogx([10,S.fc], [-3, -3], '-k');
    S.FCx = semilogx([S.fc, S.fc], [-3, -80], '-k');
    S.gradSlope = semilogx([S.fc, 20000], [-3, S.highGain]);
    hold off
    grid on
    S.ax.XLim = [20 24000];
    S.ax.YLim = [-80 20];
    set(get(gca,'XLabel'),'String','Frequency (Hz)')
    set(get(gca,'YLabel'),'String','Magnitude (dB)')

    % Create sliders
    S = SliderZ(S,length(S.z));
    S = SliderP(S,length(S.p));
    S = SliderK(S);
    
    S.update = update;
    S.update(S);
    guidata(S.fh, S);  % Store S struct in the figure
end