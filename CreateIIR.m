function CreateIIR(z, p, k, update, fs)
    
    S.fs = fs;
    S.fh = figure('units','pixels', 'position',[10 10 1900 1020],...
                  'menubar','none',...
                  'name','slider_plot',...
                  'numbertitle','off',...
                  'resize','off');    
    S.ax = axes('units', 'pixels', 'position',[100 250 1000 700]);
    
    S.fcn = @(z, p, k) IIRFilter(z, p, k, S.fs);
    S.z = z;
    S.p = p;
    S.k   = k;
    [tfmag, fvec] = IIRFilter(S.z, S.p, S.k, S.fs);
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
    % Slider for z:
    S = CreateZ(S,length(S.z));
    % Slider for p:
    S = CreateP(S,length(S.p));
    % Slider for k:
    S = CreateK(S);
    S.update = update;
    S.update(S);
    guidata(S.fh, S);  % Store S struct in the figure
end