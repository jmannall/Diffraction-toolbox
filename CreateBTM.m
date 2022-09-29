function CreateBTM(w, bA, mA, update, fs)
    
    S.fs = fs;
    S.fh = figure('units','pixels', 'position',[10 10 1900 1020],...
                  'menubar','none',...
                  'name','slider_plot',...
                  'numbertitle','off',...
                  'resize','off');    
    S.ax = axes('units', 'pixels', 'position',[100 250 1000 700]);
    
    S.fcn = @(w, bA, mA) BTM(w, bA, mA, S.fs);
    S.w = w;
    S.bA = bA;
    S.mA   = mA;
    [tfmag, fvec] = BTM(S.w, S.bA, S.mA, S.fs);
    figure(S.fh);
    S.x   = fvec;
    S.LN  = semilogx(fvec,tfmag);
    grid on
    S.ax.XLim = [20 24000];
    S.ax.YLim = [-80 20];
    set(get(gca,'XLabel'),'String','Frequency (Hz)')
    set(get(gca,'YLabel'),'String','Magnitude (dB)')
    % Slider for w:
    S = CreateSliders(S);
    % Slider for p:
    S.update = update;
    S.update(S);
    guidata(S.fh, S);  % Store S struct in the figure
end