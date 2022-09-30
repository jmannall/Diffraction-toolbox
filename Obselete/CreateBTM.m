%% Create interactive BTM plot - IS it possible to generalise these functions

function CreateBTM(S, plotFcn, sliderFcn, updateFcn, fs)

    % Format figure
    S.fh = figure('units','pixels', 'position',[10 10 1900 1020],...
                  'menubar','none',...
                  'name','slider_plot',...
                  'numbertitle','off',...
                  'resize','off');    
    S.ax = axes('units', 'pixels', 'position',[100 250 1000 700]);
    
    % Inputs    
    S.fs = fs;
    S.fcn = plotFcn;

    % Initialise first response
    [tfmag, fvec] = S.fcn(S);


    % Initialise figure
    figure(S.fh);
    S.x   = fvec;
    S.LN  = semilogx(fvec,tfmag);
    grid on
    S.ax.XLim = [20 24000];
    S.ax.YLim = [-80 20];
    set(get(gca,'XLabel'),'String','Frequency (Hz)')
    set(get(gca,'YLabel'),'String','Magnitude (dB)')

    % Create sliders
    S = sliderFcn(S);
    S.update = updateFcn;
    S.update(S);
    guidata(S.fh, S);  % Store S struct in the figure
end