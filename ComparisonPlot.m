% Create plots of BTM against variables

function ComparisonPlot(numxy,ixy,result, control, variable, flip)
    c1 = control.c1;
    c2 = control.c2;
    c1Name = control.c1Name;
    c2Name = control.c2Name;
    cIndex = control.cIndex;

    v = variable.v;
    vName = variable.vName;

    c1FileName = strrep(c1Name,' ', '_');
    c2FileName = strrep(c2Name,' ', '_');
    vFileName = lower(strrep(vName,' ', '_'));
    
    t2v = lower(vName);
    t2 = ['Frequency against magnitude for each ', t2v];
    t3 = [vName, ' against magnitude for set frequencies'];
    
    folder = ['figures/', vFileName];
    if not(isfolder(folder))
        mkdir(folder)
    end

    newcolors = [0 0.4470 0.7410
        0.8500 0.3250 0.0980
        0.9290 0.6940 0.1250
        0.4940 0.1840 0.5560
        0.4660 0.6740 0.1880
        0.3010 0.7450 0.9330
        0.6350 0.0780 0.1840
        1 0 0
        0 1 0
        0 0 1
        1 1 0
        1 0 1
        0 1 1
        0 0 0
        0.83 0.14 0.14
        1.00 0.54 0.00
        0.47 0.25 0.80
        0.25 0.80 0.54];
    
    for i = 1:numxy
        plotIndex = ixy == i;
        store = result(plotIndex);
        numVar = length(store);
        fvec = store.fvec;
        res = [store.tfmag];
        if flip
            vCurrent = v(((length(v) + 1) - numVar):end);
        else
            vCurrent = v(1:numVar);
        end
    
        %% Create figure
        figure('Units','normalized','Position',[0 0 1 1]);
        colororder(newcolors)
        sgtitle([c1Name, ': ', num2str(c1(cIndex(i,1))), ', ', c2Name, ': ', num2str(c2(cIndex(i,2)))]);
    
        % Subplot 1
        subplot(1,3,1)
        if numVar > 1
            ax = gca;
            mesh(ax,vCurrent,fvec,res);
            xlim([min(v) max(v)])
            ylim([20, 20000])
            zlim([-40 0])
            xlabel(vName)
            ylabel('Frequency (Hz)')
            zlabel('Magnitude (dB)')
            title('3D plot')
            set(gca, 'yscale', 'log')
        end
    
        % Subplot 2
        fvec = store(1).fvec;
        subplot(1,3,2)
        for j = 1:numVar
            semilogx(fvec, store(j).tfmag)
            hold on
        end
        hold off
        xlim([20, 20000])
        ylim([-40 0])
        xlabel('Frequency (Hz)')
        ylabel('Magnitude (dB)')
        title(t2);
        legend(num2str(vCurrent))
    
        % Subplot 3
        numFreq = 8;
        freq = logspace(log10(20),log10(20000),numFreq);
        tfnew = zeros(numFreq,numVar);
        for k = 1:numFreq
            indexStore = find(fvec > freq(k));
            for j = 1:numVar
                tfnew(k,j) = store(j).tfmag(indexStore(1));
            end
        end
        subplot(1,3,3)
        for k = 1:numFreq
            plot(vCurrent,tfnew(k,:))
            hold on
        end
        hold off
        title(t3)
        xlim([min(v) max(v)])
        ylim([-40 0])
        xlabel(vName)
        ylabel('Magnitude (dB)')
        legend(num2str(transpose(round(freq))))
        
        saveas(gcf,['figures/', vFileName, '/', c1FileName, '_', num2str(c1(cIndex(i,1))), '_', c2FileName, '_', num2str(c2(cIndex(i,2))),'.png'])
        disp('Save figure');

        close all
        % Plot change as wedge changes for frequencies.
    end
end
