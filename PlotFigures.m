function PlotFigures(pso, wedge, index, variable, geometry)
    poles = [pso.p];
    zeros = [pso.z];
    gain = [pso.k];
    
    for i = 1:(size(index, 1)-1)
        for j = 1:(size(index, 2)-1)
            indexing = index(i,j).i;
            figure(4)
            plot(wedge(indexing), poles(1,indexing), '-O')
            hold on
            plot(wedge(indexing), poles(2,indexing), '-x')
            title("Zeros");
            xlabel(variable);
            
            figure(5)
            plot(wedge(indexing), zeros(1,indexing), '-O')
            hold on
            plot(wedge(indexing), zeros(2,indexing), '-x')
            title("Poles");
            xlabel(variable);
            
            figure(6)
            plot(wedge(indexing), gain(indexing), '-x')
            hold on
            title("Gain");
            xlabel(variable);
        end
    end
end