function lines = CreateAnimatedLinePlot()

    figure
    C = colororder;
    lines.iteration = animatedline('Color',C(1,:));
    lines.epoch = animatedline('Color',C(2,:));
    lines.test = animatedline('Color',C(3,:));
    ylim([0 inf])
    xlabel("Iteration")
    ylabel("Loss")
    grid on
end