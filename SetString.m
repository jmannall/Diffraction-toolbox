function SetString(object, name, variable)
    numObjects = length(object);
    for i = 1:numObjects
        set(object(i),'String',[name, ' ', num2str(i), ' (', num2str(round(variable(i))), ')']);
    end
end