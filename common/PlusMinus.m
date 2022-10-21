%% Flip the sign of the input based on the bool 'plus'

function output = PlusMinus(input, plus)
    if plus
        output = input;
    else
        output = -input;
    end
end