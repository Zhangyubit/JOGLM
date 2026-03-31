function NMSE = calculate_nmse(original_signal,recovered_signal)
    % Calculate NMSE
    numerator = norm(original_signal - recovered_signal, 2)^2;
    denominator = norm(original_signal, 2)^2;

    % Prevents the denominator from being zero
    if denominator == 0
        NMSE = NaN;
    else
        NMSE = numerator / denominator;
    end
end