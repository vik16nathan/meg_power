function produceOutputCSV(region_dict, output_file)
    region_names = keys(region_dict);
    values = cellfun(@(k) region_dict(k), region_names)';
    T = table(region_names(:), values, 'VariableNames', {'region', 'value'});
    writetable(T, output_file);
end
