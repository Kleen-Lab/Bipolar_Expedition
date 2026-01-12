
% Code to generate linear frequency table used for supplementary figure 1
% Supplementary figure 1 is generated using bipolar_supplementary1.py

% Takes in 'grid_trm_data.mat' file generated in
% bipolarexpedition_Linear_2025.m

% Output is a spreadsheet 'tableoutput.xlsx', displaying the mean sqrt
% transformed power at each bipolar distance, canonical frequency band, for
% each patient.

datadir = getenv("BIPOLAR_DATA");
data_results = getenv("RESULTS");
    
save_data = true;
save_data_path = data_results; 

load(fullfile(save_data_path,'grid_trm_data.mat'));


grids = trm_data.grids;

freq_bands = struct();
freq_bands.delta = [2,4];
freq_bands.theta = [4,8];
freq_bands.alpha = [8,13];
freq_bands.beta = [13,25];
freq_bands.low_gamma = [25,50];
freq_bands.high_gamma = [50,200];
band_names = fieldnames(freq_bands);

display_names = struct();
display_names.delta = 'delta';
display_names.theta = 'theta';
display_names.alpha = 'alpha';
display_names.beta = 'beta';
display_names.low_gamma = 'gamma';
display_names.high_gamma = 'high-gamma';


frx = trm_data.frequency;

all_data = [];
[n_distances, n_patients] = size(grids);

% Processing each patient
for p = 1:n_patients
    if isfield(trm_data, 'patients') && ~isempty(trm_data.patients)
        patient_name = trm_data.patients{p};
    else
        patient_name = sprintf('P%03d', p);
    end
    
    for d = 1:n_distances
        bpd_value = d - 1;  % saved as 0, 1, 2, 3, 4, 5

        data_matrix = grids{d, p};
       
        if ~isempty(data_matrix)
            % data_matrix is electrodes × frequencies
            % Need to average across electrodes
            mean_spectrum = nanmean(data_matrix, 1);  % 1 × 100 frequencies
            
            % no negative values b/c take square root transform in previous
            % linear analysis

            for b = 1:length(band_names)
                band = band_names{b};
                freq_range = freq_bands.(band);
                freq_idx = frx >= freq_range(1) & frx <= freq_range(2);
                band_power = nanmean(mean_spectrum(freq_idx));
                band_display = display_names.(band);
                all_data = [all_data; {patient_name, band_display, bpd_value, band_power}];
            end
        end
    end
end

% PUTTING INTO TABLE FORMAT
output_table = cell2table(all_data, ...
    'VariableNames', {'patient', 'freq', 'bpd', 'power'});

% Save to excel
if save_data
    writetable(output_table, fullfile(save_data_path,'tableoutput.xlsx'));
    disp(['Table saved to ' fullfile(save_data_path,'tableoutput.xlsx')]);
end

fprintf('\nUnique patients: %d\n', length(unique(output_table.patient)));
fprintf('Unique bpd values: ');
disp(unique(output_table.bpd)');
fprintf('\nFrequency bands used:\n');

for b = 1:length(band_names)
    fprintf('%s: %d-%d Hz\n', display_names.(band_names{b}), freq_bands.(band_names{b})(1), freq_bands.(band_names{b})(2));
end

% Can now plot using bipolar_supplementary1.py file with this
% tableoutput.xlsx file
