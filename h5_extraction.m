% Extract task information from .h5 files

% Version 1.0 by Zheng Liu at the University of Rochester on 2019-10-03
	% To view a .h5 file directly, use HDFView
	% Tested with Matlab R2018b, with imagesci toolbox

%% Set up
clc
clearvars
%data_folder = pwd; % assuming .h5 file saved in the same folder of script 
 data_folder = 'E:\Data\Q_ICMS';
filename = '\QCOT_20201015_ICMS_dim4.h5';
complete_filename = strcat(data_folder, filename);
if exist(complete_filename, 'file') ~= 2
	error('h5 file not exist in the current path')
end
file_info = h5info(complete_filename);

%% Extract trial info
% clear trial_data
for trial = 1:length(file_info.Groups)
    trial_info = h5info(complete_filename, file_info.Groups(trial).Name);
	
	% Extract information saved in trial Attributes
	for attribute_ind = 1:size(trial_info.Attributes,1)
        trial_data(trial).(genvarname(trial_info.Attributes(attribute_ind).Name)) = ...
            trial_info.Attributes(attribute_ind).Value;
	end
	
    % load trial information under root trial group
    for dataset_ind = 1:length(trial_info.Datasets)
        trial_data(trial).(genvarname(trial_info.Datasets(dataset_ind).Name)) = ...
            h5read(complete_filename, [file_info.Groups(trial).Name '/' trial_info.Datasets(dataset_ind).Name]);
        % generate variable name using trial_info.Datasets(dataset_ind).Name
        % trial_data stores all datasets from BCI controlled test
    end
    
    % read data in decoder_data
    trial_info_decoder_data = h5info(complete_filename, [file_info.Groups(trial).Name '/decoder_data']);
    for dataset_ind = 1:length(trial_info_decoder_data.Datasets)
        trial_data(trial).(genvarname(trial_info_decoder_data.Datasets(dataset_ind).Name)) = ...
            h5read(complete_filename, [file_info.Groups(trial).Name '/decoder_data/' trial_info_decoder_data.Datasets(dataset_ind).Name]);
    end
    
    % load trial information in /trial_xxxx/decoder_data
    for dataset_ind = 1:length(trial_info_decoder_data.Attributes)
        trial_data(trial).(genvarname(trial_info_decoder_data.Attributes(dataset_ind).Name)) = ...
            h5readatt(complete_filename, [file_info.Groups(trial).Name '/decoder_data'], ...
            trial_info_decoder_data.Attributes(dataset_ind).Name);
	end
	
	% Extract scales and offsets from decoder strings
% 	for unit_in_decoder = 1:min(size(trial_data(trial).channels_list.filter,1),4)
% 		scale_text = regexp(trial_data(trial).channels_list.filter{unit_in_decoder}, '* [0-9]+(.[0-9]+)?', 'match');
% 		if length(scale_text)~=1
% 			error('found no or multiple scale')
% 		end
% 		trial_data(trial).scale(unit_in_decoder,1) = str2double(scale_text{1}(3:end));
% 		
% 		offset_text = regexp(trial_data(trial).channels_list.filter{unit_in_decoder}, '- [0-9]+(.[0-9]+)?', 'match');
% 		if length(offset_text)~=1
% 			error('found no or multiple offset')
% 		end
% 		trial_data(trial).offset(unit_in_decoder,1) = str2num(offset_text{1}); % do not use isdouble -- not able to process negative numbers
% 	end
end









