function [rclrs,changes] = COTReportTargetColorChange(filename, report)
%COTREPORTTARGETCOLORCHANGE Summary of this function goes here
%   Detailed explanation goes here


%% check and process inputs
narginchk(0, 2);

if nargin==0
    [file, path] = uigetfilex({...
        '*.hdf5;*.h5', 'HDF5 file (*.hdf5,*.h5)'}, ...
        'Select a File');
    if isequal(file, 0)
        disp('operation canceled');
        return;
    else
        filename = fullfile(path, file);
    end
end

if nargin == 1
    report = true;
end

trials_count = h5readatt(filename, '/', 'trials_count');
clrs = strings(1,trials_count);
for i=1:trials_count
    clrs(i) = string(h5readatt(filename, ['/trial_' pad(num2str(i), 4, 'left', '0')], 'target_instruction_color_rgb'));
end

if nargout>0
    rclrs = clrs;
end

if report
    disp(['''' filename ''' has ' num2str(trials_count) ' trials.']);
    color_changes = 0;
    changes =[];
    for i=2:trials_count
        if ~strcmp(clrs(i-1), clrs(i))
            color_changes = color_changes+1;
            disp(['color changed from ' char(clrs(i-1)) ' to ' char(clrs(i)) ' at trial ' num2str(i)]);
            changes(color_changes,1)=i;
        end
    end
    if color_changes == 0
        disp('no color changes detected');
    elseif color_changes == 1
        disp('1 color change detected');
    else
        disp([num2str(color_changes) ' color changes detected']);
    end
end
       

end
