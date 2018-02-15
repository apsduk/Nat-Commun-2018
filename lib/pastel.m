function cmap = pastel(N,ext)

%% pastel
%
% Makes matrix of default 'pastel-like' colors
%
% Color list: [ blue; orange; yellow; purple; green; light blue; maroon]
%
%   with an extension to grey and black by setting the ext option.
%
% Input/outputs
%
%       cmap = pastel(N,ext)
%
% Alexander Darlington 12/2016
%

%% set colormap
basecmap(1,:) = [0.0000	0.4431	0.7373]; % blue
basecmap(2,:) = [0.8510	0.3255	0.0980]; % orange
basecmap(3,:) = [0.9294	0.6941	0.1255]; % yellow
basecmap(4,:) = [0.4941	0.1843	0.5569]; % purple
basecmap(5,:) = [0.4667	0.6745	0.1882]; % green
basecmap(6,:) = [0.3020	0.7451	0.9333]; % light blue
basecmap(7,:) = [0.6353	0.0784	0.1843]; % maroon

% extend colormap with grey and black
if ext == 1
    basecmap(8,:) = [0.4902	0.4902	0.4902]; % grey
    basecmap(9,:) = [0.0000	0.0000	0.0000]; % black
end

%% extend colmap as needed
if N > length(basecmap(1,:))
    
    r = ceil(N./length(basecmap(1,:)));
    
    basecmap = repmat(basecmap,[r,1]);
    
end

%% make colormap
cmap = basecmap(1:N,:);

end