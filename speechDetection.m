function [Limits] = speechDetection(x, stEnergy, stZcr, win, step, fs)

% Apply median filtering in the feature sequences (twice), using 5 windows:
E = medfilt1(stEnergy, 5); E = medfilt1(E, 5);
Z = medfilt1(stZcr, 5); C = medfilt1(Z, 5);
% Get the average values of the smoothed feature sequences:
E_mean = mean(E);
Z_mean = mean(Z);

Weight = 10; % used in the threshold estimation method


% Find energy threshold:
[HistE, X_E] = hist(E, round(length(E) / 10));  % histogram computation
[MaximaE] = findMaxima(HistE, 3); % find the local maxima of the histogram
if (size(MaximaE,2)>=2) % if at least two local maxima have been found in the histogram:
    T_E = (Weight*X_E(MaximaE(1,1))+X_E(MaximaE(1,2))) / (Weight+1); % ... then compute the threshold as the weighted average between the two first histogram's local maxima.
else
    T_E = E_mean / 2;
end

% Find zcr threshold:
[HistZ, X_Z] = hist(Z, round(length(Z) / 10));
[MaximaZ] = findMaxima(HistZ, 3);
if (size(MaximaZ,2)>=2)
    T_Z = (Weight*X_Z(MaximaZ(1,1))+X_Z(MaximaZ(1,2))) / (Weight+1);
else
    T_Z = Z_mean / 2;
end

% Thresholding:
Flags1 = (E>=T_E);
Flags3 = (Z>=T_Z);
%flags = Flags1 & Flags2 & Flags3;
flags = Flags1 & Flags3;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SPEECH SEGMENTS DETECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
posWin = 1; %First flag element
Limits = [];
while (posWin < length(flags)) % while there are windows to be processed:
	% initilize:
	%curX = [];	
	countTemp = 1;
	% while flags=1:
	while ((flags(posWin)==1) && (posWin < length(flags)))
		if (countTemp==1) % if this is the first of the current speech segment:
			%Limit1 = round((posWin-win)*step*fs)+1; % set start limit:
            Limit1 = round((posWin * step)+(win-step)-(win))+1;
			if (Limit1<1)	
                Limit1 = 1; 
            end        
		end	
		posWin = posWin + 1; 		% next window
		countTemp = countTemp + 1;	% next sample
	end

	if (countTemp>1) % if at least one segment has been found in the current loop:
		%Limit2 = round((count+WIN)*step*fs);			% set end counter
        Limit2 = round((posWin * step)+(win-step));
		if (Limit2>length(x))
            Limit2 = length(x);
        end
        
        Limits(end+1, 1) = Limit1;
        Limits(end,   2) = Limit2;
    end
	posWin = posWin + 1; % next window
end

%%%%%%%%%%%%%%%%%%%%%%%
% POST - PROCESS      %
%%%%%%%%%%%%%%%%%%%%%%%

% A. MERGE OVERLAPPING SEGMENTS:
RUN = 1;
while (RUN==1)
    RUN = 0;
    for (i=1:size(Limits,1)-1) % for each segment
        if (Limits(i,2) >= Limits(i+1,1) - round(fs/2))
            RUN = 1;
            Limits(i,2) = Limits(i+1,2);
            Limits(i+1,:) = [];
            break;
        end
    end
end

% B. Get final segments:
% segments = {};
% for (i=1:size(Limits,1))
%     segments{end+1} = x(Limits(i,1):Limits(i,2)); 
% end

