% Copyright (c) 2014, Swallowing
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.
% 
% * Neither the name of AHI nor the names of its
%   contributors may be used to endorse or promote products derived from
%   this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% The data is based on LS-120 monitoring with the recorded id: 17812 
% The time frame is between 0:20 - 0:50 (30 min frame)
% Last Modified (record list below):
% ---- 11-Nov-2014 by Adam Lin
% ---- record required if edit by others
%

function output = showPlot(filedata, max_val, min_val, plot_number, times)
 if filedata ~= NaN
     org_value = csvread(filedata, 10, 5);
     z_size = [1:1:size(org_value)]';
     %----- First Plot -----%
     if strcmpi(plot_number,'original')
        plotMaker(z_size, org_value, 'k', 1/times);
        hold off
        output = [];
     end
     %----- Second Plot -----%
	 if strcmpi(plot_number,'diagnosis')
        
        if max_val > 0.1 || min_val > 0.1 % prevent for crash if max_val or min_val are not seted
            max_val = 0.01;
            min_val = 0.01;
        end
        
        x = z_size;
        ti_min = min(org_value) * max_val; % Depend on the value used (movement detection). 
        te_min = max(org_value) * min_val; % The method decode min and max of the activity of breath by 0.05
        time_frame = 125;                % 10 sec time frame for AHI --- TO DO: 70 is not right number ....
        countAHI = 0;
        new_end_loop = 1;
        max_AHI = 0;
        
        for jj = 1 : length(x)
            new_end_loop = new_end_loop + jj;
            
            if new_end_loop >= length(x)
                break;
            end
            
            for ii = new_end_loop : length(x)
                if org_value(ii) >= ti_min || org_value(ii) <= te_min
                    pstart = ii;
                    break;
                end
            end
            
            new_start_loop = pstart;
            
            for kk = new_start_loop: length(x)
                if (org_value(kk) < ti_min && org_value(kk) > 0) || org_value(kk) > te_min && org_value(kk)
                    pend = kk;
                    break;
                end
            end
            
            % Check the max of AHI continous time in second 
            gap_AHI = pend - pstart;
            if gap_AHI >= max_AHI
                max_AHI = gap_AHI;
            end
            
            if pstart ~= pend && pstart < pend
                % fprintf('pstart: %d --- pend: %d \n', pstart, pend);
                drawColorValue(org_value, pstart, pend, time_frame);
                % Print out AHI, AI and HI time per hour
                if pend - pstart >= time_frame
                    countAHI = countAHI + 1;
                end
            end
            new_end_loop = pend;
            
            if pstart > pend
                break;
            end
        end
        axis([0,numel(org_value), -2000, 2000]);
        set(gca,'ytick',[]);
        %plotMaker(z_size, org_value, 'k', 1);  
        hold off
        datacursormode on
        %set(handles.ahi, 'value', countAHI);
        %fprintf('AHI(RDI): %d T/H - AI: %d T/H \n', countAHI * 2, countAHI);
        %fprintf('Max AHI continous time: %0.2f Second \n', max_AHI/125);
        output = [countAHI*2 max_AHI*0.07]; % 30 min with 22490 frame 
     end
     %----- Third Plot -----%
     if strcmpi(plot_number, 'cancellation')
       ti_min = min(org_value) * 0.01; % Depend on the value used (movement detection). 
       te_min = max(org_value) * 0.01; % The method decode min and max of the activity of breath by 0.05
       
       compareBetweenValue(org_value, te_min, ti_min);
       plot(z_size, org_value, 'k');
       axis([0,numel(org_value), -2000, 2000]);
       set(gca, 'Color', 'None');
       set(gca,'ytick',[]);
       set(gca,'xtick',[]);
       
       hold off
       output = [];
     end
 end
end


function plotMaker(size, value, color, time_value)
    plot(size, value, color);
    set(gca, 'Color', 'None');
    hight_pick = -2000 * time_value;
    low_pick = 2000 * time_value;
    axis([0,numel(value), hight_pick, low_pick]);
    xlabel('Time (sec)');
    title('Respiratory Flow Signal');
    %grid on 
    %grid minor
end

function compareBetweenValue(value, upper_value, lower_value)
    % Compare the value infront and after to see the value of changes
    first_value = 0;
    second_value = 0;
    obtain_value = 0;
    draw_value = 0;
    y = [min(value) max(value)];
    for ii = 1:length(value)
        if first_value == 0 
           first_value = value(ii);
        else
           second_value = value(ii);
        end
        % others: 
        % detect big change s - f > upper_value || s - f < lower_value
        if second_value - first_value > upper_value || second_value - first_value < lower_value
           %fprintf('Show the number %d : %d Where ii is in %d\n', first_value, second_value, ii);
           if second_value ~= 0
               for pos = ii: ii + 1
                   x = [pos pos];
                   plot (x, y, 'm');
                   hold on
               end
           end
        else
            draw_value = draw_value + 1;
        end
    
        if second_value ~= 0
            first_value = second_value;
            second_value = 0;
        end
    end
end

function drawALine(value)
    y = [min(value) max(value)];
           for pos = ii: ii + 1
               x = [pos pos];
               plot (x, y, 'r');
           end
           hold on
end

function drawColorValue(org_value, fstart, fend, time_frame)
    % TO DO LIST:
    % The apneas (pauses in breathing) must last for at least 10 seconds!
    % The method below DID NOT check 10 seconds frame (10sec is 125(70)) in this
    % frame).
    y = [min(org_value) max(org_value)];
    if (fstart < fend)
        for pos = fstart: fend
            x = [pos pos];
            if fend - fstart >= time_frame
                if fend - fstart <= 250 
                    plot(x, y, 'b');
                else
                    plot(x, y, 'r');
                end
                hold on
            end
        end
    end
end
 