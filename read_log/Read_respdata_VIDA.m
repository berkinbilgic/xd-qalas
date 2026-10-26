% Function Read_respdata_VIDA
% Input: Path to .rep file
% Output: Structure file containing all time and resp signals
% run in command line -> Resp = Read_respdata_VIDA ('ROCARD020.resp')


function Resp = Read_respdata_VIDA(pulsefile)
% [filepath,name,ext] = fileparts(pulsefile);
[fid,emsg]= fopen(pulsefile);
if fid ==-1;
    error(emsg);
end

tline = fgetl(fid);
C = [];
while ischar(tline)
    C{end+1,1}= tline;
    tline = fgetl(fid);
end

% acquire start and stop time
Log1 = strsplit(cell2mat(C(end-4,1)));
Log2 = strsplit(cell2mat(C(end-3,1)));
Log3 = strsplit(cell2mat(C(end-2,1)));
Log4 = strsplit(cell2mat(C(end-1,1)));

StartMDHTime = str2double(Log1(2));
StopMDHTime = str2double(Log2(2));
StartMPCUTime = str2double(Log3(2));
StopMPCUTime = str2double(Log4(2));

fclose(fid);
temp = C(1,1);
temp1 = strsplit(cell2mat(temp));
temp2 = temp1(1,7:end);
temp3 = cellfun(@str2double,temp2);

Dimension = size(temp3);

j = 1;
k = 1
trigger_signal = 0;

% avoid 5003: terminate point of the actual data
for (i = 3:Dimension(2)-1)
    % trigger not working here
    % Resp.trigger(j) = trigger_signal;
    Resp.time{j} = convert_time(StartMDHTime + (j-1)*2.5);
    
    if (temp3(i)) == 5000;
        trigger_signal(k) = j;
        k = k + 1;
        %Resp.trigger(j) = trigger_signal;
    end

    if (temp3(i)> 0 && temp3(i)< 4097)
        Resp.bellows(j) = temp3(i);
    elseif (temp3(i)> 33554431 && temp3(i)< 50331649)
        Resp.Rsr_raw_HeadFirst(j) = temp3(i);
    elseif (temp3(i)> 65535 && temp3(i)< 131073)
        Resp.Rsr_preprcd_HeadFirst(j) = temp3(i);
    elseif (temp3(i)> 67108863 && temp3(i)< 83886081)
        Resp.Rsr_raw_FeetFirst(j) = temp3(i);
    elseif (temp3(i)> 131071 && temp3(i)< 196609)
        Resp.Rsr_preprcd_FeetFirst(j) = temp3(i);
        j=j+1;
    end
end

% Stop Time and Duration based on resp data point and 400Hz sampling rate
% based on MDH Start time

clearvars -except Resp

%% TO do check which signals are actively logged and delete the rest
%% Save to file in same filepath
end