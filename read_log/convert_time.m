function time_string = convert_time(MDHTime)

MDHTime_in_sec = MDHTime/1000;
hr = floor(MDHTime_in_sec/3600);
min = floor(mod(MDHTime_in_sec,3600)/60);
sec = mod(MDHTime_in_sec,60);

time_string = [num2str(hr), ':', num2str(min), ':', num2str(sec)];
end