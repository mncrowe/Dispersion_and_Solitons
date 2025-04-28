function counter(t,p,units,txt)
% Gives current progress and time remaining
% t: current time, e.g. determined using toc function
% p: current progress (normalised to 1)
% units: 0 - determine automatically, 1 - secs, 2 - mins, 3 - hours, 4 - days, default is 0
% txt: optional additional text to display after remaining time

if nargin<3; units=0; end
if units>4; units=0; end
if nargin>3; txt = [', ' txt]; else; txt = ''; end

t_rem=t/p*(1-p);

if units==0
    units=1;
    if t_rem>60; units=2; end
    if t_rem>60*60; units=3; end
    if t_rem>60*60*24; units=4; end
end    

if units==1
	units_str='secs';
end
if units==2
	units_str='mins';
	t_rem=t_rem/60;
end
if units==3
	units_str='hours';
	t_rem=t_rem/(60*60);
end
if units==4
	units_str='days';
	t_rem=t_rem/(60*60*24);
end

t_rem=round(100*t_rem)/100;

disp([num2str(round(p*1000)/10) ' % complete, remaining time ~ ' num2str(t_rem) ' ' units_str txt]);

end
