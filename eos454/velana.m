function [vel] = velana(Data,SuTraceHeaders,SuHeader,CDPnumber,vmin,vmax,pick)
% this function is aiming to got velocity vector for Data volumn for certain 
% cdp point with manual pick(cdp adapative)

% name : file name of .su file
% CDPnumber : intended CDP nummber to carry out vel analysis
% vmin : minum velocity in velocity semblance
% vmax : maxmum velocity in velocity semblance
% pick : an vector given the manual picked point on the velocity semblance
%        eg: for cdp 341, pick = [6,12,49]
%%name: Lina Maio
%%Student No: 74721119

for k=1:length(SuTraceHeaders)
    SuTraceHeaders(k).offset=SuTraceHeaders(k).SourceX - SuTraceHeaders(k).GroupX;
end
% form a midpoint gather(using CDP by assuming horizontal layers)
j = 1;
for i = 1:size(Data,2)  
    if SuTraceHeaders(i).cdp == CDPnumber
        CDPGA(:,j) = Data(:,i); 
        CDPGAHeader(:,j) = SuTraceHeaders(i);
        j = j + 1;
    end
end
if CDPGAHeader(1).offset > CDPGAHeader(end).offset
    CDPGA = fliplr(CDPGA);
    CDPGAHeader = fliplr(CDPGAHeader);
end
% extract the time, offset needed for 'nmor' function
%t0 = SuHeader.dtOrig;
t0 = 0;
dt = SuHeader.dt;
nt = size(CDPGA,1);
t = t0+[1:1:nt].*dt; % time coordinate
t = t./10^6; % time coordinate in second
x = [CDPGAHeader.offset]; % offset coordinate
% form velocity span semblance
v_low = vmin;
v_high = vmax;
n = 50;%sampling in velocity
for i = 1:n
    v = v_low+(i-1)*((v_high-v_low)/n);
    vnmo = v.*ones(size(t));
    s = nmor(CDPGA,vec(t),x,vec(vnmo),1);
    S(:,i) = sum(s.^2,2); 
end
figure;imagesc(S);hold on;
% set(gca,'XTickLabel',[linspace(v_low,v_high,10)]);
[maxm,a] = max(abs(S),[],2);
index = 1:10:size(a,1); % sample in the axe of t0
asample = a(index);
P = polyfit(index,asample',3);
% velocity picking
hold on;plot(asample,index,'*-r');
pick = [6 12 49];% default picking points
fprintf('please pick the velocities, eg: by default pick=[6 12 49], for CDP:%3.0f\n',CDPnumber);
keyboard;
index_pick = index(pick);
% handle the shallow and deep part
if 1 < pick(1)
    asample(1) = asample(pick(1));
    pick = [1 pick];
end

if pick(end) < length(index)
    asample(end) = asample(pick(end));
    pick = [pick length(index)];
end
index_pick = index(pick);
% piecewise spline interpolation
yy = interp1(index_pick,asample(pick),1:1:size(CDPGA,1),'pchip');
% piecewise linear interpolation
yy = interp1(index_pick,asample(pick),1:1:size(CDPGA,1),'linear');
% compute velocity from the indices information
v = v_low+(yy-1).*((v_high-v_low)/n);
vel = v;
% figure;imagesc(S);hold on;
hold on; plot(yy,[1:1:size(CDPGA,1)],'-w');title('veloctity analysis curve')













