%% Parabolic Radon transform
if ~exist('B','var')
     [B,SuTraceHeaders,SuHeader] = ReadSu('parab.su');
end
% [B,SuTraceHeaders,SuHeader] = ReadSu('parab.su');
for k=1:length(SuTraceHeaders)
    h(k)=SuTraceHeaders(k).SourceX - SuTraceHeaders(k).GroupX;
end
t = [0:SuHeader.ns-1]'*SuHeader.dt*1e-6;
figure;imagesc(h,t,B);colormap(gray);
xlabel('offset [m]');ylabel('time [s]');
% choosing p[m-2], using matlab polyfit function.
p = 1e-6*[-3:.05:3];
B_tp = lpradon(B,t, h, p, 2, 1);
figure;imagesc(p,t,B_tp);colormap(gray);
% demultiple
% one way of doing this: in this data, it seems the events are in
% parabolic, and all the mutiples are in same curvature with different
% shifts. So cut off those points lay in a line. 
F_tp = ones(size(B_tp));
% ptemp = round(.20*size(F_tp,1)):round(.99*size(F_tp,1));
% ttemp = round(.53*size(F_tp,2)):round(.99*size(F_tp,2));
ptemp1 = round(.23*size(F_tp,1)):round(.27*size(F_tp,1));
ptemp2 = round(.35*size(F_tp,1)):round(.39*size(F_tp,1));
ptemp3 = round(.48*size(F_tp,1)):round(.52*size(F_tp,1));
ptemp4 = round(.60*size(F_tp,1)):round(.64*size(F_tp,1));
ptemp = [ptemp1 ptemp2 ptemp3 ptemp4];
ttemp = round(.01*size(F_tp,2)):round(.99*size(F_tp,2));
F_tp(ptemp,ttemp) = 0;
F_tp = my_smooth(F_tp,5,5);
figure;imagesc(p,t,F_tp);
B_tpF = B_tp.*F_tp;
figure;imagesc(B_tpF);colormap(gray)
B_itp = lpradon(B_tpF,t, h, p, 2, -1);
figure; imagesc(fliplr(B_itp));colormap(gray);


% % demultiple
% % Another way is using NMO. After NMO with velocity between primaries and
% % multiples, events of multiple stay undercorrected--large curveture, while 
% % events for primaries are overcorrected--small curveture.
% 
% % velocity analysis, pick velocity for multiples
% vmin = 1000;
% vmax = 4000;
% % pick = [7 8 10 13 16 34 59];
% % velocity = velana(B,SuTraceHeaders,SuHeader,0,vmin,vmax,pick);
% velocity = linspace(1050,1080,length(t));
% 
% % NMO
% sout = nmor(B,vec(t),h,vec(velocity),1); 
% figure;imagesc(sout);colormap(gray)
% % filter in radon domain
% B_tp = lpradon(sout,t, h, p, 2, 1);
% figure;imagesc(p,t,B_tp);colormap(gray);
% F_tp = ones(size(B_tp));
% ptemp = round(.01*size(F_tp,1)):round(.99*size(F_tp,1));
% ttemp = round(.49*size(F_tp,2)):round(.99*size(F_tp,2));
% F_tp(ptemp,ttemp) = 0;
% F_tp = my_smooth(F_tp,5,5);
% figure;imagesc(p,t,F_tp);
% figure;imagesc(B_tp.*F_tp);colormap(gray)
% B_itp = lpradon(B_tp.*F_tp,t, h, p, 2, -1);
% sout = nmor(B_itp,vec(t),h,vec(velocity),-1);
% figure;imagesc(sout);colormap(gray)
% % in this way, we only assume we should cut off points corresponding to
% % large curveture, so I do not know whether it is proper to defilter window
% % in a horizontal way.
