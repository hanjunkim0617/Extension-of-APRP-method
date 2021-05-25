function [dswtoa dswsfc] = aprp_extension_hjkim(rsdsm1,rsusm1,rsutm1,rsdtm1,rsutcsm1,rsdscsm1,rsuscsm1,cltm1, ...
									            rsdsm2,rsusm2,rsutm2,rsdtm2,rsutcsm2,rsdscsm2,rsuscsm2,cltm2, ...
									            lon,lat, ...
									            flag_model,ar,rr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extended versions of APRP method of Taylor et al., 2007 (https://journals.ametsoc.org/view/journals/clim/20/11/jcli4143.1.xml)            %%
%% Contact : Hanjun Kim (hanjunkim0617@gmail.com)                                                                                            %%
%% All input data should be seasonal data                                                                                                    %%
%% In this version of APRP method, surface flux components are added and one-layer radiation model could be customized using the parameters  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -----------------------------------   Parameter to customize the one-layer radiation model   -----------------------------------%
%                                                                                                                                 %
% flag_model = 1 : multiple absorption model (Addition of multiple absorption to Taylor's model)                                  %
% flag_model = 2 : multiple absorption model (reflection and absorption at the same time, idea from Donohoe and Battisti, 2011)   %
%                                                                                                                                 %
% ar : ratio of upward beam absorptivity to downward beam absorptivity                                                            %
% rr : ratio of upward beam reflectivity to downward beam reflecitivity                                                           %
%                                                                                                                                 %
% for original APRP method, use flag_model = 1 / ar = 0 / rr = 1                                                                  %
% --------------------------------------------------------------------------------------------------------------------------------%

%% radiative parameters %%
% mu : atmospheric absorption coefficient
% ga : atmospheric scattering coefficient
% a : surface albedo

% consider the unit of clt
if nanmean(nanmean(nanmean(cltm1)))>1; cltm1 = cltm1/100; end;
if nanmean(nanmean(nanmean(cltm2)))>1; cltm2 = cltm2/100; end;

% calculate overcast fluxes for rsds, rsus, rsut
rsdsocm1 = (rsdsm1 - (1-cltm1).*rsdscsm1)./cltm1; rsdsocm2 = (rsdsm2 - (1-cltm2).*rsdscsm2)./cltm2;
rsusocm1 = (rsusm1 - (1-cltm1).*rsuscsm1)./cltm1; rsusocm2 = (rsusm2 - (1-cltm2).*rsuscsm2)./cltm2;
rsutocm1 = (rsutm1 - (1-cltm1).*rsutcsm1)./cltm1; rsutocm2 = (rsutm2 - (1-cltm2).*rsutcsm2)./cltm2;

% treatment for overcast-sky calculations
% erronusly large overcast-sky value for too small clt cases (dividing by small value --> too big values)
% cloud fraction less than 3% is considered as clear-sky
rsdsocm1(find(cltm1<0.03)) = rsdscsm1(find(cltm1<0.03)); rsdsocm2(find(cltm2<0.03)) = rsdscsm2(find(cltm2<0.03)); 
rsusocm1(find(cltm1<0.03)) = rsuscsm1(find(cltm1<0.03)); rsusocm2(find(cltm2<0.03)) = rsuscsm2(find(cltm2<0.03)); 
rsutocm1(find(cltm1<0.03)) = rsutcsm1(find(cltm1<0.03)); rsutocm2(find(cltm2<0.03)) = rsutcsm2(find(cltm2<0.03)); 

% treatment for albedo calculations
% NaN value for no-insolation cases (dividing by small value --> too big values)
smask1 = find(rsdtm1<3); smask2 = find(rsdtm2<3);
a_rsdtm1=rsdtm1; a_rsdtm2=rsdtm2;
a_rsdscsm1=rsdscsm1; a_rsdscsm2=rsdscsm2;
a_rsdsocm1=rsdsocm1; a_rsdsocm2=rsdsocm2;
a_rsdtm1(smask1) = NaN; a_rsdtm2(smask2) = NaN;
a_rsdscsm1(smask1) = NaN; a_rsdscsm2(smask2) = NaN;
a_rsdsocm1(smask1) = NaN; a_rsdsocm2(smask2) = NaN;

% calculate clear- and overcast- sky parameters from one-layer radiation model
aclr1=rsuscsm1./a_rsdscsm1; aoc1=rsusocm1./a_rsdsocm1;
Qclr1=a_rsdscsm1./a_rsdtm1; Qoc1=a_rsdsocm1./a_rsdtm1;
aclr2=rsuscsm2./a_rsdscsm2; aoc2=rsusocm2./a_rsdsocm2;
Qclr2=a_rsdscsm2./a_rsdtm2; Qoc2=a_rsdsocm2./a_rsdtm2;

if flag_model == 1
	muclr1=(1-rsutcsm1./a_rsdtm1-Qclr1.*(1-aclr1))./(1+ar.*Qclr1.*aclr1);
	muoc1=(1-rsutocm1./a_rsdtm1-Qoc1.*(1-aoc1))./(1+ar.*Qoc1.*aoc1);
	gaclr1=(1-muclr1-Qclr1)./(1-muclr1-(1-ar.*muclr1).*rr.*Qclr1.*aclr1);
	gaoc1=(1-muoc1-Qoc1)./(1-muoc1-(1-ar.*muoc1).*rr.*Qoc1.*aoc1);
	muclr2=(1-rsutcsm2./a_rsdtm2-Qclr2.*(1-aclr2))./(1+ar.*Qclr2.*aclr2);
	muoc2=(1-rsutocm2./a_rsdtm2-Qoc2.*(1-aoc2))./(1+ar.*Qoc2.*aoc2);
	gaclr2=(1-muclr2-Qclr2)./(1-muclr2-(1-ar.*muclr2).*rr.*Qclr2.*aclr2);
	gaoc2=(1-muoc2-Qoc2)./(1-muoc2-(1-ar.*muoc2).*rr.*Qoc2.*aoc2);
elseif flag_model == 2
	muclr1=(1-rsutcsm1./a_rsdtm1-Qclr1.*(1-aclr1))./(1+ar.*Qclr1.*aclr1);
	muoc1=(1-rsutocm1./a_rsdtm1-Qoc1.*(1-aoc1))./(1+ar.*Qoc1.*aoc1);
	gaclr1=(1-Qclr1-muclr1)./(1-rr.*Qclr1.*aclr1);
	gaoc1=(1-Qoc1-muoc1)./(1-rr.*Qoc1.*aoc1);
	muclr2=(1-rsutcsm2./a_rsdtm2-Qclr2.*(1-aclr2))./(1+ar.*Qclr2.*aclr2);
	muoc2=(1-rsutocm2./a_rsdtm2-Qoc2.*(1-aoc2))./(1+ar.*Qoc2.*aoc2);
	gaclr2=(1-Qclr2-muclr2)./(1-rr.*Qclr2.*aclr2);
	gaoc2=(1-Qoc2-muoc2)./(1-rr.*Qoc2.*aoc2);
end

% calculating cloud parameters through clear sky and overcast parameters
mucld1=(muoc1-1)./(1-muclr1)+1; gacld1=(gaoc1-1)./(1-gaclr1)+1;
mucld2=(muoc2-1)./(1-muclr2)+1; gacld2=(gaoc2-1)./(1-gaclr2)+1;

% decomposition of swtoa and swsfc change into several components (downward positive)
% dswtoa = ds*(1-A) - S*dA , dswsfc = ds*Qs*(1-a) - S*Qs*da + S*dQs*(1-a) 
% insolation
S = (rsdtm1 + rsdtm2)/2; 
dS = rsdtm2 - rsdtm1;

% surface albedo
a = ((rsusm1./rsdsm1) + (rsusm2./rsdsm2))/2; 
da = (rsusm2./rsdsm2) - (rsusm1./rsdsm1);

% calculate A (albedo) using the one-layer raidation model (use the attached functions)
A = ((rsutm1./rsdtm1) + (rsutm2./rsdtm2))/2;
A_1=albedo(cltm1,aclr1,aoc1,muclr1,mucld1,gaclr1,gacld1,flag_model,ar,rr); A_2=albedo(cltm2,aclr2,aoc2,muclr2,mucld2,gaclr2,gacld2,flag_model,ar,rr);
dA.dc= 0.5*(albedo(cltm2,aclr1,aoc1,muclr1,mucld1,gaclr1,gacld1,flag_model,ar,rr)-A_1)+0.5*(A_2-albedo(cltm1,aclr2,aoc2,muclr2,mucld2,gaclr2,gacld2,flag_model,ar,rr)) ;
dA.aclr= 0.5*(albedo(cltm1,aclr2,aoc1,muclr1,mucld1,gaclr1,gacld1,flag_model,ar,rr)-A_1)+0.5*(A_2-albedo(cltm2,aclr1,aoc2,muclr2,mucld2,gaclr2,gacld2,flag_model,ar,rr)) ;
dA.aoc= 0.5*(albedo(cltm1,aclr1,aoc2,muclr1,mucld1,gaclr1,gacld1,flag_model,ar,rr)-A_1)+0.5*(A_2-albedo(cltm2,aclr2,aoc1,muclr2,mucld2,gaclr2,gacld2,flag_model,ar,rr)) ;
dA.mucld= 0.5*(albedo(cltm1,aclr1,aoc1,muclr1,mucld2,gaclr1,gacld1,flag_model,ar,rr)-A_1)+0.5*(A_2-albedo(cltm2,aclr2,aoc2,muclr2,mucld1,gaclr2,gacld2,flag_model,ar,rr)) ;
dA.muclr= 0.5*(albedo(cltm1,aclr1,aoc1,muclr2,mucld1,gaclr1,gacld1,flag_model,ar,rr)-A_1)+0.5*(A_2-albedo(cltm2,aclr2,aoc2,muclr1,mucld2,gaclr2,gacld2,flag_model,ar,rr)) ;
dA.gacld= 0.5*(albedo(cltm1,aclr1,aoc1,muclr1,mucld1,gaclr1,gacld2,flag_model,ar,rr)-A_1)+0.5*(A_2-albedo(cltm2,aclr2,aoc2,muclr2,mucld2,gaclr2,gacld1,flag_model,ar,rr)) ;
dA.gaclr= 0.5*(albedo(cltm1,aclr1,aoc1,muclr1,mucld1,gaclr2,gacld1,flag_model,ar,rr)-A_1)+0.5*(A_2-albedo(cltm2,aclr2,aoc2,muclr2,mucld2,gaclr1,gacld2,flag_model,ar,rr)) ;
dA.sum = dA.dc + dA.aclr + dA.aoc + dA.mucld + dA.muclr +dA.gacld + dA.gaclr;
dA.total=  A_2 - A_1 ;
dA.model = (rsutm2./rsdtm2) - (rsutm1./rsdtm1);

% calculate Qs (incident sw) using one-layer raidation model (use the attached functions)
Qs = ((rsdsm1./rsdtm1) + (rsdsm2./rsdtm2))/2;
Qs_1=incident_sw_ratio(cltm1,aclr1,aoc1,muclr1,mucld1,gaclr1,gacld1,flag_model,ar,rr); Qs_2=incident_sw_ratio(cltm2,aclr2,aoc2,muclr2,mucld2,gaclr2,gacld2,flag_model,ar,rr);
dQs.dc= 0.5*(incident_sw_ratio(cltm2,aclr1,aoc1,muclr1,mucld1,gaclr1,gacld1,flag_model,ar,rr)-Qs_1)+0.5*(Qs_2-incident_sw_ratio(cltm1,aclr2,aoc2,muclr2,mucld2,gaclr2,gacld2,flag_model,ar,rr));
dQs.aclr= 0.5*(incident_sw_ratio(cltm1,aclr2,aoc1,muclr1,mucld1,gaclr1,gacld1,flag_model,ar,rr)-Qs_1)+0.5*(Qs_2-incident_sw_ratio(cltm2,aclr1,aoc2,muclr2,mucld2,gaclr2,gacld2,flag_model,ar,rr));
dQs.aoc= 0.5*(incident_sw_ratio(cltm1,aclr1,aoc2,muclr1,mucld1,gaclr1,gacld1,flag_model,ar,rr)-Qs_1)+0.5*(Qs_2-incident_sw_ratio(cltm2,aclr2,aoc1,muclr2,mucld2,gaclr2,gacld2,flag_model,ar,rr));
dQs.mucld= 0.5*(incident_sw_ratio(cltm1,aclr1,aoc1,muclr1,mucld2,gaclr1,gacld1,flag_model,ar,rr)-Qs_1)+0.5*(Qs_2-incident_sw_ratio(cltm2,aclr2,aoc2,muclr2,mucld1,gaclr2,gacld2,flag_model,ar,rr));
dQs.muclr= 0.5*(incident_sw_ratio(cltm1,aclr1,aoc1,muclr2,mucld1,gaclr1,gacld1,flag_model,ar,rr)-Qs_1)+0.5*(Qs_2-incident_sw_ratio(cltm2,aclr2,aoc2,muclr1,mucld2,gaclr2,gacld2,flag_model,ar,rr));
dQs.gacld= 0.5*(incident_sw_ratio(cltm1,aclr1,aoc1,muclr1,mucld1,gaclr1,gacld2,flag_model,ar,rr)-Qs_1)+0.5*(Qs_2-incident_sw_ratio(cltm2,aclr2,aoc2,muclr2,mucld2,gaclr2,gacld1,flag_model,ar,rr));
dQs.gaclr= 0.5*(incident_sw_ratio(cltm1,aclr1,aoc1,muclr1,mucld1,gaclr2,gacld1,flag_model,ar,rr)-Qs_1)+0.5*(Qs_2-incident_sw_ratio(cltm2,aclr2,aoc2,muclr2,mucld2,gaclr1,gacld2,flag_model,ar,rr));
dQs.sum = dQs.dc + dQs.aclr + dQs.aoc + dQs.mucld + dQs.muclr + dQs.gacld + dQs.gaclr;
dQs.total=  Qs_2 - Qs_1;
dQs.model =  (rsdsm2./rsdtm2) - (rsdsm1./rsdtm1);

% calculate sw change directly from model output
dswtoa.model = (rsdtm2-rsutm2) - (rsdtm1-rsutm1);
dswsfc.model = (rsdsm2-rsusm2) - (rsdsm1-rsusm1);

% assign each component (downward positive for both TOA and SFC)
% note that region with insufficient insolation has no change (0 values) (all NaN value become 0 here)
% TOA --> sum = ds + cld + clr + alb
dswtoa.ds = dS.*(1-A); dswtoa.ds(smask1) = 0; dswtoa.ds(smask2) = 0;
dswtoa.A_dc = -S.*dA.dc; dswtoa.A_dc(smask1) = 0; dswtoa.A_dc(smask2) = 0;
dswtoa.A_aclr = -S.*dA.aclr; dswtoa.A_aclr(smask1) = 0; dswtoa.A_aclr(smask2) = 0; 
dswtoa.A_aoc = -S.*dA.aoc; dswtoa.A_aoc(smask1) = 0; dswtoa.A_aoc(smask2) = 0;
dswtoa.A_mucld = -S.*dA.mucld; dswtoa.A_mucld(smask1) = 0; dswtoa.A_mucld(smask2) = 0;
dswtoa.A_muclr = -S.*dA.muclr; dswtoa.A_muclr(smask1) = 0; dswtoa.A_muclr(smask2) = 0;
dswtoa.A_gacld = -S.*dA.gacld; dswtoa.A_gacld(smask1) = 0; dswtoa.A_gacld(smask2) = 0;
dswtoa.A_gaclr = -S.*dA.gaclr; dswtoa.A_gaclr(smask1) = 0; dswtoa.A_gaclr(smask2) = 0;
dswtoa.cld = dswtoa.A_dc + dswtoa.A_mucld + dswtoa.A_gacld;
dswtoa.ncld = dswtoa.A_muclr + dswtoa.A_gaclr;
dswtoa.alb = dswtoa.A_aclr + dswtoa.A_aoc;
dswtoa.sum = dswtoa.ds + dswtoa.cld + dswtoa.ncld + dswtoa.alb;
% if the change is larger than 100, we will consider it as abnormal value --> (it occur hardly)
dswtoa.ds(dswtoa.ds>100|dswtoa.ds<-100)=NaN;
dswtoa.A_dc(dswtoa.A_dc>100|dswtoa.A_dc<-100)=NaN;
dswtoa.A_aclr(dswtoa.A_aclr>100|dswtoa.A_aclr<-100)=NaN;
dswtoa.A_aoc(dswtoa.A_aoc>100|dswtoa.A_aoc<-100)=NaN;
dswtoa.A_mucld(dswtoa.A_mucld>100|dswtoa.A_mucld<-100)=NaN;
dswtoa.A_muclr(dswtoa.A_muclr>100|dswtoa.A_muclr<-100)=NaN;
dswtoa.A_gacld(dswtoa.A_gacld>100|dswtoa.A_gacld<-100)=NaN;
dswtoa.A_gacld(dswtoa.A_gacld>100|dswtoa.A_gacld<-100)=NaN;
dswtoa.cld(dswtoa.cld>100|dswtoa.cld<-100)=NaN;
dswtoa.ncld(dswtoa.ncld>100|dswtoa.ncld<-100)=NaN;
dswtoa.alb(dswtoa.alb>100|dswtoa.alb<-100)=NaN;
dswtoa.sum(dswtoa.sum>100|dswtoa.sum<-100)=NaN;


% SFC --> sum = ds + cld + clr + alb
dswsfc.ds = dS.*Qs.*(1-a); dswsfc.ds(smask1) = 0; dswsfc.ds(smask2) = 0;
dswsfc.da = -S.*Qs.*da; dswsfc.da(smask1) = 0; dswsfc.da(smask2) = 0;
dswsfc.Qs_dc = S.*dQs.dc.*(1-a); dswsfc.Qs_dc(smask1) = 0; dswsfc.Qs_dc(smask2) = 0;
dswsfc.Qs_aclr = S.*dQs.aclr.*(1-a); dswsfc.Qs_aclr(smask1) = 0; dswsfc.Qs_aclr(smask2) = 0;
dswsfc.Qs_aoc = S.*dQs.aoc.*(1-a); dswsfc.Qs_aoc(smask1) = 0; dswsfc.Qs_aoc(smask2) = 0;
dswsfc.Qs_mucld = S.*dQs.mucld.*(1-a); dswsfc.Qs_mucld(smask1) = 0; dswsfc.Qs_mucld(smask2) = 0;
dswsfc.Qs_muclr = S.*dQs.muclr.*(1-a); dswsfc.Qs_muclr(smask1) = 0; dswsfc.Qs_muclr(smask2) = 0;
dswsfc.Qs_gacld = S.*dQs.gacld.*(1-a); dswsfc.Qs_gacld(smask1) = 0; dswsfc.Qs_gacld(smask2) = 0;
dswsfc.Qs_gaclr = S.*dQs.gaclr.*(1-a); dswsfc.Qs_gaclr(smask1) = 0; dswsfc.Qs_gaclr(smask2) = 0;
dswsfc.cld = dswsfc.Qs_dc + dswsfc.Qs_mucld + dswsfc.Qs_gacld;
dswsfc.ncld = dswsfc.Qs_muclr + dswsfc.Qs_gaclr;
dswsfc.alb = dswsfc.da + dswsfc.Qs_aclr + dswsfc.Qs_aoc;
dswsfc.sum = dswsfc.ds + dswsfc.cld + dswsfc.ncld + dswsfc.alb;
% if the change is larger than 100, we will consider it as abnormal value --> (it occur hardly)
dswsfc.ds(dswsfc.ds>100|dswsfc.ds<-100)=NaN;
dswsfc.da(dswsfc.da>100|dswsfc.da<-100)=NaN;
dswsfc.Qs_dc(dswsfc.Qs_dc>100|dswsfc.Qs_dc<-100)=NaN;
dswsfc.Qs_aclr(dswsfc.Qs_aclr>100|dswsfc.Qs_aclr<-100)=NaN;
dswsfc.Qs_aoc(dswsfc.Qs_aoc>100|dswsfc.Qs_aoc<-100)=NaN;
dswsfc.Qs_mucld(dswsfc.Qs_mucld>100|dswsfc.Qs_mucld<-100)=NaN;
dswsfc.Qs_muclr(dswsfc.Qs_muclr>100|dswsfc.Qs_muclr<-100)=NaN;
dswsfc.Qs_gacld(dswsfc.Qs_gacld>100|dswsfc.Qs_gacld<-100)=NaN;
dswsfc.Qs_gaclr(dswsfc.Qs_gaclr>100|dswsfc.Qs_gaclr<-100)=NaN;
dswsfc.cld(dswsfc.cld>100|dswsfc.cld<-100)=NaN;
dswsfc.ncld(dswsfc.ncld>100|dswsfc.ncld<-100)=NaN;
dswsfc.alb(dswsfc.alb>100|dswsfc.alb<-100)=NaN;
dswsfc.sum(dswsfc.sum>100|dswsfc.sum<-100)=NaN;


%% get global-mean value to check whether the global-mean model SW changes are corresponding to APRP sum or not
% area weighting 
if size(lat,2) > 10; lat = lat'; end;
if size(lon,2) > 10; lon = lon'; end;
a = 6371e3; nlat = length(lat); nlon = length(lon'); 
dx = a * repmat(cosd(lat),[1 nlon]) .* deg2rad(repmat(gradient(lon'),[nlat 1]));
dy = a * deg2rad(repmat(gradient(lat),[1 nlon]));
aw = dx.*dy; % [m^2] area weight
awtr = repmat(reshape(aw,[1 nlat nlon]),[12 1 1]);

dswtoa.sum(isnan(dswtoa.sum)) = 0; awt = awtr; awt(isnan(dswtoa.sum))=0;
dswtoa.model(isnan(dswtoa.model)) = 0; awt = awtr; awt(isnan(dswtoa.model))=0;
dswsfc.sum(isnan(dswsfc.sum)) = 0; awt = awtr; awt(isnan(dswsfc.sum))=0;
dswsfc.model(isnan(dswsfc.model)) = 0; awt = awtr; awt(isnan(dswsfc.model))=0;

% global mean value for validations
gmtoasum = mean(sum(sum(dswtoa.sum.*awt,2),3)./sum(sum(awt,2),3),1);
gmtoamodel = mean(sum(sum(dswtoa.model.*awt,2),3)./sum(sum(awt,2),3),1);
gmsfcsum = mean(sum(sum(dswsfc.sum.*awt,2),3)./sum(sum(awt,2),3),1);
gmsfcmodel = mean(sum(sum(dswsfc.model.*awt,2),3)./sum(sum(awt,2),3),1);
display([ 'TOA_{sum,model} : ' num2str(gmtoasum,'%10.3f') ',' num2str(gmtoamodel,'%10.3f') ' / SFC_{sum,model} : ' num2str(gmsfcsum,'%10.3f') ',' num2str(gmsfcmodel,'%10.3f') ])

end


% % validation figure (optional)
% load coast_hjkim.mat;

% b=10; nc=40;
% figure('Position',[10 10 1800 900]);
% subplot(2,3,1);
% contourf(lon,lat,squeeze(mean(dswtoa.sum,1)),linspace(-3*b,3*b,3*nc+1),'linestyle','none'); hold on;
% plot(lonmap+360,latmap,'k');
% plot(lonmap,latmap,'k');
% xlim([2 358]);
% colormap(french(nc)); colorbar;
% caxis([-b b]);
% title('\DeltaSW_{TOA} : APRP');

% subplot(2,3,2);
% contourf(lon,lat,squeeze(mean(dswtoa.cld,1)),linspace(-3*b,3*b,3*nc+1),'linestyle','none'); hold on;
% plot(lonmap+360,latmap,'k');
% plot(lonmap,latmap,'k');
% xlim([2 358]);
% colormap(french(nc)); colorbar;
% caxis([-b b]);
% title('\DeltaSW_{TOA} : cloud');

% subplot(2,3,3);
% contourf(lon,lat,squeeze(mean(dswtoa.ncld,1)),linspace(-3*b,3*b,3*nc+1),'linestyle','none'); hold on;
% plot(lonmap+360,latmap,'k');
% plot(lonmap,latmap,'k');
% xlim([2 358]);
% colormap(french(nc)); colorbar;
% caxis([-b b]);
% title('\DeltaSW_{TOA} : non-cloud');

% subplot(2,3,4);
% contourf(lon,lat,squeeze(mean(dswtoa.alb,1)),linspace(-3*b,3*b,3*nc+1),'linestyle','none'); hold on;
% plot(lonmap+360,latmap,'k');
% plot(lonmap,latmap,'k');
% xlim([2 358]);
% colormap(french(nc)); colorbar;
% caxis([-b b]);
% title('\DeltaSW_{TOA} : albedo');

% subplot(2,3,5);
% contourf(lon,lat,squeeze(mean(dswtoa.sum-dswtoa.model,1)),linspace(-3*b,3*b,3*nc+1),'linestyle','none'); hold on;
% plot(lonmap+360,latmap,'k');
% plot(lonmap,latmap,'k');
% xlim([2 358]);
% colormap(french(nc)); colorbar;
% caxis([-b b]);
% title('\DeltaSW_{TOA} : APRP-model');

% subplot(2,3,6);
% contourf(lon,lat,squeeze(mean(dswsfc.sum-dswsfc.model,1)),linspace(-3*b,3*b,3*nc+1),'linestyle','none'); hold on;
% plot(lonmap+360,latmap,'k');
% plot(lonmap,latmap,'k');
% xlim([2 358]);
% colormap(french(nc)); colorbar;
% caxis([-b b]);
% title('\DeltaSW_{SFC} : APRP-model');


