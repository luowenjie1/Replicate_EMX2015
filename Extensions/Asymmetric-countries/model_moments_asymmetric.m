fprintf('...calculating model moments... \n');

%%%%% A. Domestic market share moments   

omegad = omegaH;
omegad = omegaH./repmat(sum(omegaH),N,1);
omegad = omegad(:,sum(omegaH)>0); 

[NN,SS] = size(omegad);

omegads = omegaFs;
omegads = omegaFs./repmat(sum(omegaFs),N,1);
omegads = omegads(:,sum(omegaFs)>0); 

[NNs,~] = size(omegads);

omegad(omegad==0)   = NaN; %% mimic data structure (use of NaN instead of zeros)
omegads(omegads==0) = NaN;

sharedom  = omegad;       %% domestic market shares
sharedoms = omegads;

hh            = nansum(sharedom.^2);
highest_share = nanmax(sharedom);

share         = sharedom(:); share = share(share>0);

hhs           = nansum(sharedoms.^2);
highest_shares = nanmax(sharedoms);

shares         = sharedoms(:); shares = shares(shares>0);


for s = 1:SS,
    
    share90(s)          = prctile(sharedom(:,s),90);
    share75(s)          = prctile(sharedom(:,s),75);
    share50(s)          = prctile(sharedom(:,s),50);
    share25(s)          = prctile(sharedom(:,s),25);
    share10(s)          = prctile(sharedom(:,s),10);
        
    sdshare(s)          = nanstd(sharedom(:,s));
        
    inversehh(s)        = 1./nansum(sharedom(:,s).^2);
        
    topshare(s)         = nanmax(sharedom(:,s));
        
    ndom(s)            = sum((sharedom(:,s)>0));
    
    ihh(s)              = 1./nansum(sharedom(:,s).^2);
    
    share90s(s)          = prctile(sharedoms(:,s),90);
    share75s(s)          = prctile(sharedoms(:,s),75);
    share50s(s)          = prctile(sharedoms(:,s),50);
    share25s(s)          = prctile(sharedoms(:,s),25);
    share10s(s)          = prctile(sharedoms(:,s),10);
        
    sdshares(s)          = nanstd(sharedoms(:,s));
        
    inversehhs(s)        = 1./nansum(sharedoms(:,s).^2);
        
    topshares(s)         = nanmax(sharedoms(:,s));
        
    ndoms(s)            = sum((sharedoms(:,s)>0));
    
    ihhs(s)              = 1./nansum(sharedoms(:,s).^2);
              
end

%%% collate

HH            = [hh];
HHs           = [hhs];

Highest_Share  = [highest_share];
Highest_Shares = [highest_shares];

Share         = [share];
Shares        = [shares];
    
Share90       = [share90];
Share75       = [share75];
Share50       = [share50];
Share25       = [share25];
Share10       = [share10];

Share90s      = [share90s];
Share75s      = [share75s];
Share50s      = [share50s];
Share25s      = [share25s];
Share10s      = [share10s];
    
Sdshare       = [sdshare];
Sdshares      = [sdshares];
    
Ndom          = [ndom];
Ndoms         = [ndoms];
    
IndustryIHH   = [inversehh];
IndustryIHHs  = [inversehhs];

IndustryTop   = [topshare];
IndustryTops  = [topshares];

ppndom            = prctile(Ndom,[10,25,50,75,90,95]); 
ppndoms           = prctile(Ndoms,[10,25,50,75,90,95]); 



%%%%% B. Size distribution moments

sales = pH.*yH; %% domestic only, corrected Chris Jan 2014 
labor = lH;

saless = pFs.*yFs; 
labors = lFs;

sales_s = sum(sales)';
labor_s = sum(labor)';

saless_s = sum(saless)';
labors_s = sum(labors)';

labor = labor(sales>0);
sales = sales(sales>0);

labors = labors(saless>0);
saless = saless(saless>0);

data = sortrows([sales, labor],1);

fRtopY1  = sum(data(end-floor(length(data)*0.01)+1:end,1))/sum(data(:,1)); %% size distribution of establishments
fRtopY5  = sum(data(end-floor(length(data)*0.05)+1:end,1))/sum(data(:,1)); %% (data counterpart is value-added)
fRtopL1  = sum(data(end-floor(length(data)*0.01)+1:end,2))/sum(data(:,2));
fRtopL5  = sum(data(end-floor(length(data)*0.05)+1:end,2))/sum(data(:,2));

datas = sortrows([saless, labors],1);

fRtopY1s  = sum(datas(end-floor(length(datas)*0.01)+1:end,1))/sum(datas(:,1)); %% size distribution of establishments
fRtopY5s  = sum(datas(end-floor(length(datas)*0.05)+1:end,1))/sum(datas(:,1)); %% (data counterpart is value-added)
fRtopL1s  = sum(datas(end-floor(length(datas)*0.01)+1:end,2))/sum(datas(:,2));
fRtopL5s  = sum(datas(end-floor(length(datas)*0.05)+1:end,2))/sum(datas(:,2));

data  = sortrows([sales_s, labor_s], 1);

fRtopR1_s = sum(data(end-floor(length(data)*0.01)+1:end,1))/sum(data(:,1)); %% size distribution of sectors
fRtopR5_s = sum(data(end-floor(length(data)*0.05)+1:end,1))/sum(data(:,1)); %% (data counterpart is sales)
fRtopL1_s = sum(data(end-floor(length(data)*0.01)+1:end,2))/sum(data(:,2));
fRtopL5_s = sum(data(end-floor(length(data)*0.05)+1:end,2))/sum(data(:,2));

datas = sortrows([saless_s, labors_s], 1);

fRtopR1s_s = sum(datas(end-floor(length(datas)*0.01)+1:end,1))/sum(datas(:,1)); %% size distribution of sectors
fRtopR5s_s = sum(datas(end-floor(length(datas)*0.05)+1:end,1))/sum(datas(:,1)); %% (data counterpart is sales)
fRtopL1s_s = sum(datas(end-floor(length(datas)*0.01)+1:end,2))/sum(datas(:,2));
fRtopL5s_s = sum(datas(end-floor(length(datas)*0.05)+1:end,2))/sum(datas(:,2));

%%%%% C. Export/Import moments

fexporters      = sum(omegaH>0 & omegaF>0)./sum(omegaH>0);
fexporterss     = sum(omegaHs>0 & omegaFs>0)./sum(omegaHs>0);
agg_fexporters  = sum(omegaH(:)>0 & omegaF(:)>0)/sum(omegaH(:)>0); 
agg_fexporterss = sum(omegaHs(:)>0 & omegaFs(:)>0)/sum(omegaFs(:)>0); 

imports_s      = sum(pF.*yF);  
import_shares  = imports_s./sum(pF.*yF+pH.*yH); 
agg_impshare   = sum(pF(:).*yF(:))/sum(pH(:).*yH(:) + pF(:).*yF(:));

importss_s      = sum(pHs.*yHs);  
import_sharess  = importss_s./sum(pFs.*yFs+pHs.*yHs); 
agg_impshares   = sum(pHs(:).*yHs(:))/sum(pHs(:).*yHs(:) + pFs(:).*yFs(:));

impshare   = sum(omegaF)';
expshare   = sum(omegaHs)';
ttrade     = impshare + expshare;

impshares   = sum(omegaHs)';
expshares   = sum(omegaF)';
ttrades     = impshares + expshares;

%%%%% GL index

intraindex  = 1 - mean(sj(ttrade>0)'.*abs(impshare(ttrade>0)-expshare(ttrade>0))./(impshare(ttrade>0) + expshare(ttrade>0)))./mean(sj(ttrade>0)');
intraindexs = 1 - mean(sjs(ttrades>0)'.*abs(impshares(ttrades>0)-expshares(ttrades>0))./(impshares(ttrades>0) + expshares(ttrades>0)))./mean(sjs(ttrades>0)');



sector_imports       = sum(pF.*yF);
total_imports        = sum(sector_imports);
sector_share_imports = sector_imports/total_imports;

sector_exports       = sum(pHs.*yHs);
total_exports        = sum(sector_exports);
sector_share_exports = sector_exports/total_exports;

sector_sales         = sum(pH.*yH); 
total_sales          = sum(sector_sales);
sector_share_sales   = sector_sales/total_sales;



sector_imports       = sum(pF.*yF);
total_imports        = sum(sector_imports);
sector_share_imports = sector_imports/total_imports;

sector_importss       = sum(pHs.*yHs);
total_importss        = sum(sector_importss);
sector_share_importss = sector_importss/total_importss;

sector_exports       = sum(pHs.*yHs);
total_exports        = sum(sector_exports);
sector_share_exports = sector_exports/total_exports;

sector_exportss       = sum(pF.*yF);
total_exportss        = sum(sector_exportss);
sector_share_exportss = sector_exportss/total_exportss;

sector_sales         = sum(pH.*yH); 
total_sales          = sum(sector_sales);
sector_share_sales   = sector_sales/total_sales;

sector_saless         = sum(pFs.*yFs);
total_saless          = sum(sector_saless);
sector_share_saless   = sector_saless/total_saless;

temp = corrcoef([sector_share_imports',sector_share_exports',sector_share_sales']);

corr_share_exports_share_sales = temp(2,3); %% correlation of sector share exports with sector share of sales
corr_share_imports_share_sales = temp(1,3); %% correlation of sector share imports with sector share of sales

temps = corrcoef([sector_share_importss',sector_share_exportss',sector_share_saless']);

corr_share_exports_share_saless = temps(2,3); %% correlation of sector share exports with sector share of sales
corr_share_imports_share_saless = temps(1,3); %% correlation of sector share imports with sector share of sales

tempX = [ones(size(sector_share_imports))',sector_share_sales'];
tempY = [sector_share_imports'];

temp_beta = inv(tempX'*tempX)*tempX'*tempY;

beta_share_imports_share_sales = temp_beta(2); %% regression coefficient of sector share imports on sector share *sales*

tempXs = [ones(size(sector_share_importss))',sector_share_saless'];
tempYs = [sector_share_importss'];

temp_betas = inv(tempXs'*tempXs)*tempXs'*tempYs;

beta_share_imports_share_saless = temp_betas(2); %% regression coefficient of sector share imports on sector share *sales*

tempX = [ones(size(sector_share_imports))',sector_share_exports'];
tempY = [sector_share_imports'];

temp_beta = inv(tempX'*tempX)*tempX'*tempY;

beta_share_imports_share_exports = temp_beta(2); %% regression coefficient of sector share imports on sector share *exports*

tempXs = [ones(size(sector_share_importss))',sector_share_exportss'];
tempYs = [sector_share_importss'];

temp_betas = inv(tempXs'*tempXs)*tempXs'*tempYs;

beta_share_imports_share_exportss = temp_betas(2); %% regression coefficient of sector share imports on sector share *exports*

% Block #1

meaninvHH         = nanmean(1./HH); 
medianinvHH       = nanmedian(1./HH);
meanmaxshare      = nanmean(Highest_Share);
medianmaxshare    = nanmedian(Highest_Share);          

Block1            = [meaninvHH,medianinvHH,meanmaxshare,medianmaxshare];

meaninvHHs         = nanmean(1./HHs); 
medianinvHHs       = nanmedian(1./HHs);
meanmaxshares      = nanmean(Highest_Shares);
medianmaxshares    = nanmedian(Highest_Shares);          

Block1s            = [meaninvHHs,medianinvHHs,meanmaxshares,medianmaxshares];

% Block #2 

meanshare         = nanmean(Share);
sdshare           = nanstd(Share); 
medianshare       = nanmedian(Share);
ppshare           = prctile(Share,[75,95,99]);         

Block2            = [meanshare,sdshare,medianshare,ppshare];

meanshares         = nanmean(Shares);
sdshares           = nanstd(Shares); 
medianshares       = nanmedian(Shares);
ppshares           = prctile(Shares,[75,95,99]);         

Block2s            = [meanshares,sdshares,medianshares,ppshares];

% Block #3


Block3            = [fRtopR1_s,fRtopR5_s,fRtopL1_s,fRtopL5_s];
Block3s           = [fRtopR1s_s,fRtopR5s_s,fRtopL1s_s,fRtopL5s_s];

% Block #4 

ppndom            = prctile(Ndom,[10,25,50,75,90,95]); 
ppndoms           = prctile(Ndoms,[10,25,50,75,90,95]); 

Block4            = [ppndom];
Block4s           = [ppndoms];

% Block #5

Block5            = [fRtopY1,fRtopY5,fRtopL1,fRtopL5];
Block5s           = [fRtopY1s,fRtopY5s,fRtopL1s,fRtopL5s];

% Block #6

Block6            = [agg_fexporters,agg_impshare,1-aweight,intraindex,corr_share_exports_share_sales,corr_share_imports_share_sales,beta_share_imports_share_sales,beta_share_imports_share_exports];
Block6s           = [agg_fexporterss,agg_impshares,1-aweights,intraindexs,corr_share_exports_share_saless,corr_share_imports_share_saless,beta_share_imports_share_saless,beta_share_imports_share_exportss];

% Block #7

mean90to10        = nanmean(Share90-Share10);
median90to10      = nanmedian(Share90-Share10);
mean75to25        = nanmean(Share75-Share25);
median75to25      = nanmedian(Share75-Share25);
meanSDshare       = nanmean(Sdshare);
medianSDshare     = nanmedian(Sdshare);                

Block7            = [mean90to10,median90to10,mean75to25,median75to25,meanSDshare,medianSDshare];

mean90to10s        = nanmean(Share90s-Share10s);
median90to10s      = nanmedian(Share90s-Share10s);
mean75to25s        = nanmean(Share75s-Share25s);
median75to25s      = nanmedian(Share75s-Share25s);
meanSDshares       = nanmean(Sdshares);
medianSDshares     = nanmedian(Sdshares);                

Block7s            = [mean90to10s,median90to10s,mean75to25s,median75to25s,meanSDshares,medianSDshares];


% Block #8

ppindustryihh     = prctile(IndustryIHH,[10,25,50,75,90,95]); 
ppindustryihhs    = prctile(IndustryIHHs,[10,25,50,75,90,95]); 

Block8            = [ppindustryihh];
Block8s           = [ppindustryihhs];



% Block #9

ppindustrytop     = prctile(IndustryTop,[10,25,50,75,90,95]); 
ppindustrytops    = prctile(IndustryTops,[10,25,50,75,90,95]); 

Block9            = [ppindustrytop];
Block9s           = [ppindustrytops];

all_model_moments  = [Block1,Block2,Block3,Block4,Block5,Block6,Block7,Block8,Block9]';
all_model_momentss = [Block1s,Block2s,Block3s,Block4s,Block5s,Block6s,Block7s,Block8s,Block9s]'

        

%%%%% list of moments
fprintf('\n');
display('Model Moments                            Home    Foreign')
fprintf('\n');
fprintf('mean HH inverse                      = %7.3f %7.3f \n',nanmean(1./HH)             ,nanmean(1./HHs));  
fprintf('median HH inverse                    = %7.3f %7.3f \n',nanmedian(1./HH)           ,nanmedian(1./HHs));
fprintf('mean highest share                   = %7.3f %7.3f \n',nanmean(Highest_Share)     ,nanmean(Highest_Shares));  
fprintf('median highest share                 = %7.3f %7.3f \n',nanmedian(Highest_Share)   ,nanmedian(Highest_Shares));
fprintf('\n');
fprintf('mean share                           = %7.3f %7.3f \n',nanmean(Share)             ,nanmean(Shares)); 
fprintf('s.d. share                           = %7.3f %7.3f \n',nanstd(Share)              ,nanstd(Shares));
fprintf('median share                         = %7.3f %7.3f \n',nanmedian(Share)           ,nanmedian(Shares)); 
fprintf('75th share                           = %7.3f %7.3f \n',prctile(Share,75)          ,prctile(Shares,75)); 
fprintf('95th share                           = %7.3f %7.3f \n',prctile(Share,95)          ,prctile(Shares,95));
fprintf('99th share                           = %7.3f %7.3f \n',prctile(Share,99)          ,prctile(Shares,99));
fprintf('\n')
fprintf('frac. sales top 0.01 sectors         = %7.3f %7.3f \n',fRtopR1_s                  ,fRtopR1s_s);
fprintf('frac. sales top 0.05 sectors         = %7.3f %7.3f \n',fRtopR5_s                  ,fRtopR5s_s);
fprintf('frac. WL of (same) top 0.01          = %7.3f %7.3f \n',fRtopL1_s                  ,fRtopL1s_s);
fprintf('frac. WL of (same) top 0.05          = %7.3f %7.3f \n',fRtopL5_s                  ,fRtopL5s_s);
fprintf('\n');
fprintf('number producers p10                 = %7.0f %7.0f \n',prctile(Ndom,10)           ,prctile(Ndoms,10));
fprintf('number producers p25                 = %7.0f %7.0f \n',prctile(Ndom,25)           ,prctile(Ndoms,25));
fprintf('number producers p50                 = %7.0f %7.0f \n',prctile(Ndom,50)           ,prctile(Ndoms,50));
fprintf('number producers p75                 = %7.0f %7.0f \n',prctile(Ndom,75)           ,prctile(Ndoms,75));
fprintf('number producers p90                 = %7.0f %7.0f \n',prctile(Ndom,90)           ,prctile(Ndoms,90));
fprintf('number producers p95                 = %7.0f %7.0f \n',prctile(Ndom,95)           ,prctile(Ndoms,95));
fprintf('\n');
fprintf('frac. value added top 0.01 establ.   = %7.3f %7.3f \n',fRtopY1                    ,fRtopY1s);
fprintf('frac. value added top 0.05 establ.   = %7.3f %7.3f \n',fRtopY5                    ,fRtopY5s);
fprintf('frac. WL of (same) top 0.01          = %7.3f %7.3f \n',fRtopL1                    ,fRtopL1s);
fprintf('frac. WL of (same) top 0.05          = %7.3f %7.3f \n',fRtopL5                    ,fRtopL5s);
fprintf('\n');
fprintf('agg. fraction exporters              = %7.3f %7.3f \n',agg_fexporters             ,agg_fexporterss);
fprintf('agg. import share                    = %7.3f %7.3f \n',agg_impshare               ,agg_impshares);
fprintf('\n');
fprintf('index import share dispers.          = %7.3f %7.3f \n',1-aweight                  ,1-aweights); % switched definition to 1- what we had before
fprintf('Grubel-Lloyd (GL) index              = %7.3f %7.3f \n',intraindex                 ,intraindexs);
fprintf('\n');
fprintf('corr(share exports , share sales)    = %7.3f %7.3f \n',corr_share_exports_share_sales  ,corr_share_exports_share_saless);  
fprintf('corr(share imports , share sales)    = %7.3f %7.3f \n',corr_share_imports_share_sales  ,corr_share_imports_share_saless);  
fprintf('\n');
fprintf('coeff. share imports wrt sh. sales   = %7.3f %7.3f \n',beta_share_imports_share_sales  ,beta_share_imports_share_saless);   
fprintf('coeff. share imports wrt sh, exports = %7.3f %7.3f \n',beta_share_imports_share_exports,beta_share_imports_share_exportss);   
fprintf('\n');
fprintf('mean   p90-p10                       = %7.3f %7.3f \n',nanmean(Share90-Share10)   ,nanmean(Share90s-Share10s)); 
fprintf('median p90-p10                       = %7.3f %7.3f \n',nanmedian(Share90-Share10) ,nanmedian(Share90s-Share10s)); 
fprintf('mean   p75-p25                       = %7.3f %7.3f \n',nanmean(Share75-Share25)   ,nanmean(Share75s-Share25s)); 
fprintf('median p75-p25                       = %7.3f %7.3f \n',nanmedian(Share75-Share25) ,nanmedian(Share75s-Share25s)); 
fprintf('mean   sd share                      = %7.3f %7.3f \n',nanmean(Sdshare)           ,nanmean(Sdshares)); 
fprintf('median sd share                      = %7.3f %7.3f \n',nanmedian(Sdshare)         ,nanmedian(Sdshares));

fprintf('\n');
fprintf('HH inverse p10                       = %7.3f %7.3f \n',prctile(IndustryIHH,10)    ,prctile(IndustryIHHs,10));
fprintf('HH inverse p25                       = %7.3f %7.3f \n',prctile(IndustryIHH,25)    ,prctile(IndustryIHHs,25));
fprintf('HH inverse p50                       = %7.3f %7.3f \n',prctile(IndustryIHH,50)    ,prctile(IndustryIHHs,50));
fprintf('HH inverse p75                       = %7.3f %7.3f \n',prctile(IndustryIHH,75)    ,prctile(IndustryIHHs,75));
fprintf('HH inverse p90                       = %7.3f %7.3f \n',prctile(IndustryIHH,90)    ,prctile(IndustryIHHs,90));
fprintf('HH inverse p95                       = %7.3f %7.3f \n',prctile(IndustryIHH,95)    ,prctile(IndustryIHHs,95));

fprintf('\n');
fprintf('top share  p10                       = %7.3f %7.3f \n',prctile(IndustryTop,10)    ,prctile(IndustryTops,10)); 
fprintf('top share  p25                       = %7.3f %7.3f \n',prctile(IndustryTop,25)    ,prctile(IndustryTops,25));
fprintf('top share  p50                       = %7.3f %7.3f \n',prctile(IndustryTop,50)    ,prctile(IndustryTops,50));
fprintf('top share  p75                       = %7.3f %7.3f \n',prctile(IndustryTop,75)    ,prctile(IndustryTops,75));
fprintf('top share  p90                       = %7.3f %7.3f \n',prctile(IndustryTop,90)    ,prctile(IndustryTops,90));
fprintf('top share  p95                       = %7.3f %7.3f \n',prctile(IndustryTop,95)    ,prctile(IndustryTops,95));


if 1
    
disp('Home aggregate labor share')
agg_lshare=[sum(labor_s)/sum(sales_s)]

disp('Home average labor share')
avg_lshare = mean(labor./sales)

disp('Home real wage')
real_wage = 1/P

disp('Foreign aggregate labor share')
agg_lshares=[Ws*sum(labors_s)/sum(saless_s)]

disp('Foreign average labor share')
avg_lshares = mean(Ws*labors./saless)

disp('Foreign real wage')
real_wages = Ws/Ps

disp('Relative real wage')
relative_real_wage = real_wages/real_wage

end 


