clear all;
clc;

fprintf('...calculating data moments... \n');

%%%%% load in daniel's data

T = 4;                             %% four years of data

data = cell(T,1);                  %% cell to store all the raw data

load data_89.mat;                  %% matrix where each row is a sector

data{1,1} = datam'; clear datam;   %% matrix where each column is a sector

load data_91.mat; 

data{2,1} = datam'; clear datam;

load data_92.mat; 

data{3,1} = datam'; clear datam;

load data_93.mat; 

data{4,1} = datam'; clear datam;


HH = []; Highest_Share = []; 

Share = []; Share90 = []; Share75 = []; Share50 = []; Share25 = []; Share10 = [];
Sdshare = []; Ndom = []; IndustryIHH = []; IndustryTop = []; Sales = []; Sales_s = [];

for t=1:T,
    
    
    data_t = data{t,1};
    
    [NN,SS] = size(data_t);  

    data_t(isinf(data_t))=NaN; %% replace all "inf" with "nan"

    sales        = data_t;
    
    sales_s      = nansum(sales); 
    sector_sales = repmat(sales_s,NN,1);

    sharedom       = sales./sector_sales; %% domestic market shares 
     
    hh            = nansum(sharedom.^2);
    highest_share = max(sharedom);
 
    share = sharedom(:)';
    
    for s = 1:SS,
        
        share90(s)          = prctile(sharedom(:,s),90);
        share75(s)          = prctile(sharedom(:,s),75);
        share50(s)          = prctile(sharedom(:,s),50);
        share25(s)          = prctile(sharedom(:,s),25);
        share10(s)          = prctile(sharedom(:,s),10);
        
        sdshare(s)          = nanstd(sharedom(:,s));
        
        inversehh(s)        = 1./nansum(sharedom(:,s).^2);
        
        topshare(s)         = max(sharedom(:,s));
        
        ndom(s)             = sum((sharedom(:,s)>0)); %sum(~isnan(sharedom(:,s)));
        
        
        
          
    end
  
    %%%%% pool all 4 years into one big data set
    
    HH            = [HH,hh];
    Highest_Share = [Highest_Share,highest_share];
    Share         = [Share, share];
    
    Share90       = [Share90, share90];
    Share75       = [Share75, share75];
    Share50       = [Share50, share50];
    Share25       = [Share25, share25];
    Share10       = [Share10, share10];
    
    Sdshare       = [Sdshare, sdshare];
    
    Ndom          = [Ndom,ndom];
    
    IndustryIHH   = [IndustryIHH,inversehh];
    IndustryTop   = [IndustryTop,topshare];
    
    
    Sales         = [Sales  , sales(sales>0)'];
    Sales_s       = [Sales_s, sales_s];
    
    
end



Sales = sortrows(Sales',1); 
    
fRtopR1 = sum(Sales(end-floor(length(Sales)*0.01)+1:end,1))/sum(Sales(:,1));
fRtopR5 = sum(Sales(end-floor(length(Sales)*0.05)+1:end,1))/sum(Sales(:,1));

Sales_s = sortrows(Sales_s',1); 
    
fRtopR1_s = sum(Sales_s(end-floor(length(Sales_s)*0.01)+1:end,1))/sum(Sales_s(:,1));
fRtopR5_s = sum(Sales_s(end-floor(length(Sales_s)*0.05)+1:end,1))/sum(Sales_s(:,1));

%alternate calculation
%x=quantile(Sales_s,0.99,3);
%[~,xi] = min(abs(x-Sales_s)); 
%fRtopR1_s = 1-sum(Sales_s(1:xi))/sum(Sales_s); % gives same answer
   

%%%%% compute moments

fprintf('\n');
fprintf('Pooled Moments Calculated in Matlab \n');
fprintf('\n');
fprintf('mean HH inverse               = %7.3f \n',nanmean(1./HH));  
fprintf('median HH inverse             = %7.3f \n',nanmedian(1./HH));
fprintf('mean highest share            = %7.3f \n',nanmean(Highest_Share));  
fprintf('median highest share          = %7.3f \n',nanmedian(Highest_Share));
fprintf('\n');
fprintf('mean share                    = %7.3f \n',nanmean(Share)); 
fprintf('s.d. share                    = %7.3f \n',nanstd(Share));
fprintf('median share                  = %7.3f \n',nanmedian(Share)); 
fprintf('75th share                    = %7.3f \n',prctile(Share,75)); 
fprintf('95th share                    = %7.3f \n',prctile(Share,95));
fprintf('99th share                    = %7.3f \n',prctile(Share,99));
fprintf('\n');
fprintf('frac. sales top 0.01 sectors  = %7.3f \n',fRtopR1_s); 
fprintf('frac. sales top 0.05 sectors  = %7.3f \n',fRtopR5_s); 
fprintf('\n');
fprintf('number producers p10          = %7.0f \n',prctile(Ndom,10));
fprintf('number producers p25          = %7.0f \n',prctile(Ndom,25));
fprintf('number producers p50          = %7.0f \n',prctile(Ndom,50));
fprintf('number producers p75          = %7.0f \n',prctile(Ndom,75));
fprintf('number producers p90          = %7.0f \n',prctile(Ndom,90));
fprintf('number producers p95          = %7.0f \n',prctile(Ndom,95));
fprintf('\n');
%fprintf('frac. sales top 0.01 establ.  = %7.3f \n',fRtopR1); % we do value-added in the data
%fprintf('frac. sales top 0.05 establ.  = %7.3f \n',fRtopR5);
fprintf('\n');
fprintf('mean   p90-p10                = %7.3f \n',nanmean(Share90-Share10)); 
fprintf('median p90-p10                = %7.3f \n',nanmedian(Share90-Share10)); 
fprintf('mean   p75-p25                = %7.3f \n',nanmean(Share75-Share25)); 
fprintf('median p75-p25                = %7.3f \n',nanmedian(Share75-Share25)); 
fprintf('mean   sd share               = %7.3f \n',nanmean(Sdshare)); 
fprintf('median sd share               = %7.3f \n',nanmedian(Sdshare));
fprintf('\n');
fprintf('HH inverse p10                = %7.3f \n',prctile(IndustryIHH,10));
fprintf('HH inverse p25                = %7.3f \n',prctile(IndustryIHH,25));
fprintf('HH inverse p50                = %7.3f \n',prctile(IndustryIHH,50));
fprintf('HH inverse p75                = %7.3f \n',prctile(IndustryIHH,75));
fprintf('HH inverse p90                = %7.3f \n',prctile(IndustryIHH,90));
fprintf('HH inverse p95                = %7.3f \n',prctile(IndustryIHH,95));
fprintf('\n');
fprintf('top share  p10                = %7.3f \n',prctile(IndustryTop,10)); 
fprintf('top share  p25                = %7.3f \n',prctile(IndustryTop,25));
fprintf('top share  p50                = %7.3f \n',prctile(IndustryTop,50));
fprintf('top share  p75                = %7.3f \n',prctile(IndustryTop,75));
fprintf('top share  p90                = %7.3f \n',prctile(IndustryTop,90));
fprintf('top share  p95                = %7.3f \n',prctile(IndustryTop,95));
fprintf('\n');

