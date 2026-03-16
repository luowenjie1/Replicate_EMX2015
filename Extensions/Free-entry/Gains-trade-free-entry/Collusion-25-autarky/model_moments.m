fprintf('...calculating model moments... \n');

%%%%% A. Domestic market share moments   

omega1H = (p1H./Pi).^(1-gamma); 
omega2H = (p2H./Pi).^(1-gamma); 
omega1F = (p1F./Pi).^(1-gamma);
omega2F = (p2F./Pi).^(1-gamma); 

dshare1 = omega1H./(n1.*omega1H + n2.*omega2H); 
dshare2 = omega2H./(n1.*omega1H + n2.*omega2H); 

hshare = max(dshare1, dshare2); 
lshare = min(dshare1, dshare2); 

cutoff = n2./(n1+n2).*(dshare1>dshare2) + n1./(n1+n2).*(dshare1<=dshare2); 

share90 = hshare.*(cutoff<=0.90) + lshare.*(cutoff>0.90); 
share75 = hshare.*(cutoff<=0.75) + lshare.*(cutoff>0.75); 
share50 = hshare.*(cutoff<=0.50) + lshare.*(cutoff>0.50); 
share25 = hshare.*(cutoff<=0.25) + lshare.*(cutoff>0.25); 
share10 = hshare.*(cutoff<=0.10) + lshare.*(cutoff>0.10); 

share90(n1==0&n2>0) = dshare2(n1==0&n2>0); 
share75(n1==0&n2>0) = dshare2(n1==0&n2>0); 
share50(n1==0&n2>0) = dshare2(n1==0&n2>0); 
share25(n1==0&n2>0) = dshare2(n1==0&n2>0); 
share10(n1==0&n2>0) = dshare2(n1==0&n2>0);

share90(n1>0&n2==0) = dshare1(n1>0&n2==0); 
share75(n1>0&n2==0) = dshare1(n1>0&n2==0); 
share50(n1>0&n2==0) = dshare1(n1>0&n2==0); 
share25(n1>0&n2==0) = dshare1(n1>0&n2==0); 
share10(n1>0&n2==0) = dshare1(n1>0&n2==0);

mshare = n1./(n1+n2).*dshare1 + n2./(n1+n2).*dshare2; 
sdshare = sqrt(n1./(n1+n2).*(dshare1-mshare).^2 + n2./(n1+n2).*(dshare2-mshare).^2);

hh = n1.*dshare1.^2 + n2.*dshare2.^2; 

topshare  = hshare.*(n1>0&n2>0) + dshare1.*(n1>0&n2==0) + dshare2.*(n1==0&n2>0);

ndom      = n1 + n2; 

load saved_data_moments_jan2014.mat;  % vector of data moments, ordered same as model moments above
                                      % particular care should be taken if there is any
                                      % reordering of moments or change in
                                      % moment selection
                                      
data      = sortrows([wi, 1./hh], 2); 
data(:,1) = cumsum(data(:,1)); 

meanihh   = wi'*(1./hh); 
medianihh = min(data(data(:,1)>=0.5,2));

data      = sortrows([wi, topshare], 2); 
data(:,1) = cumsum(data(:,1)); 

meantop   = wi'*topshare;
mediantop = min(data(data(:,1)>=0.5,2));

display('                                         Model    Data')
fprintf('\n');
fprintf('mean HH inverse                      = %7.3f %7.3f \n', meanihh              ,all_data_moments(1));  
fprintf('median HH inverse                    = %7.3f %7.3f \n', medianihh            ,all_data_moments(2));
fprintf('mean highest share                   = %7.3f %7.3f \n', meantop              ,all_data_moments(3));  
fprintf('median highest share                 = %7.3f %7.3f \n', mediantop            ,all_data_moments(4));
fprintf('\n');

Block1 = [meanihh, medianihh, meantop, mediantop]; 

% pool market shares of all firms in all sectors and report key
% statistics

share = [wi.*n1, dshare1; wi.*n2, dshare2]; 
share(share(:,1) == 0,:) = [];

weight = share(:,1)/sum(share(:,1)); 
share  = share(:,2); 

psharemean = weight'*share;
psharesd   = sqrt(weight'*(share - psharemean).^2); 

data      = [weight, share]; 
data      = sortrows(data, 2); 
data(:,1) = cumsum(data(:,1)); 

pshare50   = min(data(data(:,1)>=.50,2));
pshare75   = min(data(data(:,1)>=.75,2));
pshare95   = min(data(data(:,1)>=.95,2));
pshare99   = min(data(data(:,1)>=.99,2));


fprintf('mean share                           = %7.3f %7.3f \n', psharemean        ,all_data_moments(5)); 
fprintf('s.d. share                           = %7.3f %7.3f \n', psharesd          ,all_data_moments(6));
fprintf('median share                         = %7.3f %7.3f \n', pshare50          ,all_data_moments(7)); 
fprintf('75th share                           = %7.3f %7.3f \n', pshare75          ,all_data_moments(8)); 
fprintf('95th share                           = %7.3f %7.3f \n', pshare95          ,all_data_moments(9));
fprintf('99th share                           = %7.3f %7.3f \n', pshare99          ,all_data_moments(10));
fprintf('\n');

Block2 = [psharemean, psharesd, pshare50, pshare75, pshare95, pshare99];

%%%%% B. Size distribution moments (ARE THESE TOTAL FOR A GIVEN FIRM?
%%%%% DOMESTIC + EXPORT SALES? IF SO, OLD CODE WAS WRONG)

% Do this for a bunch of different values of mu, the Pareto shape
% parameter, and pick whichever fits the data the best

nmu = 16; % number of different values of mu we search over (this is actually mu/(theta-1) where mu is the tail of the productivity distribution)

muz = nodeunif(nmu, 1.1, 5); 
Kz  = 21; % number of nodes to integrate

sizemoments = zeros(8, nmu); 
gap         = zeros(nmu, 1); 

fprintf('\n');
fprintf('Search over sectoral Pareto shape that best fits size distribution \n');
fprintf('\n');

for itt = 1 : nmu
    
 [zz, wz] = qnwunif(Kz, 0, 1); 
 zz       = (1-zz).^(-1/muz(itt)); 
   
dsales1 = p1H.*y1H; 
dsales2 = p2H.*y2H;
xsales1 = tau*p1Hs.*y1Hs;
xsales2 = tau*p2Hs.*y2Hs;

tsales1 = dsales1 + xsales1; 
tsales2 = dsales2 + xsales2; 

dlabor1 = y1H./a1H;
dlabor2 = y2H./a2H;
xlabor1 = tau*y1Hs; 
xlabor2 = tau*y2Hs/ebar;

tlabor1 = dlabor1 + xlabor1; 
tlabor2 = dlabor2 + xlabor2; 

data    = [wi.*n1, dsales1, dlabor1; wi.*n2, dsales2, dlabor2]; 

data(data(:,1) == 0,:) = [];
data(:,1) = data(:,1)/sum(data(:,1)); 

ns        = size(data, 1); 

data      = repmat(data, Kz, 1); 

data(:,1) = data(:,1).*kron(wz, ones(ns,1));
data(:,2) = data(:,2).*kron(zz, ones(ns,1)); 
data(:,3) = data(:,3).*kron(zz, ones(ns,1)); 

data      = sortrows(data, 2); 

cshare  = cumsum(data(:,1)); 

fRtopY1 = data(cshare>=0.99,1)'*data(cshare>=0.99,2)/(data(:,1)'*data(:,2));
fRtopY5 = data(cshare>=0.95,1)'*data(cshare>=0.95,2)/(data(:,1)'*data(:,2));

fRtopL1 = data(cshare>=0.99,1)'*data(cshare>=0.99,3)/(data(:,1)'*data(:,3));
fRtopL5 = data(cshare>=0.95,1)'*data(cshare>=0.95,3)/(data(:,1)'*data(:,3));

sizemoments(1, itt) = fRtopY1; 
sizemoments(2, itt) = fRtopY5; 
sizemoments(3, itt) = fRtopL1; 
sizemoments(4, itt) = fRtopL5; 

% Sectoral moments

sales_s = n1.*dsales1 + n2.*dsales2;
labor_s = n1.*dlabor1 + n2.*dlabor2;

data    = [wi, sales_s, labor_s];

data(data(:,1) == 0,:) = [];
data(:,1) = data(:,1)/sum(data(:,1)); 

ns = size(data, 1); 

data      = repmat(data, Kz, 1); 

data(:,1) = data(:,1).*kron(wz, ones(ns,1));
data(:,2) = data(:,2).*kron(zz, ones(ns,1)); 
data(:,3) = data(:,3).*kron(zz, ones(ns,1)); 

data      = sortrows(data, 2); 

cshare  = cumsum(data(:,1)); 

fRtopY1_s = data(cshare>=0.99,1)'*data(cshare>=0.99,2)/(data(:,1)'*data(:,2));
fRtopY5_s = data(cshare>=0.95,1)'*data(cshare>=0.95,2)/(data(:,1)'*data(:,2));

fRtopL1_s = data(cshare>=0.99,1)'*data(cshare>=0.99,3)/(data(:,1)'*data(:,3));
fRtopL5_s = data(cshare>=0.95,1)'*data(cshare>=0.95,3)/(data(:,1)'*data(:,3));


sizemoments(5, itt) = fRtopY1_s; 
sizemoments(6, itt) = fRtopY5_s; 
sizemoments(7, itt) = fRtopL1_s; 
sizemoments(8, itt) = fRtopL5_s; 

gap(itt) = norm(abs(sizemoments(:,itt) - [all_data_moments(21:24); all_data_moments(11:14)])); 

end


[ii, jj]  = min(gap); 


fprintf('\n');
fprintf('\n');
fprintf('Pareto shape sectoral draws          = %7.3f  \n',muz(jj)                  );
fprintf('\n');
fprintf('\n');
fprintf('frac. sales top 0.01 sectors         = %7.3f %7.3f \n',sizemoments(5,jj)                  ,all_data_moments(11));
fprintf('frac. sales top 0.05 sectors         = %7.3f %7.3f \n',sizemoments(6,jj)                  ,all_data_moments(12));
fprintf('frac. WL of (same) top 0.01          = %7.3f %7.3f \n',sizemoments(7,jj)                  ,all_data_moments(13));
fprintf('frac. WL of (same) top 0.05          = %7.3f %7.3f \n',sizemoments(8,jj)                  ,all_data_moments(14));
fprintf('\n');

data      = sortrows([wi, ndom], 2); 
data(:,1) = cumsum(data(:,1)); 

p10ndom = min(data(data(:,1)>=0.10,2));
p25ndom = min(data(data(:,1)>=0.25,2));
p50ndom = min(data(data(:,1)>=0.50,2));
p75ndom = min(data(data(:,1)>=0.75,2));
p90ndom = min(data(data(:,1)>=0.90,2));
p95ndom = min(data(data(:,1)>=0.95,2));


fprintf('number producers p10                 = %7.0f %7.0f \n',p10ndom                             ,all_data_moments(15));
fprintf('number producers p25                 = %7.0f %7.0f \n',p25ndom                             ,all_data_moments(16));
fprintf('number producers p50                 = %7.0f %7.0f \n',p50ndom                             ,all_data_moments(17));
fprintf('number producers p75                 = %7.0f %7.0f \n',p75ndom                             ,all_data_moments(18));
fprintf('number producers p90                 = %7.0f %7.0f \n',p90ndom                             ,all_data_moments(19));
fprintf('number producers p95                 = %7.0f %7.0f \n',p95ndom                             ,all_data_moments(20));
fprintf('\n');
fprintf('frac. value added top 0.01 establ.   = %7.3f %7.3f \n',sizemoments(1,jj)                    ,all_data_moments(21));
fprintf('frac. value added top 0.05 establ.   = %7.3f %7.3f \n',sizemoments(2,jj)                    ,all_data_moments(22));
fprintf('frac. WL of (same) top 0.01          = %7.3f %7.3f \n',sizemoments(3,jj)                    ,all_data_moments(23));
fprintf('frac. WL of (same) top 0.05          = %7.3f %7.3f \n',sizemoments(4,jj)                    ,all_data_moments(24));

Block3 = [sizemoments(5:8,jj)', p10ndom, p25ndom, p50ndom, p75ndom, p90ndom, p95ndom, ...
          sizemoments(1:4,jj)'];


%%%%% C. Export/Import moments

impshare = 1 - lambdai; 
expshare = (n1.*p1Hs.*y1Hs.*tau + n2.*p2Hs.*y2Hs.*tau)./(n1.*p1H.*y1H + n2.*p2H.*y2H + n1.*p1Hs.*y1Hs.*tau + n2.*p2Hs.*y2Hs.*tau);
ttrade   = impshare + expshare;
GLi      = 1 - abs(impshare - expshare)./ttrade;  
GL       = wi'*(si.*GLi);          % 1 means all trade is within industry

impshare = wi'*(si.*impshare); 

fprintf('\n');
fprintf('agg. import share                    = %7.3f %7.3f \n', impshare           ,all_data_moments(26));
fprintf('\n');
fprintf('index import share dispers.          = %7.3f %7.3f \n', aweight            ,all_data_moments(27));
fprintf('Grubel-Lloyd (GL) index              = %7.3f %7.3f \n', GL                 ,all_data_moments(28));
fprintf('\n');

Block4 = [impshare, aweight, GL]; 

%%%%% D. Across indistry moments

data      = sortrows([wi, share90 - share10], 2); 
data(:,1) = cumsum(data(:,1)); 

mean9010   = wi'*(share90 - share10);
median9010 = min(data(data(:,1)>=0.5,2));

data      = sortrows([wi, share75 - share25], 2); 
data(:,1) = cumsum(data(:,1)); 

mean7525   = wi'*(share75 - share25);
median7525 = min(data(data(:,1)>=0.5,2));

data      = sortrows([wi, sdshare], 2); 
data(:,1) = cumsum(data(:,1)); 

meansdshare   = wi'*sdshare;
mediansdshare = min(data(data(:,1)>=0.5,2));


data      = sortrows([wi, 1./hh], 2); 
data(:,1) = cumsum(data(:,1)); 

p10ihh = min(data(data(:,1)>=0.10,2));
p25ihh = min(data(data(:,1)>=0.25,2));
p50ihh = min(data(data(:,1)>=0.50,2));
p75ihh = min(data(data(:,1)>=0.75,2));
p90ihh = min(data(data(:,1)>=0.90,2));
p95ihh = min(data(data(:,1)>=0.95,2));

data      = sortrows([wi, topshare], 2); 
data(:,1) = cumsum(data(:,1)); 

p10topshare = min(data(data(:,1)>=0.10,2));
p25topshare = min(data(data(:,1)>=0.25,2));
p50topshare = min(data(data(:,1)>=0.50,2));
p75topshare = min(data(data(:,1)>=0.75,2));
p90topshare = min(data(data(:,1)>=0.90,2));
p95topshare = min(data(data(:,1)>=0.95,2));

smarkup     = ((n1.*y1H./a1H + n2.*y2H./a2H)./(n1.*p1H.*y1H + n2.*p2H.*y2H)).^(-1);

data      = sortrows([wi, smarkup], 2); 
data(:,1) = cumsum(data(:,1)); 

p25smarkup = min(data(data(:,1)>=0.25,2));
p50smarkup = min(data(data(:,1)>=0.50,2));
p75smarkup = min(data(data(:,1)>=0.75,2));
p90smarkup = min(data(data(:,1)>=0.90,2));
p95smarkup = min(data(data(:,1)>=0.95,2));
p99smarkup = min(data(data(:,1)>=0.99,2));


fprintf('\n');
fprintf('mean   p90-p10                       = %7.3f %7.3f \n',mean9010                   ,all_data_moments(33)); 
fprintf('median p90-p10                       = %7.3f %7.3f \n',median9010                 ,all_data_moments(34)); 
fprintf('mean   p75-p25                       = %7.3f %7.3f \n',mean7525                   ,all_data_moments(35)); 
fprintf('median p75-p25                       = %7.3f %7.3f \n',median7525                 ,all_data_moments(36)); 
fprintf('mean   sd share                      = %7.3f %7.3f \n',meansdshare                ,all_data_moments(37)); 
fprintf('median sd share                      = %7.3f %7.3f \n',mediansdshare              ,all_data_moments(38));
fprintf('\n');
fprintf('HH inverse p10                       = %7.3f %7.3f \n',p10ihh                     ,all_data_moments(39));
fprintf('HH inverse p25                       = %7.3f %7.3f \n',p25ihh                     ,all_data_moments(40));
fprintf('HH inverse p50                       = %7.3f %7.3f \n',p50ihh                     ,all_data_moments(41));
fprintf('HH inverse p75                       = %7.3f %7.3f \n',p75ihh                     ,all_data_moments(42));
fprintf('HH inverse p90                       = %7.3f %7.3f \n',p90ihh                     ,all_data_moments(43));
fprintf('HH inverse p95                       = %7.3f %7.3f \n',p95ihh                     ,all_data_moments(44));
fprintf('\n');
fprintf('top share  p10                       = %7.3f %7.3f \n',p10topshare                ,all_data_moments(45)); 
fprintf('top share  p25                       = %7.3f %7.3f \n',p25topshare                ,all_data_moments(46));
fprintf('top share  p50                       = %7.3f %7.3f \n',p50topshare                ,all_data_moments(47));
fprintf('top share  p75                       = %7.3f %7.3f \n',p75topshare                ,all_data_moments(48));
fprintf('top share  p90                       = %7.3f %7.3f \n',p90topshare                ,all_data_moments(49));
fprintf('top share  p95                       = %7.3f %7.3f \n',p95topshare                ,all_data_moments(50));
fprintf('\n');

all_data_moments = [all_data_moments; 1.61; 1.14; 1.43; 1.82; 3.16];

fprintf('sectoral markup  p50                 = %7.3f %7.3f \n',p50smarkup                 ,all_data_moments(51)); 
fprintf('sectoral markup  p75 to p50          = %7.3f %7.3f \n',p75smarkup/p50smarkup      ,all_data_moments(52)); 
fprintf('sectoral markup  p90 to p50          = %7.3f %7.3f \n',p90smarkup/p50smarkup      ,all_data_moments(53)); 
fprintf('sectoral markup  p95 to p50          = %7.3f %7.3f \n',p95smarkup/p50smarkup      ,all_data_moments(54)); 
fprintf('sectoral markup  p99 to p50          = %7.3f %7.3f \n',p99smarkup/p50smarkup      ,all_data_moments(55)); 

Block5 = [mean9010, median9010, mean7525, median7525, meansdshare, mediansdshare];

Block6 = [p10ihh, p25ihh, p50ihh, p75ihh, p90ihh, p95ihh];

Block7 = [p10topshare, p25topshare, p50topshare, p75topshare, p90topshare, p95topshare];

Block8 = [p50smarkup, p75smarkup/p50smarkup, p90smarkup/p50smarkup, p95smarkup/p50smarkup, p99smarkup/p50smarkup];

all_model_moments = [Block1, Block2, Block3, Block4, Block5, Block6, Block7, Block8]';


%%%%% load data moments from saved file

%load saved_data_moments_sept2013.mat;   % vector of data moments, ordered same as model moments above
                                        % particular care should be taken if there is any
                                        % reordering of moments or change in
                                        % moment selection
                                  

all_data_moments([25; 29; 30; 31; 32]) = [];                                  
                                     
                                     
if length(all_model_moments)~=length(all_data_moments);
    
    display('problem with number of model and data moments!')
    break
end

loss = (all_data_moments-all_model_moments)./all_data_moments; 

fprintf('\n');
fprintf('\n');
fprintf('mean absolute percentage loss        = %7.3f   \n'    ,mean(abs(loss)));
fprintf('\n');   
fprintf('\n');        
   



