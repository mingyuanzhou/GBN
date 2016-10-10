% Plots of PRG & Bessel distributions with different parameters
% Yulai Cong
% 2016 07 13
clear,clc,close all
LineSettings    =  {':*',':o',':s',':d',':>',':h',':x',':p'}    ;
xlabelsize  =   16  ;   ylabelsize  =   16  ;   gcafontsize     =   18  ;
linewidthsize   =   1.5   ;   Markersize = 10 ;

ColorSet    =   [0 0 0 ; 0 0 1 ; 0 1 0 ; 0 1 1 ; 1 0 0 ; 1 0 1 ; 1 1 0 ; 1 1 1]     ;
ColorSet    =   [ColorSet ; ColorSet(2:end,:) * 0.5]    ;
LineStyleCellSet    =   {'-' , '--' , ':' , '-.' }  ;
MarkerCellSet       =   {'o' , 'x' , 's' , 'd' , '^' , 'v' , '<' , '>' , 'p' , 'h' , '.' , '+' , '*'}   ;


%% PRG distributions: x ~ PRG( lambda , c )
xmax    =   20  ;

x   =   0 : 0.01 : xmax     ;
y   =   x(2:end)    ;

lambdaSet   =       [1.0    2.0     3.0     5.0     9.0     7.5     0.5]        ;   
cSet        =   1./ [2.0    2.0     2.0     1.0     0.5     1.0     1.0]        ;
StrLegend   =   cell(1,length(lambdaSet))  ;

for i = 1:length(lambdaSet)
    
    lambda  =   lambdaSet(i)    ;   c       =   cSet(i)         ;
    tmp     =   ['\lambda=',num2str(lambda),' , c=',num2str(c)]   ;
    StrLegend{i}   =   tmp  ; 
    
    figure(1),hold on,box on,
    plot(0, exp(-lambda) ,['-',MarkerCellSet{i}],'color',ColorSet(i,:),'LineWidth',linewidthsize,'MarkerSize',Markersize),
    
end
legend(StrLegend,'FontSize',ylabelsize),

for i = 1:length(lambdaSet)
    
    lambda  =   lambdaSet(i)    ;   c       =   cSet(i)         ;
    tmp     =   ['\lambda=',num2str(lambda),' , c=',num2str(c)]   ;
    StrLegend{i}   =   tmp  ;
    
    PRGy    =   exp(-lambda - c * y) .* sqrt(lambda*c ./ y) .* besseli(-1 , 2 * sqrt(lambda*c*y))   ;
    
    figure(1),hold on,box on,
    plot(y,PRGy,'color',ColorSet(i,:),'LineWidth',linewidthsize),
  
end
set(gca,'fontsize',gcafontsize)
ylim([0 0.61])

filename    =   'PRGFigure'    ;
saveas(gcf,filename,'pdf')     ;




%% PRG distributions: x ~ PRG( lambda , c )
xmax    =   20  ;

x   =   0 : 0.01 : xmax     ;
y   =   x(2:end)    ;

lambdaSet   =       [1.0    2.0     3.0     5.0     9.0     7.5     0.5]        ;   
cSet        =   1./ [2.0    2.0     2.0     1.0     0.5     1.0     1.0]        ;
StrLegend   =   cell(1,length(lambdaSet))  ;

for i = 1:length(lambdaSet)
    
    lambda  =   lambdaSet(i)    ;   c       =   cSet(i)         ;
    tmp     =   ['\lambda=',num2str(lambda),' , c=',num2str(c)]   ;
    StrLegend{i}   =   tmp  ;
    
    PRGy    =   exp(-lambda - c * y) .* sqrt(lambda*c ./ y) .* besseli(-1 , 2 * sqrt(lambda*c*y))   ;
    
    figure(2),hold on,box on,
    plot(y,PRGy,'color',ColorSet(i,:),'LineWidth',linewidthsize),
  
end
legend(StrLegend,'FontSize',ylabelsize),
for i = 1:length(lambdaSet)
    
    lambda  =   lambdaSet(i)    ;   c       =   cSet(i)         ;
    tmp     =   ['\lambda=',num2str(lambda),' , c=',num2str(c)]   ;
    StrLegend{i}   =   tmp  ; 
    
    figure(2),hold on,box on,
    plot(0, exp(-lambda) ,['-',MarkerCellSet{i}],'color',ColorSet(i,:),'LineWidth',linewidthsize,'MarkerSize',Markersize),
    
end
set(gca,'fontsize',gcafontsize)
ylim([0 0.61])

filename    =   'PRGFigure1'    ;
saveas(gcf,filename,'pdf')     ;




%% bessel distributions : n ~ bessel_{-1} (alpha)
nmax    =   10  ;

n   =   1 : nmax     ;

alphaSet   =       [0.1     1.0    5.0     10.0     15.0    20.0]        ;   
StrLegend   =   cell(1,length(alphaSet))  ;

for i = 1:length(alphaSet)
    
    alpha   =   alphaSet(i)    ; 
    tmp     =   ['\alpha=',num2str(alpha)]   ;
    StrLegend{i}   =   tmp  ;
    
    BESSELn     =   (alpha / 2).^(2*n-1) ./ besseli(-1 , alpha) ./ gamma(n+1) ./ gamma(n)    ;
    
    figure(3);hold on,box on,
    plot(n,BESSELn,LineSettings{i},'LineWidth',linewidthsize/1.2,'MarkerSize',Markersize),
  
end
legend(StrLegend,'FontSize',ylabelsize)
set(gca,'fontsize',gcafontsize),set(gca, 'XTick', 1:1:nmax)
xlim([1 nmax])

filename    =   'BesselFigure'    ;
saveas(gcf,filename,'pdf')     ;



