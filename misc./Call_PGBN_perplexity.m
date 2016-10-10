function Call_PGBN_perplexity(i1,i2,i3,i4,i5,i6)

K0_all = [25,50,100,200,400,600,800]
eta_all=[0.01,0.05,0.1];

switch i1
    case 1
        SystemInRuningLinux=1;
        ToBeAnalized    =   9   ;
        IsBinaryClassificaiton  =   false   ;
        dataset=5;
        % IsBinaryClassificaiton  =  true  ;
        % dataset=1;
        LoadData_GBN    ;
        dataname = '20newsTop2000'
    case 2
        load data/nips12raw_str602
        X_all = counts;
        dataname = 'NIPSTop2000'
        prepar=1:size(X_all,2);
        X_all = X_all(1:2000,:);
end

K0 = K0_all(i2) %[50,100,200,400,600]
eta0=eta_all(i3)
trial = i4
Settings.TrimStrategy  =   i5
Settings.PreTrainYes = i6

T   =   5
Settings.TrainBurnin  =   [1000,500*ones(1,T-1)]     ;       Settings.TrainCollection  =   [500,500*ones(1,T-1)];            Settings.TrainSampleSpace   =   1   ;

Settings.percentage=0.3;

dataname = [dataname,'_',num2str(Settings.percentage*100)]

Demo_PGBN_Perplexity



if 0
    serverpath = 'server/';
    %     if ~isdir(serverpath)
    %         mkdir(serverpath)
    %     end
    addpath(serverpath)
    core = 'PGBN_perplexity';
    
    submit=[];
    %for i1=[0,1,3,4]
    for i1=1:2
        RunTime = 24;
        for i2=1:7
            corejob=[core,'_30_dataset',num2str(i1),'_Kinit',num2str(i2)];
            fid = fopen([corejob,'.q'],'W');
            TaskNum=0;
            %for i1 = dataset
            for i3=[1,2,3]
                for i4 = 1:5
                    for i5=[1,4,5]
                        for i6=0
                            TaskNum = TaskNum+1;
                            fprintf(fid,'matlab  -nodisplay -nosplash -nodesktop -r "Call_%s(%d,%d,%d,%d,%d,%d)" -logfile %s_%d_%d_%d_%d_%d_%d.txt\n', core,i1,i2,i3,i4,i5,i6, core,i1,i2,i3,i4,i5,i6);
                        end
                    end
                end
            end
            fclose(fid);
            switch i2
                case {1,2}
                    TaskPerNode = 12;
                case {3,4}
                    TaskPerNode = 8;
                case {5}
                    TaskPerNode = 6;
                otherwise
                    TaskPerNode = 4;
            end
            
            filename=jobs_Lonestar(corejob,TaskNum,RunTime,TaskPerNode);
            submit=[submit,'qsub ',filename,'; '];
            
            %TaksPerNode = 12;
            %filename=jobs_Stampede(corejob,TaskNum,RunTime,TaksPerNode);
            %submit=[submit,'sbatch ',filename,'; '];
        end
        %end
        
        %
    end
    
end


%

