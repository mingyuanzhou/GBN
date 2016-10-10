function Call_PRG_GBN(i1,i2,i3,i4)
if nargin<4
    i4=1;
end


K0_all = [25,50,100,200,400,600,800,1000]
eta_all=[0.01,0.05,0.1];
K0 = K0_all(i1) %[50,100,200,400,600]
eta0=eta_all(i2)
trial = i3
Settings.TrimStrategy  =   5
Settings.PreTrainYes = 0

T   =   5

%eta0 = [eta0,0.1*ones(1,T-1)];
Settings.TrainBurnin  =   [750,750*ones(1,T-1)]     ;
Settings.TrainCollection  =  [250,250*ones(1,T-1)];
Settings.TrainSampleSpace   =   1   ;

Settings.TestBurnin  =   500     ;
Settings.TestCollection  =   500     ;
Settings.TestSampleSpace   =   1   ;


Settings.ParallelProcessing      =   i4

POOL.NumWorkers = 24;
Settings.NumBlockParallel   =   POOL.NumWorkers*1;

% if Settings.ParallelProcessing==1
%     %if K0<=400
%         matlabpool;
%     %else
%     %    matlabpool(12);
%     %end
% end

Demo_PRG_GBN_FeatureExtraction






if 0
    
    
    core = 'PRG_GBN';
    
    submit=[];
    for i1=6:-1:3
        %for i1=[1,3,4]
        
        RunTime = 48;
        
        
        for i2=[1,2]
            %for i1 = dataset
            for i3=1:5
                TaskNum=0;
                corejob=[core,'_K0',num2str(i1),'_eta0',num2str(i2),'_trial',num2str(i3)];
                fid = fopen([corejob,'.q'],'W');
                
                i4=1;
                TaskNum = TaskNum+1;
                fprintf(fid,'matlab  -nodisplay -nosplash -nodesktop -r "Call_%s(%d,%d,%d,%d)" -logfile %s_%d_%d_%d_%d.txt\n', core,i1,i2,i3,i4, core,i1,i2,i3,i4);
                fclose(fid);
                %filename=jobs_Lonestar(corejob,TaskNum,RunTime);
                TaskPerNode = 1;
                filename=jobs_Stampede(corejob,TaskNum,RunTime,TaskPerNode);
                submit=[submit,'sbatch ',filename,'; '];
                %filename=jobs_Lonestar(corejob,TaskNum,RunTime,TaskPerNode)
                
                %submit=[submit,'qsub ',filename,'; '];
            end
            
        end
        %end
        
        %
    end
    
    for i1=2:-1:1
        %for i1=[1,3,4]
        
        RunTime = 48;
        
        for i2=[1,2]
            corejob=[core,'_K0',num2str(i1),'_eta0',num2str(i2)];
            fid = fopen([corejob,'.q'],'W');
            TaskNum=0;
            %for i1 = dataset
            for i3=1:5
                i4=0;
                TaskNum = TaskNum+1;
                fprintf(fid,'matlab  -nodisplay -nosplash -nodesktop -r "Call_%s(%d,%d,%d,%d)" -logfile %s_%d_%d_%d_%d.txt\n', core,i1,i2,i3,i4, core,i1,i2,i3,i4);
            end
            fclose(fid);
            %filename=jobs_Lonestar(corejob,TaskNum,RunTime);
            TaskPerNode = 5;
            filename=jobs_Stampede(corejob,TaskNum,RunTime,TaskPerNode);
            submit=[submit,'sbatch ',filename,'; '];
        end
        %end
        
        %
    end
    
end
%



if 0
    
    
    core = 'PRG_GBN';
    
    submit=[];
    for i1=5:-1:4
        %for i1=[1,3,4]
        
        RunTime = 24;
        
        
        for i2=1
            
            %for i1 = dataset
            for i3=1:5
                TaskNum=0;
                corejob=[core,'Lone_K0',num2str(i1),'_eta0',num2str(i2),'_trial',num2str(i3)];
                fid = fopen([corejob,'.q'],'W');
                
                i4=1;
                TaskNum = TaskNum+1;
                fprintf(fid,'matlab  -nodisplay -nosplash -nodesktop -r "Call_%s(%d,%d,%d,%d)" -logfile %s_%d_%d_%d_%d.txt\n', core,i1,i2,i3,i4, core,i1,i2,i3,i4);
                fclose(fid);
                %filename=jobs_Lonestar(corejob,TaskNum,RunTime);
                TaskPerNode = 1;
                %filename=jobs_Stampede(corejob,TaskNum,RunTime,TaskPerNode);
                % submit=[submit,'sbatch ',filename,'; '];
                filename=jobs_Lonestar(corejob,TaskNum,RunTime,TaskPerNode)
                
                submit=[submit,'qsub ',filename,'; '];
            end
            
            
        end
        %end
        
        %
    end
    
    for i1=3:-1:1
        %for i1=[1,3,4]
        
        RunTime = 24;
        
        for i2=1
            corejob=[core,'Lone_K0',num2str(i1),'_eta0',num2str(i2)];
            fid = fopen([corejob,'.q'],'W');
            TaskNum=0;
            %for i1 = dataset
            for i3=1:5
                i4=0;
                TaskNum = TaskNum+1;
                fprintf(fid,'matlab  -nodisplay -nosplash -nodesktop -r "Call_%s(%d,%d,%d,%d)" -logfile %s_%d_%d_%d_%d.txt\n', core,i1,i2,i3,i4, core,i1,i2,i3,i4);
            end
            fclose(fid);
            %filename=jobs_Stampede(corejob,TaskNum,RunTime,TaskPerNode);
            %submit=[submit,'sbatch ',filename,'; '];
            TaskPerNode = 6;
            filename=jobs_Lonestar(corejob,TaskNum,RunTime,TaskPerNode)
            
            submit=[submit,'qsub ',filename,'; '];
        end
        %end
        
        %
    end
    
end






if 0
    
    
    core = 'PRG_GBN';
    
    submit=[];
    for i1=[3,6]
        %for i1=[1,3,4]
        
        RunTime = 48;
        
        
        for i2=[2,3]
            %for i1 = dataset
            for i3=1:5
                TaskNum=0;
                corejob=[core,'_K0',num2str(i1),'_eta0',num2str(i2),'_trial',num2str(i3)];
                fid = fopen([corejob,'.q'],'W');
                
                i4=1;
                TaskNum = TaskNum+1;
                fprintf(fid,'matlab  -nodisplay -nosplash -nodesktop -r "Call_%s(%d,%d,%d,%d)" -logfile %s_%d_%d_%d_%d.txt\n', core,i1,i2,i3,i4, core,i1,i2,i3,i4);
                fclose(fid);
                %filename=jobs_Lonestar(corejob,TaskNum,RunTime);
                TaskPerNode = 1;
                filename=jobs_Stampede(corejob,TaskNum,RunTime,TaskPerNode);
                submit=[submit,'sbatch ',filename,'; '];
                %filename=jobs_Lonestar(corejob,TaskNum,RunTime,TaskPerNode)
                
                %submit=[submit,'qsub ',filename,'; '];
            end
            
        end
        %end
        
        %
    end
   
    
end
%
