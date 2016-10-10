function Call_PGBN(i1,i2,i3,i4,i5,i6)

ToBeAnalized    =   9   ;


K0_all = [5,10,25,50,100,200,400,600,800,1000]
eta_all=[0.01,0.05,0.1];

dataset = i1 %0,1,3,4

if dataset==5
    K0_all = [16,32,64,128,256,512,1024]
end

K0 = K0_all(i2) %[50,100,200,400,600]
eta0=eta_all(i3)
trial = i4
Settings.TrimStrategy  =   i5
Settings.PreTrainYes =i6

if i1==0 || i1==5
    T   =   5
    Settings.TrainBurnin  =   [1000,1000*ones(1,T-1)]     ;       Settings.TrainCollection  =   [500,500*ones(1,T-1)];            Settings.TrainSampleSpace   =   1   ;
else
    T   =   8
    Settings.TrainBurnin  =   [1000,1000*ones(1,T-1)]     ;       Settings.TrainCollection  =   [1000,1000*ones(1,T-1)];            Settings.TrainSampleSpace   =   1   ;
end

if i1==0 && K0>500
    T   =   5
    Settings.TrainBurnin  =   [1000,500*ones(1,T-1)]     ;       Settings.TrainCollection  =   [500,500*ones(1,T-1)];            Settings.TrainSampleSpace   =   1   ;
end


Demo_GBN_FeatureExtraction_v1






if 0
    
    core = 'PGBN';
    
    submit=[];
    %for i1=[0,1,3,4]
    for i1=[0,5]
        
        K0_all = [5,10,25,50,100,200,400,600,800,1000]
        
        if i1==5
            K0_all = [16,32,64,128,256,512,1024]
        end
        
        
        RunTime = 24;
        I2max=9;
        if i1==0
            I2dex = 4:7;
        else
            I2dex = 2:5;
        end
        for i2=I2dex
            K=K0_all(i2);
            corejob=[core,'_dataset',num2str(i1),'_Kinit',num2str(i2)];
            fid = fopen([corejob,'.q'],'W');
            TaskNum=0;
            %for i1 = dataset
            for i3=2
                for i4 = 1:12
                    for i5=5
                        for i6=0
                            TaskNum = TaskNum+1;
                            fprintf(fid,'matlab  -nodisplay -nosplash -nodesktop -r "Call_%s(%d,%d,%d,%d,%d,%d)" -logfile %s_%d_%d_%d_%d_%d_%d.txt\n', core,i1,i2,i3,i4,i5,i6, core,i1,i2,i3,i4,i5,i6);
                        end
                    end
                end
            end
            fclose(fid);
            
            TaskPerNode = 12;
            if K>=200
                TaskPerNode = 6;
            end
            %
            %             if K>200
            %                 TaksPerNode = 4;
            %             end
            
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

if 0
    serverpath = 'server/';
    %     if ~isdir(serverpath)
    %         mkdir(serverpath)
    %     end
    addpath(serverpath)
    core = 'PGBN';
    
    submit=[];
    %for i1=[0,1,3,4]
    for i1=[3,4]
        if i1==0
            RunTime = 48;
            I2max=9;
        else
            RunTime = 24;
            I2max=9;
        end
        for i2=1:I2max
            corejob=[core,'_dataset',num2str(i1),'_Kinit',num2str(i2)];
            fid = fopen([corejob,'.q'],'W');
            TaskNum=0;
            %for i1 = dataset
            for i3=[1,2,3]
                for i4 = 1:12
                    for i5=[1,4,5]
                        for i6=0
                            TaskNum = TaskNum+1;
                            fprintf(fid,'matlab  -nodisplay -nosplash -nodesktop -r "Call_%s(%d,%d,%d,%d,%d,%d)" -logfile %s_%d_%d_%d_%d_%d_%d.txt\n', core,i1,i2,i3,i4,i5,i6, core,i1,i2,i3,i4,i5,i6);
                        end
                    end
                end
            end
            fclose(fid);
            
            filename=jobs_Lonestar(corejob,TaskNum,RunTime);
            submit=[submit,'qsub ',filename,'; '];
            
            %TaksPerNode = 12;
            %filename=jobs_Stampede(corejob,TaskNum,RunTime,TaksPerNode);
            %submit=[submit,'sbatch ',filename,'; '];
        end
        %end
        
        %
    end
    
end


if 0
    serverpath = 'server/';
    %     if ~isdir(serverpath)
    %         mkdir(serverpath)
    %     end
    addpath(serverpath)
    core = 'PGBN';
    
    submit=[];
    for i1=0
        %for i1=[1,3,4]
        if i1==0
            RunTime = 48;
            I2max=9;
            T   =   5  ;
        else
            RunTime = 24;
            I2max=9;
            T   =   8  ;
        end
        for i2=1:I2max
            corejob=[core,'_dataset',num2str(i1),'_Kinit',num2str(i2)];
            fid = fopen([corejob,'.q'],'W');
            TaskNum=0;
            %for i1 = dataset
            for i3=[1,2,3]
                for i4 = 1:12
                    for i5=[1,4,5]
                        for i6=0
                            %:1
                            TaskNum = TaskNum+1;
                            fprintf(fid,'matlab  -nodisplay -nosplash -nodesktop -r "Call_%s(%d,%d,%d,%d,%d,%d)" -logfile %s_%d_%d_%d_%d_%d_%d.txt\n', core,i1,i2,i3,i4,i5,i6, core,i1,i2,i3,i4,i5,i6);
                        end
                    end
                end
            end
            fclose(fid);
            %filename=jobs_Lonestar(corejob,TaskNum,RunTime);
            TaksPerNode = 12;
            filename=jobs_Stampede(corejob,TaskNum,RunTime,TaksPerNode);
            submit=[submit,'sbatch ',filename,'; '];
        end
        %end
        
        %
    end
    
end
%



if 0
    serverpath = 'server/';
    %     if ~isdir(serverpath)
    %         mkdir(serverpath)
    %     end
    addpath(serverpath)
    core = 'PGBN';
    
    submit=[];
    %for i1=[0,1,3,4]
    for i1=[1]
        if i1==0
            RunTime = 48;
            I2max=9;
        else
            RunTime = 24;
            I2max=9;
        end
        for i2=1:I2max
            corejob=[core,'_dataset',num2str(i1),'_Kinit',num2str(i2)];
            fid = fopen([corejob,'.q'],'W');
            TaskNum=0;
            %for i1 = dataset
            for i3=[1,2,3]
                for i4 = 1:12
                    for i5=[4,5]
                        for i6=0
                            TaskNum = TaskNum+1;
                            fprintf(fid,'matlab  -nodisplay -nosplash -nodesktop -r "Call_%s(%d,%d,%d,%d,%d,%d)" -logfile %s_%d_%d_%d_%d_%d_%d.txt\n', core,i1,i2,i3,i4,i5,i6, core,i1,i2,i3,i4,i5,i6);
                        end
                    end
                end
            end
            fclose(fid);
            
            filename=jobs_Lonestar(corejob,TaskNum,RunTime);
            submit=[submit,'qsub ',filename,'; '];
            
            %TaksPerNode = 12;
            %filename=jobs_Stampede(corejob,TaskNum,RunTime,TaksPerNode);
            %submit=[submit,'sbatch ',filename,'; '];
        end
        %end
        
        %
    end
    
end





if 0
    serverpath = 'server/';
    %     if ~isdir(serverpath)
    %         mkdir(serverpath)
    %     end
    addpath(serverpath)
    core = 'PGBN';
    
    submit=[];
    %for i1=[0,1,3,4]
    for i1=5
        
        RunTime = 48;
        I2max=7;
        
        for i2=1:I2max
            corejob=[core,'_dataset',num2str(i1),'_Kinit',num2str(i2)];
            fid = fopen([corejob,'.q'],'W');
            TaskNum=0;
            %for i1 = dataset
            for i3=[1,2,3]
                for i4 = 1:12
                    for i5=[1,4,5]
                        for i6=0
                            TaskNum = TaskNum+1;
                            fprintf(fid,'matlab  -nodisplay -nosplash -nodesktop -r "Call_%s(%d,%d,%d,%d,%d,%d)" -logfile %s_%d_%d_%d_%d_%d_%d.txt\n', core,i1,i2,i3,i4,i5,i6, core,i1,i2,i3,i4,i5,i6);
                        end
                    end
                end
            end
            fclose(fid);
            
            %filename=jobs_Lonestar(corejob,TaskNum,RunTime);
            %submit=[submit,'qsub ',filename,'; '];
            
            TaksPerNode = 12;
            filename=jobs_Stampede(corejob,TaskNum,RunTime,TaksPerNode);
            submit=[submit,'sbatch ',filename,'; '];
        end
        %end
        
        %
    end
    
end




if 0
    serverpath = 'server/';
    %     if ~isdir(serverpath)
    %         mkdir(serverpath)
    %     end
    addpath(serverpath)
    core = 'PGBN';
    
    submit=[];
    %for i1=[0,1,3,4]
    for i1=0
        
        RunTime = 48;
        I2max=7;
        
        for i2=8:9
            
            %for i1 = dataset
            for i3=[1,2,3]
                corejob=[core,'_dataset',num2str(i1),'_Kinit',num2str(i2),'_eta',num2str(i3)];
                fid = fopen([corejob,'.q'],'W');
                TaskNum=0;
                for i4 = 1:12
                    for i5=[1,4,5]
                        for i6=0
                            TaskNum = TaskNum+1;
                            fprintf(fid,'matlab  -nodisplay -nosplash -nodesktop -r "Call_%s(%d,%d,%d,%d,%d,%d)" -logfile %s_%d_%d_%d_%d_%d_%d.txt\n', core,i1,i2,i3,i4,i5,i6, core,i1,i2,i3,i4,i5,i6);
                        end
                    end
                end
                fclose(fid);
                
                %filename=jobs_Lonestar(corejob,TaskNum,RunTime);
                %submit=[submit,'qsub ',filename,'; '];
                
                TaksPerNode = 6;
                filename=jobs_Stampede(corejob,TaskNum,RunTime,TaksPerNode);
                submit=[submit,'sbatch ',filename,'; '];
            end
            
        end
        %end
        
        %
    end
    
end
