%function AccuracyTmp = TestingInside(Theta,Para,Parellel,BlockSize)
function test_accuracy = TestingInside(Theta,Para,Parellel,BlockSize)
%Logistic regression implemented in liblinear
if ~exist('Parellel','var')
    Parellel = false;
end
if ~exist('BlockSize','var')
    BlockSize = 2;
end
if ~isfield(Para,'AddBiasTerm')
    Para.AddBiasTerm = false;
end

%% Layer One Testing

%%================= ProcessMethod End ===================
CCstart = -10;
CCend = 15;
CCstep = 1;


ThetaTmp = bsxfun(@rdivide,Theta, max(realmin,sum(Theta,1)) );

sparse_TrainX_Transpose = sparse(ThetaTmp(:,Para.train_idx)');  
TrainY = Para.Y(Para.train_idx);

sparse_TestX_Transpose = sparse(ThetaTmp(:,Para.test_idx)');  
TestY = Para.Y(Para.test_idx);

CC=2.^(CCstart : CCstep : CCend);
ModelOut=zeros(1,length(CC));
if Parellel==false
    for ij=1:length(CC)
        if Para.AddBiasTerm == false
            option = ['-s 0 -c ', num2str(CC(ij)), ' -v 5 -q '];
        else
            option =['-B 1 -s 0 -c ', num2str(CC(ij)), ' -v 5 -q ']
        end
        ModelOut(ij) = train(double(TrainY), sparse_TrainX_Transpose, option);
    end
else
    ModelOut = train_par(double(TrainY), sparse_TrainX_Transpose, CC,BlockSize);
end
[~,maxdex]=max(ModelOut);
if Para.AddBiasTerm == false
    option = ['-s 0 -c ', num2str(CC(maxdex)), ' -q'];
else
    option = ['-B 1 -s 0 -c ', num2str(CC(maxdex)), ' -q'];
end
model = train(double(TrainY), sparse_TrainX_Transpose, option);
[predicted_label, accuracy, prob_estimates] = predict(double(TestY), sparse_TestX_Transpose, model, ' -b 1');

test_accuracy = accuracy(1);
% AccAver = 0;
% for iii = 1:max(TestY)
%     ind = find(TestY==iii);
%     AccAver = AccAver + 1/max(TestY) * sum(predicted_label(ind)==iii) / length(ind);
% end
% AccuracyTmp.Accuracy_ThetaFreq_AccAver = AccAver;
% AccuracyTmp.Accuracy_ThetaFreq_predicted_label = predicted_label;
%AccuracyTmp.Accuracy_ThetaFreq_accuracy = accuracy;
%AccuracyTmp.Accuracy_ThetaFreq_prob_estimates = prob_estimates;





