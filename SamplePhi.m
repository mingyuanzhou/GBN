function  Phi = SamplePhi(WSZS,Eta,IsNoSample)
if nargin<3
    IsNoSample=false;
end
if ~IsNoSample
    Phi = randg(bsxfun(@plus,Eta,WSZS));
    temp = sum(Phi,1);
    tempdex=temp>0;
    Phi(:,tempdex) = bsxfun(@rdivide, Phi(:,tempdex),temp(tempdex));
    %Phi{t}(:,~tempdex) = (R{t}(:)+eps)/sum(R{t}+eps)*ones(1,nnz(~tempdex));
    Phi(:,~tempdex) = 0;
   % Phi = bsxfun(@rdivide, Phi,temp);
    if nnz(isnan(Phi))
        warning('Phi Nan');
        tempdex=temp>0;
        Phi(:,~tempdex) = 0;
    end
else
    Phi = bsxfun(@plus,Eta,WSZS);
    temp = sum(Phi,1);
    Phi = bsxfun(@rdivide, Phi,temp);
    if nnz(isnan(Phi))
        warning('Phi Nan');
        tempdex=temp>0;
        Phi(:,~tempdex) = 0;
    end
end
