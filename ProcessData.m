function [Est,Data]=ProcessData(c,s,reAttach,p,Href,SiteName,Method)

%% process centerline
z1 = utmzone([64.78 -148]); %6W
[ellipsoid,estr] = utmgeoid(z1);
utmstruct = defaultm('utm'); 
utmstruct.zone = z1; 
utmstruct.geoid = ellipsoid; 
utmstruct = defaultm(utmstruct);

cutm=c;
for i=1:length(cutm),
    [cutm(i).X,cutm(i).Y] = mfwdtran(utmstruct,c(i).Y,c(i).X);
end

%%
Nc=length(cutm.X);
FD=0;
for i=2:Nc,
    FD(i)=FD(i-1)+ sqrt( (cutm.X(i)-cutm.X(i-1)).^2 + (cutm.Y(i)-cutm.Y(i-1)).^2 );
end
FDkm=FD./1000;

%% process height observation
NumberScenes=length(s);
for i=1:NumberScenes,
    
    disp(['Processing scene ' num2str(i) '/' num2str(NumberScenes) '.'])
    
    Data{i}.lat=s{i}(:,1);
    Data{i}.lon=s{i}(:,2);
    Data{i}.h=s{i}(:,3);
    [Data{i}.X,Data{i}.Y]=mfwdtran(utmstruct,Data{i}.lat,Data{i}.lon);
    
    %% attach each height measurement to centerline
    Data{i}.N=length(Data{i}.h);
    if reAttach,
        for j=1:Data{i}.N,
            dist=sqrt( (Data{i}.Y(j)-cutm.Y).^2 + (Data{i}.X(j) - cutm.X).^2);
            [distmin(j),k(j)]=min(dist);    
            FDh(j)=FDkm(k(j));
        end
        name=['Elevations/' SiteName 'distsi' num2str(i) '.mat'];
        save(name,'distmin', 'k', 'FDh');
    else
        name=['Elevations/' SiteName 'distsi' num2str(i) '.mat'];
        load(name)
    end

    %% estimate heights for each bin
    Data{i}.iuse=distmin<p.distmax;    clear distmin;
    Data{i}.FDh=FDh; clear FDh
    Est{i}.Hhat=nan(1,p.N);
    Est{i}.HhatStd=nan(1,p.N);
        
    for j=1:p.N,
        if p.x(j)==59,
            stop=1;
        end
        dists=abs(Data{i}.FDh-p.x(j));
        m=dists<p.dx/2;
        if ~isnan(Href{i}(j)),
            iref= abs(Data{i}.h'-Href{i}(j))<p.HrefCut;
        else
            iref=ones(size(Data{i}.FDh));
        end
        
        nhat(j)=sum(Data{i}.iuse&m&iref);
        
        if j==p.N-4,
            stop=1;
        end
        
        if nhat(j)>p.ncut,    
            Est{i}.Hhat(j)=prctile(Data{i}.h(Data{i}.iuse&m&iref),p.pct);
            Est{i}.HhatStd(j)=std(Data{i}.h(Data{i}.iuse&m&iref));
        else
            Est{i}.Hhat(j)=nan;
            Est{i}.HhatStd(j)=nan;
        end
    end
        
    Est{i}.Use=Est{i}.HhatStd<p.StdMax & ~isnan(Est{i}.Hhat);
    iObs=find(Est{i}.Use);
    Nobs=sum(Est{i}.Use);
        
    %% Constrain to go downhill    
    if strcmpi(Method,'LP'),            
        Est{i}.Hc=ConstrainHeights(p.N,Nobs,Est{i}.Hhat(Est{i}.Use)',iObs);         
        xMax=max(p.x(Est{i}.Use));
        xMin=min(p.x(Est{i}.Use));
        Est{i}.iSolUse=p.x > xMin & p.x<xMax;
    elseif strcmpi(Method,'SLM')     
        Est{i}.InBounds=false(size(Est{i}.Use));
        Est{i}.InBounds(find(Est{i}.Use,1,'first'):find(Est{i}.Use,1,'last'))=true;
        x=p.x(Est{i}.Use); y=Est{i}.Hhat(Est{i}.Use)';
        Ltot=range(x);
        ReachMin=6;
        nReach=floor(Ltot/ReachMin);
        [Est{i}.slm,Est{i}.xp,Est{i}.yp] = slmengine(x,y,'increasing','on','result','slm','knots',nReach+1,'degree',3);
        Est{i}.Slope=slmeval(Est{i}.xp,Est{i}.slm,1).*100;
        for j=1:length(p.x),
            [~,k]=min( abs(Est{i}.xp-p.x(j)) );
            if Est{i}.InBounds(j) && abs(p.x(j)-x(1)) > 1 && abs(p.x(j)-x(end)) > 1,
                Est{i}.Hc(j)=Est{i}.yp(k);
                Est{i}.iSolUse(j)=true;
            else
                Est{i}.Hc(j)=nan;
                Est{i}.iSolUse(j)=false;                
            end            
        end    
        Est{i}.res=Est{i}.Hc(Est{i}.Use)-y';
        Est{i}.stdres=nanstd(Est{i}.res);        
    end
    
end %loop over 

return

