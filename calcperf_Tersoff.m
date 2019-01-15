function [perf,Ex,Vhat]=calcperf_Tersoff(net,rs,V,Q,param,clust_size,type)

perf=0.0;

%looping over no. of data pts. Q
for iQ=1:1:Q
    %converting rs to my format
    for i=1:1:clust_size
         for j=1:1:clust_size
             if i~=j
                 r(i,j)=rs(i,j,iQ);
             end
         end
    end
    
    %calculating Tersoff energy for iQ configuration
    Vhat_iQ=0.0;
    for i=1:1:clust_size
        for j=1:1:clust_size
            if i==j
                continue;
            end
            bi=net.b{1};  %sim(net,r(i,j));
            w=net.IW{1};
            A=w(1)*r(i,j)+bi(1);
            B=w(2)*r(i,j)+bi(2);%param(type(i),type(j),type(j),2);
            lambda1=param(type(i),type(j),type(j),3);
            lambda2=param(type(i),type(j),type(j),4);
    %         lambda3=param(type(i),type(j),type(j),5);
    %         alpha=param(type(i),type(j),type(j),6);
            beta=param(type(i),type(j),type(j),7);
            eta=param(type(i),type(j),type(j),8);
    %         c=param(type(i),type(j),type(j),9);
    %         d=param(type(i),type(j),type(j),10);
    %         h=param(type(i),type(j),type(j),11);
            R=param(type(i),type(j),type(j),12);
            D=param(type(i),type(j),type(j),13);        
            [fC]=fc(i,j,r,R,D);
            [fR]=fr(i,j,r,A,lambda1);
            [fA]=fa(i,j,r,B,lambda2);
            zetaij=0.0;
            for k=1:1:clust_size
                if i==j || i==k
                    continue;
                end
    %             A=param(type(i),type(j),type(k),1);
    %             B=param(type(i),type(j),type(k),2);
    %             lambda1=param(type(i),type(j),type(k),3);
    %             lambda2=param(type(i),type(j),type(k),4);
    %             lambda3=param(type(i),type(j),type(k),5);
    %             alpha=param(type(i),type(j),type(k),6);
    %             beta=param(type(i),type(j),type(k),7);
    %             eta=param(type(i),type(j),type(k),8);
                c=param(type(i),type(j),type(k),9);
                d=param(type(i),type(j),type(k),10);
                h=param(type(i),type(j),type(k),11);
                R=param(type(i),type(j),type(k),12);
                D=param(type(i),type(j),type(k),13); 
                [zetaij] = zetaij+zeta(i, j, k, r, c, d, h, R,D);		
            end

            beta_eta=beta^eta;
            zeta_eta=zetaij^eta;
            [bij]=b(zetaij,beta_eta, zeta_eta, eta);
            Vhat_iQ=Vhat_iQ+fC*(fR+bij*fA);
        end
    end
    %error vector
    Ex(iQ,1)=V(iQ)-0.5*Vhat_iQ;
    Vhat(iQ)=0.5*Vhat_iQ;
    %SSE
    perf=perf+Ex(iQ,1)*Ex(iQ,1);
   
end %iQ