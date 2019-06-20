function [ time_total, energy_total,power_total,path1,succ,flag_pred ] = power_routing_paths_pred_new_pathloss( dist,R,f,t,q,p)
no=size(dist,1);
dist_0=dist;
energy_total=0;
time_total=0;
path1=0;
p_old=p;
p_0=p;
q_0=q;
flag_pred=0;
flaga_1=0;
flaga_2=0;
power_total=0;
val1=0;
i=1;
c_vel=3e8;

%add for path loss
PL_d0 = 55; 
n_pl = 4; 
d0 = 1; 
X_sigma = normrnd(0,4); 
PRR = 0.95;  
phi = 30; 
Pn = -145; 
rho = 50e-9; 
epsilon = 100e-12; 
Tb = 52.08e-6; 

while q~=p || i~=no
    
    if isinf(dist(q,p))==true %%need retransmission
        counter=0;
        for j=1:1:no
            
            
            r=dist(j,p);
            if j==p
                r=inf;
            end
            s=dist(j,q);
            %power_direct(j,1)=(a*(r^alfa))+c+(s*c*((((a*(alfa-1)/c)))^(1/alfa)))+(s*a*(((a*(alfa-1))/c))^((1-alfa)/alfa));
            
            PL = PL_d0+10*n_pl*log10(r/d0)+X_sigma;
            psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
            Pt = psi+PL+Pn;
            Pt = 10^(Pt*0.1);
            Et = rho+Pt*epsilon*Tb;
            power_direct(j,1)=Et/Tb;
            %power_direct(j,1)=(a*(r^alfa))+c;%
            
            if j==q
                d=dist(q,p);
                power_direct(j,1)=inf;
                
                PL = PL_d0+10*n_pl*log10(d/d0)+X_sigma;
                psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
                Pt = psi+PL+Pn;
                Pt = 10^(Pt*0.1);
                Et = rho+Pt*epsilon*Tb;
                power_dir=Et/Tb;
                %power_dir=(a*(d^alfa))+c;%
            end
            helper2(j,:)=[r,power_direct(end)];
            
        end
        
        [val,indeksik]=min(power_direct);
        if val<power_dir && isinf(power_dir)==0 && isinf(val)==0
            path1=1;
            succ=true;
            power_total=power_dir;
            time_total=(d/(c_vel));
            break;
        end
                
        
        power=power_direct;
       
        p_old=p_old(p_old>0);
        if p_old(end)==q
            break;
        end
        
        flaga_1=0;
       
        for ii=1:1:numel(power)
            if isinf(power(ii))==true || isnan(power(ii))==true
                counter=counter+1;
            end
        end
        p=find(power==min(power));
        if isempty(intersect(p_old,p))==false
            
            flaga_1=1;
            
        end
        
        if counter==no
            power_total=0;
            path1=0;
            succ=false;
            time_total=0;
            flaga_2=1;
        end
        
        if isempty(intersect(p_old,p))==false && flaga_2==1
            succ=false;
            power_total=0;
            path1=0;
            flag_pred=1;
            time_total=0;
            break;
        end
        p_old(i+1,1)=p;
        if flaga_1==0
            power_total(i,1)=power(p);
            time_total=(helper2(p,1)/(c_vel));
            path1=path1+1;
        end
            
        
        
    else %%retransmission is not needed
       
        d=dist(q,p);
        if size(dist,1)<=2
            
            PL = PL_d0+10*n_pl*log10(d/d0)+X_sigma;
            psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
            Pt = psi+PL+Pn;
            Pt = 10^(Pt*0.1);
            Et = rho+Pt*epsilon*Tb;
            power_total=Et/Tb;
            %power_total=(a*(d^alfa))+c;%
            
            time_total=(d/(c_vel));
            path=1;
            succ=true;
            break;
            
        else
            
            for u=1:1:size(dist,1)
                r=dist(u,p);             
                s=dist(u,q);
                 
                PL = PL_d0+10*n_pl*log10(r/d0)+X_sigma;
                psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
                Pt = psi+PL+Pn;
                Pt = 10^(Pt*0.1);
                Et = rho+Pt*epsilon*Tb;
                power_direct(u,1)=Et/Tb;
                %power_direct(u,1)=(a*(r^alfa))+c;%
                
                helper(u,:)=[r,power_direct(u,1)];
                if u==q
                    d=dist(q,p);
                    power_direct(u,1)=inf;
                    
                    PL = PL_d0+10*n_pl*log10(d/d0)+X_sigma;
                    psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
                    Pt = psi+PL+Pn;
                    Pt = 10^(Pt*0.1);
                    Et = rho+Pt*epsilon*Tb;
                    power_dir=Et/Tb;
                    %power_dir=(a*(d^alfa))+c;%
                    
                    d_power_dir=d;
                end
            end
            [~,indeks]=min(power_direct);
            if power_direct(indeks,1)<power_dir
                path1=path1+1;
                p=indeks;
                power_total(i,1)=power_direct(p);
                time_total(i,1)=(helper(p,1)/c_vel);
                p_old(i+1,1)=p;
               
            else
                succ=true;
                % forcing part takes part of it hopefully %
%                 power_total(i,1)=power_dir;
%                 time_total(i,1)=((d_power_dir)/c_vel);
                break;
            end
            
        end
               
    end
    
    for jj=1:1:numel(p_old)
        for kk=1:1:numel(p_old)
            dist(p_old(jj),p_old(kk))=inf;
            dist(p_old(kk),p_old(jj))=inf;
        end
    end 
    i=i+1;
end
%%%%%%%%%% forcing part %%%%%%%%%%%%%%%
d=dist_0(q_0,p(end));
%%%%%%%%% end of forcing part %%%%%%%%%
time_total(i,1)=(d/c_vel);


PL = PL_d0+10*n_pl*log10(d/d0)+X_sigma;
psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
Pt = psi+PL+Pn;
Pt = 10^(Pt*0.1);
Et = rho+Pt*epsilon*Tb;
power_total(i,1)=Et/Tb;
%power_total(i,1)=(a*(d^alfa))+c;%

path1=path1+1;

energy=power_total.*time_total;
energy_total=sum(energy,1);
power_total=sum(power_total,1);
time_total=sum(time_total,1);
%%%%%%%%% safety feature %%%%%%%%%%%
if isinf(power_total)==true
    power_total=0;
    time_total=0;
    energy_total=0;
    succ=false;
    path1=0;
end
%%%%%% end of safety feature %%%%%%%

if power_total~=0 && path1~=0 && isinf(power_total)==false
    val1=1;
end

if isinf(d)==false && val1==1
    succ=true;
else
    succ=false;
end
if path1==0
    flag_pred=1;
end

end

