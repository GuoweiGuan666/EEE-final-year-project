% function [ cost_total,path1,succ,flag_pred ] = power_cost_routing_0_paths_pred( dist,R,f,f_prim,t,q,p,a,alfa,c )
function [time_total, energy_total, power_total,power_cost_total,path1,succ,flag_pred ] = power_cost_routing_0_paths_pred_new_pathloss(dist,R,f,f_prim,t,q,p)
no=size(dist,1);

dist_0=dist;
path1=0;
p_old=p;
q_0=q;
p_0=p;
flag_pred=0;
flaga_1=0;
flaga_2=0;
power_cost_total=0;
i=1;
vall=0;
c_vel=3e8;
time_total=0;
power_total=0;
energy_total=0;
power_cost_total=0;

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


while p~=q || i~=no
%         if isempty(intersect(p_old(1:end-1),p))==false
%         break;
%     end
    if isinf(dist(q,p))==true
        counter=0;
        for j=1:1:no
            
            r=dist(j,p);
            if j==p
                r=inf;
            end
            s=dist(j,q);
            
            %every parameters at here is on r%
            PL = PL_d0+10*n_pl*log10(r/d0)+X_sigma;
            psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
            Pt = psi+PL+Pn;
            Pt = 10^(Pt*0.1);
            Et = rho+Pt*epsilon*Tb;
            power_cost(j,1)=f*(Et/Tb); %Et/Tb is the direct transmission power replace (a*(r^alfa))+c 
            
            %power_cost(j,1)=f*((a*(r^alfa))+c);% %on r%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            power_direct(j,1)= Et/Tb;
            %power_direct(j,1)=((a*(r^alfa))+c);% %on r%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %v(s) of path loss is on s!!!, need to recalculate psi and Pt on s%
            PL = PL_d0+10*n_pl*log10(s/d0)+X_sigma;
            psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
            Pt = psi+PL+Pn;
            Pt = 10^(Pt*0.1);
            
            %v(s) of path loss is on s!!!% %the equation follows the f_prim* is the power transmission via neighouring node of path loss version designed by Guowei Guan%
            power_cost_tot(j,1)=power_cost(j,1)+f_prim*(s*d0^(-1)*(10^(0.1*n_pl^(-1)*(((10*n_pl)/log(10))-(rho/(epsilon*Tb))-psi-PL_d0-X_sigma-Pn)))^(-1)*(rho+epsilon*Tb*Pt));
            %v(s) of old version is on s!!!%
            %power_cost_tot(j,1)=power_cost(j,1)+f_prim*((s*c*((((a*(alfa-1)/c)))^(1/alfa)))+(s*a*(((a*(alfa-1))/c))^((1-alfa)/alfa)));%
            
            
            if isinf(power_cost_tot(j,1))==true
                power_cost_tot(j,1)=power_cost(j,1);
            end
            if j==q
                d=dist(q,p);
                power_cost_tot(j,1)=inf;
                
                %every parameters at here is on d%
                PL = PL_d0+10*n_pl*log10(d/d0)+X_sigma;
                psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
                Pt = psi+PL+Pn;
                Pt = 10^(Pt*0.1);
                Et = rho+Pt*epsilon*Tb;
                power_cost_dir=f*(Et/Tb); %Et/Tb is direct transmission power on d to replace (a*(d^alfa))+c %
                %on d%
                %power_cost_dir=f*((a*(d^alfa))+c);%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                power_direct(j,1)=inf;
                
                power_dir=Et/Tb;
                %power_dir=(a*(d^alfa))+c;%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            helper2(j,:)=[r,power_direct(end)];
            
        end
        
        [val,~]=min(power_cost_tot);
        if val<power_cost_dir && isinf(power_cost_dir)==0 && isinf(val)==0
            path1=1;
            succ=true;
            power_total=power_dir;
            time_total=(d/(c_vel));
            break;
        end
        
        p_old=p_old(p_old>0);
        if p_old(end)==q
            break;
        end
        
        flaga_1=0;
        
        for ii=1:1:numel(power_cost_tot)
            if isinf(power_cost_tot(ii))==true || isnan(power_cost_tot(ii))==true
                counter=counter+1;
            end
        end
        
        p=find(power_cost_tot==min(power_cost_tot));
        
        if isempty(intersect(p_old,p))==false
            
            flaga_1=1;
            break;
            
        end
        
        if counter==no
            power_cost_total=0;
            power_total=0;
            time_total=0;
            path1=0;
            succ=false;
            flaga_2=1;
        end
        
        if isempty(intersect(p_old,p))==false && flaga_2==1
            succ=false;
            power_cost_total=0;
            power_total=0;
            time_total=0;
            path1=0;
            flag_pred=1;
            break;
        end
        
        p_old(i+1,1)=p;
        
        if flaga_1==0
            power_cost_total(i,1)=power_cost_tot(p);
            power_total(i,1)=helper2(p,2);
            time_total(i,1)=((helper2(p,1))/(c_vel));
            path1=path1+1;
        end
        
    else
        d=dist(q,p);
        if size(dist,1)<=2
            
            %every parameters at here is on d%
            PL = PL_d0+10*n_pl*log10(d/d0)+X_sigma;
            psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
            Pt = psi+PL+Pn;
            Pt = 10^(Pt*0.1);
            Et = rho+Pt*epsilon*Tb;
            power_cost_total=f*(Et/Tb);%Et/Tb is direct transmission power on d to replace (a*(d^alfa))+c %
            %on d%
            %power_cost_total=f*((a*(d^alfa))+c);%
            
            power_total=Et/Tb;
            %power_total=((a*(d^alfa))+c);%
            
            time_total=(d/(c_vel));
            path1=1;
            succ=true;
            break;
        else
            
            for u=1:1:size(dist,1)
                r=dist(u,p);
                s=dist(u,q);
                
                %every parameters at here is on r%
                PL = PL_d0+10*n_pl*log10(r/d0)+X_sigma;
                psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
                Pt = psi+PL+Pn;
                Pt = 10^(Pt*0.1);
                Et = rho+Pt*epsilon*Tb;
                power_cost(u,1)=f*(Et/Tb);%Et/Tb is the direct transmission power replace (a*(r^alfa))+c
                %on r%
                %power_cost(u,1)=f*((a*(r^alfa))+c);%
                
                
                
                %v(s) of path loss is on s!!!, need to recalculate psi and Pt on s%
                PL = PL_d0+10*n_pl*log10(s/d0)+X_sigma;
                psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
                Pt = psi+PL+Pn;
                Pt = 10^(Pt*0.1);
                
                %v(s) of path loss is on s!!!% %the equation follows the f_prim* is the power transmission via neighouring node of path loss version designed by Guowei Guan%
                power_cost_tot(u,1)=power_cost(u,1)+f_prim*(s*d0^(-1)*(10^(0.1*n_pl^(-1)*(((10*n_pl)/log(10))-(rho/(epsilon*Tb))-psi-PL_d0-X_sigma-Pn)))^(-1)*(rho+epsilon*Tb*Pt));
                %v(s) of old version is on s!!!%
                %power_cost_tot(u,1)=power_cost(u,1)+f_prim*((s*c*((((a*(alfa-1)/c)))^(1/alfa)))+(s*a*(((a*(alfa-1))/c))^((1-alfa)/alfa)));%
                
                
                if isinf(power_cost_tot(u,1))==true
                    power_cost_tot(u,1)=power_cost(u,1);
                end
                
                %every parameters at here is on r%
                PL = PL_d0+10*n_pl*log10(r/d0)+X_sigma;
                psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
                Pt = psi+PL+Pn;
                Pt = 10^(Pt*0.1);
                Et = rho+Pt*epsilon*Tb;
                power_direct(u,1)=Et/Tb;
                %on r%
                %power_direct(u,1)=(a*(r^alfa))+c;%
                helper(u,:)=[r,power_direct(u,1)];
                if u==q
                    d=dist(q,p);
                    power_cost_tot(u,1)=inf;
                    power_direct(u,1)=inf;
                    
                    %every parameters at here is on d%
                    PL = PL_d0+10*n_pl*log10(d/d0)+X_sigma;
                    psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
                    Pt = psi+PL+Pn;
                    Pt = 10^(Pt*0.1);
                    Et = rho+Pt*epsilon*Tb;
                    power_dir=Et/Tb;
                    %on d%
                    %power_dir=(a*(d^alfa))+c;%
                    d_power_dir=d;
                end
            end
            [~,indeks]=min(power_cost_tot);
            
            
            path1=path1+1;
            p=indeks;
            power_total(i,1)=power_direct(p);
            time_total(i,1)=(helper(p,1)/c_vel);
            power_cost_total(i,1)=power_cost_tot(p);
            p_old(i+1,1)=p;
            
            
        end
    end
    
    for jj=1:1:numel(p_old)
        for kk=1:1:numel(p_old)
            dist(p_old(jj),p_old(kk))=inf;
            dist(p_old(kk),p_old(jj))=inf;
        end
    end
    
%%%%%%%%%%%%%
    if isempty(intersect(p_old(1:end-1),p))==false
        break;
    end
%%%%%%%%%%%%
    i=i+1;
end

d=dist_0(q_0,p(end));
if isnan(dist_0(q_0,p(end)))==true
   power_cost_total(i,1)=0;
else
   %every parameters at here is on d%
   PL = PL_d0+10*n_pl*log10(d/d0)+X_sigma;
   psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
   Pt = psi+PL+Pn;
   Pt = 10^(Pt*0.1);
   Et = rho+Pt*epsilon*Tb;
   power_cost_total(i,1)=Et/Tb;
   %on d% 
   %power_cost_total(i,1)=(a*(d^alfa))+c;%
end

if p(end)==q
    cost_temp=power_cost_total(end);
    if nnz(isinf(cost_temp))>0
       power_cost_total(isinf(power_cost_total))=0; 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_total(i,1)=(d/c_vel);

%every parameters at here is on d%
PL = PL_d0+10*n_pl*log10(d/d0)+X_sigma;
psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
Pt = psi+PL+Pn;
Pt = 10^(Pt*0.1);
Et = rho+Pt*epsilon*Tb;
power_total(i,1)=Et/Tb;
%on d%
%power_total(i,1)=(a*(d^alfa))+c;%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path1=path1+1;
energy=power_total.*time_total;
energy_total=sum(energy,1);
power_total=sum(power_total,1);
time_total=sum(time_total,1);

if isinf(power_cost_total(i,1))==true
    power_cost_total=0;
    succ=false;
    path1=0;
end

power_cost_total=sum(power_cost_total,1);
if isinf(power_cost_total)==true
    power_total=0;
    energy_total=0;
    time_total=0;
    power_cost_total=0;
end

if power_cost_total~=0 && path1~=0 && isinf(power_cost_total)==false
    vall=1;
end
if isinf(d)==false && vall==1
    succ=true;
else 
    succ=false;
    path1=0;
    power_cost_total=0;
end

if path1==0
    flag_pred=1;
end
end

