function [ time_total,energy_total,power_total,cost_total,path1,succ,flag_pred ] = cost_routing_paths_pred_new_pathloss( dist,R,f,t,q,p )
no=size(dist,1);

dist_0=dist;
path1=0;
p_old=p;
flag_pred=0;
flaga_1=0;
flaga_2=0;
i=1;
val1=0;
time_total=0;
energy_total=0;
power_total=0;
cost_total=0;
c_vel=3e8;
succ=true;

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
    % node with packet check if it is not in close proximity of destination
    if isinf(dist(q,p))==true % not in close proximity, so alternative route via neighbouring node is chosen
        counter=0;
        
        for j=1:1:no % j is the neighbouring node
            s=dist(j,q);
            cost(j,1)=f+((t*s)/R);
            
            PL = PL_d0+10*n_pl*log10(s/d0)+X_sigma;
            psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
            Pt = psi+PL+Pn;
            Pt = 10^(Pt*0.1);
            Et = rho+Pt*epsilon*Tb;
            power_direct(j,q) = Et/Tb;
            
            %power_direct(j,q)=(a*(s^alfa))+c;
            
            helper2(j,:)=[s,power_direct(end),cost(j,1)];
            
        end
        
        p_old=p_old(p_old>0);
        cost(p_old)=inf;
        flaga_1=0;
        
        for ii=1:1:numel(cost)
            if isinf(cost(ii))==true || isnan(cost(ii))==true
                counter=counter+1;
            end
        end
        
        p=find(cost==min(cost));
        
        if isempty(intersect(p_old,p))==false   
            flaga_1=1;
        end
        
        if counter==no %all of the element in cost is either inf or NaN
            cost_total=0;
            energy_total=0;
            time_total=0;
            power_total=0;
            path1=0;
            succ=false;
            flaga_2=1;
        end
        
        if isempty(intersect(p_old,p))==false && flaga_2==1
            succ=false;
            cost_total=0;
            power_total=0;
            time_total=0;
            energy_total=0;
            path1=0;
            flag_pred=1;
            break;
        end
        %%%%%%%%%%%%%%%%%%%%%%
        if size(p,1)==2 || size(p,2)==2
            p=p(1);
        end
        %%%%%%%%%%%%%%%%%%%%%%
        
        p_old(i+1,1)=p;
        
        if flaga_1==0
            cost_total(i,1)=cost(p);
            time_total(i,1)=((helper2(p,1))/(c_vel));
            power_total(i,1)=((helper2(p,2)));
            path1=path1+1;
        end
        
    else % the node with packet is in close proximity of destination, packet directly send to destination
        path1=path1+1;
        flag=1;
        s=dist_0(q,p);
        cost=f+((t*s)/R);
        cost_total(i,1)=cost;
        time_total(i,1)=(s/c_vel);
        
        PL = PL_d0+10*n_pl*log10(s/d0)+X_sigma;
        psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));
        Pt = psi+PL+Pn;
        Pt = 10^(Pt*0.1);
        Et = rho+Pt*epsilon*Tb;
        power_total(i,1) = Et/Tb;
        
        %power_total(i,1)=(a*(s^alfa))+c;
        
        break;
    end
    i=i+1;
end
energy=power_total.*time_total;
energy_total=sum(energy,1);
power_total=sum(power_total,1);
time_total=sum(time_total,1);
cost_total=sum(cost_total,1);
%%%%%%%%%%%%%%%%%%% small modification %%%%%%%%%%%%%%%%%%%%%%
if isinf(cost_total)==true
   cost_total=0;
   power_total=0;
   time_total=0;
   energy_total=0;
   succ=false;
end
%%%%%%%%%%%%%%%% end of small modification %%%%%%%%%%%%%%%%%%
if cost_total~=0 && path1~=0 && isinf(cost_total)==false
    val1=1;
end

if (isinf(s)==false && isnan(s)==false) && val1==1
    succ=true;
else
    succ=false;
end
if path1==0
    flag_pred=1;
end
% if succ==true
%     disp('aaa');
% end
end

