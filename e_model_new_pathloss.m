
function [n_optimum,distance,power,Pr] = e_model_new_pathloss(x,y,r) %modified HCB model to account for path loss is used in here
k=size(x,1);



%for step 1
PL_d0 = 55; %path loss at reference distance in dB
n_pl = 4; %path loss exponent, rate at which the signal decays
d0 = 1; %reference distance
X_sigma = normrnd(0,4); % random number from normal distribution (0,4)

%for step 2
PRR = 0.95;  %given target PRR value should be set in the range(0.9,1)
phi = 30; %packet size in bytes

%for step 3
Pn = -145; %noise floor

%for step 4
rho = 50e-9; %enerfy dissipation on circuit
epsilon = 100e-12; %transmitter efficiency
%alpha=2;
Tb = 52.08e-6; %duration of one bit

%for step 5
%Pt and PL calculated from step 3 and step 1





const=0.031; % n_optimum = 0.031d
iterator = 0;

x_1=zeros(k);
y_1=zeros(k);

while iterator~=k
    x_1(:,iterator+1)=x-x(iterator+1);
    y_1(:,iterator+1)=y-y(iterator+1);
    iterator=iterator+1;
end

x_1=x_1.^2;
y_1=y_1.^2;
% distance=unique(x_1+y_1);
distance = x_1+y_1;

% distance(distance==0)=[];
distance = sqrt(distance); % this vector distance contains the distance dij of all node pair in network



%step 1: obtian path loss PL(dij), dB between all node pairs in network
PL = PL_d0+10*n_pl*log10(distance./d0)+X_sigma;

%step 2: calculate the SNR value, psi(dij), by a given target PRR(dij) between 0.9 and 1  
psi = -2*0.64*log(2*(1-PRR^(1/(8*phi))));

%step 3: calculate the transmission power values for all node pairs in network, Pt
Pt = psi+PL+Pn;
Pt = 10^(Pt.*0.1); %convert Pt from dB to Watt to incorporate the code in the old version

%step 4: calculate the energy required to transmit one bit(1 dB) for all node pair in the network, Et;
% and then find the power needed to transmit a signal for all node pair in the network,P
Et = rho+Pt.*epsilon*Tb; %final equation for HCB in paper
power = Et./Tb; %convert energy per bit to power for signal in order to incorporrate the code of the old version

%step 5: calculate power recieved for all node pairs in network, Pr
Pr = Pt-PL;
Pr = 10^(Pr.*0.1);
%a new parameter, Pr, to test the performance added in this step


vect=distance>r;
power(vect)=0;
%distance(vect)=0;

n_optimum=const*distance; %for HCB, n_optmum = 0.031d
n_optimum=round(n_optimum); %n is rounded to the nearest integer
n_optimum=(n_optimum);
end

