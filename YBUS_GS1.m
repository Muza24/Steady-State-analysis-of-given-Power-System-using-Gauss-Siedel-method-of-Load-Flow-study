clc
clear
format compact
[Z,LD,Y_Bus,nbus,Y_b] = Ybusfunction2();
BaseMVA=100;

%1-Slack bus,2-PQ bus,3-PV bus
%          1     2         3      4   5    6       7    8 
%      bus_no bus_type   vmag  vphase Pi   Qi    Qmax Qmin 
BD=[     1       1       1.05     0    0    0      10  -10;
         2       2       1.00     0   -4  -2.5   10  -10;
         3       3       1.04     0    2    0      10  -10];


V=(BD(:,3).*exp(1i*BD(:,4)))';
V_prev=V;
Tolerance=10^-10;
iteration=0;
delV=1;
while delV>Tolerance
    iteration=iteration+1;
    
    for b=1:nbus
        if BD(b,2)==2 %PQ Bus
            V(b)=((BD(b,5)-1i*BD(b,6))/(conj(V(b)))-(sum(Y_Bus(b,:).*V))+Y_Bus(b,b)*V(b))/Y_Bus(b,b);
                %Voltage and magnitude are both updated
                %abs(V(b))
                %angle(V(b))
                BD(b,4)=angle(V(b));
                BD(b,3)=abs(V(b));
        elseif BD(b,2)==3 %PV bus
            BD(b,6)=-imag(conj(V(b))*sum(Y_Bus(b,:).*V));%Qi calculation
            if BD(b,6)>BD(b,7)
                BD(b,6)=BD(b,7);
                BD(b,2)=2;
                V(b)=((BD(b,5)-1i*BD(b,6))/(conj(V(b)))-(sum(Y_Bus(b,:).*V))+Y_Bus(b,b)*V(b))/Y_Bus(b,b);
                %Voltage and magnitude are both updated
                BD(b,4)=angle(V(b));
                BD(b,3)=abs(V(b));
            elseif BD(b,6)<BD(b,8)
                BD(b,6)=BD(b,8);
                BD(b,2)=2;
                V(b)=((BD(b,5)-1i*BD(b,6))/(conj(V(b)))-(sum(Y_Bus(b,:).*V))+Y_Bus(b,b)*V(b))/Y_Bus(b,b);
                %Voltage and magnitude are both updated
                BD(b,4)=angle(V(b));
                BD(b,3)=abs(V(b));
            else
                BD(b,4)=angle(((BD(b,5)-1i*BD(b,6))/(conj(V(b)))-(sum(Y_Bus(b,:).*V))+Y_Bus(b,b)*V(b))/Y_Bus(b,b));
                V(b)=BD(b,3)*exp(1i*BD(b,4));%Magnitude remains same, only the angle gets updated
            end
        end
    end
    delV=max(abs(V-V_prev));
    V_prev=V;
end
% %We have the load flow solution, now we need to find the Slack bus power
% %and line losses

G=real(Y_Bus);
B=imag(Y_Bus);
theta=zeros(size(G,1)); 
R=zeros(size(G,1)); 
for i=1:size(G,1)
    for j=1:size(G,2)
        [theta(i,j),R(i,j)]=cart2pol(G(i,j),B(i,j));
    end
end
V_mag=abs(V);P_g=zeros(nbus); 
del=angle(V);Q_g=zeros(nbus);
%Calculating the power flow between buses i and k and total bus generation
for i=1:nbus
    for k=1:nbus
        P_g(i,k)=V_mag(i)*V_mag(k)*R(i,k)*cos(theta(i,k)-del(i)+del(k));
        Q_g(i,k)=-V_mag(i)*V_mag(k)*R(i,k)*sin(theta(i,k)-del(i)+del(k));
    end
end
P_g=sum(P_g,2);
Q_g=sum(Q_g,2);
% %Now we will caculate the line losses
P=zeros(nbus);
Q=zeros(nbus);
for i=1:nbus
    for k=1:nbus
        P(i,k)=-(V_mag(i)^2)*R(i,k)*cos(theta(i,k))+ V_mag(i)*V_mag(k)*R(i,k)*cos(theta(i,k)-del(i)+del(k));
        Q(i,k)=(V_mag(i)^2)*R(i,k)*sin(theta(i,k))- V_mag(i)*V_mag(k)*R(i,k)*sin(theta(i,k)-del(i)+del(k))-(V_mag(i)^2)*imag(Y_b(i,k));
    end
end
S=P+1i*Q;
% P=P*BaseMVA;
% Q= Q*BaseMVA;
%For loss in each line
P_loss=zeros(nbus);
Q_loss=zeros(nbus);
for i=1:nbus
    for k=1:nbus
        P_loss(i,k)=(P(i,k)+P(k,i));
        Q_loss(i,k)=Q(i,k)+Q(k,i);
    end
end
BD(:,4)=BD(:,4)*180/pi;
P_loss=P_loss*BaseMVA;
Q_loss=Q_loss*BaseMVA;
P_g=P_g*BaseMVA;
Q_g=Q_g*BaseMVA;
Total_loss=(sum(P_loss(:))+1i*sum(Q_loss(:)))/2;
%Total Generation = Total Load+ Total Losses

fprintf('\nY-Bus matrix for the 3-bus system:\n');
Y_Bus

% %Exporting data to excel
% Lineloss=zeros(3,4);j=1;l=1;
% for i=1:nbus
%     for k=1:nbus
%         if i~=k && i<k && Y_Bus(i,k)~=0
%          Lineloss(j,l)=i;
%          Lineloss(j,l+1)=k;
%          Lineloss(j,l+2)=round(P_loss(i,k),5);
%          Lineloss(j,l+3)=round(Q_loss(i,k),5);
%          j=j+1;
%         end
%     end
% end
% A=Lineloss(2,:);
% Lineloss(2,:)=Lineloss(3,:);
% Lineloss(3,:)=A;
% writematrix(Lineloss,'testdata1.xlsx','Range','A2');
% Powerflow=zeros(3,5);
% for i=1:nbus
%     Powerflow(i,1)= BD(i,1);
%     Powerflow(i,2)=round(BD(i,3),5);
%     Powerflow(i,3)=round(BD(i,4),5);
%     Powerflow(i,4)=round(P_g(i),5);
%     Powerflow(i,5)=round(Q_g(i),5);
% end
% writematrix(Powerflow,'testdata1.xlsx','Range','A8');

fprintf('Complex Power Flow in lines between ith and jth bus(in p.u.):\n');

S
fprintf('Line losses:  Line   P_Loss(MW)  Q_Loss(MVar)\n');
for i=1:nbus
    for k=1:nbus
        if i~=k && i<k
         fprintf('\t\t\t %i - %i   %f    %f\n',i,k,P_loss(i,k),Q_loss(i,k));
        end
    end
end

fprintf('Total loss: %f + i%f MVA\n', real(Total_loss),imag(Total_loss));
fprintf('\nFinal Soln: Bus   Voltage  Angle(degree)    Pi(MW)         Qi(MVar) \n');
for i=1:nbus
    fprintf(' \t\t\t %i   %f  %f       %f      %f\n', BD(i,1),BD(i,3),BD(i,4),P_g(i),Q_g(i));
end
iteration