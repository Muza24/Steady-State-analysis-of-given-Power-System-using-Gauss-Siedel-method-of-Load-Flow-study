clc
clear
Base_mva=100;

        % Busno.   type of bus   Vmag               Vt.ph      Pi     Qi       Pd     Qd         Qmax    Qmin  
BD=[      1          1           1.06                0.0       0       0       0      0          1.05   -10.0 
          2          2           1.00                0.0       0       0       0      0          0.00    0.00 
          3          2           1.00                0.0     -0.2    -0.1     0.2     0.1        0.00    0.00 
          4          2           1.00                0.0     -0.4    -0.05    0.4     0.05       0.00    0.00 
          5          2           1.00                0.0     -0.45   -0.1     0.45     0.1       0.00    0.00 
          6          2           1.00                0.0     -0.6    -0.1     0.6     0.1        0.00    0.00 
          7          1           1.00                0.0      0.4       0       0       0        10.0   -10.0 ];
[z,Y]=Ybus_7bus();

Vo=(BD(:,3).*exp(1i*BD(:,4)))';
tolerance=1e-6;
iteration=0;
err=10;
Nbus=7;
alpha=1.5;


 while err>tolerance
    V0=Vo;
   iteration=iteration+1;
for b=1:Nbus
    if BD(b,2)==2  %pq bus
       
         Vo(b)= ((BD(b,5)-1i*BD(b,6))/conj(Vo(b)) - sum(Y(b,:).*Vo) + Y(b,b)*Vo(b)) / Y(b,b);
       
        
        
       BD(b,3)=abs(Vo(b));
       BD(b,4)=angle(Vo(b));
    
    elseif BD(b,2)==3  %pv bus
           BD(b,6)=-imag(conj(Vo(b))*(sum(Y(b,:).*Vo)));
    
        if BD(b,6)>BD(b,9)            %Q>Qmax
            BD(b,6)=BD(b,9);
            BD(b,2)=2;
            BD(b,3)=1.00;
             Vo(b)= ((BD(b,5)-1i*BD(b,6))/conj(Vo(b))- sum(Y(b,:).*Vo) + Y(b,b)*Vo(b)) / Y(b,b);
        elseif BD(b,6)<BD(b,10)         %Q<Qmin
             BD(b,6)=BD(b,10);
              BD(b,2)=2;
              BD(b,3)=1.00;
               Vo(b)= ((BD(b,5)-1i*BD(b,6))/conj(Vo(b))- sum(Y(b,:).*Vo) + Y(b,b)*Vo(b)) / Y(b,b);
        end

        BD(b,4)=angle((((BD(b,5)-1i*BD(b,6))/conj(Vo(b)))-(sum(Y(b,:).*Vo)-(Y(b,b).*Vo(b))))/Y(b,b));
        Vo(b)=BD(b,3)*exp(1i*BD(b,4));
    end
end

p=Vo-V0;
err=max(abs(p));
Vo=V0+alpha*(Vo-V0); % accelaration factor
end



Y_mag=abs(Y);
V_mag=abs(Vo);
del=angle(Vo);

G=real(Y);
B=imag(Y);

for i=1:size(G,1)
    for k=1:size(B,1)
        [theta(i,k),R(i,k)]=cart2pol(G(i,k),B(i,k));
    end
end

%slack bus power

PG=zeros(size(Y));
QG=zeros(size(Y));

for i=1:Nbus
    for k=1:Nbus

        PG(i,k)=PG(i,k)+(V_mag(i)*V_mag(k)*Y_mag(i,k)*cos(theta(i,k)-del(i)+del(k)));
        QG(i,k)=QG(i,k)-(V_mag(i)*V_mag(k)*Y_mag(i,k)*sin(theta(i,k)-del(i)+del(k)));
    end
end

PG;
QG;

PG=sum(PG,2);
QG=sum(QG,2);

PG=PG*Base_mva;
QG=QG*Base_mva;

%computation of line flows

for i=1:Nbus
    for k=1:Nbus

        if i~=k
            P(i,k)= -(V_mag(i)^2)*Y_mag(i,k)*cos(theta(i,k) )+ (Y_mag(i,k)*V_mag(i)*V_mag(k)*cos(theta(i,k)-del(i)+del(k)));

             Q(i,k)= (V_mag(i)^2)*Y_mag(i,k)*sin(theta(i,k) )- (Y_mag(i,k)*V_mag(i)*V_mag(k)*sin(theta(i,k)-del(i)+del(k)));
        end
    end
end



%computation of losses
 for i=1:Nbus
    for k=1:Nbus
        ploss(i,k)=P(i,k)+P(k,i);
        qloss(i,k)=Q(i,k)+Q(k,i);
    end
 end
 ploss=ploss*Base_mva
 qloss=qloss*Base_mva

 totalploss=sum(ploss(:))
 totalqloss=sum(qloss(:))

