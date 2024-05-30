function [Z,YBUS]=Ybus_7bus()

%               from bus   to bus                impedance         line charging
BD = [            1        2                  0.0000+1i*0.1000       1i*0.000 
                  6        2                  0.0175+1i*0.0628       1i*0.030 
                  6        5                  0.0777+1i*0.2013       1i*0.025
                  2        5                  0.0573+1i*0.158        1i*0.020 
                  2        4                  0.0607+1i*0.171        1i*0.020 
                  2        3                  0.0431+1i*0.14         1i*0.015 
                  5        4                  0.011+1i*0.028         1i*0.010 
                  4        3                  0.0810+1i*0.2065       1i*0.025 
                  6        7                  0.0000+1i*0.1000       1i*0.010];
Z = zeros(7);
ych=zeros(7);
for i=1:size(BD,1)
    Z(BD(i,1),BD(i,2))=BD(i,3);    
    Z(BD(i,2),BD(i,1))=BD(i,3);  
    ych(BD(i,1),BD(i,2))=BD(i,4);    
    ych(BD(i,2),BD(i,1))=BD(i,4);
end

ybus=zeros(7);
for i=1:7
    for j=1:7
        
            if Z(i,j)~=0
            ybus(i,j)=1/Z(i,j);
                      
        end
    end
end

YBUS=zeros(7);
for i=1:7
    for j=1:7
        if i==j
            YBUS(i,j)=sum(ybus(i,:))+sum(ych(i,:));
        else
            YBUS(i,j)=-ybus(i,j);
        end
    end
end


           