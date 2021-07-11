%% Transient Analysis for a Multimachine System
clc
clear
%% Reading original bus data
bdata.busdata='busdata4';
bdata=busdata4(bdata);
businfo=bdata.bustype;
lbusinfo=length(businfo);

%% Sorting bus numbers with slack bus first followed by the generator and load bus

ctr=1;
for x=1:3
    for y=1:lbusinfo
        if(businfo(y,2)==x&&x==1)
            temp=businfo(y,:);
            businfo(y,:)=[];
            businfo=[temp;businfo];
            temp=[];
        end
        if(businfo(y,2)==x&&x==2)
            ctr=ctr+1;
            temp=businfo(y,:);
            businfo(y,:)=[];
            temp2=businfo((ctr):(lbusinfo-1),:);
            businfo=[businfo(1:ctr-1,:);temp;temp2];
            temp2=[];
        end
    end
end
%%  Reading line data and renaming bus data and line data with the new bus number
linedata=bdata.line(:,1:5);
r=length(linedata(:,1));
c=length(linedata(1,:));
linerenamed=linedata;
dummy=ones(length(linedata(:,1)),2);
origbus=bdata.bus;
%Sorting line data and bus data. Ensure busdata rows<=line data rows for
%the code to work
for x=1:lbusinfo
    for y=1:length(linedata(:,1))
 % Code for Sorting bus data
         if y<=lbusinfo
             if businfo(x,1)==origbus(y,1)
               modifiedbus(x,:)=origbus(y,:);
               modifiedbus(x,1)=x;
             end
         end

 %Code for Sorting of bus data ends here 
 %sorting line data       
       if linerenamed(y,1)==businfo(x,1)&&dummy(y,1)==1
           linerenamed(y,1)=x;
           dummy(y,1)=0;
       end
       if linerenamed(y,2)==businfo(x,1)&&dummy(y,2)==1
           linerenamed(y,2)=x;
           dummy(y,2)=0;
       end
  % Sorting line data ends     
    end
end
%% Read the sorted line data values and compute the Ybus matrix
%Reading sorted line data values
linedata=linerenamed;
frombus=linedata(:,1);
tobus=linedata(:,2);
R=linedata(:,3);
X=linedata(:,4);
B=linedata(:,5);
G=1./(R+X*i);
buslength=length(R);


for x=1:buslength
    Y1(frombus(x),tobus(x))=G(x)+B(x)*i;
end

length=max(max(frombus),max(tobus));
l=length;
Y=zeros(l,l);
mf=max(frombus);
mt=max(tobus);


for x=1:mf
    for y=1:mt
        
        if x~=y
            if(Y1(x,y)~=0)
            Y(x,y)=-Y1(x,y);
            Y(y,x)=-Y1(x,y);
            end
           
        end
    end
end
Y3=Y;
for x=1:l
    for y=1:l
        
        if x~=y
           
            Y(x,x)=Y(x,x)-Y(x,y);
        end
    end
end
spy(Y)
title('spy of sorted Ybus');

n=lbusinfo;
%% Powerflow calculation

%1=slack bus 2=PV bus 3=PQ bus
bust=modifiedbus(:,10)';
for x=1:lbusinfo
    if bust(x)==1||bust(x)==2
        V(x)=modifiedbus(x,2);
        delta(x)=modifiedbus(x,3)*pi/180;
        Pload(x)=modifiedbus(x,6);
        Pgen(x)=modifiedbus(x,4);

        
    else
        V(x)=1;
        delta(x)=0;
        Pgen(x)=modifiedbus(x,4);
        Qgen(x)=modifiedbus(x,5);
        Pload(x)=modifiedbus(x,6);
        Qload(x)=modifiedbus(x,7);
    end
end
Ymag=abs(Y);
Yang=angle(Y);

error=100;% use while loop to run until error is less than specifeid 
itcount=0; %iteration count
while  error > 0.001 % specified error upto which the iterations are executed
itcount = itcount+1;
cnt = 1;
% Finding the mismatch equations, 1 equation each for PV bus and 2 equations each for PQ bus
% Mismatch of P 
for p = 1:n
   if ((bust(p) == 3) || (bust(p)==2))
        sum1 = 0;
        for q= 1:n
            sum1 = sum1 + V(q) * Ymag(p,q) * cos(delta(p)-delta(q)-Yang(p,q));
        end
        Mismt(cnt) = Pgen(p) - Pload(p) - V(p) * sum1;
        cnt = cnt + 1;
    end  
end
% Mismatch of Q
for p = 1:n
    if (bust(p) == 3)
        sum1 = 0;
        for q = 1:n
            sum1 = sum1 + V(q) * Ymag(p,q) * sin(delta(p)-delta(q)-Yang(p,q));
        end
        Mismt(cnt) = Qgen(p) - Qload(p) - V(p) * sum1;
        cnt = cnt + 1;
    end   
end


% Calculating the Jacobian Matrix 
 p1=0; p2=1;p3=1;p4=0;p5=1;p6=1;
for p=1:n
    if bust(p)==3||bust(p)==2
        p1=p1+1;
        p2=1;
        p3=1;
        for q=1:n
            if bust(q)==3||bust(q)==2
               sum2 = 0;
                if (p == q)
                    for r = 1:n
                        sum2 = sum2 + V(r) * Ymag(p,r) * sin(delta(p) - delta(r) - Yang(p,r)); 
                    end
                    J1(p1,p2) = V(p)* sum2 + V(p) * V(p) * Ymag(p,p) * sin(Yang(p,p));
                 else
                    J1(p1,p2) = -V(p) * V(q) * Ymag(p,q) * sin(delta(p) - delta(q) - Yang(p,q));
                end  
                p2=p2+1;
            end
           if bust(q)==3
               sum2 = 0;
                if (p == q)
                    for r = 1:n
                        sum2 = sum2 + V(r) * Ymag(p,r) * cos(delta(p) - delta(r) - Yang(p,r)); 
                    end
                    J2(p1,p3) = -sum2 - V(p) * Ymag(p,p) * cos(Yang(p,p));
                 else
                    J2(p1,p3) = -V(p) * Ymag(p,q) * cos(delta(p) - delta(q) - Yang(p,q));
                end  
                p3=p3+1;
               
           end
        end
    end
        
        
   if bust(p)==3
        p4=p4+1;
        p5=1;
        p6=1;
        for q=1:n
            if bust(q)==3||bust(q)==2
               sum2 = 0;
                if (p == q)
                    for r = 1:n
                         sum2 = sum2 + V(p) * V(r) * Ymag(p,r) * cos(delta(p) - delta(r) - Yang(p,r)); 
                    end
                        J3(p4,p5)=-sum2 + V(p) * V(p) * Ymag(p,p) * cos(Yang(p,p));
    
                else
                    J3(p4,p5) = V(p) * V(q) * Ymag(p,q) * cos(delta(p) - delta(q) - Yang(p,q));
                end  
                p5=p5+1;
            end
           if bust(q)==3
               sum2 = 0;
                if (p == q)
                    for r = 1:n
                    sum2 = sum2 + V(r) * Ymag(p,r) * sin(delta(p) - delta(r) - Yang(p,r)); 
                    end
                    J4(p4,p6) = -sum2 + V(p) * Ymag(p,p) * sin(Yang(p,p));
                    
                else
                    J4(p4,p6) = -V(p) * Ymag(p,q) * sin(delta(p) - delta(q) - Yang(p,q));
                end
                p6=p6+1;
           end   
       end   
   end
end

% matrix inversion to calculate  the increments
J5=[J1 J2;J3 J4];
Inc = -1 * inv(J5) * Mismt';
inccount=1;
%Updating delta and voltage
for p=1:n
    if bust(p)==3||bust(p)==2
        delta(p)=delta(p)+Inc(inccount);
        inccount=inccount+1;
    end
end
for p=1:n
    if bust(p)==3
        V(p)=V(p)+Inc(inccount);
        inccount=inccount+1;
    end
end    
error = max(abs(Mismt));

end % end of main loo

% Load flow iteration ends here
% Load flow post-calculations

for p=1:n
   sum1=0;
   sum2=0;
    for q=1:n
    sum1=sum1+V(q)*Ymag(p,q)*cos(delta(p)-delta(q)-Yang(p,q));
    sum2=sum2+V(q)*Ymag(p,q)*sin(delta(p)-delta(q)-Yang(p,q));
    end
    Pinj(p)=V(p)*sum1; 
    Qinj(p)=V(p)*sum2; 
end
    Pgen=Pinj+Pload;
    Qgen=Qinj+Qload;

    m=bdata.generatorno;   % first 10 generators buses
    xd=bdata.xd;
%% For each generator, calculating internal voltage and initial rotor angle
    
    
for p=1:m;
    I(p)=conj((Pgen(p)+i*Qgen(p))/(V(p)*exp(i*delta(p))));
    E(p)=V(p)*exp(i*delta(p))+i*xd(p)*I(p);
end

clear delta;
% now the variable delta has angle of gen internal  voltage not terminal voltage 
delta=angle(E);
Eabs=abs(E);
%% For each load in system, Pload and Qload are converted to admittances

% Coverting load to the equivalent impedance, GB=1/(r+jx);
for p=m+1:n;
    GB(p)=(Pload(p)-1i*Qload(p))/V(p)^2;  
end

%% Load admittances are added to the diagonal of the reordered matrix 

Ybus=Ymag.*exp(1i*Yang);
for p=1:n;
    Ybus(p,p)=Ybus(p,p)+GB(p);
end
%% Each generator  in the system is augmented by adding an internal bus connected to the terminal bus with xd'

Y1=zeros(m,m);
Y2=zeros (m,n);
Y3=zeros(n,m);
for p=1:m;
    Y1(p,p)=1./(1i*xd(p));
    Y2(p,p)=-1./(1i*xd(p));
    Y3(p,p)=-1./(1i*xd(p));
end

% Adding xd' to the diagonal of generator Y bus
    Y4=Ybus+[Y1 zeros(m,n-m);
                zeros(n-m,m) zeros(n-m,n-m)];
            
%% Kron reduction for pre-fault condition
 
% Y-reduced before the faults
Yred=Y1-Y2*inv(Y4)*Y3;

%% Reading Fault bus data
nf=bdata.faultyBus;   %faulty node
node1=bdata.nodeOpen(1);  %fault cleared by opening node 1 & node 2
node2=bdata.nodeOpen(2);

%mapping the original faulty bus to reordered bus

for z=1:lbusinfo
    if nf==businfo(z,1)
        nfTemp=z;
    end
    if node1==businfo(z,1)
        node1Temp=z;
    end
    if node2==businfo(z,1)
        node2Temp=z;
    end
    
end
nf=nfTemp;
node1=node1Temp;
node2=node2Temp;


%% Performing kron reduction  during fault condition
Ybus=[Y1 Y2; Y3 Y4]; % The enlarged Y bus before kron reduction where Y1,Y2,Y3 holds only generator bus information and Y4 holds the actual Y bus information
count1=0;
%removing row and column of the faulty node in Y4 where Y4 is the original sorted bus 
for p=1:m+n
    if p~=nf+m
    count1=count1+1;
    count2=0;
    for q=1:m+n
        if q~=nf+m 
        count2=count2+1;
        Ybusf(count1,count2)=Ybus(p,q);
        end
    end
    end
end
%Performing kron reduction
Y1f=Ybusf(1:m,1:m);
Y2f=Ybusf(1:m,m+1:m+n-1);
Y3f=Ybusf(m+1:m+n-1,1:m);
Y4f=Ybusf(m+1:m+n-1,m+1:m+n-1);
Yredf=Y1f-Y2f*inv(Y4f)*Y3f;
%% Kron reduction for post fault condition by removing connection between node1 and node2

%  fault cleared by opening line between node1 and node2


Ybus(node1+m,node1+m)=Ybus(node1+m,node1+m)+Ybus(node1+m,node2+m);
Ybus(node2+m,node2+m)=Ybus(node2+m,node2+m)+Ybus(node1+m,node2+m);
Ybus(node1+m,node2+m)=0;
Ybus(node2+m,node1+m)=0;
Ybuspf=Ybus;
Y1pf=Ybuspf(1:m,1:m);
Y2pf=Ybuspf(1:m,m+1:m+n);
Y3pf=Ybuspf(m+1:m+n,1:m);
Y4pf=Ybuspf(m+1:m+n,m+1:m+n);
Yredpf=Y1pf-Y2pf*inv(Y4pf)*Y3pf;
%% Runge kutta starts here by calculationg inital Pm
delT=0.001;
% Fault inception time
tfi=bdata.faultInceptionTime;
% Fault clearance time tfi plus fault dration
tfc=bdata.faultClearanceTime;
% staudy time horizon
ttot=bdata.runtime;

M=bdata.H*2/(2*pi*60);
%Initial power angle
del(1,:)=delta;
omega(1,:)=zeros(1,m);
t=0;
count=1;
noit=ttot/delT+1;
%Mechanical input
Pe=loadflow(t,tfi,tfc,Yred,Yredf,Yredpf,Eabs,delta); %initial Pe value which is same as Pm
Pm=Pe;
D=bdata.damping;
while t<=ttot  
time(count)=t; 
 
 k1(count,:)=omega(count,:)*delT;
 Pe=loadflow(time(count),tfi,tfc,Yred,Yredf,Yredpf,Eabs,del(count,:));
 l1(count,:)=delT*(Pm-Pe-(D.*k1(count,:)/delT))./M;
 k2(count,:)=(omega(count,:)+l1(count,:)/2)*delT;
 Pe=loadflow(time(count),tfi,tfc,Yred,Yredf,Yredpf,Eabs,del(count,:)+k1(count,:)/2);
 l2(count,:)=delT*(Pm-Pe-(D.*k2(count,:)/delT))./M;
 k3(count,:)=(omega(count,:)+l2(count,:)/2)*delT;
 Pe=loadflow(time(count),tfi,tfc,Yred,Yredf,Yredpf,Eabs,del(count,:)+k2(count,:)/2);
 l3(count,:)=delT*(Pm-Pe-(D.*k3(count,:)/delT))./M;
 k4(count,:)=(omega(count,:)+l3(count,:))*delT;
 Pe=loadflow(time(count),tfi,tfc,Yred,Yredf,Yredpf,Eabs,del(count,:)+k3(count,:));
 l4(count,:)=delT*(Pm-Pe-(D.*k4(count,:)/delT))./M;
 count=count+1; 
 del(count,:)=del(count-1,:)+(1/6)*(k1(count-1,:)+2*k2(count-1,:)+2*k3(count-1,:)+k4(count-1,:));
 omega(count,:)=omega(count-1,:)+(1/6)*(l1(count-1,:)+2*l2(count-1,:)+2*l3(count-1,:)+l4(count-1,:));
 t=t+delT;
 
end 
time(count)=t; 
%% Plot of individual bus angle
% figure
% plot(time,omega(:,1)*180/pi,'r')
% hold on;
% plot(time, omega(:,2)*180/pi,'b')
% hold on;
% plot(time, omega(:,3)*180/pi,'k')
% hold on;
% plot(time, omega(:,4)*180/pi,'y')
% hold on;
% plot(time, omega(:,5)*180/pi,'m')
% hold on;
% plot(time, omega(:,6)*180/pi,'c')
% hold on;
% plot(time, omega(:,7)*180/pi,'g')
% hold on;
% plot(time, omega(:,8)*180/pi,'k')
% hold on;
% plot(time, omega(:,9)*180/pi,'r')
% hold on;
% plot(time, omega(:,10)*180/pi,'m')
% hold on;
% xlabel('Time, s');
% ylabel('Delta, degree');
% grid on;

%% Plotting difference in bus angle with slack bus

set(0,'defaultlinelinewidth',2);
figure
plot(time, (del(:,2)-del(:,1))*180/pi,'b')
hold on;
title('Generator connected to Bus#: 30, 31,32 33 34 35 36 37 38');
plot(time, (del(:,3)-del(:,1))*180/pi,'k')
hold on;
plot(time, (del(:,4)-del(:,1))*180/pi,'y')
hold on
plot(time, (del(:,5)-del(:,1))*180/pi,'m')
hold on;
plot(time, (del(:,6)-del(:,1))*180/pi,'c')
hold on;
plot(time, (del(:,7)-del(:,1))*180/pi,'g')
hold on
plot(time, (del(:,8)-del(:,1))*180/pi,'Color',[0.6 0.3 0.2])
hold on
plot(time, (del(:,9)-del(:,1))*180/pi,'r')
hold on;
plot(time, (del(:,10)-del(:,1))*180/pi,'Color',[0.5 0.3 0.7])
hold off;
legend('Bus 30','Bus 31','Bus 32', 'Bus 33', 'Bus 34', 'Bus 35', 'Bus 36', 'Bus 37', 'Bus 38');
xlabel('Time, s');
ylabel('Delta, degree');
grid on;

%% 





