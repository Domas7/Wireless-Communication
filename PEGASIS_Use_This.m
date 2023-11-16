%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              PEGASIS                                 %
%        Power-Efficient Gathering in Sensor Information Systems       %
%                              Protocol                                %
%                           Part of Thesis:                            %
%       Energy-Efficient Protocols In Wireless Sensor Networks         %
%  First Test (Consumed per trans.) -> 0.00101                                                                 %
% (c) Alexandros Nikolaos Zattas                                       %
% University of Peloponnese                                            %
% Department of Informatics and Telecommunications                     %
% For any related questions or additional information                  %
% send an e-mail to:                                                   %
% alexzattas@gmail.com                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;
%%%%%%%%%%%%%%%%%%%% Network Establishment Parameters %%%%%%%%%%%%%%%%%%%%
%%% Area of Operation %%%
% Field Dimensions in meters %
xm=1000;
ym=1000;
maxRounds = 100;
packets_received = 0;
%x=0; % added for better display results of the plot
%y=0; % added for better display results of the plot
% Number of Nodes in the field %
n=100;
% Number of Dead Nodes in the beggining %
dead_nodes=0;
% Coordinates of the Sink (location is predetermined in this simulation) %
sinkx=xm/2;
sinky=ym/2;
%%% Energy Values %%%
% Initial Energy of a Node (in Joules) % 
Eo=0.02; % units in Joules
% Energy required to run circuity (both for transmitter and receiver) %
Eelec=50*0.000000001; % units in Joules/bit
ETX=50*0.000000001; % units in Joules/bit
ERX=50*0.000000001;
% Transmit Amplifier Types %
Efs=0.0013*0.000000000001; %10*0.000000000001 units in Joules/bit/m^2 (amount of energy spent by the amplifier to transmit the bits)
% Data Aggregation Energy %
EDA=5*0.000000001; % units in Joules/bit
tx_power=30;
% Size of data package %
packetLength=6400; % units in bits 
% Round of Operation %
rnd=0;
ctrPacketLength = 200;
% Current Number of operating Nodes %
operating_nodes=n;
transmissions=0;
d(n,n)=0;
temp_dead=0;
dead_nodes=0;
selected=0;
flag1stdead=0;
received_packets_array = zeros(1, flag1stdead + 1);
count=0;
turn=0;
temp_val=0;
min_power=-80;
mountain = [1, 1; 1, 500; 500, 500; 500, 1; 1,1];
forest = [600, 500; 600, 900; 900, 900; 900, 500; 600, 500];
lake = [1, 600; 1,1000; 400, 1000; 400, 600; 1,600];
%%%%%%%%%%%%%%%%%%%%%%%%%%% End of Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Creation of the Wireless Sensor Network %%%
% Plotting the WSN %
for i=1:n
    
    SN(i).id=i;	% sensor's ID number
    SN(i).x=rand(1,1)*xm;	% X-axis coordinates of sensor node
    XR(i)=SN(i).x;
    SN(i).y=rand(1,1)*ym;	% Y-axis coordinates of sensor node
    YR(i)=SN(i).y;
    SN(i).E=Eo;     % nodes energy levels (initially set to be equal to "Eo"
    SN(i).cond=1;   % States the current condition of the node. when the node is operational its value is =1 and when dead =0
    SN(i).dts=0;    % nodes distance from the sink
    SN(i).role=0;   % node acts as normal if the value is '0', if elected as a cluster head it  gets the value '1' (initially all nodes are normal)
    SN(i).pos=0;
    SN(i).closest=0;
    SN(i).prev=0;
    S(i).packToSend=1; % number of data to send
    SN(i).dis=0;	% distance between two nodes headin towards to the cluster head from position 1
    SN(i).dis2=0;   % distance between two nodes headin towards to the cluster head from position 2
    SN(i).order=0;
    SN(i).sel=0;    % states if the node has already operated for this round or not (if 0 then no, if 1 then yes) 
    SN(i).rop=0;    % number of rounds node was operational
    SN(i).tel=0;    % states how many times the node was elected as a Cluster Head
    order(i)=0;
    hold on;
    figure(1)
    plot(xm,ym,SN(i).x,SN(i).y,'ob',sinkx,sinky,'*r');
    plot(mountain(:, 1), mountain(:, 2), 'r-', 'LineWidth', 2);
    plot(forest(:, 1), forest(:, 2), 'g-', 'LineWidth', 2);
    plot(lake(:,1), lake(:,2), 'b-', 'LineWidth', 2);
    title 'Wireless Sensor Network';
    xlabel '(m)';
    ylabel '(m)';
    
end
SN(n+1).x = sinkx;
SN(n+1).y = sinky;
 
    % Calculates Distance Between Each Node and the Sink (Base Station) %
 for i=1:n
    SN(i).dts=sqrt((sinkx-SN(i).x)^2 + (sinky-SN(i).y)^2);
    SN(i).Esink=Eelec*packetLength + Efs*packetLength*(SN(i).dts)^2;
    T(i)=SN(i).dts;
 end
 
 
    A=sort(T,'descend'); % Creates array A containing the distance between each node and the sink,
                % sorted in an asceding order
     
     A_id(1:n)=0;
     % Creates array A_id which is sorted in a way that it's elements are
     % aligned with those of A. Contains the node ID
     for i=1:n
         for j=1:n
            if A(i)==SN(j).dts
               A_id(i)=SN(j).id;
            end
         end
     end
     
     % Creation of d Array with shortest distances %
            for i=1:n
                SN(i).closest=0;
                for j=1:n
                    d(j,i)=sqrt((SN(i).x-SN(j).x)^2 + (SN(i).y-SN(j).y)^2);
                    if d(j,i)==0
                        d(j,i)=9999;
                    end
                end
            end
         
        for i=1:n     
            [M,I]=min(d(:,i)); % finds the minimum distance of node to CH
            [Row, Col] = ind2sub(size(d),I); % displays the Cluster Number in which this node belongs too
            SN(i).closest=Row; % assigns node to the cluster
            SN(i).dis= d(Row,i); % assigns the distance of node to CH
        end
     
        % Choosing furthest node from sink %
    for i=1:n
         if SN(A_id(i)).E>0 && SN(A_id(i)).sel==0 && SN(A_id(i)).cond==1
            set= A_id(i);
            SN(set).sel=1;
            SN(set).pos=1;
            break;
         end
    end
     order(1)=set;
     temp=1;   

    while selected<n
        min_dis=9999;
        for i=1:n
            if  SN(i).sel==0 
                d=sqrt((SN(i).x-SN(set).x)^2 + (SN(i).y-SN(set).y)^2);
                if d<min_dis
                    min_dis=d;
                    next=i; 
                end
            end
        end
        selected=selected+1;
        SN(set).closest=next;
        SN(set).dis=min_dis;
        SN(next).sel=1;
        SN(next).prev=set;
        SN(next).dis2=sqrt((SN(set).x-SN(next).x)^2 + (SN(set).y-SN(next).y)^2);
        plot([SN(set).x SN(next).x], [SN(set).y SN(next).y])
        hold on;
        set=next;
        temp=temp+1;
        order(temp)=set;
    end
        order(n+1)=[];
        SN(set).pos=2;
        SN(set).dis=0;
        SN(set).closest=0;
        for i=1:n
            if SN(i).closest==set && SN(i).pos==0;
               SN(set).prev=i;
               SN(set).dis2=sqrt((SN(i).x-SN(set).x)^2 + (SN(i).y-SN(set).y)^2);
            end
        end
    
        
        
        
for rnd = 1:maxRounds
                    % Data Transmission %
    if operating_nodes == 1
        fprintf("Transmittion Stopped, only one node is alive.\n");
        break;
    end
    energy=0;
    for i=1:n
        SN(i).role=0;
    end
     
        % Cluster Head Election %
    cluster_head=mod(turn,n) +1;
    if SN(cluster_head).cond==0
        while SN(cluster_head).cond==0
            turn=turn+1;
            cluster_head=mod(turn,n) +1;
        end
    end
    if SN(cluster_head).cond==1
        SN(cluster_head).role=1;
        SN(cluster_head).tel=SN(cluster_head).tel+1;
        figure(1)
        plot(SN(cluster_head).x,SN(cluster_head).y,'+r')
    end
    for i=1:n
        if order(i)==cluster_head
            cl_pos=i;
            break;
        end
    end
        
    for i=1:n
        if SN(order(i)).E>0 && SN(order(i)).cond==1
            if i<cl_pos       
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              if SN(order(i)).pos==1 && SN(order(i)).role==0 
                  SN(i).E = SN(i).E -(ETX*(ctrPacketLength) + Efs*ctrPacketLength); 
                  SN(i).E = SN(i).E -(ETX*(packetLength) + Efs*packetLength);
              end
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              if SN(order(i)).pos==0 && SN(order(i)).role==0
                  SN(i).E = SN(i).E -(ETX*(ctrPacketLength) + Efs*ctrPacketLength); 
                  SN(i).E = SN(i).E -(ETX*(packetLength) + Efs*packetLength);
                  
                  power=tx_power;
                  source = [SN(i).x, SN(i).y];
                  dest = [SN(n+1).x, SN(n+1).y];
                  distTemp=pdist2(source,dest);
                  %distTemp = norm(source - dest); % Euclidean distance calculation
    
                  [in1, dist1]=distInRegion(source,dest,mountain);
                  [in2, dist2]=distInRegion(source,dest,forest);
                  [in3, dist3]=distInRegion(source,dest,lake);
                  % Is the packet received?
                  
                  power=calcLossEnv(in1,dist1,3,power);
                  distTemp=distTemp-dist1;
        
                  power=calcLossEnv(in2,dist2,2,power);
                  distTemp=distTemp-dist2;
        
                  power=calcLossEnv(in3,dist3,4,power);
                  distTemp=distTemp-dist3;
        
                  power_lost=SendSignal(1,distTemp);
                  power=power-power_lost;
                  if((power)>min_power)
                      SN(order(i)).E =SN(order(i)).E - ((ERX + EDA)*packetLength);
                      SN(order(i)).E =SN(order(i)).E - ERX * ctrPacketLength;
                      %ERX1=(EDA+Eelec)*(packetLength);
                      %ERX2=(Eelec*ctrPacketLength);
                      %SN(order(i)).E =SN(order(i)).E - ERX1;
                      %SN(order(i)).E =SN(order(i)).E - ERX2;
                     % SN(order(i)).E =SN(order(i)).E - ERX *ctrPacketLength;
                  end
              end     
          end
               
           if i>cl_pos       
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              if SN(order(i)).pos==2 && SN(order(i)).role==0 
                  SN(i).E = SN(i).E -(ETX*(ctrPacketLength) + Efs*ctrPacketLength); 
                  SN(i).E = SN(i).E -(ETX*(packetLength) + Efs*packetLength);
              end
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              if SN(order(i)).pos==0 && SN(order(i)).role==0
                  SN(i).E = SN(i).E -(ETX*(ctrPacketLength) + Efs*ctrPacketLength); 
                  SN(i).E = SN(i).E -(ETX*(packetLength) + Efs*packetLength);

                  power=tx_power;
                  %is the paquet going through hard terrain?
                  source=[SN(i).x,SN(i).y];
                  dest=[(SN(n+1).x),(SN(n+1).y)];
                  distTemp=pdist2(source,dest);
                  [in1,dist1]=distInRegion(source,dest,mountain);
                  [in2,dist2]=distInRegion(source,dest,forest);
                  [in3, dist3]=distInRegion(source,dest,lake);
                  % Is the packet received?
                  
                  power=calcLossEnv(in1,dist1,3,power);
                  distTemp=distTemp-dist1;
        
                  power=calcLossEnv(in2,dist2,2,power);
                  distTemp=distTemp-dist2;
        
                  power=calcLossEnv(in3,dist3,4,power);
                  distTemp=distTemp-dist3;
        
                  power_lost=SendSignal(1,distTemp);
                  power=power-power_lost;
                  if((power)>min_power)
                      SN(order(i)).E =SN(order(i)).E - ((ERX + EDA)*packetLength);
                      SN(order(i)).E =SN(order(i)).E - ERX * ctrPacketLength;
                  end
               end    
            end
            
           if i==cl_pos       
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              if SN(order(i)).pos==0  && SN(order(i)).role==1
                  SN(i).E = SN(i).E -(ETX*(ctrPacketLength) + Efs*ctrPacketLength); 
                  SN(i).E = SN(i).E -(ETX*(packetLength) + Efs*packetLength);

                  power=tx_power;
                  %is the paquet going through hard terrain?
                  source=[SN(i).x,SN(i).y];
                  dest=[(SN(n+1).x),(SN(n+1).y)];
                  distTemp=pdist2(source,dest);
                  [in1,dist1]=distInRegion(source,dest,mountain);
                  [in2,dist2]=distInRegion(source,dest,forest);
                  [in3, dist3]=distInRegion(source,dest,lake);
                  % Is the packet received?
                  
                  power=calcLossEnv(in1,dist1,3,power);
                  distTemp=distTemp-dist1;
        
                  power=calcLossEnv(in2,dist2,2,power);
                  distTemp=distTemp-dist2;
        
                  power=calcLossEnv(in3,dist3,4,power);
                  distTemp=distTemp-dist3;
        
                  power_lost=SendSignal(1,distTemp);
                  power=power-power_lost;
                  if((power)>min_power)
                      SN(order(i)).E =SN(order(i)).E - ((ERX + EDA)*packetLength);
                      SN(order(i)).E =SN(order(i)).E - ERX * ctrPacketLength;
                  end
              end
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              if SN(order(i)).pos==1  && SN(order(i)).role==1
                  SN(i).E = SN(i).E -(ETX*(ctrPacketLength) + Efs*ctrPacketLength); 
                  SN(i).E = SN(i).E -(ETX*(packetLength) + Efs*packetLength);

                  power=tx_power;
                  %is the paquet going through hard terrain?
                  source=[SN(i).x,SN(i).y];
                  dest=[(SN(n+1).x),(SN(n+1).y)];
                  distTemp=pdist2(source,dest);
                  [in1,dist1]=distInRegion(source,dest,mountain);
                  [in2,dist2]=distInRegion(source,dest,forest);
                  [in3, dist3]=distInRegion(source,dest,lake);
                  % Is the packet received?
                  
                  power=calcLossEnv(in1,dist1,3,power);
                  distTemp=distTemp-dist1;
        
                  power=calcLossEnv(in2,dist2,2,power);
                  distTemp=distTemp-dist2;
        
                  power=calcLossEnv(in3,dist3,4,power);
                  distTemp=distTemp-dist3;
        
                  power_lost=SendSignal(1,distTemp);
                  power=power-power_lost;
                  if((power)>min_power)
                      SN(order(i)).E =SN(order(i)).E - ((ERX + EDA)*packetLength);
                      SN(order(i)).E =SN(order(i)).E - ERX * ctrPacketLength;
                  end
              end
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              if SN(order(i)).pos==2  && SN(order(i)).role==1
                  SN(i).E = SN(i).E -(ETX*(ctrPacketLength) + Efs*ctrPacketLength); 
                  SN(i).E = SN(i).E -(ETX*(packetLength) + Efs*packetLength);

                  power=tx_power;
                  %is the paquet going through hard terrain?
                  source=[SN(i).x,SN(i).y];
                  dest=[(SN(n+1).x),(SN(n+1).y)];
                  distTemp=pdist2(source,dest);
                  [in1,dist1]=distInRegion(source,dest,mountain);
                  [in2,dist2]=distInRegion(source,dest,forest);
                  [in3, dist3]=distInRegion(source,dest,lake);
                  % Is the packet received?
                  
                  power=calcLossEnv(in1,dist1,3,power);
                  distTemp=distTemp-dist1;
        
                  power=calcLossEnv(in2,dist2,2,power);
                  distTemp=distTemp-dist2;
        
                  power=calcLossEnv(in3,dist3,4,power);
                  distTemp=distTemp-dist3;
        
                  power_lost=SendSignal(1,distTemp);
                  power=power-power_lost;
                  if((power)>min_power)
                      SN(order(i)).E =SN(order(i)).E - ((ERX + EDA)*packetLength);
                      SN(order(i)).E =SN(order(i)).E - ERX * ctrPacketLength;
                  end
               end 
           end
           packets_received = packets_received + 1;        
           received_packets_array(rnd+1) = packets_received;
        end
        if SN(order(i)).E<=0 && SN(order(i)).cond==1
             SN(order(i)).cond=0;
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             if SN(order(i)).pos==1 && SN(order(i)).role==0
                if operating_nodes==1
                done='OK'
                else
                t=i;
                while SN(order(t)).cond==0 && t<n
                    t=t+1;
                end
                SN(order(t)).pos=1;
                SN(order(t)).prev=0;
                end
             end
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             if SN(order(i)).pos==2 && SN(order(i)).role==0
                if operating_nodes==1
                done='OK';
                else
                t=i;
                while SN(order(t)).cond==0 && t>1
                    t=t-1;
                end
                SN(order(t)).pos=2; 
                SN(order(t)).closest=0;
                end
             end
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             if SN(order(i)).pos==0 && SN(order(i)).role==1
                SN(order(i)).role=0;
                after=i;
                for l=after:n 
                    if SN(order(l)).cond==1
                    break;
                    end
                end
                bef=i;
                for j=bef:-1:1 
                    if SN(order(j)).cond==1
                    break;
                    end
                end
                SN(order(j)).closest=order(l);
                SN(order(l)).prev=order(j);
                SN(order(j)).dis=sqrt((SN(order(l)).x-SN(order(j)).x)^2 + (SN(order(l)).y-SN(order(j)).y)^2);
                SN(order(l)).dis2=SN(order(j)).closest; 
             end
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             if SN(order(i)).pos==1 && SN(order(i)).role==1
                SN(order(i)).role=0;
                t=i;
                while SN(order(t)).cond==0 && t<n
                    t=t+1;
                end
                SN(order(t)).pos=1;
             end
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             if SN(order(i)).pos==2 && SN(order(i)).role==1
                SN(order(i)).role=0;
                t=i;
                while SN(order(t)).cond==0 && t>1
                    t=t-1;
                end
                SN(order(t)).pos=2; 
             end
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             if SN(order(i)).pos==0 && SN(order(i)).role==0
                after=i;
                for l=after:n 
                    if SN(order(l)).cond==1
                    break;
                    end
                end
                bef=i;
                for j=bef:-1:1 
                    if SN(order(j)).cond==1
                    break;
                    end
                end
                SN(order(j)).closest=order(l);
                SN(order(l)).prev=order(j);
                SN(order(j)).dis=sqrt((SN(order(l)).x-SN(order(j)).x)^2 + (SN(order(l)).y-SN(order(j)).y)^2);
                SN(order(l)).dis2=SN(order(j)).dis;       
             end
             operating_nodes=operating_nodes-1;
             dead_nodes=dead_nodes+1;
             SN(order(i)).closest=0; 
             SN(order(i)).prev=0;
             SN(order(i)).dis=0;
             SN(order(i)).dis2=0;
             SN(order(i)).pos=101;
             SN(order(i)).rop=rnd;
         end
     end
     if operating_nodes<n && temp_val==0
        temp_val=1;
        flag1stdead=rnd;
     end
    % if energy>0
    %     nrg(rnd)=energy;
    % end
     rnd=rnd+1;
     turn=turn+1;
     op(rnd)=operating_nodes;
end


%figure(1); % Make sure you're plotting on the correct figure
%plot([SN(i).x sinkx], [SN(i).y sinky], 'g-'); % Change color and style as needed

%plot(x, y, xm, ym, SN(i).x, SN(i).y, 'ob', sinkx, sinky, '*r');


fprintf('Number of packets Received by Sink: %d\n', packets_received);
%Plotting Figure

%plot(0:rnd, received_packets_array, '-o');
%title('Number of Received Packets Over Time');
%xlabel("Round");
%ylabel("Number of Packets Received");

sum=0;
%for i=1:flag1stdead
%    sum = nrg(i) + sum;
%end
temp1=sum/flag1stdead;
temp2=temp1/n;
for i=1:flag1stdead
    avg_node(i)=temp2;
end        
   
    % Plotting Simulation Results "Operating Nodes per Transmission" %
    figure(2)
    plot(1:rnd,op(1:rnd),'-r','Linewidth',2);
    title ({'PEGASIS'; 'Operational Nodes per Transmission';})
    xlabel 'Transmissions';
    ylabel 'Operational Nodes';
    hold on;
    
    % Check Energy Consumption Data
%    disp(nrg);
    %%
    % Plotting Simulation Results "Energy consumed per Transmission" %
    %figure(3)
    %plot(1:flag1stdead, nrg(1:flag1stdead), '-r', 'Linewidth', 2);
    %title({'PEGASIS'; 'Energy consumed per Transmission';})
    %xlabel 'Transmission';
    %ylabel 'Energy ( J )';
    %hold on;

    % Plotting Number of Packets Received by the Base Station %
    figure(4)
    plot(1:rnd, received_packets_array, '-b', 'Linewidth', 2);
    title({'PEGASIS'; 'Number of Packets Received by the Base Station';})
    xlabel 'Transmissions';
    ylabel 'Number of Packets Received';
    hold on;

    % Plotting Number of Packets Received per Transmission
    figure;
    plot(1:rnd-1, diff(received_packets_array), '-o', 'LineWidth', 2);
    title('Packets Received per Transmission');
    xlabel('Transmission');
    ylabel('Packets Received');
    grid on;

    % Check if avg_node is not empty before plotting
    if exist('avg_node', 'var')
        % Plotting Simulation Results "Average Energy consumed by a Node per Transmission"
        figure(5)
        plot(1:flag1stdead, avg_node(1:flag1stdead), '-r', 'Linewidth', 2);
        title({'PEGASIS'; 'Average Energy consumed by a Node per Transmission';})
        xlabel 'Transmissions';
        ylabel 'Energy ( J )';
        hold on;
    else
        disp('avg_node is empty. Check the code for errors.');
    end

    
    % Check Average Node Energy Data
    disp(avg_node);
    

  