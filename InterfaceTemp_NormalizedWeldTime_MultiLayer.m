%% Predicting temperature of the interface between two printed layers
% Predictions based on thermal models derived by Rebecca Kurfess
% v2 uses the correct time in the recursive formulas (v1 does not)
% Assumes instantaneous extrusion of a layer
% Assumes both layers have the same dimensions
% Last updated: 6/10/19 by Rebecca Kurfess

close all
clc
clear all



%% Inputs
t = 0:5:1750;                   % total time to calculate cooling of both layer 1 and layer 2 [s]
t_print = 150;                 % print time of a layer (does not include non-printing motions) [s]
t_extra = 0;                   % time to complete non-print motions in each layer (tip wipe, spiral in, etc.) [s]
v_print_in_s = 5;             % print velocity [in/s]
d_nozzle = 0.00762;                      % nozzle diameter [m]
w = d_nozzle*1.5;                % bead width (here calculated to be 1.5 times the nozzle diameter of 0.3in) [m]
z = d_nozzle*0.5;                % bead height (here calculated to be 0.5 times the nozzle diameter of 0.3in) [m]

% Temperatures
T_b = 70+273;               % Bed temperature [K]
T_amb = 35+273;             % Ambient temperature [K]
T_i = 180+273;              % Initial extrusion temperature [K}

% WLF equation parameters
T_r = 250;        %250          % reference temperature for WLF equation (from Christine A) [C]
C_1_wlf = 2.12;              % constant for WLF equation (from Christine A) []
C_2_wlf = 141;              % constant for WLF equation (from Christine A) [K]
 
% Read file 
% filepath = '/Users/gregorydreifus/Downloads/abaqus_rpt.txt';
% list = fileread(filepath);
delta = 1;
numnodes=1638;

testingtest=15;
testingtest
% fid_node1 = fopen('/Users/gregorydreifus/Downloads/shortPrintSimLonger-nodes.txt');
% fid_node2 = fopen('/Users/gregorydreifus/Downloads/complexpath-nodes.txt');
% 
% fid1 = fopen('/Volumes/VERBATIM HD/abaqus.rpt'); %short print
% fid2 = fopen('/Volumes/VERBATIM HD/testing_20190919/abaqus_complex1.rpt'); %complex print

longprintid = fopen('/Users/gregorydreifus/Dropbox (Personal)/GregDreifus_Pum-copy/block03in10bead16layer.rpt');
testingtest1=16;
testingtest1
longprintnodes = fopen('/Users/gregorydreifus/Dropbox (Personal)/GregDreifus_Pum-copy/Job-2.inp');
testingtest2=16;
testingtest2
numnodes_longprint = 23744;

count = 1;
temp_count = count;

weldtimes = [];
% fid = fopen('/Users/gregorydreifus/Downloads/abaqus_rpt.txt');
for w=1:1
    
%     if w==1
%         fid = fid1;
%         numnodes=1638;
%         nid = fid_node1;
%     end
%     if w==2
%         fid = fid2;
%         numnodes=1404;
%         nid = fid_node2;
%     end
    
    nid = longprintid;
    fid = longprintnodes;
    numnodes = numnodes_longprint;
    
    nodes = textscan(nid,'%s','Delimiter','','endofline','\n');
    testingtest3=17;
    testingtest3
    fclose(nid)
    nodesCol = nodes{1};
    nodeList = regexp(nodesCol,'\d+,\s+(-)?\d+(\.)?(\d+)?,\s+(-)?\d+(\.)?(\d+)?,\s+\d+(\.)?(\d+)?','match');
    testingtest4=18;
    testingtest4
    fclose(fid)
    %organize nodes based on layer height 
    list2 = zeros(numnodes,1);
    for a = 1:numnodes
        nodeList1 = nodeList{a};
        nodeListNum = strsplit(nodeList1{1},',');
        nodeListNum1 = regexp(nodeListNum{1},'\d+(\.)?(\d+)?','match'); % node number
        nodeListNum2 = regexp(nodeListNum{2},'\d+(\.)?(\d+)?','match'); % x value
        nodeListNum3 = regexp(nodeListNum{3},'\d+(\.)?(\d+)?','match'); % y value
        nodeListNum4 = regexp(nodeListNum{4},'\d+(\.)?(\d+)?','match'); % z value
        
        list2(a,1) = str2num(nodeListNum4{1}); % layer height = list 2
    end
    testingtest5=19;
    testingtest5
    % sort the layer heights into an ordered list, where b is the 
    maximal = max(list2);
    layerheights = unique(list2);
    
    [amount_of_layers,waste] = size(layerheights);
    
    list3 = ones(amount_of_layers,1);
    testing1=1;
    testing1
    text = textscan(fid,'%s','Delimiter','','endofline','\n');
    text1=text{1};
%     % text2=text1{1};
%     % test=regexp(text1,'\s*\d*\s*\d*\.(\d*)?','match'); %prior version
    test=regexp(text1,'\d+\s+\d+\.(\d*)?','match'); %new test versions
    test3=[test{:}];
    [one two]=size(test3);
%     y=str2num(test3{numnodes});
    count=1;
    list = zeros(numnodes,2);
    tempcount=0;
    twonew=two-1;
    m=1;

    temperatures_unsorted = ones(numnodes,1);
    temperatures_sorted = ones(numnodes,1);
    testing2=2;
    testing2
    for c = 1:amount_of_layers
        for b = 1:numnodes
            x = str2num(test3{b}); %node number
            %if the layer height in a list of all layer heights is equal to
            %the value in a list of unique layer heights add the node value
            %to list 3 in the same row as the layer number 
            % i.e. list3(1,:) == all the node numbers of layer 1
            if list2(b,1) == layerheights(c,1) 
                d = find(list3(c,:),1,'last');
                
                list3(c,d+1) = d+1;
                list3(c,d) = x(1);
            end
        end
    end
    
%     [b,edges] = histcounts(list2,4); 
  

    
    % temperatures(1,15)=15;
    % temprow=temperatures.';
    % size(temperatures,2)
    % size(temprow,1)
    % temprow(1,:)
    % % arrayfun(@(x)find(temperatures(:,x),1,'last'),1:size(temperatures,2))
    %  find(temperatures(1,:),1,'last')
    %  temperatures(1,2)=2;
    %  find(temperatures(1,:),1,'last')
    %  find(temperatures(2,:),1,'last')
    % arrayfun(@(x)find(temprow(x,:),1,'last'),1:size(temprow,1))

    temp_sort_by_layer = ones(amount_of_layers,1);
%     for n = 1:twonew
%     
%         x = str2num(test3{n});
%         count = mod(n-1,numnodes)+1;
%     
%         % x(1)==node number
%         % x(2)==temperature
%         if x(2)<=250 && x(2)>110
%             %list has the node value in column 1 and the count of
%             %temperature values within a desired range in column two
%             list(count,1) = x(1);
%             list(count,2) = list(count,2)+1;
%             
%             
%             
%             %temperature is an arrary of arrays of temperature values for each
%             %node greater than glass transition by less than or equal to 250C
%             % sort temperatures into bins corresponding to each node; count
%             % tracks which bin you're in and m allows tracking of the last
%             % nonzero temperature value of each bin
%             m = find(temperatures_unsorted(count,:),1,'last');
%             temperatures_unsorted(count,m+1) = m+1;
%             temperatures_unsorted(count,m) = x(2); %sort temp by node
%             
% %             for k = 1:amount_of_layers
% %                 if ismember(x(1),list3(k,:)) == 1
% %                     o = find(temperatures_sorted(k,:),1,'last');
% %                     temperatures_sorted(k,o+1) = o+1;
% %                     temperatures_sorted(k,o) = x(2); % temperature sorted by layer
% %                 end
% %             end
% 
%         
%         end
%   
%     end
    
    list4 = zeros(numnodes,2);
    temp_sort_by_layer = {};
    temperary_temperatures_sorted_by_nodes = ones(numnodes,1);
    for i = 1:amount_of_layers
        
        for j = 1:twonew
            
            xx = str2num(test3{j});
            count1 = mod(j-1,numnodes)+1;
            
            if xx(2)<=250 && xx(2)>110
                
                list4(count1,1) = xx(1);
                list4(count1,2) = list(count1,2)+1;
                
                if ismember(xx(1),list3(i,:)) == 1
                    m = find(temperary_temperatures_sorted_by_nodes(count1,:),1,'last');
                    j
                    temperary_temperatures_sorted_by_nodes(count1,m+1) = m+1;
                    temperary_temperatures_sorted_by_nodes(count1,m) = xx(2); %sort temp by node
                end
            end
        end
        temp_sort_by_layer = [temp_sort_by_layer, temperary_temperatures_sorted_by_nodes];
        temperary_temperatures_sorted_by_nodes = ones(numnodes,1);
    end

    %list 2  == layers heights
    %list3(n,:) == all nodes in layer n
    %temperatures_sorted(n,:) == all temperatures in layer n
    %temperatures_unsorted(n,:) == all temperatures in at node n
    % sum up all temperature within a range corresponding to each node
    % number
    
    size1= size(temp_sort_by_layer);
    size2=size(temp_sort_by_layer(1));
    size2=size(temp_sort_by_layer{1});
    
    %size1
    testA=temp_sort_by_layer{1};
    testB=temp_sort_by_layer{2};
    testC=temp_sort_by_layer{3};
    testD=temp_sort_by_layer{4};
    testE=temp_sort_by_layer{5};
    testF=temp_sort_by_layer{6};
    testG=temp_sort_by_layer{7};
    testA(1)
    size1
    size2
    
    a_T = 0;
    temp_lay_list = zeros(amount_of_layers);
    for q = 1:amount_of_layers
        [sizex,sizey] = size(temp_sort_by_layer{q});
        temp_sort = temp_sort_by_layer{q};
        temp_list = zeros(sizex);
        for w = 1:sizex
            for u = 1:sizey
                if tempsort(w,u) < 250
                    a_T = 1/(10^((-C_1_wlf*(temp_sort(w,u)-T_r))./(C_2_wlf+(temp_sort(w,u)-T_r))));
                    
                    t_weld = t_weld+delta/2*(a_T_prev+a_T);
                    a_T_prev = a_T;
                end
                
            end
            temp_list(w) = t_weld;
        end
        temp_lay_list(q) = temp_list;
    end
    
    
    %temp_list = zeros(numnodes,1);
    %UNDO COMMENT HERE
%     temp_list_layers = zeros(numnodes,amount_of_layers);
%     a_T_prev = 0;
%     a_T_prev_sort = 0;
%     for m = 1:numnodes
%         t_weld_sort = 0;
%         for n = 1:amount_of_layers 
%             t_weld = 0;
%             
%             for q = 1:list(m,2)
%                 if temperatures_unsorted(m,q) < 250
%                     if ismember(temperatures_unsorted(m,q),temperatures_sorted(n,:)) == 1  
%                         a_T_sort = 1/(10^((-C_1_wlf*(temperatures_sorted(m,q)-T_r))./(C_2_wlf+(temperatures_sorted(m,q)-T_r))));
%                         t_weld_sort = t_weld+delta/2*(a_T_prev_sort+a_T_sort);
%                         a_T_prev_sort = a_T_sort;
%                     end
%                     %a_T = 1/(10^((-C_1_wlf*(temperatures_unsorted(m,q)-T_r))./(C_2_wlf+(temperatures_unsorted(m,q)-T_r))));
%                     %t_weld = t_weld+delta/2*(a_T_prev+a_T);
%                     %a_T_prev = a_T;
%                 end
%             end
%             %temp_list(m) = t_weld;
%             lay = n;
%         end
%         temp_list_layers(m,lay) = t_weld_sort;
%     end
%     
%     if w==1
%         %temp_list1 = temp_list;
%         temp_list_sort1 = temp_list_layers;
%     end
%     
%     if w==2
%         %temp_list2 = temp_list;
%         temp_list_sort2 = temp_list_layers;
%     end
    
end

% bar(temp_list1)
% title('Normalized Weld Time of All Nodes in a Simulation')
% xlabel('Nodes')
% ylabel('Normalized Weld Time [t]')
% bar(temp_list2)
% title('Normalized Weld Time of All Nodes in a Simulation')
% xlabel('Nodes')
% ylabel('Normalized Weld Time [t]')
% tiledlayout(1,2);

% % UNDO COMMENT HERE
% bar(temp_list_sort2)
% title('Normalized Weld Time of All Nodes in a Simulation 2')
% xlabel('Nodes')
% ylabel('Normalized Weld Time [t]')
% ylim([7.28425*10^4 7.28426*10^4])

% nexttile
% bar(temp_list_sort1)
% title('Normalized Weld Time of All Nodes in a Simulation 1')
% xlabel('Nodes')
% ylabel('Normalized Weld Time [t]')
% ylim([7*10^4 8*10^4])
% bar(temp_list_sort2)

% grp = [zeros(1,1638),ones(1,1404)];
% maxtemplist1=max(temp_list1)
% mintemplist1=min(temp_list1)
% avgtemplist1=mean(temp_list1)
% 
% maxtemplist2=max(temp_list2)
% mintemplist2=min(temp_list2)
% avgtemplist2=mean(temp_list2)
% boxplot(vertcat(temp_list1,temp_list2),grp,'Labels',{'Short Print','Complex Print'})
% % boxplot([temp_list1 temp_list2],[zeros(1638,1),ones(1404,1)],'Labels',{'Short Print','Complex Print'})
% title('Normalized Weld Time of All Nodes in a Simulation')
% xlabel('Nodes')
% ylabel('Normalized Weld Time [t]')


% fid1 = fopen('/Volumes/VERBATIM HD/abaqus.rpt');
% fid2 = fopen('/Volumes/VERBATIM HD/testing_20190919/abaqus_complex1.rpt');
% 
% 
% text = textscan(fid,'%s','Delimiter','','endofline','\n');
% text1=text{1};
% % text2=text1{1};
% % test=regexp(text1,'\s*\d*\s*\d*\.(\d*)?','match'); %prior version
% test=regexp(text1,'\d+\s+\d+\.(\d*)?','match'); %new test versions
% test3=[test{:}];
% [one two]=size(test3);
% y=str2num(test3{1638});
% 
% 
% count=1;
% list = zeros(numnodes,2);
% tempcount=0;
% twonew=two-1;
% m=1;
% 
% 
% temperatures = ones(numnodes,1);
% % temperatures(1,15)=15;
% % temprow=temperatures.';
% % size(temperatures,2)
% % size(temprow,1)
% % temprow(1,:)
% % % arrayfun(@(x)find(temperatures(:,x),1,'last'),1:size(temperatures,2))
% %  find(temperatures(1,:),1,'last')
% %  temperatures(1,2)=2;
% %  find(temperatures(1,:),1,'last')
% %  find(temperatures(2,:),1,'last')
% % arrayfun(@(x)find(temprow(x,:),1,'last'),1:size(temprow,1))
% 
% 
% for n = 1:twonew
%     
%     x = str2num(test3{n});
%     count = mod(n-1,numnodes)+1;
%     
%     % x(1)==node number
%     % x(2)==temperature
%     if x(2)<=250 && x(2)>110
%         list(count,1) = x(1);
%         list(count,2) = list(count,2)+1;
%         
%         %temperature is an arrary of arrays of temperature values for each
%         %node greater than glass transition by less than or equal to 250C
%         % m = arrayfun(@(x)find(temperatures(:,x),1,'last'),1:size(temperatures,2))
%         m = find(temperatures(count,:),1,'last');
%         temperatures(count,m+1) = m+1;
%         temperatures(count,m) = x(2);
%         
%     end
%     
%     
%     
%     
% end
% 
% 
% 
% 
% temp_list = zeros(numnodes,1);
% a_T_prev = 0;
% for m = 1:numnodes
%     t_weld = 0;
%     for q = 1:list(m,2)
%         if temperatures(m,q) < 250
%             a_T = 1/(10^((-C_1_wlf*(temperatures(m,q)-T_r))./(C_2_wlf+(temperatures(m,q)-T_r))));
%             t_weld = t_weld+delta/2*(a_T_prev+a_T);
%             a_T_prev = a_T;
%         end
%     end
%     temp_list(m) = t_weld;
% end

% max(temp_list)
% min(temp_list)
% 
% bar(temp_list)
% ylim([7.6*10^4 7.8*10^4])


% test3 = textscan(test{1},'%s %s','Delimiter',' ');
% 
% 
% test3{1}


% iter=regexp(text1,'---------------------------------','token');
% nodeAndtemperature = regexp(text, '\s*\d*\s*\d*\.(\d*)?\s*','match');
% while ~feof(fid)
%     tline = fgetl(fid);
%     strfind(list, 'Step Time =');
%     U=strfind(tline, 'Step Time =')+1;
%     tline(U)
% end

% fid = fopen(filepath,'r');
% while 1
%     tline = fgetl(fid);
%     if strfind(tline, 'Step Time =')>0
%         U=strfind(tline, 'Step Time =')+1
%         tline(U)
%         sprintf('%s   %f',tline)
%     end    
%     if ~ischar(tline)
%        break
%     end
% % end
% % 
% for row in list:
%     num = row[0];
%     if num == count:
%     
%         temp = row[1];
%         a_T = 1/(10^((-C_1_wlf*(temp-T_r))./(C_2_wlf+(temp_int_34-T_r))));
%     
%         if temp >= T_g
%             t_weld = t_weld+delta/2*(a_T_prev+a_T);
%         end
%     
%         weldtimes = append(weldtimes,t_weld);
%         a_T_prev = a_T;
%         
%     end
%     
%     count = count + 1;
%     
% end
%                                                              
% % Material properties
% rho_cfabs = 1076;           % Density of 20CF ABS [kg/m3]
% c_p_cfabs=1584;             % specific heat of CF ABS [J/kg-K]
% epsilon_cfabs=0.9;          % emissivity of CF ABS []
% sigma=0.0000000567;         % Boltzmann constant [W/m2-K4]
% k_cfabs=0.262;              % thermal conductivity of CF ABS [W/m-K]
% h_air=4;                    % thermal coefficient of CFair [W/m2-K]
% T_g = 100;                  % glass transition temperature [C]
% 
% 
% 
% %% Conversions and Calculations
% v_print = 0.0254*v_print_in_s;  % print velocity [m_print/s]
% x = v_print*t_print;            % length of the print based on velocity and time to print [m_print]
% A_print = x*w;                  % top surface area of print (=bottom surface area of print) [m2]
% m_print = x*w*z*rho_cfabs;      % layer mass [kg]
% delta = t(2)-t(1);              % time interval [s]
% % t_L = t_print+t_extra;          % total time layer has to cool before next layer is extruded [s]
% % t_f = t(end);                   % end time [s]
% t_1 = 0:delta:t_L;              % time vector for layer 1 cooling time [s]
% t_2 = t_L:delta:2*t_L;          % time vector for layer 2 cooling time [s]
% t_3 = 2*t_L:delta:3*t_L;        % time vector for layer 3 cooling time [s]
% t_4 = 3*t_L:delta:4*t_L;        % time vector for layer 4 cooling time [s]
% t_5 = 4*t_L:delta:t(end);       % time vector for layer 5 cooling time [s]    

% [t_weld_12 t_weld_23 t_weld_34 t_weld_45]
% 
% %% Plotting
% T_glass = (T_g+273)*ones(length(t),1);
% figure
% plot(t_1,temp_layer1_0,'b','LineWidth',3,'HandleVisibility','off')
% hold on
% plot(t_L:delta:t(end),temp_layer1, '-b','LineWidth',3)
% plot(t_L:delta:t(end),temp_layer2, '-r','LineWidth',3)
% plot(t_L:delta:t(end), temp_int_12,'-','Color',[.8 0.1 .7],'LineWidth',3)
% plot(2*t_L:delta:t(end), temp_layer3, '-','Color',[0.9 0.9 0.1],'LineWidth',3)
% plot(2*t_L:delta:t(end), temp_int_23, '-','Color',[0.96 0.64 0.22],'LineWidth',3)
% plot(3*t_L:delta:t(end), temp_layer4, '-','Color',[0.2 0.8 0.8],'LineWidth',3)
% plot(3*t_L:delta:t(end), temp_int_34, '-','Color',[0.8 0.9 0.5],'LineWidth',3)
% plot(4*t_L:delta:t(end), temp_layer5, '-','Color',[0.84 0.37 0.81],'LineWidth',3)
% plot(4*t_L:delta:t(end), temp_int_45, '-','Color',[0.13 0.78 0.09],'LineWidth',3)
% plot(t,T_glass-273,'--k','LineWidth',3,'HandleVisibility','off')
% legend('Layer 1','Layer 2','1-2 Interface', 'Layer 3', '2-3 Interface', 'Layer 4', '3-4 Interface', 'Layer 5', '4-5 Interface')
% xlabel('Time (s)')
% ylabel(['Temperature (', char(0176),'C)'])
% xlim([0 t(end)])
% pbaspect([16 9 1])
% %title(['Interface Temperature for Layer Time of ',num2str(t_L),'s'])
% %text(550,103,'Glass Transition Temperature','FontSize',12,'FontName','Times New Roman')
% set(gca,'FontSize',20,'FontName','Times New Roman')

idx=find(test3{k},'last');
prev=test3{k}(idx);
if idx<249.9 && b<249.9 && isempty(test3{k})==0
    test3{k}=zeros(idx)
end

                                                        
                                                             
