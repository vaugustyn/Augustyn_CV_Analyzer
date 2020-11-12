% INTRODUCTION
% 1. Please put the '.mpt' files you want to analyze in the same folder as
% this file.
% 2. Update the mass and molecular weight
% 3. Click 'Run'
% 4. You should expect a text file that contains the analyzed outcome

% mass = 0.49*(105.53*0.933/(105.53*0.933+20.36+1056.28*0.05))/1000; %mass of the electrode (unit: g)
mass = 1;
MW = 265.81; % molecular weight of your sample 

% Load files
allFiles = dir( '*.mpt' );% get the information all .mpt file in the current folder
allNames = {allFiles(~[allFiles.isdir]).name};% get the names of the .mpt files
Num_El = numel(allNames);% count the number of files in the folder
Ewe = NaN(Num_El,1);% create a matrix to contain working electrode potential
I = NaN(Num_El,1);% create a matrix to contain current
Cath_Cap = NaN(Num_El,1);% create a matrix to contain cathodic capacity
Anod_Cap = NaN(Num_El,1);% create a matrix to contain anodic capacity
sweep_rates = NaN(Num_El,1);% create a matrix to sweep rates
 for i = 1:1:Num_El% a for loop function to extract data from each .mpt file
     % To find out headerlines and sweep rate
    textFileName = char(allNames(i)); % convert the name of the i-th file into a string format readable by Matlab
    fileID = fopen(textFileName,'rt');% open the i-th file for getting headerline and sweep rate infos
    D = textscan(fileID,'%s','Delimiter','\n');% scan the text in the i-th file and save to D
    fclose(fileID);% close the file
    Headerlines_Num = str2num(D{1}{2}(19))*10+str2num(D{1}{2}(20));% Find the number of headerlines in the 2nd row of D{1} which are stored in the 19th and 20th positions and convert the string format to a number format
        % find the sweep rate
    for sweep_rate_finder = 1:Headerlines_Num
        Sweep_Rate_Line = D{1}{sweep_rate_finder};
        Sweep_Rate_Key = 'dE/dt               ';
        Sweep_Rate_Index = strfind(Sweep_Rate_Line, Sweep_Rate_Key);
        if Sweep_Rate_Index == 1
            sweep_rate(i) = sscanf(Sweep_Rate_Line(Sweep_Rate_Index(1) + length(Sweep_Rate_Key):end), '%g', 1);% scan the line that stores the sweep rate data from the position that the Sweep_Rate_Key ends to where the line ends. The number is stored in the i-th position in the sweep_rate matrix.
        end
    end
    % to read the columns
    fileID = fopen(textFileName,'rt');% open the i-th file for getting actual data
    C = textscan(fileID,'%d %d %d %d %d %f %f %f %f %d %f %f','Headerlines',Headerlines_Num);% obtain data
    fclose(fileID);% close the file
    CC = [double(C{8}),double(C{9}),double(C{10}),(double(C{6})-12000)/3600];% save contained working electrode potential, current, cycle number, and time
    if i == 1
        CC(any(CC(:,3)==1 | CC(:,3)==2 | CC(:,3)==4, 2),:) = [];
    else
        CC(any(CC(:,3)==1 | CC(:,3)==3,2),:) = [];% delete cycle number 1 and 3
    end
    if i == 1% for the first time running the loop
        Ewe = CC(:,1);% obtain working electrode potential
        I = CC(:,2)/mass;% obtain specific current
    else% other runs
        % the following lines are able to concatenate matrices even if they
        % have different sizes
        Size_File = numel(CC(:,1));% get the size of the file
        Size_Mat = numel(Ewe(:,1));% get the size of the current matrix containing the previous data
        Max = max(Size_File, Size_Mat);% get the larger one of the two sizes above
        Ewe = [[Ewe;NaN(abs(Max - Size_Mat),i-1)],[CC(:,1);NaN(abs(Max - Size_File),1)]];% Concatenate matrices
        I = [[I;NaN(abs(Max - Size_Mat),i-1)],[CC(:,2)/mass;NaN(abs(Max - Size_File),1)]];% Concatenate matrices
    end% end if
    % for calculating cathodic and anodic capacities
    It_neg = [CC(CC(:,2)<0,4),CC(CC(:,2)<0,2)];% Get time and current that is less than 0
    It_pos = [CC(CC(:,2)>=0,4),CC(CC(:,2)>=0,2)];% Get time and current that is greater than or equal to 0
    Cath_Cap(i) = abs(trapz(It_neg(:,1), It_neg(:,2))/mass);% integration to find cathodic capacity
    Anod_Cap(i) = abs(trapz(It_pos(:,1), It_pos(:,2))/mass);% integration to find anodic capacity
    sweep_rates(i) = sweep_rate(i);% this line is stupid.
 end
 
for i = 0:length(sweep_rates)% sort capacities from small sweep rate to large
    j = 2;
    while j <= length(sweep_rates)
        if sweep_rates(j-1,1) > sweep_rates(j,1)
            temp = sweep_rates(j-1,1);
            sweep_rates(j-1,1) = sweep_rates(j,1);
            sweep_rates(j,1) = temp;
            temp_Cath_Cap = Cath_Cap(j-1);
            temp_Anod_Cap = Anod_Cap(j-1);
            Cath_Cap(j-1) = Cath_Cap(j);
            Anod_Cap(j-1) = Anod_Cap(j);
            Cath_Cap(j) = temp_Cath_Cap;
            Anod_Cap(j) = temp_Anod_Cap;
        end
        j=j+1;
    end
end
 
% calculate Coulombic efficiency, number of electrons stored, and rate
% capability
Coulombic_Efficiency = Anod_Cap./Cath_Cap;
e_stored_Cath = Cath_Cap * MW / 96485.3329 * 3.6;
e_stored_Anod = Anod_Cap * MW / 96485.3329 * 3.6;
rate_capability_Cath = Cath_Cap / max(Cath_Cap);
rate_capability_Anod = Anod_Cap / max(Anod_Cap);

%text file output
Name = char(allNames(1));
New_fileName = strrep(Name, Name(end-13:end), '_analysis.txt');% create new file names
fileID = fopen(char(New_fileName), 'w');% create a file
A = [sweep_rates,Cath_Cap,Anod_Cap,Coulombic_Efficiency, e_stored_Cath,e_stored_Anod,rate_capability_Cath,rate_capability_Anod];% create a matrix to contain the data to be exported
fprintf(fileID,'%10s\t%18s\t%24s\t%20s\t%14s\t%12s\t%27s\t%25s\r\n','Sweep Rate','Cathodic Specific Capacity', 'Anodic Specific Capacity', 'Coulombic Efficiency', '#e in Cathodic', '#e in Anodic', 'Rate Capability in Cathodic', 'Rate Capability in Anodic');% Name of the columns
fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n','mV/s', 'mAh/g', 'mAh/g', '', '','','','');% units of the columns
fprintf(fileID,'%10.1f\t%19.5f\t%25.5f\t%21.5f\t%22.5f\t%20.5f\t%28.5f\t%26.5f\r\n',A');% data exported
fclose(fileID);% close the file
