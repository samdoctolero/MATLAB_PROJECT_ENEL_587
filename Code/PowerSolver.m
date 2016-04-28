classdef PowerSolver < handle
%PowerSolver solves power flow equations by reading in data from an IEEE
%power data format. It then creates an array of bus and branch data. From
%bus and branch data it creates an admittance matrix. When these two
%variables are filled, it creates a function that has all the power flow
%equations. The power solver then calls a newton-raphson solver which reads
%values from the function that has the power flow equations. 
    
    properties
        branch;     %array of branches
        bus;        %array of busses
        Ymat;       %Admittance matrix (complex)
        Sbase;      %Base power MVA
        unknownVars %array of unknow variables
    end
    
    properties(Access = private)
        ready;      %Indicates if the PowerSolver is ready to be solved
    end
    
    methods (Access = private)
        function [branch, bus, sbase] = ReadIEEEPowerData(this, filename)
            %Description:
            %Reads the IEEE formatted power flow data from a text file.
            %Input:
            %filename = string file path and its name (i.e. \...\folder\data.txt)
            %Output:
            %bus = array of Bus objects
            %branch = array of Branch objects
            %Not used yet....
            %LossZoneData = structured loss zone data
            %InterchangeData = structured interchange data
            %TieLineData = structured tie line data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fID = fopen(filename, 'r');
            line = fgetl(fID); %TAPE line ?
            line = fgetl(fID);
            sbase = str2double(line(1,32:37));
            while(ischar(line))
                if(regexp(line,'END') > 0)
                    %If end is read then break out of the loop
                    %fprintf('File read ended...\n');
                    break;
                end
                if(regexp(line,'BUS DATA FOLLOWS') > 0)
                    %fprintf('Parsing bus data information...\n');
                    bus = this.BusInfo(fID);
                    
                elseif(regexp(line,'BRANCH DATA FOLLOWS') > 0)
                    %fprintf('Parsing branch data information...\n');
                    branch = this.BranchInfo(fID);
                    
                elseif(regexp(line,'LOSS ZONES FOLLOW') > 0)
                    %fprintf('Parsing loss zones information...\n');
                    %LossZoneData = LossZoneInfo(fID);
                    
                elseif(regexp(line,'INTERCHANGE DATA FOLLOWS') > 0)
                    %fprintf('Parsing interchange information...\n');
                    %InterchangeData = InterchangeInfo(fID);
                    
                elseif(regexp(line,'TIE LINES FOLLOW') > 0)
                    %fprintf('Parsing tie lines information...\n');
                    %TieLineData = TieLineInfo(fID);
                end
                line = fgetl(fID);
            end
            fclose(fID);
        end
        
        function [ arrBus ] = BusInfo(this,fID)
            %Description:
            %Read lines that have the Bus Data as string and formats them in a
            %array of Bus objects.
            %Input:
            %fID = file handler
            %Output:
            %arrBus = array of Bus Objects
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            arrBus = [];
            line = fgetl(fID);
            while(~strcmp(line,'-999'))
                %Bus(num,v,deg,mw,mvar,type)
                num = str2double(line(1,1:4));
                v = str2double(line(1,28:33));
                deg = str2double(line(1,34:40));
                type = 0;
                mw = 0;
                mvar = 0;
                switch str2double(line(1,25:26))
                    case {0, 1}
                        type = BusType.PQ;
                        mw = -str2double(line(1,41:49)); %Load -ve
                        mvar = -str2double(line(1,50:59)); %Load -ve
                    case 2
                        type = BusType.PV;
                        mw = str2double(line(1,60:67));
                    case 3
                        type = BusType.Slack;
                end
                data = Bus(num,v,deg,mw,mvar,type);
                if(numel(arrBus) > 0)
                    arrBus = [arrBus ; data];
                else
                    arrBus = data;
                end
                line = fgetl(fID);
                if(regexp(line,'-999')>0)
                    break;
                end
            end
        end
        
        function [ arrBranch ] = BranchInfo(this,fID)
            %Description:
            %Read lines that have the Branch data as string and formats them in an
            %array of Branch objects.
            %Input:
            %fID = file handler
            %Output:
            %BranchData = array of Branch objects
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            arrBranch = [];
            line = fgetl(fID);
            while(~strcmp(line,'-999'))
                fromBus = str2double(line(1,1:4));
                toBus = str2double(line(1,6:9));
                r = str2double(line(1,20:29));
                x = str2double(line(1,30:40));
                b = str2double(line(1,41:50));
                temp = Branch(fromBus, toBus, r, x,b);
                if (numel(arrBranch) > 0)
                    arrBranch = [arrBranch ; temp];
                else
                    arrBranch = temp;
                end
                line = fgetl(fID);
                if(regexp(line,'-999')>0)
                    break;
                end
            end
        end
        
        function [ x,vars ] = CreatePowerFlowEqs(this, funcName, initVoltage, initAngle)
            if ~exist('funcName','var')
                funcName = 'flowEquations';
            end
            if ~exist('initVoltage','var')
                initVoltage = 1;
            end
            if~exist('initVoltage','var')
                initAngle = 0;
            end
            
            tBus = this.bus;
            tY = this.Ymat;
            tSBase = this.Sbase;
            
            fID = fopen([funcName '.m'],'w');
            comment = '%';
            %Setup function keywords
            fprintf(fID,'function [ F ] = %s(x)\n',funcName);
            fprintf(fID,'%sThis is a temporary function meant to list all the\n',comment);
            fprintf(fID,'%spower flow equations.\n\n',comment);
            
            %Find how many x variables there are and list it in the function
            %Also set the initial values for the x vector
            %Set the output to what each x is supposed to be (ie. V1..VN, or
            %theta1...thetaN
            fprintf(fID,'%sList all the x variables and what they correspond to.\n',comment);
            pfEqCount = 0;
            vars = [];
            x = [];
            for i = 1:numel(tBus)
                temp = [];
                
                busNum = tBus(i).number;
                switch tBus(i).type
                    case BusType.PQ
                        pfEqCount = pfEqCount + 1;
                        temp = UnknownVariable('V',initVoltage,busNum,pfEqCount);
                        fprintf(fID,'V%g = x(%g);\n',busNum,pfEqCount);
                        if ~isempty(x)
                            x = [x; initVoltage];
                        else
                            x = initVoltage;
                        end
                        pfEqCount = pfEqCount + 1;
                        temp = [temp ; ...
                            UnknownVariable('Theta',initAngle,busNum,pfEqCount)];
                        fprintf(fID,'theta%g = x(%g);\n',busNum,pfEqCount);
                        x = [x;initAngle];
                    case BusType.PV
                        temp = UnknownVariable('Q',-1,busNum,-1);
                        pfEqCount = pfEqCount + 1;
                        temp = [temp ; UnknownVariable('Theta',initAngle,busNum,pfEqCount)];
                        fprintf(fID,'theta%g = x(%g);\n',busNum,pfEqCount);
                        if ~isempty(x)
                            x = [x;initAngle];
                        else
                            x = initAngle;
                        end
                    case BusType.Slack
                        temp = UnknownVariable('P',-1,busNum,-1);
                        temp = [temp ; UnknownVariable('Q',-1,busNum,-1)];
                end
                if ~isempty(vars)
                    vars = [vars ; temp];
                else
                    vars = temp;
                end
            end
            
            %Initialize a size for the output:
            fprintf(fID,'\nF = zeros(1,%g);\n',pfEqCount);
            
            %Now, create the power flow equations:
            fprintf(fID,'\n%sList all the equations: \n',comment);
            strActive = []; %Store all the active power eqs to print in this variable
            %And write it all at the end
            strReactive = []; %Store all the reactive power eqs to print in this variable
            %and write it all at the end
            for i = 1:numel(tBus)
                if (tBus(i).type ~= BusType.Slack) %If bus is slack then skip it
                    %Take into account the power values for PQ and PV
                    %busses
                    if(~isempty(strActive))
                        strActive = [strActive  sprintf('F(%g) = (%g) ', pfEqCount,(tBus(i).powerMW/tSBase))];
                        pfEqCount = pfEqCount - 1;
                        if(tBus(i).type == BusType.PQ)
                            strReactive = [strReactive  sprintf('F(%g) = (%g) ', pfEqCount,(tBus(i).powerMVar/tSBase))];
                            pfEqCount = pfEqCount - 1;
                        end
                    else
                        strActive = sprintf('F(%g) = (%g) ', pfEqCount,(tBus(i).powerMW/tSBase));
                        pfEqCount = pfEqCount - 1;
                        if(tBus(i).type == BusType.PQ)
                            strReactive = sprintf('F(%g) = (%g) ', pfEqCount,(tBus(i).powerMVar/tSBase));
                            pfEqCount = pfEqCount - 1;
                        end
                    end
                    %Now include the rest of the power equations:
                    % P = Vi*Vj*Yij*cos(thetaij - phi_i + phi_j)
                    % Q = Vi*Vj*Yij*sin(thetaij - phi_i + phi_j)
                    for j = 1:numel(tBus)
                        magY = abs(tY(i,j));
                        theta = angle(tY(i,j));
                        if(tBus(i).type == BusType.PV) %Voltage given
                            if(tBus(j).type == BusType.Slack) %V and Theta given
                                strActive = [strActive  sprintf('+ (%g)*(%g)*(%g)*cos(%g - theta%g + %g)...\n',tBus(i).voltage, tBus(j).voltage, magY, theta, tBus(i).number, tBus(i).angleRad)];
                            elseif(tBus(j).type == BusType.PV) %Voltage given
                                strActive = [strActive  sprintf('+ (%g)*(%g)*(%g)*cos(%g - theta%g + theta%g)...\n',tBus(i).voltage, tBus(j).voltage,magY, theta, tBus(i).number, tBus(j).number)];
                            else %V and theta not given
                                strActive = [strActive  sprintf('+ (%g)*V%g*(%g)*cos(%g - theta%g + theta%g)...\n',tBus(i).voltage, tBus(j).number, magY, theta, tBus(i).number, tBus(j).number)];
                            end
                        else %Neither voltage and phi given
                            %Only PQ has Q equations
                            if(tBus(j).type == BusType.Slack) %V and THeta given
                                strActive = [strActive  sprintf('+ V%g*(%g)*(%g)*cos(%g - theta%g + %g)...\n',tBus(i).number,tBus(j).voltage,magY, theta,tBus(i).number, tBus(j).angleRad)];
                                strReactive = [strReactive  sprintf('- V%g*(%g)*(%g)*sin(%g - theta%g + %g)...\n',tBus(i).number,tBus(j).voltage,magY, theta,tBus(i).number, tBus(j).angleRad)];
                            elseif(tBus(j).type == BusType.PV)
                                strActive = [strActive  sprintf('+ V%g*(%g)*(%g)*cos(%g - theta%g + theta%g)...\n',tBus(i).number, tBus(j).voltage,magY, theta, tBus(i).number, tBus(j).number)];
                                strReactive = [strReactive  sprintf('- V%g*(%g)*(%g)*sin(%g - theta%g + theta%g)...\n',tBus(i).number, tBus(j).voltage,magY, theta, tBus(i).number, tBus(j).number)];
                            else
                                strActive = [strActive  sprintf('+ V%g*V%g*(%g)*cos(%g - theta%g + theta%g)...\n',tBus(i).number, tBus(j).number, magY, theta, tBus(i).number, tBus(j).number)];
                                strReactive = [strReactive  sprintf('- V%g*V%g*(%g)*sin(%g - theta%g + theta%g)...\n',tBus(i).number, tBus(j).number, magY, theta, tBus(i).number, tBus(j).number)];
                            end
                        end
                    end
                    strActive = [strActive  ';\n'];
                    if(tBus(i).type == BusType.PQ)
                        strReactive = [strReactive  ';\n'];
                    end
                end
            end
            fprintf(fID,strActive);
            fprintf(fID,strReactive);
            fprintf(fID,'end\n');
            fclose(fID);
        end
    end
    
    methods (Static)
        function [x,dError] = NewtonRaphson(x,Fparam,e,h)
            %Summary:
            %Solve a nonlinear set of equations using the newton raphson method
            %It incorporates the LU function for a faster way of getting the inverse
            %Jacobian. And it uses the sparse function to compress data and save
            %storage.
            %x = initial set of values
            %Fparam = name of the function that contains the set of nonlinear equations
            %e = (optional) is the error threshold. Default is 1e-10
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            N = numel(x);
            maxIteration = 100;
            if ~exist('h','var')
                h = 1e-10; %default step size
            end
            if ~exist('e','var')
                e = 1e-10; %default threshold
            end
            done = false;
            F = feval(Fparam,x)';
            for iter = 1:maxIteration
                F = feval(Fparam,x)';
                %Check to see if each F is less than the error threshold 'e'
                for c = 1:N
                    %if at least one of them is above e then exit and update x
                    if(abs(F(c)) > e)
                        done = false;
                        break;
                    end
                    done = true;
                end
                if(done)
                    break;%if done then exit the loop
                end
                J = zeros(N) ;%Initialize jacobian
                %Fill jacobian with the proper values using numerical differentiation
                %J = F' = (F(x) - F(x+h))/h
                xh = x;
                for i = 1 :N
                    xh(i) = x(i) + h;
                    Fplus = feval(Fparam,xh)';
                    xh(i) = x(i) - 2*h;
                    Fminus = feval(Fparam,xh)';
                    %Numerical central difference formula
                    J(:,i) = (Fplus - Fminus)/(2*h);
                    xh(i) = x(i);
                end
                
                J = sparse(J); %to save memory
                
                %Update the values of x by doing x = F/J; (Jx = F)
                %Get inverse of J using LU factorization method and
                %forward and backward substitution
                [L,U,P] = lu(J); %L and U matrices for the Jacobian
                dx = x; %initialize a vector dx
                F = P*F; % Include the permutation matrix from the LU function
                d = dx; %initilize a vector d
                %Forward substitution (Ld = F)
                d(1) = F(1);
                for i = 2:N
                    
                    d(i) = F(i);
                    for j = 1: i-1
                        d(i) = d(i) -  L(i,j)*F(j);
                    end
                end
                
                %Backward substitution (Ux = d)
                dx(N) = d(N)/U(N,N);
                for i = N-1:-1:1
                    dx(i) = d(i);
                    for j = i+1:N
                        dx(i) = dx(i) - U(i,j)*dx(j);
                    end
                    dx(i) = dx(i)/U(i,i);
                    %If at least one of dx is NaN or Inf then the whole process
                    %is compromised
                    if(isnan(dx(i)) || dx(i) == Inf)
                        fprintf('Solution does not converge. Change initial values\n');
                        x = [];
                        if(nargout > 1)
                            dError = inf;
                        end
                        return;
                    end
                end
                %Update x
                x = x - dx;
            end
            
            if(~done)
                fprintf('Max iteration reached.\n');
            else
                fprintf('Method successful with %d iterations\n',iter);
            end
            if(nargout > 1)
                dError = F;
            end
        end
    end
    
    methods (Access = public)
        function [this] = PowerSolver(dataPath)
            this.ready = false;
            [this.branch, this.bus, this.Sbase] = ...
                this.ReadIEEEPowerData(dataPath);
            this.Ymat = Branch.AdmittanceMatrixGen(this.branch, numel(this.bus));
            if ~isempty(this.branch) && ~isempty(this.bus) && ~isempty(this.Ymat) &&  this.Sbase ~= -1
                this.ready = true;
            end
        end
        
        function [ output ] = Start(this,eqFuncName,initVoltage,initAngle,threshold)
            if ~exist('eqFuncName','var') || isempty(eqFuncName)
                eqFuncName = 'flowEquations';
            end
            if ~exist('initVoltage','var') || isempty(initVoltage)
                initVoltage = 1;
            end
            if ~exist('initAngle','var') || isempty(initAngle)
                initAngle = 0;
            end
            
            if ~this.ready
                output = -1;
                fprintf('Not ready to solve the power flow problem. Make sure the bus, branch, and admittance matrix are filled.\n');
                return;
            end
            %Create the power flow equations
            [initX,this.unknownVars] = this.CreatePowerFlowEqs(eqFuncName,initVoltage,initAngle);
            
            %Call the newton-raphson function to solve the set of
            %non-linear equations
            x = [];
            if ~exist('threshold','var') || isempty(threshold)
                [x] = this.NewtonRaphson(initX,eqFuncName);
            else
                [x] = this.NewtonRaphson(initX,eqFuncName,threshold);
            end
            
            if isempty(x)
                output = [];
                return;
            end
            varTemp = this.unknownVars;
            
            %Organize the output and store them in their proper bus data
            for i = 1:numel(varTemp)
                switch varTemp(i).type
                    case VariableType.ANGLE
                        this.bus(varTemp(i).busNumber).angleRad = x(varTemp(i).number);
                        this.bus(varTemp(i).busNumber).angleDeg = x(varTemp(i).number)*180/pi;
                    case VariableType.VOLTAGE
                        this.bus(varTemp(i).busNumber).voltage = x(varTemp(i).number);
                end
            end
            %Now update the power in each bus
            for i = 1:numel(varTemp)
                switch varTemp(i).type
                    case VariableType.ACTIVE_POWER
                        this.bus(varTemp(i).busNumber).powerMW = Bus.MWPowerSolver(varTemp(i).busNumber,this.bus,this.Ymat)*this.Sbase;
                    case VariableType.REACTIVE_POWER
                        this.bus(varTemp(i).busNumber).powerMVar = Bus.MVarPowerSolver(varTemp(i).busNumber,this.bus,this.Ymat)*this.Sbase;
                end
            end
            
            %Store the output in a formatted table
            out = [];
            for i = 1:numel(this.bus)
                b = this.bus(i);
                temp = struct(...
                    'BusNumber',b.number,...
                    'Voltage',b.voltage,...
                    'AngleRad',b.angleRad,...
                    'AngleDeg',b.angleDeg,...
                    'PowerMW',b.powerMW,...
                    'PowerMVar',b.powerMVar);
                if ~isempty(out)
                    out = [out;temp];
                else
                    out = temp;
                end
            end
            output = struct2table(out);
            
        end
    end
    
end

