classdef Bus < handle
    properties
        number;         %Bus number
        voltage;        %Voltage of bus in PU
        angleDeg;       %Angle of bus in degrees
        angleRad;       %Angle of bus in radians
        powerMW;        %Active power in bus in MW (negative if injected and positive if leaving)
        powerMVar;      %Reactive Power in bus in MVar (negative if injected and positive if leaving)
        type;           %Bus type (Slack, PQ, or PV)
    end
    
    methods (Access = public)
        function [this] = Bus(num, v, deg, mw, mvar, type)
            this.number = num;
            this.voltage = v;
            this.angleDeg = deg;
            this.angleRad = deg*pi/180;
            this.powerMW = mw;
            this.powerMVar = mvar;
            switch type
                case {0, 1}
                    this.type = BusType.PQ;
                case 2
                    this.type = BusType.PV;
                case 3
                    this.type = BusType.Slack;
            end
        end
    end
    
    methods (Static)
        function [P] = MWPowerSolver(i, bus, Y)
            P = 0;
            for j = 1:numel(bus)
                temp = (bus(i).voltage)*(bus(j).voltage) ...
                    *(abs(Y(i,j)))*cos(angle((Y(i,j))) ...
                    - (bus(i).angleRad) + (bus(j).angleRad));
                P = P + temp;
            end
            if(P <= 1e-10)
                P = 0;
            end
        end
        
        function [Q] = MVarPowerSolver(i, bus, Y)
            Q = 0;
            for j = 1:numel(bus)
                temp = (bus(i).voltage)*(bus(j).voltage) ...
                    *(abs(Y(i,j)))*sin(angle(Y(i,j)) ...
                    - (bus(i).angleRad) + (bus(j).angleRad));
                Q = Q - temp;
            end
            if (Q <= 1e-10)
                Q = 0;
            end
        end
    end
end

