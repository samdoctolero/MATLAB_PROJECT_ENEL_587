classdef UnknownVariable
    properties
        type;       %Type of variable (Power, Reactive Power, Volts, Angle)
        value;      %Value (radians, volts, MW, or MVAR)
        busNumber;  %Number of the bus at which this variable belongs
        number;     %Index of the variable in the X unknown variables
    end
    
    methods
        function [this] = UnknownVariable(type, value, busNumber, varNum)
            switch type
                case 'P'
                    this.type = VariableType.ACTIVE_POWER;
                case 'Q'
                    this.type = VariableType.REACTIVE_POWER;
                case 'V'
                    this.type = VariableType.VOLTAGE;
                case 'Theta'
                    this.type = VariableType.ANGLE;
            end
            this.value = value;
            this.busNumber = busNumber;
            this.number = varNum;
        end
    end
    
end

