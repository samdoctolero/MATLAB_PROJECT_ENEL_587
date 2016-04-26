classdef Impedance
    
    properties
        magnitude;              %Polar magnitude
        angleDeg;               %Polar angle in degrees
        angleRad;               %Polar angle in radians
        real;                   %Real part of the complex impedance
        imag;                   %Imaginary part of the complex impedance
        z;                      %Complex impedance (a + ib)
    end
    
    methods 
        function [this] = Impedance(real,imag)
            this.z = real + 1i*imag;
            this.magnitude = abs(this.z);
            this.angleRad = angle(this.z);
            this.angleDeg = this.angleRad*180/pi;
            this.real = real;
            this.imag = imag;
        end
    end
end

