classdef Admittance

    properties
        magnitude;          %Magnitude of the admittance (polar)
        angleDeg;           %Angle in degrees (polar)
        angleRad;           %Angle in radians (polar)
        real;               %Real part of the complex admittance
        imag;               %Imaginary part of the complex admittance
        y;                  %Admittance in rectangular format (a + ib)
    end
    
    methods
        function [this] = Admittance(real,imag)
            this.y = real + 1i*imag;
            this.magnitude = abs(this.y);
            this.angleRad = angle(this.y);
            this.angleDeg = this.angleRad*180/pi;
            this.real = real;
            this.imag = imag;
        end
    end
    
end

