classdef Branch
    
    properties
        busFrom;        %Starting point (bus) of the branch    
        busTo;          %Ending point (bus) of the branch
        Z;              %Impedance of the branch
        Y;              %Admittance of the branch (1/X)
    end
    
    methods (Access = public)
        function [this] = Branch(from, to, r, x)
            %Constructor
            this.busFrom = from;
            this.busTo = to;
            this.Z = Impedance(r,x);
            y = 1/this.Z.z;
            this.Y = Admittance(real(y), imag(y));
        end
    end
    
    methods(Static)
            function [Ymat] = AdmittanceMatrixGen(branch, N)
            %Description:
            %To generate a Ybus matrix using the branch data.
            %Input:
            %branch = branch data that includes all the admittances and impedances on
            %each branch
            %N = number of busses in the system
            %Output:
            %Ybus matrix
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            Ymat = zeros(N);
            for j = 1: numel(branch)
                %Off diagonal Yij = -1/z;
                Ymat(branch(j).busFrom, branch(j).busTo) = -branch(j).Y.y;
                Ymat(branch(j).busTo, branch(j).busFrom) = -branch(j).Y.y;
                %Diagonal Yii = sum(1/z);
                Ymat(branch(j).busFrom, branch(j).busFrom) = ...
                    Ymat(branch(j).busFrom, branch(j).busFrom) + branch(j).Y.y;
                Ymat(branch(j).busTo, branch(j).busTo) = ...
                    Ymat(branch(j).busTo, branch(j).busTo) + branch(j).Y.y;
            end
        end
    end
    
end

