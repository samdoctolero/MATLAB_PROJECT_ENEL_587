# MATLAB_PROJECT_ENEL_587
Newton-Raphson power flow solver using MATLAB.

Classes:
- PowerSolver
- Bus
- Branch
- Impedance
- Admittance
- UnknownVariable

Enumerations:
- BusType
- VariableType

How to:
Create a PowerSolver object and give it an path file as an input. The inputted file must be an IEEE power systems data formatted. 
Once the object is created, it will try to initialize its Bus and Branch arrays, and Admittance matrix. If that is successful,
then run Start.

Example Script:
P = PowerSolver('16bus.txt');
P.Start;
