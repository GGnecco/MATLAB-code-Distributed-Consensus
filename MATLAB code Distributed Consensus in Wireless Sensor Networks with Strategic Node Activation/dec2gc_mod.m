function gcode = dec2gc_mod(dec,N)

% Description: This function converts a decimal number to its equivalent
% gray code representation.
%
% Call Syntax: [output_variables] = function_name(input_variables)
%
% Input Arguments:
%	Name: dec
%	Type: vector
%	Description: real positive numeric column vector

%	Name: N
%	Type: scalar
%	Description: precision for the gray code

% Output Arguments:
%	Name: gcode
%	Type: matrix
%	Description: gray code of size length(dec) x N as an integer matrix
%
% Creation Date: 9/17/2005
% Author : Arjun Srinivasan Rangamani

%*************************************************************************

%------------------
% Check valid input
%------------------
if (nargin ~= 2)
   error('Error (dec2gc): must have 2 input arguments.');
end   

if (isnumeric(dec) ~= 1 | isreal(dec) ~= 1 | min(dec) < 0)
    error('Error (dec2gc): input parameter must be a real positive numeric column vector ');
end

bin = dec2bin(round(dec),N); %change the real value to binary representation

len = size(bin,2);
gcode = zeros(size(bin));


for row = 1:size(bin,1)
    temp = bin(row,:);
    
    if len > 1
        
       for col  = len:-1:2            
            if (temp(col-1)-48) == 1
                gcode(row,col) = 1 - (temp(col)-48);
            else
                gcode(row,col) = (temp(col)-48);
            end            
            gcode(row,1) = (temp(1)-48);
        end
        
    else       
        gcode(row,len) = (temp-48);        
    end
    
end