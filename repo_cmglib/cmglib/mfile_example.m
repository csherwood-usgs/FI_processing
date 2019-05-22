function y = MFILE_EXAMPLE( x )
% MFILE_EXAMPLE - Example of basic m-file documentation standards
% y = MFILE_EXAMPLE( x )
%
% Input:
%  x - The function arguments (units)
%
% Returned:
%  y - What is returns (ideally a structure)
%
% Expanded discussion goes here...anything that should appear after the 
% help command. The first line in the file should be the function call; the
% second line should list the m-file name and one-line description. The
% third line should be an example call...often exactly line the first line,
% but not an executable.
% Things that might be included are more info on inputs and returned 
% values. If the input and return value are structures,
% addtional return info can be added later without breaking calling
% programs.
%
% Optional; related routines and dependencies can be listed.

% Blank line means that this and further
% documentation will not appear after the help command.
%
% Any copyrights, licensing statements and disclaimers go here.
%
% Authors and revision history
% csherwood@usgs.gov 
%
% TODO list...plans for future improvements
rev = '05-Oct-2006';

% handle input arguments
if(nargin>0),fprintf(1,'This mfile ignores all input.\n'),end

%return a structure
y.mfile = mfilename('fullpath');
y.rev = rev;
y.text=sprintf('This mfile contains basic examples of programming and documentation standards.\n');