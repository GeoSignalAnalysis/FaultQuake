
clear; clc

% Sample structure
myStruct.field1 = [1; 2; 3; 5];
myStruct.field2 = [4; 5; 6; 8];
myStruct.field3 = [7; 8; 9; 10];

% Initialize the matrix
numFields = numel(fieldnames(myStruct));
numElements = numel(myStruct.field1); % Assuming all fields have the same number of elements
outputMatrix = zeros(numElements, numFields);

% Loop through the fields and save their values in the matrix
fields = fieldnames(myStruct);
for i = 1:numFields
    outputMatrix(:, i) = myStruct.(fields{i});
end

% Save the matrix to a .mat file
%save('outputMatrix.mat', 'outputMatrix');
