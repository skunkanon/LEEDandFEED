%% General Variables
fileID = fopen('test_data.c',"w");
fileID_H = fopen('test_data.h',"w");
M_tests = [1,7,12,40,77,100,150,200];
num_values = 300;
testnum = 1;

% Write general variables to c file
fprintf(fileID, "int M_count = %d;\n", length(M_tests));
fprintf(fileID, "int num_values = %d;\n", num_values);
fprintf(fileID, "int M_values[%d] = {", length(M_tests));
for i = 1:(length(M_tests) - 1)
    fprintf(fileID, "%d,", M_tests(i));
end
fprintf(fileID,'%d',M_tests(length(M_tests)));
fprintf(fileID,'};\n');

% Write general variables to h file
fprintf(fileID_H, '#ifndef TEST_DATA_H\n#define TEST_DATA_H\n\n\n');
fprintf(fileID_H, 'extern int M_values[];\n');
fprintf(fileID_H, 'extern int M_count;\n');
fprintf(fileID_H, 'extern int num_values;\n');
fprintf(fileID_H, 'extern int num_tests;\n');

%% Test Case 1: SIN
test_case_name = "SIN";
fprintf(fileID, "\n// Test Case %d: %s\n\n", testnum, test_case_name);
t = linspace(0,10*pi,num_values);
x = sin(t);

% Write test case variables to h file
fprintf(fileID_H, "extern float test%d_x[];\n", testnum);
fprintf(fileID_H, "extern float test%d_y[][%d];\n", testnum, num_values);

% Write x array
fprintf(fileID, "float test%d_x[%d] = {", testnum, num_values);
for i = 1:(length(x) - 1)
    fprintf(fileID, "%.5f,", x(i));
end
fprintf(fileID,'%.5f',x(length(x)));
fprintf(fileID,'};\n');

% Write y array
fprintf(fileID, "float test%d_y[%d][%d] = {", testnum, length(M_tests), num_values);
for j = 1:(length(M_tests))
    M = M_tests(j);
    y = conv( ones(1,M)/M , x); % Perform Rolling Mean

    % Write Y Array to File
    if (j == 1)
        fprintf(fileID, "{");
    else
        fprintf(fileID, "\t\t\t\t\t\t\t{");
    end
    for i = 1:(num_values - 1)
        fprintf(fileID,'%.5f,',y(i));
    end
    fprintf(fileID,'%.5f',y(num_values));
    if (j ~= length(M_tests)) 
        fprintf(fileID,'},\n');
    else
        fprintf(fileID,'}');
    end
end
fprintf(fileID,'};\n\n');

testnum = testnum + 1;

%% Test Case 2: SQUARE
test_case_name = "SQUARE";
fprintf(fileID, "\n// Test Case %d: %s\n\n", testnum, test_case_name);
t = linspace(0,10,num_values);
x = t.^2;

% Write test case variables to h file
fprintf(fileID_H, "extern float test%d_x[];\n", testnum);
fprintf(fileID_H, "extern float test%d_y[][%d];\n", testnum, num_values);

% Write x array
fprintf(fileID, "float test%d_x[%d] = {", testnum, num_values);
for i = 1:(length(x) - 1)
    fprintf(fileID, "%.5f,", x(i));
end
fprintf(fileID,'%.5f',x(length(x)));
fprintf(fileID,'};\n');

% Write y array
fprintf(fileID, "float test%d_y[%d][%d] = {", testnum, length(M_tests), num_values);
for j = 1:(length(M_tests))
    M = M_tests(j);
    y = conv( ones(1,M)/M , x); % Perform Rolling Mean

    % Write Y Array to File
    if (j == 1)
        fprintf(fileID, "{");
    else
        fprintf(fileID, "\t\t\t\t\t\t\t{");
    end
    for i = 1:(num_values - 1)
        fprintf(fileID,'%.5f,',y(i));
    end
    fprintf(fileID,'%.5f',y(num_values));
    if (j ~= length(M_tests)) 
        fprintf(fileID,'},\n');
    else
        fprintf(fileID,'}');
    end
end
fprintf(fileID,'};\n\n');

testnum = testnum + 1;


%% Test Case 3: RAND
test_case_name = "RAND";
fprintf(fileID, "\n// Test Case %d: %s\n\n", testnum, test_case_name);
%t = linspace(0,10*pi,num_values);
x = rand(num_values,1);

% Write test case variables to h file
fprintf(fileID_H, "extern float test%d_x[];\n", testnum);
fprintf(fileID_H, "extern float test%d_y[][%d];\n", testnum, num_values);

% Write x array
fprintf(fileID, "float test%d_x[%d] = {", testnum, num_values);
for i = 1:(length(x) - 1)
    fprintf(fileID, "%.5f,", x(i));
end
fprintf(fileID,'%.5f',x(length(x)));
fprintf(fileID,'};\n');

% Write y array
fprintf(fileID, "float test%d_y[%d][%d] = {", testnum, length(M_tests), num_values);
for j = 1:(length(M_tests))
    M = M_tests(j);
    y = conv( ones(1,M)/M , x); % Perform Rolling Mean

    % Write Y Array to File
    if (j == 1)
        fprintf(fileID, "{");
    else
        fprintf(fileID, "\t\t\t\t\t\t\t{");
    end
    for i = 1:(num_values - 1)
        fprintf(fileID,'%.5f,',y(i));
    end
    fprintf(fileID,'%.5f',y(num_values));
    if (j ~= length(M_tests)) 
        fprintf(fileID,'},\n');
    else
        fprintf(fileID,'}');
    end
end
fprintf(fileID,'};\n\n');

testnum = testnum + 1;

%Test Case 2

%t = linspace(0,10*pi,100);
%x = t^2;
%M = 77;
%y = conv( ones(1,M)/M , x);

%Test Case 3

%t = linspace(0,10*pi,100);
%x = t^2;
%M = 40;
%y = conv( ones(1,M)/M , x);

% Create array of all tests
fprintf(fileID, "float* all_x_tests[%d] = {", testnum - 1);
for i = 1:testnum-2
    fprintf(fileID, "test%d_x,", i);
end
fprintf(fileID, "test%d_x};\n", testnum-1);

fprintf(fileID, "float (*all_y_tests[%d])[%d][%d] = {", testnum - 1, length(M_tests), num_values);
for i = 1:testnum-2
    fprintf(fileID, "&test%d_y,", i);
end
fprintf(fileID, "&test%d_y};\n", testnum-1);

fprintf(fileID_H, "extern float* all_x_tests[%d];\n", testnum - 1);
fprintf(fileID_H, "extern float (*all_y_tests[%d])[%d][%d];\n", testnum - 1, length(M_tests), num_values);

fprintf(fileID, "int num_tests = %d;\n", testnum - 1);
fprintf(fileID_H, "\n\n#endif");
fclose(fileID);
fclose(fileID_H);