% 3/19 - Figure 3, Chlorine-saturated Pd(111) and Pt(111) Desorption 
fprintf('############ NEW INSTANCE ############ \n');

Ea_Pd = 60.5; %kcal/mol
Ea_Pt = 47.5; %kcal/mol
kcal_to_J = 4184;

[t_Pd,x_Pd,y_Pd] = polyani_wigner(8,550, Ea_Pd * kcal_to_J,0,10^13,1,1050);
[t_Pt,x_Pt,y_Pt] = polyani_wigner(8,550, Ea_Pt * kcal_to_J,0,10^13,1,1050);


%% 3/19 - Making Plots

figure;
hold on;
plot(x_Pd, y_Pd);
plot(x_Pt,y_Pt);

[M, I] = max(y_Pd);
plot(x_Pd(I), M, 'o', 'MarkerSize',10);





hold off;




%% 3/19 - Make function that plots desorption and labels them nicely? Making desorption spectra 
% as cells, not linear arrays

x_Erley = cell(1,2);
y_Erley = cell(1,2);

[~,x_Erley{1}, y_Erley{1}] = polyani_wigner(8,550, Ea_Pd * kcal_to_J,0,10^13,1,1050);
[~,x_Erley{2}, y_Erley{2}] = polyani_wigner(8,550, Ea_Pt * kcal_to_J,0,10^13,1,1050);




%% Testing out the plot() function 


labels = ["Saturated Pd (111)", "Saturated Pt (111)"];
tight = 'Cl Desorption at 8 K/s';


polyani_plot(x_Erley,y_Erley, labels, tight);

