%6/5. Organized background reduction, and integration/normalization for coverage as distinct functions at this point. Testing them out, along with
%functions that implement the Polyani-Wigner equation. 

%From erley80_test: 
Ea_Pt = 47.5; %kcal/mol
kcal_to_J = 4184;
init_K = 300; %Kelvin 
pre_exp = 10^15; %s^-1 

[t_Pt_01, x_Pt_01, y_Pt_01,c_Pt_01] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,0,pre_exp,N_01,1050); %Choosing 0.1 = N_0 as the starting value, same as Cox '81
[t_Pt_015, x_Pt_015, y_Pt_015,c_Pt_015] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,0,pre_exp,N_015,1050);
[t_Pt_02, x_Pt_02, y_Pt_02,c_Pt_02] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,0,pre_exp,N_02,1050);
[t_Pt_025, x_Pt_025, y_Pt_025,c_Pt_025] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,0,pre_exp,N_025,1050);
[t_Pt_03, x_Pt_03, y_Pt_03,c_Pt_03] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,0,pre_exp,N_03,1050);
[t_Pt_035, x_Pt_035, y_Pt_035,c_Pt_035] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,0,pre_exp,N_035,1050);

figure(1); clf(1); 
cmap = parula(6);
hold on;
plot(t_Pt_01, y_Pt_01, '--', t_Pt_01, c_Pt_01, 'Color', cmap(1,:));
plot(t_Pt_015, y_Pt_015, '--', t_Pt_015, c_Pt_015,'Color', cmap(2,:));
plot(t_Pt_02, y_Pt_02, '--', t_Pt_02, c_Pt_02, 'Color', cmap(3,:));
plot(t_Pt_025, y_Pt_025, '--' , t_Pt_025, c_Pt_025,'Color', cmap(4,:));
plot(t_Pt_03, y_Pt_03, '--' , t_Pt_03, c_Pt_03,'Color', cmap(5,:));
plot(t_Pt_035, y_Pt_035, '--', t_Pt_035, c_Pt_035, 'Color', cmap(6,:));
legend('0.1 rate','0.1 coverage', '0.15 rate','0.15 coverage', '0.2 rate', '0.2 coverage', '0.25 rate', '0.25 coverage' ...
    , '0.3 rate', '0.3 coverage', '0.35 rate', '0.35 coverage');
hold off;



%%

