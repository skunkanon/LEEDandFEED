%OCT 28, RUN HORVARTHSCAN FIRST 
Mg_timeSPAN = Mg_scan(:,1);
Mg_rateSPAN = Mg_scan(:,2);
Mg_tempSPAN = Mg_scan(:,3);
figure(9); clf;
plot(Mg_tempSPAN, Mg_rateSPAN);

%%


Mg_beta = (max(Mg_tempSPAN)-min(Mg_tempSPAN)) / (max(Mg_timeSPAN));
Mg_rateNORM = trapz(Mg_timeSPAN, Mg_rateSPAN);
Mg_rateSPANACTUAL = Mg_rateSPAN / Mg_rateNORM;


Mg_initparams = [145 * 1000, 10 * 1000 , 10^13, 20000];
Mg_N0 = 1;

figure(10); clf; hold on;
plot(Mg_timeSPAN, Mg_rateSPANACTUAL); %assumes initial coverage of 1 
[sim_Mg_timeSPAN, ~, sim_Mg_rateSPAN, ~] = polyani_wigner_niemant(Mg_beta, 300, ...
    Mg_initparams(1), Mg_initparams(2), Mg_initparams(3), Mg_N0, 1000, Mg_initparams(4));
plot(sim_Mg_timeSPAN, sim_Mg_rateSPAN);

hold off;

%%
%fit_polyani_wigner_niemant(Mg_timeSPAN, Mg_rateSPAN, Mg_initparams,Mg_beta, 300, 1000, Mg_beta);


fprintf('\n NEW INSTANCE \n');

Mg_exp_data = {Mg_timeSPAN, Mg_rateSPAN, Mg_N0};
 
[~,~] = fit_polyani_wigner_niemant_tripletest(Mg_exp_data, Mg_initparams, Mg_beta, ...
    min(Mg_tempSPAN), max(Mg_timeSPAN));

