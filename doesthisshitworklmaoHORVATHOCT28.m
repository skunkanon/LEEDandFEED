%OCT 28, RUN HORVARTHSCAN FIRST
Mg_timeSPAN = Mg_scan(:,1);
Mg_rateSPAN = Mg_scan(:,2);
Mg_tempSPAN = Mg_scan(:,3);
figure(9); clf;
plot(Mg_tempSPAN, Mg_rateSPAN);

[post_Mg_tempSPAN_row, post_Mg_rateSPAN_row, ~] = reduceBGgetspectra(Mg_tempSPAN.', Mg_rateSPAN.', 420, 760, 'r');
post_Mg_tempSPAN = post_Mg_tempSPAN_row(:);
post_Mg_rateSPAN = post_Mg_rateSPAN_row(:);

%%


Mg_beta = (max(Mg_tempSPAN) - min(Mg_tempSPAN)) / (Mg_timeSPAN(end) - Mg_timeSPAN(1));
post_Mg_timeSPAN = (post_Mg_tempSPAN - post_Mg_tempSPAN(1)) / Mg_beta;
Mg_rateNORM = trapz(post_Mg_timeSPAN, post_Mg_rateSPAN);
Mg_rateSPANACTUAL = post_Mg_rateSPAN / Mg_rateNORM;


Mg_initparams = [160 * 1000, 20 * 1000 , 10^13, 2000];
Mg_N0 = 1;

figure(10); clf; hold on;
plot(post_Mg_timeSPAN, Mg_rateSPANACTUAL); %assumes initial coverage of 1
[sim_Mg_timeSPAN, ~, sim_Mg_rateSPAN, ~] = polyani_wigner_niemant(Mg_beta, post_Mg_tempSPAN(1), ...
    Mg_initparams(1), Mg_initparams(2), Mg_initparams(3), Mg_N0, post_Mg_tempSPAN(end), Mg_initparams(4));
plot(sim_Mg_timeSPAN, sim_Mg_rateSPAN);

hold off;

%%
%fit_polyani_wigner_niemant(Mg_timeSPAN, Mg_rateSPAN, Mg_initparams,Mg_beta, 300, 1000, Mg_beta);


fprintf('\n NEW INSTANCE \n');

Mg_unc = ones(size(post_Mg_rateSPAN));
Mg_exp_data = {{post_Mg_timeSPAN, Mg_rateSPANACTUAL, Mg_N0, Mg_unc}};

[~,~] = fit_polyani_wigner_niemant_tripletest(Mg_exp_data, Mg_initparams, Mg_beta, ...
    post_Mg_tempSPAN(1), post_Mg_tempSPAN(end));
