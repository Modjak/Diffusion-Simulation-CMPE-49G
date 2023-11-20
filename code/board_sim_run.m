%% Set Parameters
sim_params.rx_center = [0, 0, 0];
sim_params.rx_r_inMicroMeters = 5;
sim_params.rx_tx_distance = 5;
sim_params.tx_emission_pt = [10, 0, 0];
sim_params.D_inMicroMeterSqrPerSecond = 75;
sim_params.reflecting_line_eqn_A = 1;
sim_params.reflecting_line_eqn_B = 1;
sim_params.reflecting_line_eqn_C = -15;
sim_params.tend = 0.6;
sim_params.delta_t = 0.0001;
sim_params.num_molecules = 50000;
sim_params.reflection_strategy = 'reflection_wrt_line';

%% SIMULATE

fprintf('\nSimulation <sim_gaussianRW_Point2Spherical_FFP_3D> \t\t[START]')
tstart = tic;
[nrx_sim_timeline_noreflect, time] = no_reflect_2d(sim_params);
fprintf('\nSimulation <sim_gaussianRW_Point2Spherical_FFP_3D> \t\t[End] \tDuration = %f\n', toc(tstart))

[nrx_sim_timeline_reflect, ~] = sim_gaussianRW_Point2Line_FFP_2D(sim_params);

%% THEORETICAL NRX

%fprintf('\nTheoretical Formula \t\t[START]')
%tstart = tic;
%[nrx_theory_timeline] = eval_theoretical_nrx_3d_Point2Spherical_FFP_3D(sim_params, time);
%fprintf('\nTheoretical Formula  \t\t[End] \tDuration = %f\n', toc(tstart))

%% PRE PLOT
nrx_cumulative_reflect = cumsum(nrx_sim_timeline_reflect);
nrx_cumulative_noreflect = cumsum(nrx_sim_timeline_noreflect);

%nrx_theory_cumulative = cumsum(nrx_theory_timeline .* sim_params.num_molecules);
%% PLOT

hFig = figure;
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [0 101 600 400])

%plot(time, nrx_sim_timeline_merged/sim_params.num_molecules, '-', 'LineWidth', 2)

plot(time, nrx_cumulative_reflect, '-', 'LineWidth', 2)

hold on
%plot(time_merged, nrx_theory_timeline_merged, '--', 'LineWidth', 2)
%plot(time, nrx_theory_cumulative, '--', 'LineWidth', 2)
plot(time, nrx_cumulative_noreflect, '--', 'LineWidth', 2)

grid on
xlabel('Time - (s)')
ylabel('Total recieved molecules until time t (Cumulative)')

if strcmp(sim_params.reflection_strategy, 'reflection_wrt_line')
    legend('With Reflection (with respect to line)', 'Without Reflection');
elseif strcmp(sim_params.reflection_strategy, 'roll-back')
    legend('With Reflection (roll back)', 'Without Reflection');
end
%title(['\Deltat=', num2str(merge_cnt*sim_params.delta_t), '; r_{rx}=', num2str(sim_params.rx_r_inMicroMeters), '; dist=', num2str(sim_params.rx_tx_distance), '; D=', num2str(sim_params.D_inMicroMeterSqrPerSecond)])

title(['\Deltat=', num2str(sim_params.delta_t), '; r_{rx}=', num2str(sim_params.rx_r_inMicroMeters), '; dist=', num2str(sim_params.rx_tx_distance), '; D=', num2str(sim_params.D_inMicroMeterSqrPerSecond)])



