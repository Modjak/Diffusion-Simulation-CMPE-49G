function [nRx_timeline, time] = sim_gaussianRW_Point2Line_FFP_2D(sim_params)

    rx_center = sim_params.rx_center(1:2);
    rx_r_inMicroMeters = sim_params.rx_r_inMicroMeters;

    tx_emission_pt = sim_params.tx_emission_pt(1:2);
    D = sim_params.D_inMicroMeterSqrPerSecond;

    tend = sim_params.tend;
    delta_t = sim_params.delta_t;
    num_molecules = sim_params.num_molecules;

    A = sim_params.reflecting_line_eqn_A;
    B = sim_params.reflecting_line_eqn_B;
    C = sim_params.reflecting_line_eqn_C;

    strategy = sim_params.reflection_strategy;
    % Standard deviation of step size of movement N(0, sigma)
    sigma = (2 * D * delta_t)^0.5;

    % Square of the Rx Radius is useful for checking the hit action, it doesn't change so evaluating once is enough
    rx_membrane_sq = rx_r_inMicroMeters^2;

    sim_step_count = round(tend / delta_t);
    nRx_timeline = zeros(1, sim_step_count);

    % Create molecules with INITIAL Coords: replicate num_molecules times
    mol_coords_BEFORE_movement = repmat(tx_emission_pt, num_molecules, 1);

    for ii = 1:sim_step_count

        mol_displacement = normrnd(0, sigma, size(mol_coords_BEFORE_movement, 1), 2);  % 2D displacement
        
        % Evaluate signs before movement
        sign_of_mol_coords_BEFORE_movement = eval_line_sign(mol_coords_BEFORE_movement, A, B, C);
        
        % Update coordinates after movement
        mol_cords_AFTER_movement = mol_coords_BEFORE_movement + mol_displacement;
        
        % Evaluate signs after movement

        sign_of_mol_cords_AFTER_movement = eval_line_sign(mol_cords_AFTER_movement, A, B, C);
        
        hitting_to_line_signs_mult = (sign_of_mol_coords_BEFORE_movement .* sign_of_mol_cords_AFTER_movement);
        
        hitting_to_line_mask = hitting_to_line_signs_mult < 0; % 1 means hit or cross
        
        % calculates the square of the distances between molecules to the Rx Center
        dist_sq_2_rcv_center = sum((mol_cords_AFTER_movement - rx_center(1:2)).^2, 2);

        % checks if the molecules are outside of the receiver (NOT HIT)
        outside_rx_membrane_mask = dist_sq_2_rcv_center > rx_membrane_sq;

        % calculate the hit of molecules to the receiver (FOR THIS SIM STEP, NOT TOTAL)
        nRx_timeline(ii) = nRx_timeline(ii) + nnz(~outside_rx_membrane_mask);

        if strcmp(strategy, 'roll-back')
            mol_cords_AFTER_movement(hitting_to_line_mask, :) = mol_coords_BEFORE_movement(hitting_to_line_mask, :);
        elseif strcmp(strategy, 'reflection_wrt_line')
            
            x_prime = mol_cords_AFTER_movement(hitting_to_line_mask, 1);
            y_prime = mol_cords_AFTER_movement(hitting_to_line_mask, 2);
            
            x_reflect = ((B^2 - A^2) * x_prime - 2 * A * (B * y_prime + C)) / (A^2 + B^2);
            y_reflect = ((A^2 - B^2) * y_prime - 2 * B * (A * x_prime + C)) / (A^2 + B^2);


            mol_cords_AFTER_movement(hitting_to_line_mask, :) = 2 .* [x_reflect, y_reflect] - mol_coords_BEFORE_movement(hitting_to_line_mask, :);
        end
     
        % Just keep the ones that did NOT hit and do the same things again
        mol_coords_BEFORE_movement = mol_cords_AFTER_movement(outside_rx_membrane_mask, :);

    end

    time = delta_t:delta_t:tend;
end

function signs = eval_line_sign(coords, A, B, C)
    % Evaluate the sign of coordinates with respect to the line Ax + By + C = 0
    % Returns -1 for points below the line, 0 for points on the line, and 1 for points above the line
    
    distances = A * coords(:, 1) + B * coords(:, 2) + C;
    signs = sign(distances);
end

