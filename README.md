# B_Under_GRH
The parameters and outcome for the calculation of B is stored in the B_result.jld2 file. To perform a similar computation with different intervals for the integrals, run create_B_parameters(midpoints, Start_z) with a new list of midpoints and appropriate starting z (for m_1 = 4*10^18, this is approximately 650). To calculate B(x,y) for some values x and y, run calculate_B_from_result(x,y) to receive an upper bound on B(x,y). This function uses the parameters and calculations from the output of the create_B_parameters function (it uses the existing B_result.jld2 file by default).
