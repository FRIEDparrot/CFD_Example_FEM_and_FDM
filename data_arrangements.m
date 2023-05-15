addpath("FDM_solution\CFD_Project3");
addpath("FEM_soution\CFD_Project1");

FDM_velocity = load("FDM_solution\CFD_Project3\velocity_data.csv");
FDM_error = load("FDM_solution\CFD_Project3\error_data.csv");
FEM_velocity = readmatrix("FEM_soution\CFD_Project1\data_fem.csv");

figure("Name", "velocity graph");

hold on
x = FDM_velocity(:,2);
plot (x, FEM_velocity(:,2), 'ro-' ,'LineWidth', 1.5);
plot (x, FDM_velocity(:,3), 'b*-' ,'LineWidth', 1.5);
axis([0 1 4.8 7])
legend("FEM-result", "FDM-result");
title("Compare of the solution of the FEM and FDM result");
xlabel("location");
ylabel("Velocity");

figure("Name","error graph");
plot (FDM_error(:, 1), FDM_error(:,2),'k-','LineWidth', 1.5);
legend("Error");
title("The solution error change with the iteration times");
xlabel("iteration times")
ylabel("The sum of error")
