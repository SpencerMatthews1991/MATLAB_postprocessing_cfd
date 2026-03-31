%% GCI script, version 1
% Initialization
clear all
close all
clc

%% Find directory
baseDir = fullfile(fileparts(mfilename('fullpath')), 'Code_Testing');

%% Mesh information from user inputs
Mesh_coarse = 127097;
Mesh_medium = 206108;
Mesh_fine   = 824432;

%% Import GCI mesh data from files
coarse = readtable(fullfile(baseDir, '20260213_AoA_16_GCI_coarse.csv'));
medium = readtable(fullfile(baseDir, '20260213_AoA_16_GCI_medium.csv'));
fine   = readtable(fullfile(baseDir, '20260213_AoA_16_GCI_fine.csv'));

% Extract Cl, Cd, Cl/Cd (columns 5, 6, 7)
Grid1_cl   = coarse{1, 5};
Grid1_cd   = coarse{1, 6};
Grid1_clcd = coarse{1, 7};

Grid2_cl   = medium{1, 5};
Grid2_cd   = medium{1, 6};
Grid2_clcd = medium{1, 7};

Grid3_cl   = fine{1, 5};
Grid3_cd   = fine{1, 6};
Grid3_clcd = fine{1, 7};

%% Mesh spacing (2D)
h_coarse = 1/sqrt(Mesh_coarse);
h_medium = 1/sqrt(Mesh_medium);
h_fine   = 1/sqrt(Mesh_fine);
h_array  = [h_coarse, h_medium, h_fine];

%% Refinement ratios
r_cm = sqrt(Mesh_medium / Mesh_coarse);
r_mf = sqrt(Mesh_fine   / Mesh_medium);
fprintf('=== Refinement Ratios ===\n');
fprintf('Coarse to Medium: %.4f\n', r_cm);
fprintf('Medium to Fine:   %.4f\n', r_mf);

%% Observed order of accuracy for each variable
p_cd   = log((Grid1_cd   - Grid2_cd)  /(Grid2_cd   - Grid3_cd))   / log(r_cm);
p_cl   = log((Grid1_cl   - Grid2_cl)  /(Grid2_cl   - Grid3_cl))   / log(r_cm);
p_clcd = log((Grid1_clcd - Grid2_clcd)/(Grid2_clcd - Grid3_clcd)) / log(r_cm);

fprintf('\n=== Observed Order of Accuracy ===\n');
fprintf('p (Cd):   %.4f\n', p_cd);
fprintf('p (Cl):   %.4f\n', p_cl);
fprintf('p (Cl/Cd): %.4f\n', p_clcd);

%% GCI calculations for each variable
% Using observed p for each variable
% Cd
epsilon_cd_12 = (Grid2_cd - Grid1_cd) / Grid1_cd;
epsilon_cd_23 = (Grid3_cd - Grid2_cd) / Grid2_cd;
GCI_cd_coarse = (3 * abs(epsilon_cd_12)) / (r_cm^p_cd - 1);
GCI_cd_medium = (3 * abs(epsilon_cd_23)) / (r_mf^p_cd - 1);
f_exact_cd    = Grid1_cd + (Grid1_cd - Grid2_cd) / (r_cm^p_cd - 1);
asymptotic_cd = GCI_cd_medium / (r_cm^p_cd * GCI_cd_coarse);

% Cl
epsilon_cl_12 = (Grid2_cl - Grid1_cl) / Grid1_cl;
epsilon_cl_23 = (Grid3_cl - Grid2_cl) / Grid2_cl;
GCI_cl_coarse = (3 * abs(epsilon_cl_12)) / (r_cm^p_cl - 1);
GCI_cl_medium = (3 * abs(epsilon_cl_23)) / (r_mf^p_cl - 1);
f_exact_cl    = Grid1_cl + (Grid1_cl - Grid2_cl) / (r_cm^p_cl - 1);
asymptotic_cl = GCI_cl_medium / (r_cm^p_cl * GCI_cl_coarse);

% Cl/Cd
epsilon_clcd_12 = (Grid2_clcd - Grid1_clcd) / Grid1_clcd;
epsilon_clcd_23 = (Grid3_clcd - Grid2_clcd) / Grid2_clcd;
GCI_clcd_coarse = (3 * abs(epsilon_clcd_12)) / (r_cm^p_clcd - 1);
GCI_clcd_medium = (3 * abs(epsilon_clcd_23)) / (r_mf^p_clcd - 1);
f_exact_clcd    = Grid1_clcd + (Grid1_clcd - Grid2_clcd) / (r_cm^p_clcd - 1);
asymptotic_clcd = GCI_clcd_medium / (r_cm^p_clcd * GCI_clcd_coarse);

%% Print results
fprintf('\n=== GCI Results: Cd ===\n');
fprintf('Grid1 Cd:          %.6f\n', Grid1_cd);
fprintf('Grid2 Cd:          %.6f\n', Grid2_cd);
fprintf('Grid3 Cd:          %.6f\n', Grid3_cd);
fprintf('f_exact Cd:        %.6f\n', f_exact_cd);
fprintf('GCI_coarse:        %.6f (%.4f%%)\n', GCI_cd_coarse, GCI_cd_coarse*100);
fprintf('GCI_medium:        %.6f (%.4f%%)\n', GCI_cd_medium, GCI_cd_medium*100);
fprintf('Asymptotic ratio:  %.4f (should be ~1.0)\n', asymptotic_cd);

fprintf('\n=== GCI Results: Cl ===\n');
fprintf('Grid1 Cl:          %.6f\n', Grid1_cl);
fprintf('Grid2 Cl:          %.6f\n', Grid2_cl);
fprintf('Grid3 Cl:          %.6f\n', Grid3_cl);
fprintf('f_exact Cl:        %.6f\n', f_exact_cl);
fprintf('GCI_coarse:        %.6f (%.4f%%)\n', GCI_cl_coarse, GCI_cl_coarse*100);
fprintf('GCI_medium:        %.6f (%.4f%%)\n', GCI_cl_medium, GCI_cl_medium*100);
fprintf('Asymptotic ratio:  %.4f (should be ~1.0)\n', asymptotic_cl);

fprintf('\n=== GCI Results: Cl/Cd ===\n');
fprintf('Grid1 Cl/Cd:       %.6f\n', Grid1_clcd);
fprintf('Grid2 Cl/Cd:       %.6f\n', Grid2_clcd);
fprintf('Grid3 Cl/Cd:       %.6f\n', Grid3_clcd);
fprintf('f_exact Cl/Cd:     %.6f\n', f_exact_clcd);
fprintf('GCI_coarse:        %.6f (%.4f%%)\n', GCI_clcd_coarse, GCI_clcd_coarse*100);
fprintf('GCI_medium:        %.6f (%.4f%%)\n', GCI_clcd_medium, GCI_clcd_medium*100);
fprintf('Asymptotic ratio:  %.4f (should be ~1.0)\n', asymptotic_clcd);

%% Convergence plots
figure
subplot(2,1,1)
plot(h_array, [Grid1_cl, Grid2_cl, Grid3_cl], 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8)
xlabel('Mesh Spacing h = 1/\surdN')
ylabel('C_L')
title('Lift Coefficient vs Mesh Spacing')
grid on

subplot(2,1,2)
plot(h_array, [Grid1_cd, Grid2_cd, Grid3_cd], 'ro-', 'LineWidth', 1.5, 'MarkerSize', 8)
xlabel('Mesh Spacing h = 1/\surdN')
ylabel('C_D')
title('Drag Coefficient vs Mesh Spacing')
grid on