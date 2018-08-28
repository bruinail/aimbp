function [U] = make_U(model, M)
U_L = stroke_perturbation(model.U.Bmax,model.U.r_B,size(M,2));
U_M = drug_emax_model(model.U.Emax,model.U.EC50,M);
U = [U_L; U_M];