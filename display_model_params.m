function [] = display_model_params(model)

xlim([0,14]);
ylim([0,4]);

if(isfield(model,'mu_0'))
    text(1,2,['$\mu_0 = ' matrix_to_latex(model.mu_0,2) '$'],'Interpreter','latex');
end

if(isfield(model,'Sigma_0'))
    text(1,1,['$\Sigma_0 = ' matrix_to_latex(model.Sigma_0,2) '$'],'Interpreter','latex');
end

if(isfield(model,'Phi'))
    text(5,3,['$\Phi = ' matrix_to_latex(model.Phi,2) '$'],'Interpreter','latex');
end

if(isfield(model,'Upsilon'))
    text(5,2,['$\Upsilon = ' matrix_to_latex(model.Upsilon,2) '$'],'Interpreter','latex');
end

if(isfield(model,'Q'))
    text(5,1,['$Q = ' matrix_to_latex(model.Q,2) '$'],'Interpreter','latex');
end

if(isfield(model,'A'))
    text(9,2,['$A = ' matrix_to_latex(model.A,2) '$'],'Interpreter','latex');
end

if(isfield(model,'R'))
    text(9,1,['$R = ' matrix_to_latex(model.R,2) '$'],'Interpreter','latex');
end