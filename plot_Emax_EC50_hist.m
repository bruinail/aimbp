function [] = plot_Emax_EC50_hist(Emax_model_params, Emax_range, EC50_range)
hist3(Emax_model_params,'Edges',{Emax_range(1):((Emax_range(2)-Emax_range(1))/100):Emax_range(2),EC50_range(1):((EC50_range(2)-EC50_range(1))/100):EC50_range(2)},'CDataMode','auto');
xlabel('Emax');
ylabel('EC50');
colorbar();
view(2);