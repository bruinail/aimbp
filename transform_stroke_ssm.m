function [Y, u] = transform_stroke_ssm(obs, input)
%TRANSFORM_STROKE_SSM Transforms stroke variable input into standard state
%space model form.
%   Take observed variables and input variables and convert them to a
%   standard set of observed variables and inputs for a state space model.
%   Input variables:
%   obs: A 1xn struct containing fields:
%      inr: Measured INR
%      ldl: LDL cholesterol
%      bp: Systolic/diastolic blood pressure
%      glucose: Blood glucose
%      anti_xa: Anti-Xa
%      temp: Body temperature
%      fluids: Fluid status
%      mobility: Mobility therapy
%   input: A 1xn struct containing fields:
%      aspirin:
%      statins:
%      bp_meds:
%      glucose_meds:
%      heparin:
%      fluids_io:
%      mob_therapy:

Y = [obs.inr; obs.ldl; obs.bp; obs.glucose; obs.anti_xa; obs.temp; obs.fluids; obs.mobility]';
u = [input.aspirin; input.statins; input.bp; input.glucose; input.heparin; input.fluids_io; input.mob_therapy]';
constraints = [];

end

