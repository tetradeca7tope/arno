function [next_pt, next_logjoint, sample, sample_logjoint] = ...
  MCMCPickNextSample(curr_pt, curr_logjoint, evalLogJoint, proposal_std)
% Function to pick the next point using MCMC.

  sample = logitinv( logit(curr_pt) + proposal_std*randn() );
  sample_logjoint = evalLogJoint(sample);

  sample_joint = exp(sample_logjoint);
  curr_joint = exp(curr_logjoint);

  if (rand() < sample_joint/curr_joint)
    next_pt = sample;
    next_logjoint = sample_logjoint;
  else
    next_pt = curr_pt;
    next_logjoint = curr_logjoint;
  end
end
