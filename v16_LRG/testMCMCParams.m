% A script to test the number of samples needed

numSampleRange = [1e2 1e3 1e4 1e5 1e6 1e7];
evalMCMCProposalStd = 0.05;
evalMCMCInitPt = logit( ...
  lrgExp.getNormCoords([0 0.762 0.1045 0.02233 0.951 0.6845 0 1.908])  )';
evalMCMCLogJoint = @(t) evalLogJoint(logitinv(t));

for i = 1:numel(numSampleRange)

  numSamples = numSampleRange(i);
  numBurnin = numSamples/10;
  numTotal = numBurnin + numSamples;

  samples1 = CustomMCMC( numTotal, evalMCMCProposalStd, evalMCMCInitPt, ...
    evalMCMCLogJoint);
  samples1 = logitinv(samples1);
  samples1 = samples1( (numBurnin+1):end, :);
  samples2 = CustomMCMC( numTotal, evalMCMCProposalStd, evalMCMCInitPt, ...
    evalMCMCLogJoint);
  samples2 = logitinv(samples2);
  samples2 = samples2( (numBurnin+1):end, :);
  
  % Now print out the L2
  l2 = l2Divergence(samples1, samples2);
  % Print out
  fprintf('Num Samples: %d; L2 = %0.5f\n', numSamples, l2);
  fprintf('   Num unique: (%d, %d) / (%d, %d)\n', ...
    numel(unique(samples1(:,1))), numel(unique(samples2(:,1))), ...
    numel(samples1(:,1)), numel(samples2(:,2)) );
  m1 = mean(samples1);
  m2 = mean(samples2);
  fprintf('   Mean1: %s\n', mat2str(m1));
  fprintf('   Mean2: %s\n', mat2str(m2));
  fprintf('   ||m1 - m2||/sqrt(m1*m2) = %0.4f\n', ...
    norm(m1 - m2) / sqrt( norm(m1) * norm(m2) ) );

end

