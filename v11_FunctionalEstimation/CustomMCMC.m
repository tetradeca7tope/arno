function [samples] = CustomMCMC(numSamples, proposalStd, initPt, evalLogProb)

  num_dims = size(initPt, 1);
  
  samples = zeros(numSamples, num_dims);
  currPt = initPt;
  currLogProb = evalLogProb(currPt');

  for sample_iter = 1:numSamples
    nextPt = currPt + proposalStd * randn(num_dims, 1);
    nextLogProb = evalLogProb(nextPt');

    log_prob_ratio = nextLogProb - currLogProb;
    if (log(rand()) < log_prob_ratio) % then accept the proposal
      currPt = nextPt;
      currLogProb = nextLogProb;
    end
%     fprintf('Curr: %s, Next: %s\n', mat2str(currPt), mat2str(nextPt));
%     size(samples), size(currPt),
    samples(sample_iter, :) = currPt';
  end

end
