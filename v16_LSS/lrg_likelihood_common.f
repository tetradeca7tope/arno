      integer maxbands, bands, kbands
      parameter(maxbands=500)
      real*8 kvec(maxbands), obs(maxbands), W(maxbands,maxbands)
      common/lss/bands, kbands, kvec, obs, W 
