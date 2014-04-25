# Dummy likelihood computer

import sys

def f(x):
  if len(x) == 1 :
    ret =  - (2*x[0]**2) - 25;
  if len(x) == 2 :
    ret =  - (2*x[0]**2 + 3*x[1]**2) - 25;
  if len(x) == 3 :
    ret =  - (2*x[0]**2 + 3*x[1]**2 + 4*(x[2] + 0.5)**2) - 25;
  return ret

if __name__ == '__main__':

  # The first argument is the input file and the second is the output.
  inFile = sys.argv[1];
  outFile = sys.argv[2];

  # read args
  fin = open(inFile, 'r');
  inArgs = fin.read().strip();
  inArgs = inArgs.split();
  x = [float(elem) for elem in inArgs];

  # Now write to a file
  fout = open(outFile, 'w');
  result = str(f(x));
  fout.write(result);

  # close files
  fin.close();
  fout.close();

