# FSKQ

This repository contains the MATLAB codes for the article  "Fully symmetric kernel quadrature" by [Toni Karvonen](https://users.aalto.fi/~karvont2/) and [Simo Särkkä](https://users.aalto.fi/~ssarkka/).

T. Karvonen and S. Särkkä (2018). SIAM Journal on Scientific Computing, 40(2):A697–A720.

[https://doi.org/10.1137/17M1121779](https://doi.org/10.1137/17M1121779)

[https://arxiv.org/abs/1703.06359](https://arxiv.org/abs/1703.06359)

CONTACT: toni.karvonen@aalto.fi

# Usage

Clone or download the repository and add `fskq` directory to your MATLAB search path. The three numerical experiments of the article (as well as some additional illustrations) are in the `examples` directory. The file `demo.m` contains some simple examples that demonstrate how to use the MATLAB functions.

# Contributors

Most of the code is written by the authors. The function `fskq/fss.m` is partly based on a MATLAB Central [answer](https://se.mathworks.com/matlabcentral/newsreader/view_thread/164470) by Bruno Luong and the function `fskq/spgetseq.m` is taken from the [Sparse Grid Interpolation Toolbox](http://www.ians.uni-stuttgart.de/spinterp/) by Andreas Klimke.
