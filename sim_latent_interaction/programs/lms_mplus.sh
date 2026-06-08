echo "
TITLE:
  LMS Test - Latent by Categorical Interactions

DATA:
  FILE IS ${1}/mplus_data.txt;

VARIABLE:
  NAMES ARE G12 X1-X${2} Y1-Y${2} G;
  USEVARIABLES ARE G X1-X${2} Y1-Y${2};
  CATEGORICAL IS G;

MODEL:
  ! measurement models
  X BY X1* X2-X${2};    ! free first loading (intercepts free by default)
  X@1;                  ! fix factor variance to 1 (mean defaults to 0)

  Y BY Y1-Y${2};        ! first loading fixed at 1 (intercepts free by default)
  [Y1@0];               ! fix first measurement intercept to 0

  ! structural model
  XG | X XWITH G;       ! latent interactions
  Y ON X G XG;          ! structural model
  [Y*];                 ! estimate structural intercept

  G ON X;               ! link latent X to G via probit regression

ANALYSIS:
  TYPE = RANDOM;
  ALGORITHM = EM;
  INTEGRATION = MONTECARLO;
  LINK = PROBIT;

OUTPUT:
  SAMPSTAT TECH1 TECH8 STDYX;

SAVEDATA:
  results = ${1}/estimates.dat;
"
