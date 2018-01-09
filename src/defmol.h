/* Definition of the parameters to define the molecule (defmol) */

typedef struct {
  double  x;
  double  y;
  double  z;
} Txyz;

typedef struct {
  double  dist;
  double  angle;
  double  dihedral;
} Tzmat;

typedef struct {
  int     cdist;
  int     cangle;
  int     cdihedral;
} Tcnct;

typedef struct mzmat_ {
  Tcnct   *cnct;
  Tzmat   *zmat;
} mzmat;

typedef struct {
  int     atom;
} Tatom;

typedef struct {
  int     natmrot;
  int     center;
  int     plane;
} Trtt;

typedef struct defmol_ {
  int     atmn;
  char    symb[2];
  int     np;
  float   mass;
  Txyz    *xyz;
} defmol;


