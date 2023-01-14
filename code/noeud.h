class Noeud {
  // Coordonnées.
  int i;
  int j;

  double x;
  double y;

public:
  // Constructeurs:
  Noeud(int i, int j);
  Noeud(int i, int j, double x, double y);

  // Question 11:
  // Retourne le numéro global d'un noeud.
  int numgb(int N, int M);

  // Question 13:
  // Détermine les indices i et j à partir du numéro global du noeud.
  void invnumgb(int N, int M, int s);

  // Question 14:
  // Retourne le numéro intérieur d'un noeud.
  int numint(int N, int M);

  // Question 16:
  // Determine les indices i et j à partir du numéro intérieur du noeud.
  void invumint(int N, int M, int k);

  // Question 17:
  // Retourne le numéro global d'un noeud à partir de son numéro intérieur
  // Met à jour les coordonnées en fonction du numéro intérieur fourni.
  int num_int_gb(int N, int M, int k);

  // Question 18:
  // Retourne le numéro intérieur d'un noeud à partir de son numéro global.
  // Met à jour les coordonnées en fonction du numéro global fourni.
  int num_gb_int(int N, int M, int s);
};