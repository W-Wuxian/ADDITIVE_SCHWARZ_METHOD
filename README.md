# Parallélisation via MPI de l'équation de la chaleur (#Parallélisation-via-MPI-de-l'équation-de-la-chaleur)
Dans le cadre du cours CHP de Madame Beaugendre à l'ENSEIRB-MATMECA.
Contributeur: V.Lederer

## Objectifs
Mise en équation par la méthode des différences finies sur une grille cartésienne 2D.
Parallélisation du code via MPI.
## Résultats
Pour le rapport voir --> [Rapport MPI](./MPI_EQ_CHALEUR.pdf)
Pour le code voir le répertoire --> [Code MPI](./MPI_EQ_CHALEUR)
Pour la compilation voir --> [Compilation MPI](./MPI_EQ_CHALEUR/README)

# ADDITIVE_SCHWARZ_METHOD
Dans le cadre du cours CHP de Madame Beaugendre à l'ENSEIRB-MATMECA.
Contributeurs: D.Sans, B.Mistral, V.Lederer

## Objectifs
ADDITIVE_SCHWARZ_DD FOR HEAT EQUATION (2D) FINITE DIFFERENCE. Il s'agit d'appliquer une stratégie de 
décomposition de domaine au problème initialement posé en [MPI_EQ](./README#Parallélisation-via-MPI-de-l'équation-de-la-chaleur)

## Résultats
### Codes
Code_Dirichlet_LAST [DIRICHLET_DD](./Code_Dirichlet_LAST) --> Pour La Décomposition de Domaine avec conditions de Dirichlet aux frontieres immergees

                    --> exécution en séquentiel ou en parallèle. Pour l'execution en sequentiel toujours mettre overlap=1
                    
DD_Neuman_TRUE_LAST [NEUMANN_DD](./DD_Neuman_TRUE_LAST) --> Pour La Décomposition de Domaine avec conditions MIXTES(Robin) aux frontieres immergees

                    --> Uniquement en Parallèle.
                    
(Par ailleurs, choisir overlap t.q le débordement du proc rang reste sur le domaine des procs rang-1 et rang+1)
### Raport
[rapport DD](./Projet_D_composition_Domaine(1).pdf)

## COMPILATION EXECUTION CHANGEMENT DES PARAMETRES:

 Pour changer les parametres tel que le nombre de points de discrétisation de l'espace
 ou le temps final et l'overlap, aller dans mod_parametres.f90
 (dépendance compilateur: mpif90, Visualisation: gnuplot)
 
    1/POUR COMPILER FAIRE: make clean;make cleanREP;make

2/Pour executer le code faire: mpirun -n X --mca pml ob1 ./run

avec X le nombre de processus (Attention: X doit etre <= Nx_g le nombre de point de calcul suivant l'axe Ox)

    3/Pour relancer le code sans changer les parametres mais sur un autre
    cas test, faire make cleanREP puis relancer l'execution.

4/la modification des parametres, tel que l'overlap, le nombre

de point de calcul (Nx_g, Ny_g), le temps final, se fait dans:
mod_parametres.f90

    5/Si vous changer les parametres, il faut alors faire:

    make clean;make cleanREP;make avant de lancer l'execution.

make cleanREP supprime le contenu des repertoires.

    #Parametres par défaut (dans mod_parametres):

Lx = Ly = 1   Longueur et Largeur

D = 1   Coef de Diffusion

Nx_g = 200  Ny_g = 100  (Nbr de points de discrétisation du domaine intérieur suivant Ox et Oy

Nt = 30 (ite temps) Tf = 2sec (temps final)

overlap = 1
