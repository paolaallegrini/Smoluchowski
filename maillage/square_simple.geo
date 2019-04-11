Mesh.MshFileVersion = 2.2;
L=1000;
h = L/(64); //Taille caractéristique des éléménts, précision du maillage
/*
Point(1) = {-L/2, 0, 0, h};   // Construction des points ext
Point(2) = {L/2, 0, 0, h};
Point(3) = {L/2, L/10, 0, h};
Point(4) = {-L/2, L/10, 0, h};
*/
Point(1) = {0, 0, 0, h};   // Construction des points ext
Point(2) = {L, 0, 0, h};
Point(3) = {L, L/10, 0, h};
Point(4) = {0, L/10, 0, h};

Point(5) = {0,0,0,h}; //ligne centre
Point(6) = {0,L/10,h};

Line(12) = {1,2};   //Carre ext
Line(23) = {2,3};
Line(34) = {3,4};
Line(41) = {4,1};
Line(15) = {1,5};
Line(52) = {5,2};
Line(36) = {3,6};
Line(64) = {6,4};
Line(56) = {5,6}; //centre
Line(65) ={6,5};
/*
Line Loop(1) = {15,56,64,41};//Définition d'un pourtour (d'un bord)
Line Loop(2) = {52,23,36,65}; 
Line Loop(3) ={15,52,23,36,64};
Plane Surface(1) = {3};//Définition d'une surface
*/
Line Loop(1) = {12,23,34,41}; 
Plane Surface(1) ={1};
Physical Line("Mur")={12,34};
//Physical Line("Centre")={56};
Physical Line("Gauche")={41}; 
Physical Line("Droit")={23};

Physical Surface("Carre") = {1};  //A sauvegarder dans le fichier de maillage

// ctrl shift s pour fichier maillage à lire avec python