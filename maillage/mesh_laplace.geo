Mesh.MshFileVersion = 2.2;
L=1;
h = L/(10); //Taille caractéristique des éléménts, précision du maillage

Point(1) = {0, -L/2, 0, h};   // Construction des points ext
Point(2) = {L, -L/2, 0, h};
Point(3) = {L, L/2, 0, h};
Point(4) = {0, L/2, 0, h};

Line(12) = {1,2};   //Carre ext
Line(23) = {2,3};
Line(34) = {3,4};
Line(41) = {4,1};

Line Loop(1) = {12,23,34,41}; 
Plane Surface(1) ={1};
Physical Line("Bord")={12,23,34,41};

Physical Surface("Carre") = {1};  //A sauvegarder dans le fichier de maillage

// ctrl shift s pour fichier maillage à lire avec python