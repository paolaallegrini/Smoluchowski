Mesh.MshFileVersion = 2.2;
h = 10; //Taille caractéristique des éléménts, précision du maillage
Point(1) = {0, 0, 0, h};   // Construction des points ext
Point(2) = {100, 0, 0, h};
Point(3) = {100, 100, 0, h};
Point(4) = {0, 100, 0, h};

Line(12) = {1,2};   //Carre ext
Line(23) = {2,3};
Line(34) = {3,4};
Line(41)={4,1};

Rh=25;
xc=50;
yc=50;
Point(5)={xc,yc,0,h}; //centre
Point(6)={xc,yc+Rh,0,h}; //up
Point(7)={xc-Rh,yc,0,h}; //left
Point(8)={xc,yc-Rh,0,h}; //down
Point(9)={xc+Rh,yc,0,h}; //right

Circle(1)={6,5,7};
Circle(2)={7,5,8};
Circle(3)={8,5,9};
Circle(4)={9,5,6};


Line Loop(1) = {12,23,34,41};//Définition d'un pourtour (d'un bord)
Line Loop(2) = {1,2,3,4}; //Bord cercle

Plane Surface(1) = {1,2};//Définition d'une surface

Physical Line("Mur")={12,34};
Physical Line("Gauche")={41};
Physical Line("Droit")={23};
Physical Line("Trou")={1,2,3,4};

Physical Surface("Carre") = {1};  //A sauvegarder dans le fichier de maillage

// ctrl shift s pour fichier maillage à lire avec python