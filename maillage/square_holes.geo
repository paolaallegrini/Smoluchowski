Mesh.MshFileVersion = 2.2;
h = 0.1; //Taille caractéristique des éléménts, précision du maillage
Point(1) = {0, 0, 0, h};   // Construction des points ext
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};

Line(12) = {1,2};   //Carre ext
Line(23) = {2,3};
Line(34) = {3,4};
Line(41)={4,1};

Rh=0.10;
xc=0.2;
yc=0.2;
Point(5)={xc,yc,0,h}; //centre
Point(6)={xc,yc+Rh,0,h}; //up
Point(7)={xc-Rh,yc,0,h}; //left
Point(8)={xc,yc-Rh,0,h}; //down
Point(9)={xc+Rh,yc,0,h}; //right

Circle(1)={6,5,7};
Circle(2)={7,5,8};
Circle(3)={8,5,9};
Circle(4)={9,5,6};

xc=0.2;
yc=0.5;
Point(10)={xc,yc,0,h}; //centre
Point(11)={xc,yc+Rh,0,h}; //up
Point(12)={xc-Rh,yc,0,h}; //left
Point(13)={xc,yc-Rh,0,h}; //down
Point(14)={xc+Rh,yc,0,h}; //right

Circle(5)={11,10,12};
Circle(6)={12,10,13};
Circle(7)={13,10,14};
Circle(8)={14,10,11};


xc=0.2;
yc=0.8;
Point(15)={xc,yc,0,h}; //centre
Point(16)={xc,yc+Rh,0,h}; //up
Point(17)={xc-Rh,yc,0,h}; //left
Point(18)={xc,yc-Rh,0,h}; //down
Point(19)={xc+Rh,yc,0,h}; //right

Circle(9)={16,15,17};
Circle(10)={17,15,18};
Circle(11)={18,15,19};
Circle(122)={19,15,16};



Line Loop(1) = {12,23,34,41};//Définition d'un pourtour (d'un bord)
Line Loop(2) = {1,2,3,4}; // Cercle intérieur
Line Loop(3) = {5,6,7,8};
Line Loop(4) = {9,10,11,122};

Plane Surface(1) = {1,-2,-3,-4};//Définition d'une surface
Physical Line("BordExt")={12,23,34,41};
Physical Line("BordInt")={1,2,3,4,5,6,7,8,9,10,11,122};

Physical Surface("Carre") = {1};  //A sauvegarder dans le fichier de maillage

// ctrl shift s pour fichier maillage à lire avec python