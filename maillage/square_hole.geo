Mesh.MshFileVersion = 2.2;
L = 1;
//h = L/8; //Taille caractéristique des éléménts, précision du maillage
hc = L/8;//L/32; // Taille caractéristique du cercle
Point(1) = {0, 0, 0, hc};   // Construction des points ext
Point(2) = {1, 0, 0, hc};
Point(3) = {1, 1, 0, hc};
Point(4) = {0, 1, 0, hc};

Line(12) = {1,2};   //Carre ext
Line(23) = {2,3};
Line(34) = {3,4};
Line(41)={4,1};

Rh=L/5;
xc=0.5;
yc=0.5;
Point(5)={xc,yc,0,hc}; //centre
Point(6)={xc,yc+Rh,0,hc}; //up
Point(7)={xc-Rh,yc,0,hc}; //left
Point(8)={xc,yc-Rh,0,hc}; //down
Point(9)={xc+Rh,yc,0,hc}; //right

Circle(1)={6,5,7};
Circle(2)={7,5,8};
Circle(3)={8,5,9};
Circle(4)={9,5,6};


Line Loop(1) = {12,23,34,41};//Définition d'un pourtour (d'un bord)
Line Loop(2) = {1,2,3,4}; // Cercle intérieur

Plane Surface(1) = {1,-2};//Définition d'une surface
Physical Line("BordExt")={12,23,34,41};
Physical Line("BordInt")={1,2,3,4};
Physical Surface("Carre") = {1};  //A sauvegarder dans le fichier de maillage

// ctrl shift s pour fichier maillage à lire avec python
/*
Field[1] = BoundaryLayer;
Field[1].EdgesList = {1,2,3,4};
Field[1].AnisoMax = 1.0;
Field[1].hfar = hc;
Field[1].hwall_n = hc/8;
Field[1].thickness = 0.1;
Field[1].IntersectMetrics = 0;
BoundaryLayer Field = 1;
*/