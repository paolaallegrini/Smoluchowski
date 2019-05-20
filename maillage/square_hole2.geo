Mesh.MshFileVersion = 2.2;
L = 1;
h=0.2;
Mesh.CharacteristicLengthFromPoints=h;
Point(1) = {0, 0, 0};   // Construction des points ext
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};

Line(12) = {1,2};   //Carre ext
Line(23) = {2,3};
Line(34) = {3,4};
Line(41)={4,1};

Rh=L/12;
xc=0.5;
yc=0.5;
Point(5)={xc,yc,0}; //centre
Point(6)={xc,yc+Rh,0}; //up
Point(7)={xc-Rh,yc,0}; //left
Point(8)={xc,yc-Rh,0}; //down
Point(9)={xc+Rh,yc,0}; //right

Circle(1)={6,5,7};
Circle(2)={7,5,8};
Circle(3)={8,5,9};
Circle(4)={9,5,6};

Line Loop(1) = {12,23,34,41};//Définition d'un pourtour (d'un bord)
Line Loop(2) = {1,2,3,4}; // Cercle intérieur
Physical Line("BordExt")={12,23,34,41};
Physical Line("BordInt")={1,2,3,4};
Plane Surface(1) = {1,-2};// Définition d'une surface
Physical Surface("Carre")=1;

Field[1] = Box;
Field[1].VIn = h/100;
Field[1].VOut = h;
Field[1].XMax = -0.5;
Field[1].XMin = 0.5;
Field[1].YMax = -0.5;
Field[1].YMin = 0.5;
//Field[1].ZMax = 275;
//Field[1].ZMin = 350;

Background Field = 1;





/*
R2=L/5;
xc=0.5;
yc=0.5;
Point(10)={xc,yc,0}; //centre
Point(11)={xc,yc+R2,0}; //up
Point(12)={xc-R2,yc,0}; //left
Point(13)={xc,yc-R2,0}; //down
Point(14)={xc+R2,yc,0}; //right

Circle(5)={11,10,12};
Circle(6)={12,10,13};
Circle(7)={13,10,14};
Circle(8)={14,10,11};


Line Loop(1) = {12,23,34,41};//Définition d'un pourtour (d'un bord)
Line Loop(2) = {1,2,3,4}; // Cercle intérieur
Line Loop(3) = {5,6,7,8}; // Cercle grand


Plane Surface(1) = {1,-3};// Définition d'une surface
Plane Surface(2) = {3,-2};// Définition donut
Physical Line("BordExt")={12,23,34,41};
Physical Line("BordInt")={1,2,3,4};
Physical Surface("Carre") = {1};  //A sauvegarder dans le fichier de maillage

// ctrl shift s pour fichier maillage à lire avec python

h = 0.3; 
h2= 0.3;

Field[1] = MathEval;
Field[1].F = "h";  // Sets the element size in region (2)

Field[2] = MathEval;
Field[2].F = "h2";  // Sets the element size in region (1)

Field[3] = Restrict;
Field[3].IField = 1;        
Field[3].FacesList = {2};

Field[4] = Restrict;
Field[4].IField = 2;
Field[4].FacesList = {1};

/*
Field[1] = BoundaryLayer;
Field[1].EdgesList = {5,6,7,8};
Field[1].AnisoMax = 1.0;
Field[1].hfar = h;
Field[1].hwall_n = h2;
Field[1].thickness = R2-Rh;
//Field[1].IntersectMetrics = 0;
BoundaryLayer Field = 1;
*/

