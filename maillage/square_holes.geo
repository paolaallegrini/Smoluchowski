Mesh.MshFileVersion = 2.2;
L=1;
h = L/8; //Taille caractéristique des éléménts, précision du maillage
hc = L/32;
Point(1) = {0, 0, 0, h};   // Construction des points ext
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};

Line(12) = {1,2};   //Carre ext
Line(23) = {2,3};
Line(34) = {3,4};
Line(41)={4,1};

Rh=0.05;
xc=0.2;
yc=0.2;

Point(5)={xc,yc,0,hc}; //centre
Point(6)={xc,yc+Rh,0,hc}; //up
Point(7)={xc-Rh,yc,0,hc}; //left
Point(8)={xc,yc-Rh,0,hc}; //down
Point(9)={xc+Rh,yc,0,hc}; //right

Circle(1)={6,5,7};
Circle(2)={7,5,8};
Circle(3)={8,5,9};
Circle(4)={9,5,6};

xc=0.2;
yc=0.5;
Point(10)={xc,yc,0,hc}; //centre
Point(11)={xc,yc+Rh,0,hc}; //up
Point(12)={xc-Rh,yc,0,hc}; //left
Point(13)={xc,yc-Rh,0,hc}; //down
Point(14)={xc+Rh,yc,0,hc}; //right

Circle(5)={11,10,12};
Circle(6)={12,10,13};
Circle(7)={13,10,14};
Circle(8)={14,10,11};


xc=0.2; yc=0.8;
Point(15)={xc,yc,0,hc}; //centre
Point(16)={xc,yc+Rh,0,hc}; //up
Point(17)={xc-Rh,yc,0,hc}; //left
Point(18)={xc,yc-Rh,0,hc}; //down
Point(19)={xc+Rh,yc,0,hc}; //right

Circle(9)={16,15,17};
Circle(10)={17,15,18};
Circle(11)={18,15,19};
Circle(122)={19,15,16};

xc=0.5; yc=0.2;
Point(20)={xc,yc,0,hc}; //centre
Point(21)={xc,yc+Rh,0,hc}; //up
Point(22)={xc-Rh,yc,0,hc}; //left
Point(23)={xc,yc-Rh,0,hc}; //down
Point(24)={xc+Rh,yc,0,hc}; //right

Circle(13)={21,20,22};
Circle(14)={22,20,23};
Circle(15)={23,20,24};
Circle(16)={24,20,21};


xc=0.5; yc=0.5;
Point(25)={xc,yc,0,hc}; //centre
Point(26)={xc,yc+Rh,0,hc}; //up
Point(27)={xc-Rh,yc,0,hc}; //left
Point(28)={xc,yc-Rh,0,hc}; //down
Point(29)={xc+Rh,yc,0,hc}; //right

Circle(17)={26,25,27};
Circle(18)={27,25,28};
Circle(19)={28,25,29};
Circle(20)={29,25,26};

xc=0.5; yc=0.8;
Point(30)={xc,yc,0,hc}; //centre
Point(31)={xc,yc+Rh,0,hc}; //up
Point(32)={xc-Rh,yc,0,hc}; //left
Point(33)={xc,yc-Rh,0,hc}; //down
Point(34)={xc+Rh,yc,0,hc}; //right

Circle(21)={31,30,32};
Circle(22)={32,30,33};
Circle(233)={33,30,34};
Circle(24)={34,30,31};

xc=0.8; yc=0.2;
Point(35)={xc,yc,0,hc}; //centre
Point(36)={xc,yc+Rh,0,hc}; //up
Point(37)={xc-Rh,yc,0,hc}; //left
Point(38)={xc,yc-Rh,0,hc}; //down
Point(39)={xc+Rh,yc,0,hc}; //right

Circle(25)={36,35,37};
Circle(26)={37,35,38};
Circle(27)={38,35,39};
Circle(28)={39,35,36};

xc=0.8; yc=0.5;
Point(40)={xc,yc,0,hc}; //centre
Point(41)={xc,yc+Rh,0,hc}; //up
Point(42)={xc-Rh,yc,0,hc}; //left
Point(43)={xc,yc-Rh,0,hc}; //down
Point(44)={xc+Rh,yc,0,hc}; //right

Circle(29)={41,40,42};
Circle(30)={42,40,43};
Circle(31)={43,40,44};
Circle(32)={44,40,41};

xc=0.8; yc=0.8;
Point(45)={xc,yc,0,hc}; //centre
Point(46)={xc,yc+Rh,0,hc}; //up
Point(47)={xc-Rh,yc,0,hc}; //left
Point(48)={xc,yc-Rh,0,hc}; //down
Point(49)={xc+Rh,yc,0,hc}; //right

Circle(33)={46,45,47};
Circle(344)={47,45,48};
Circle(35)={48,45,49};
Circle(36)={49,45,46};






Line Loop(1) = {12,23,34,41};//Définition d'un pourtour (d'un bord)

//cercles
Line Loop(2) = {1,2,3,4}; // Cercles 1
Line Loop(3) = {5,6,7,8};
Line Loop(4) = {9,10,11,122};

Line Loop(5) = {13,14,15,16}; // Cercles 2
Line Loop(6) = {17,18,19,20};
Line Loop(7) = {21,22,233,24};

Line Loop(8) = {25,26,27,28}; // Cercles 3
Line Loop(9) = {29,30,31,32};
Line Loop(10) = {33,344,35,36};


//Plane Surface(1) = {1,-2,-3,-4,-5,-6,-7,-8,-9,-10};//Définition d'une surface
//Plane Surface(1) = {1};//Définition d'une surface

//For t In {2:10}

  //  Plane Surface(t) ={t};
//EndFor



Plane Surface(1) = {1,2,3,4,5,6,7,8,9,10};//Définition d'une surface

// Bord ext=1
Physical Line("BordExt")={12,23,34,41};

// Bords cercles
Physical Line("Cercle4")={1,2,3,4,5,6,7,8,9,10,11,122,13,14,15,16,17,18,19,20,21,22,233,24,25,26,27,28,29,30,31,32,33,344,35,36};
Physical Line("C11")=2;//{1,2,3,4}; 
Physical Line("C12")={5,6,7,8};
Physical Line("C13")={9,10,11,122};

Physical Line("C21")={13,14,15,16};
Physical Line("C22")={17,18,19,20};
Physical Line("C23")={21,22,233,24};

Physical Line("C31")={25,26,27,28};
Physical Line("C32")={29,30,31,32};
Physical Line("C33")={33,344,35,36};


//Physical Surface("Carre") = {2};  //A sauvegarder dans le fichier de maillage

// ctrl shift s pour fichier maillage à lire avec python