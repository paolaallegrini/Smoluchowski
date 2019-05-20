Mesh.MshFileVersion = 2.2;
L = 1; //GetValue("Size of square's border L : ", 1);
h = GetValue("Characteristic element size h : ", L/8);
hc=2*h;//hc = GetValue("Characteristic element size around circle hc : ", L/8);   //Taille caractéristique des éléménts, précision du maillage
Rh = GetValue("Radius of the holes Rh : ", 0.05);
bool = GetValue("Boundary Layer ? ", 0);

Point(1) = {0, 0, 0, h};   // Construction des points ext
Point(2) = {L, 0, 0, h};
Point(3) = {L, L, 0, h};
Point(4) = {0, L, 0, h};

Line(12) = {1,2};   // Carre ext
Line(23) = {2,3};
Line(34) = {3,4};
Line(41) = {4,1};
theloops[0] = newreg;
Line Loop(theloops[0]) = {12,23,34,41};
Physical Line("BordExt") = {12,23,34,41};//theloops[0];

// Given the xc, yc et Rh : creates a circle of center (xc,yc) and Rh as radius
Macro Holes

    p1 = newp; Point(p1) = {xc,yc,0,hc}; //centre
    p2 = newp; Point(p2) = {xc,yc+Rh,0,hc}; //up
    p3 = newp; Point(p3) = {xc-Rh,yc,0,hc}; //left
    p4 = newp; Point(p4) = {xc,yc-Rh,0,hc}; //down
    p5 = newp; Point(p5) = {xc+Rh,yc,0,hc}; //right
    
    c1 = newreg; Circle(c1) = {p2,p1,p3};
    c2 = newreg; Circle(c2) = {p3,p1,p4};
    c3 = newreg; Circle(c3) = {p4,p1,p5};
    c4 = newreg; Circle(c4) = {p5,p1,p2};
    Printf("Cercle tag %g : LineLoop {%g, %g, %g, %g}",nbc,c1, c2, c3, c4);
    
    //Add boundary layer if bool=True
    If (bool)
        Field[nbc] = BoundaryLayer;
        Field[nbc].EdgesList = {c1,c2,c3,c4};
        Field[nbc].AnisoMax = 1.0;
        Field[nbc].hfar = h;
        Field[nbc].hwall_n = hc/32;
        Field[nbc].thickness = Rtot - Rh;
        Field[nbc].IntersectMetrics = 0;
        BoundaryLayer Field = nbc;
    EndIf
    
    l1=newreg; Line Loop(l1)={c1,c2,c3,c4};
    theloops[nbc]=-l1;
    Physical Line(nbc+1)={c1,c2,c3,c4}; //To identify if circle border in .msh
    Printf("theloops[%g] = %g",nbc,theloops[nbc]); 
Return

Ctot = 9; //nb of circles
nbc = 1;
xc = -0.1; yc = -0.1; Rtot=0.15;

//Construction of Ctot circles
For idx In {1:3}
    xc += 0.3;
    
    For idy In {1:3}
        yc += 0.3;
        Call Holes;
        nbc += 1;
    EndFor
    yc = -0.1;
EndFor
 

Plane Surface(1) = {theloops[]}; // not the same result if we replace theloops[] by its value
Physical Surface("Carre")=1;


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
