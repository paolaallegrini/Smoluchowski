Mesh.MshFileVersion = 2.2;
L=1;
h=L/8;
hc = L/32;   //Taille caractéristique des éléménts, précision du maillage

Point(1) = {0, 0, 0, h};   // Construction des points ext
Point(2) = {L, 0, 0, h};
Point(3) = {L, L, 0, h};
Point(4) = {0, L, 0, h};

Line(12) = {1,2};   // Carre ext
Line(23) = {2,3};
Line(34) = {3,4};
Line(41)={4,1};

// Given the xc, yc et Rh : creates a circle of center (xc,yc) et Rh as radius
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
    Printf("Cercle tag %g :%g, %g, %g, %g",nbc,c1, c2, c3, c4);
    
    //Creation boundary layer
    Field[nbc] = BoundaryLayer;
    Field[nbc].EdgesList = {c1,c2,c3,c4};
    Field[nbc].AnisoMax = 1.0;
    Field[nbc].hfar = h;
    Field[nbc].hwall_n = hc/32;
    Field[nbc].thickness = Rtot - Rh;
    Field[nbc].IntersectMetrics = 0;
    BoundaryLayer Field = nbc;
    
    l1=newreg; Line Loop(l1)={c1,c2,c3,c4};
    theloops[nbc]=-l1;
    Physical Line(nbc)=l1;//- theloops[nbc]; // Ctot =Max nb_holes

    Printf("l1 %g , theloops[%g] = %g",l1, nbc,theloops[nbc]);
    
Return

Ctot = 9;
nbc = 1;
xc = -0.1; yc = -0.1; Rh = 0.05; Rtot=0.1;

For idx In {1:3}
    xc += 0.3;
    
    For idy In {1:3}
        yc += 0.3;
        Call Holes;
        nbc += 1;
    EndFor
    yc = -0.1;

EndFor
 
 //Exterior border
theloops[0] = newreg;
Line Loop(theloops[0]) = {12,23,34,41};
Physical Line("BordExt") = theloops[0];
Plane Surface(1) = {theloops[]};
