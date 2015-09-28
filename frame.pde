FR F0 = F(), F1 = F(), F1rF0, E0, E0rF0, F, E, Ft, Et, Fk; // defined all frame variables
int k=7;


FR F() {return new FR();}
FR F(pt A, pt B) {return new FR(A,B);}
FR F(vec I, vec J, pt O) {return new FR(I,J, O);}
FR F(vec I, pt O) {vec J = R(I); return new FR(I,J, O);}
FR F(FR F) {return F(F.I,F.O);}
FR F(FR F0, float t, FR F1) {
  float a = angle(F0.I,F1.I);   // rotation angle
  float s = n(F1.I)/n(F0.I);    // spiral scale 
  pt G = spiralCenter(a,s,F0.O,F1.O); show(G,5); // spiral center
  vec I = S(pow(s,t),R(F0.I,t*a));               // rotated frame I axis
  pt O = P(G,pow(s,t),R(V(G,F0.O),t*a));         // rotated frame O center
  return F(I,O);                                 // the new frame
}

//compute all of this before calling similarity function
float sc =0,dis =0,angle=0; vec normalN; pt G;
void preCompute(){
  sc = Fk.scaleBetweenPoints(F0);
  vec normalN =Fk.normalCompute(F0,sc);
  //Fk.NormalCompute(F0);
  dis = Fk.displacement(F0,normalN);
  angle = Fk.angleCompute(F0);
  G = Fk.centerPoint(F0, normalN, dis, angle, sc);
  showZ(G,10);
};


void showArrow(FR F) {F.showArrow();}
void showArrows(FR F) {F.showArrows();}

class FR {
  // let 2x2 matrix M = [I, J]
  pt O; vec I; vec J; vec K;
  //ways to make new frame
  FR () {O=P(); I=V(1,0,0); J=V(0,1,0);K=V(0,0,1);}
  FR(vec II, vec JJ, pt OO) {I=V(II); J=V(JJ); O=P(OO);K = cross(I,J);}
  FR(pt A, pt B) {O=P(A); I=V(A,B); J=R(I); K = cross(I,J);
}
  
  //multiply components
  vec of(vec V) {return W(V.x,I,V.y,J,V.z,K);}    // M * V
   // vec of(vec V) {return W(V.x,I,V.y,J,V.z,cross(I,J));}    // M * V
  //pt of(pt P) {return P(O,W(P.x,I,P.y,J));} // M * V + O
  //pt of(pt P) {return P(O,W(P.x,I,P.y,J,P.z,cross(I,J)));} // M * V + O
  pt of(pt P,float zz) {  
    pt calPt = P(O,W(P.x,I,P.y,J,0,cross(I,J)));
    calPt.z += zz;
  return calPt;

} // M * V + O
  
  
  
  FR of(FR F,float zV) {  
  return F(of(F.I),of(F.J),of(F.O,zV));} // M*I, M*J, M*O
  
  //inverse function
  vec invertedOf(vec V) {    
    float a = det(V,J)/det(I,J);
    float b = det(V,I)/det(J,I);
    float c = 0;                  
  return V(a,b,c);} // M^{-1} * V // z component is doubtful det(V,K)/det(K,I)    det(V,K)/det(I,K)
  
  
  
  pt invertedOf(pt P) {
  vec V = V(O,P);
return P(det(V,J)/det(I,J),det(V,I)/det(J,I),0);

} // M^{-1} * OP
  FR invertedOf(FR F) {
  return F(invertedOf(F.I),invertedOf(F.J),invertedOf(F.O));} // M^{-1}* F.I  // M^{-1}* F.J  // M^{-1}* F.O
  
// custom functions added +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FR FSIM(FR F0, float t) {
  //float a = angle(F0.I,I);   // rotation angle
  //float s = n(I)/n(F0.I);    // spiral scale 
  //pt G = spiralCenter(a,s,F0.O,O); show(G,5); // spiral center  
  vec I = S(pow(sc,t),R(F0.I,t*angle));               // rotated frame I axis
  pt O = P(G,pow(sc,t),R(V(G,F0.O),t*angle));         // rotated frame O center
  
  return F(I,O);                                 // the new frame
}
float scaleBetweenVectors(FR F){  
  float scale = (n(J)-n(I))/(n(F.J) - n(F.I)); 
  println("scale sc = " + n(F.J));
  println("scale sc = " + n(F.I));
  println("scale sc = " + n(J));
  println("scale sc = " + n(I));
  
  if (scale != Float.NaN) {scale = 1;};
  println("scale sc = " + scale);
 return scale;  
}

float scaleBetweenPoints(FR F){  
  float a = sqrt(sq(J.x-I.x)+sq(J.y-I.y)+sq(J.z-I.z));
  float b = sqrt(sq(F.J.x-F.I.x)+sq(F.J.y-F.I.y)+sq(F.J.z-F.I.z));
  float scale= a/b;
  println("scale sc = " + scale);
 return scale;  
}

vec normalCompute(FR F, float z) {
    vec u1 = I; 
    vec u2 = J;
    vec u3 = K; 
    vec v1 = F.I; 
    vec v2 = F.J; 
    vec v3 = F.K;
    vec m1 = M(V(z, u1), v1); 
    vec m2 = M(V(z, u2), v2); 
    vec m3 = M(V(z, u3), v3); 
    vec normal = A(N(m1, m2), N(m2, m3)); 
    normal.normalize(); 
     println("type 1 : " + normal.z+ ","+normal.y+","+normal.x);
    return normal;
  }

vec NormalCompute(FR F){  
  vec delI = V(F.I.x-I.x,F.I.y-I.y,F.I.z-I.z);
  vec delJ = V(F.J.x-J.x,F.J.y-J.y,F.J.z-J.z);
  vec delK = V(F.K.x-K.x,F.K.y-K.y,F.K.z-K.z);
  
  vec A = N(delI,delJ);
  vec B = N(delJ,delK);
  
 vec normal = V(A.x+B.x,A.y+B.y,A.z+B.z);  
 println("type 2 : " + normal.z+ ","+normal.y+","+normal.x);
 normal.normalize();
 return normal;  
}


float angleCompute(FR F){
  vec A = M(J,I);
  vec B = M(F.J,F.I);  
  float dotP = dot(A,B);
  float nA = n(A);
  float nB = n(B);  
  float angle = acos(dotP/(nA*nB));
  println("Angle is : " + angle);
  return angle;
}
  /** Identity matrix minus matrix A*/
  vec[] ImA(vec[] A){
    vec[] IA = new vec[3];
    for(int i = 0; i < 3; i++) IA[i] = new vec();
     
    IA[0].setTo(A[0].mul(-1));
    IA[0].x += 1;
  
    IA[1].setTo(A[1].mul(-1));
    IA[1].y += 1;

    IA[2].setTo(A[2].mul(-1));
    IA[2].z += 1;

    return IA;
  }
  
//spiralCenter(p1, q1, N, d, Angle, Scale);
/** P_k = F + FP_0(a) + dn */ 
//pt spiralCenter(pt p0, pt pk, vec n, float d, float a, float z){
  pt centerPoint(FR F, vec N, float d, float a, float z){
    vec[] R = RotateMatrix(N, a);
    R[0].mul(z);
    R[1].mul(z);
    R[2].mul(z);
    
    vec zRp0 = Ax(R, F.O);
    vec p = new vec();
    p.x = O.x - d * N.x - zRp0.x;
    p.y = O.y - d * N.y - zRp0.y;
    p.z = O.z - d * N.z - zRp0.z;

    vec[] IR = ImA(R);
    vec po = solveAxb(IR, p);
    println("points are : " +po.x+","+po.y+","+po.z);
    return new pt(po.x, po.y, po.z);
  }//end of centerPoint
  
  /** 3x3 matrix A times 3x1 vector x */ 
  vec Ax(vec[] A, vec x){
    float a = A[0].x * x.x + A[1].x * x.y + A[2].x * x.z;
    float b = A[0].y * x.x + A[1].y * x.y + A[2].y * x.z;
    float c = A[0].z * x.x + A[1].z * x.y + A[2].z * x.z;

    return new vec(a, b, c);
  }

  vec Ax(vec[] A, pt x){
    vec y = new vec(x.x, x.y, x.z);
    return Ax(A, y);
  }

 /** construct rotation matrix from rotaion axis u and angle a */
  vec[] RotateMatrix(vec u, float a){
    float c = cos(a);
    float s = sin(a);
    vec[] R = new vec[3];
    for(int i = 0; i < 3; i++) R[i] = new vec();
    R[0].set(c + u.x*u.x*(1-c), u.y*u.x*(1-c) + u.z*s, u.z*u.x*(1-c) - u.y*s);
    R[1].set(u.x*u.y*(1-c) - u.z*s, c + u.y*u.y*(1-c), u.z*u.y*(1-c) + u.x*s);
    R[2].set(u.x*u.z*(1-c) + u.y*s, u.y*u.z*(1-c) - u.x*s, c + u.z*u.z*(1-c));
    return R;
  }


  /** determinant of a 3x3 matrix formed by 3x1 vectors v1, v2, v3 */ 
  float detM33(vec v1, vec v2, vec v3){
    return v1.x * (v2.y*v3.z - v2.z*v3.y) - v2.x * (v1.y*v3.z - v3.y*v1.z)
      + v3.x * (v1.y*v2.z - v2.y*v1.z);    
  }

  /** using cramer's rule to solve Ax = b equation */ 
  vec solveAxb(vec[] A, vec b){
    float d = detM33(A[0], A[1], A[2]);
    float x = detM33(b, A[1], A[2]) / d;
    float y = detM33(A[0], b, A[2]) / d;
    float z = detM33(A[0], A[1], b) / d;

    return new vec(x, y, z);
  }

float displacement(FR F,vec N){  
  vec C = V(O, F.O);
  float dis = dot(C,N); 
  //println("displacement is : " + dis);
  return dis;
}



// custom functions added +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                                          
                                                                                   
  FR showArrow() {
    //show(O,4); 
arrow(O,I); return this;}
  FR showArrows() {
    //show(O,4); 
arrow(O,I); arrow(O,J);
arrow(O,K); 
return this; }
}// end of frame class
