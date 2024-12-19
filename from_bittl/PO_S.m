(* ::Package:: *)

(* ::Input:: *)
(*sx={{0,1},{1,0}}/2*)

(* ::Input:: *)
(*sy={{0,-I},{I,0}}/2*)

(* ::Input:: *)
(*sz={{1,0},{0,-1}}/2*)

(* ::Input:: *)
(*e= {{1,0},{0,1}}*)

(* ::Input:: *)
(*{{1,0},{0,1}}*)

(* ::Input:: *)
(*z={{0,0},{0,0}}*)

(* ::Input:: *)
(*Sza = TensorProduct[sz,e]//ArrayFlatten*)

(* ::Input:: *)
(*Szb= TensorProduct[e, sz]//ArrayFlatten*)

(* ::Input:: *)
(*Sxa=TensorProduct[sx,e]//ArrayFlatten*)

(* ::Input:: *)
(*Sxb= TensorProduct[e, sx]//ArrayFlatten*)

(* ::Input:: *)
(*Sya=TensorProduct[sy,e]//ArrayFlatten*)

(* ::Input:: *)
(*Syb= TensorProduct[e, sy]//ArrayFlatten*)

(* ::Input:: *)
(*EE= TensorProduct[e, e]/2//ArrayFlatten*)

(* ::Input:: *)
(*ZZ= TensorProduct[z, z]//ArrayFlatten*)

(* ::Input:: *)
(*SxSx=2TensorProduct[sx,sx]//ArrayFlatten*)

(* ::Input:: *)
(*SxSy=2TensorProduct[sx,sy]//ArrayFlatten*)

(* ::Input:: *)
(*SxSz=2TensorProduct[sx,sz]//ArrayFlatten*)

(* ::Input:: *)
(*SySx=2TensorProduct[sy,sx]//ArrayFlatten*)

(* ::Input:: *)
(*SySy=2TensorProduct[sy,sy]//ArrayFlatten*)

(* ::Input:: *)
(*SySz=2TensorProduct[sy,sz]//ArrayFlatten*)

(* ::Input:: *)
(*SzSx=2TensorProduct[sz,sx]//ArrayFlatten*)

(* ::Input:: *)
(*SzSy=2TensorProduct[sz,sy]//ArrayFlatten*)

(* ::Input:: *)
(*SzSz=2TensorProduct[sz,sz]//ArrayFlatten*)

(* ::Input:: *)
(*(* S1= {Sxa, Sya,Sza} *)*)

(* ::Input:: *)
(*(* S2= {Sxb, Syb,Szb} *)*)

(* ::Input:: *)
(*ProductBasis = {EE, Sxa,Sya,Sza,Sxb,Syb,Szb,SxSx,SxSy,SxSz,SySx,SySy,SySz,SzSx,SzSy,SzSz}*)

(* ::Input:: *)
(*(* ProductBasisAlpha = {"EE", "S1x","S1y","S1z","S2x","S2y","S2z","S1xS2x","S1xS2y","S1xS2z","S1yS2x","S1yS2y","S1yS2z","S1zS2x","S1zS2y","S1zS2z"} *)*)

(* ::Input:: *)
(*ProductBasisAlpha = {"EE", "Sxa","Sya","Sza","Sxb","Syb","Szb","SxSx","SxSy","SxSz","SySx","SySy","SySz","SzSx","SzSy","SzSz"}*)

(* ::Input:: *)
(*ProductBasisIndex=Table[i, {i,Length[ProductBasis]}]*)

(* ::Input:: *)
(*MatrixTrace[a_]:= Module[{sum=0,i}, For[i=1,i<=Length[a[[1]]], i++, sum=sum+a[[i,i]]]; sum]*)

(* ::Input:: *)
(*BasisProject[a_]:=Module[{i,out=Table[0, {i,Length[ProductBasis]}]}, For[i=1,i<= Length[ProductBasis],i++, out[[i]]=MatrixTrace[ProductBasis[[i]] . a]]; out]*)

(* ::Input:: *)
(*BasisProjectAlpha[a_] :=Simplify[ BasisProject[a] . ProductBasisAlpha] *)

(* ::Input:: *)
(*comm[a_,b_]:= a . b-b . a*)

(* ::Input:: *)
(*cComm[a_,b_]:= comm[a,comm[a,b]]*)

(* ::Input:: *)
(*(* PropagateBaseOp[a_,b_, \[Xi]_]:= Module[{index=BasisProject[b].ProductBasisIndex, out= ZZ},If[Simplify[comm[a,b]]==ZZ,out=a,If[index<8, out=MatrixExp[-I*\[Xi]*b].a.MatrixExp[I*\[Xi]*b],out=MatrixExp[-I*\[Xi]/2*b].a.MatrixExp[I*\[Xi]/2*b]]]; out= ExpToTrig[out];TrigReduce[out]] *)*)

(* ::Input:: *)
(*PropagateBaseOp[a_,b_, \[Xi]_]:= Module[{index=BasisProject[b] . ProductBasisIndex, out= ZZ},If[Simplify[comm[a,b]]==ZZ,out=a,out=MatrixExp[-I*\[Xi]*b] . a . MatrixExp[I*\[Xi]*b]; out= ExpToTrig[out]];TrigReduce[out]]*)

(* ::Input:: *)
(*PropagateOp[a_,b_, \[Xi]_]:= Module[{prop=ZZ, coff=0},For[i=1,i<=16,i++, coff=BasisProject[a][[i]];If[ToString[coff]!="0",  prop= prop+coff*PropagateBaseOp[ProductBasis[[i]],b,\[Xi]]]]; prop]*)

(* ::Input:: *)
(*rho0=-(SxSx+SySy+SzSz)/2 (* Q_S = 1/4 - S_a.S_b; EE/2 omitted! *) *)

(* ::Input:: *)
(*rho0A=BasisProjectAlpha[rho0]*)

(* ::Input:: *)
(*(* H0a= Sza*Omegaa+ Szb*Omegab+ dmJ*SzSz- dp2J*(SxSx+SySy) *)*)

(* ::Input:: *)
(*(* Eigensystem[H0a] *)*)

(* ::Input:: *)
(*(* Reorder= {{0,1,0,0}, {0,0,0,1},{0,0,1,0},{1,0,0,0}} *)*)

(* ::Input:: *)
(*(* H0aEVa//FullSimplify *)*)

(* ::Input:: *)
(*(* H0aEVaR=Reorder.H0aEVa//FullSimplify *)*)

(* ::Input:: *)
(*(* H0aEVeR=Reorder.H0aEVe *)*)

(* ::Input:: *)
(*H0=q(Sza-Szb)+ dmJ*SzSz (* \[Omega]_s(Sza+Szb) omitted *)*)

(* ::Input:: *)
(*rho0e=PropagateOp[rho0,(SxSy-SySx), \[Xi]/2]*)

(* ::Input:: *)
(*rho0eA=BasisProjectAlpha[rho0e]//FullSimplify*)

(* ::Input:: *)
(*rho1e=PropagateOp[rho0e,Sza,q*T]//Simplify*)

(* ::Input:: *)
(*rho1e=PropagateOp[rho1e,Szb,-q*T]//Simplify*)

(* ::Input:: *)
(*rho1e=PropagateOp[rho1e,SzSz,dmJ*T]//FullSimplify*)

(* ::Input:: *)
(*rho1woZQCe = T/(2\[Pi])Integrate[ rho1e, {q,0,2\[Pi]/T}]- MatrixTrace[rho1e]/MatrixTrace[EE]EE//FullSimplify (* ZQCs and trace omitted *)*)

(* ::Input:: *)
(* rho1woZQC=PropagateOp[rho1woZQCe, (SxSy-SySx), -\[Xi]/2]//FullSimplify *)

(* ::Input:: *)
(*rho2=PropagateOp[rho1woZQC, (Sxa+Sxb), \[Beta]]//FullSimplify*)

(* ::Input:: *)
(*rho2e=PropagateOp[rho2, (SxSy-SySx), \[Xi]/2]//FullSimplify*)

(* ::Input:: *)
(*coeffs= BasisProject[rho2e]//FullSimplify*)

(* ::Input:: *)
(*Y1=coeffs[[3]]//TrigExpand//FullSimplify*)

(* ::Input:: *)
(*Y2=coeffs[[6]]// TrigExpand//FullSimplify*)

(* ::Input:: *)
(*Y2==-Y1//FullSimplify*)

(* ::Input:: *)
(*Y1 ==( Cos[\[Xi]/2]Sin[2\[Xi]]Sin[\[Beta]] - Cos[\[Xi]]^2 Sin[\[Xi]/2]Sin[2\[Beta]])/4//FullSimplify*)

(* ::Input:: *)
(*Y1 =( Cos[\[Xi]/2]Sin[2\[Xi]]Sin[\[Beta]] - Cos[\[Xi]]^2 Sin[\[Xi]/2]Sin[2\[Beta]])/4*)

(* ::Input:: *)
(*Y2 == -( Cos[\[Xi]/2]Sin[2\[Xi]]Sin[\[Beta]] - Cos[\[Xi]]^2 Sin[\[Xi]/2]Sin[2\[Beta]])/4//TrigExpand//Simplify*)

(* ::Input:: *)
(*Y2 = -( Cos[\[Xi]/2]Sin[2\[Xi]]Sin[\[Beta]] - Cos[\[Xi]]^2 Sin[\[Xi]/2]Sin[2\[Beta]])/4*)

(* ::Input:: *)
(*YZ= coeffs[[13]]*)

(* ::Input:: *)
(*ZY== ZY //FullSimplify*)

(* ::Input:: *)
(*YZ ==(Cos[\[Xi]/2]Cos[\[Xi]]^2 Sin[2\[Beta]]/4+Sin[\[Xi]/2]Sin[2\[Xi]]Sin[\[Beta]]/4)//FullSimplify//Expand//Simplify*)

(* ::Input:: *)
(*YZ =(Cos[\[Xi]/2]Cos[\[Xi]]^2 Sin[2\[Beta]]/4+Sin[\[Xi]/2]Sin[2\[Xi]]Sin[\[Beta]]/4)*)

(* ::Input:: *)
(*ZY=YZ*)

(* ::Input:: *)
(*(* Sya-Syb Path *)*)

(* ::Input:: *)
(*rho3e=PropagateOp[(Sya-Syb), SzSz, dmJ]//FullSimplify*)

(* ::Input:: *)
(*rho3e=PropagateOp[rho3e, (Sza-Szb), q]//FullSimplify*)

(* ::Input:: *)
(*rho3= PropagateOp[rho3e,(SxSy-SySx),-\[Xi]/2]//ExpToTrig//FullSimplify*)

(* ::Input:: *)
(*rho4=PropagateOp[rho3 ,(Sxa+Sxb), \[Pi]]//FullSimplify*)

(* ::Input:: *)
(*rho4e= PropagateOp[rho4,(SxSy-SySx),\[Xi]/2]//ExpToTrig//FullSimplify*)

(* ::Input:: *)
(*rho5e= PropagateOp[rho4e, SzSz,dmJ]//FullSimplify*)

(* ::Input:: *)
(*rho5e= PropagateOp[rho5e,(Sza-Szb),q]//FullSimplify*)

(* ::Input:: *)
(*rho5= PropagateOp[rho5e,(SxSy-SySx),-\[Xi]/2]//ExpToTrig// FullSimplify*)

(* ::Input:: *)
(*rho5Integ=Integrate[rho5, {q,0,2\[Pi]}]/(2\[Pi])//FullSimplify*)

(* ::Input:: *)
(*coeffs=BasisProject[rho5Integ]//ExpToTrig //FullSimplify*)

(* ::Input:: *)
(*coeffs[[2]]==coeffs[[5]]*)

(* ::Input:: *)
(*coeffs[[3]]==-coeffs[[6]]*)

(* ::Input:: *)
(*XpXrefYmYcoeff= coeffs[[2]]//ExpToTrig //FullSimplify*)

(* ::Input:: *)
(*XpXrefYmYcoeff== -Sin[2 dmJ]Cos[\[Xi]]Sin[\[Xi]/2]//FullSimplify*)

(* ::Input:: *)
(*XpXrefYmYcoeff= -Sin[2 dmJ]Cos[\[Xi]]Sin[\[Xi]/2]*)

(* ::Input:: *)
(*XpXYmYpathcoeff= Y1*YmYrefXpXcoeff*)

(* ::Input:: *)
(*XpXYmYpathcoeff==-1/4(Cos[\[Xi]]^2Sin[\[Xi]]^2Sin[2dmJ] Sin[\[Beta]]-Cos[\[Xi]]^3 Sin[\[Xi]/2]^2 Sin[2dmJ] Sin[2\[Beta]])//FullSimplify*)

(* ::Input:: *)
(*XpXYmYpathcoeff=-1/4(Cos[\[Xi]]^2Sin[\[Xi]]^2Sin[2dmJ] Sin[\[Beta]]-Cos[\[Xi]]^3 Sin[\[Xi]/2]^2 Sin[2dmJ] Sin[2\[Beta]])*)

(* ::Input:: *)
(*(* SySz+ SzSy path *)*)

(* ::Input:: *)
(*\[Rho]=PropagateOp[SySz+SzSy, SzSz,dmJ]*)

(* ::Input:: *)
(*\[Rho] = PropagateOp[\[Rho], (Sza-Szb),q]//FullSimplify*)

(* ::Input:: *)
(*\[Rho] = PropagateOp[\[Rho], (SxSy-SySx), -\[Xi]/2]//FullSimplify*)

(* ::Input:: *)
(*\[Rho] = PropagateOp[\[Rho], (Sxa+Sxb), \[Pi]]//FullSimplify*)

(* ::Input:: *)
(*\[Rho] = PropagateOp[\[Rho], (SxSy-SySx), \[Xi]/2]//FullSimplify*)

(* ::Input:: *)
(*\[Rho]=PropagateOp[\[Rho], SzSz,dmJ]//FullSimplify*)

(* ::Input:: *)
(*\[Rho]=PropagateOp[\[Rho], (Sza-Szb),q]//FullSimplify*)

(* ::Input:: *)
(*\[Rho]=PropagateOp[\[Rho],(SxSy-SySx), -\[Xi]/2]//FullSimplify*)

(* ::Input:: *)
(*\[Rho]Int=Integrate[\[Rho], {q,0,2\[Pi]}]/(2\[Pi])//ExpToTrig// Simplify*)

(* ::Input:: *)
(*coeffs= BasisProject[\[Rho]Int]//FullSimplify*)

(* ::Input:: *)
(*coeffs[[2]]==coeffs[[5]]*)

(* ::Input:: *)
(*coeffs[[3]]==-coeffs[[6]]*)

(* ::Input:: *)
(*XpXYZpZYcoeff=coeffs[[2]]*)

(* ::Input:: *)
(*XpXYZpZYpathcoeff=XpXYZpZYcoeff*YZ*)

(* ::Input:: *)
(*XpXYZpZYpathcoeff==-1/4(Cos[\[Xi]]^2Sin[\[Xi]]^2Sin[\[Beta]]+  Cos[\[Xi]/2]^2 Cos[\[Xi]]^3 Sin[2\[Beta]])Sin[2dmJ]//FullSimplify*)

(* ::Input:: *)
(*XpXYZpZYpathcoeff=-1/4(Cos[\[Xi]]^2Sin[\[Xi]]^2Sin[\[Beta]]+  Cos[\[Xi]/2]^2 Cos[\[Xi]]^3 Sin[2\[Beta]])Sin[2dmJ]*)

(* ::Input:: *)
(*XpXcombinedcoeff=XpXYmYpathcoeff+XpXYZpZYpathcoeff*)

(* ::Input:: *)
(*XpXcombinedcoeff==-Sin[2 dmJ](Sin[\[Beta]]Sin[2\[Xi]]^2/2+Sin[2\[Beta]] Cos[\[Xi]]^4)/4//FullSimplify*)

(* ::Input:: *)
(*XpXcombinedcoeff=-Sin[2 dmJ](Sin[\[Beta]]Sin[2\[Xi]]^2/2+Sin[2\[Beta]] Cos[\[Xi]]^4)/4*)

(* ::Input:: *)
(*signalX= MatrixTrace[(Sxa+Sxb) . (Sxa+Sxb)]XpXcombinedcoeff*)
