function [dmdbasis y0 omega Atilde A] = dmd_comp_Q(X,X1,r,dt)


 [U1 S1 V1]=svd(X,'econ');
 U=U1(:,1:r);
 S=S1(1:r,1:r);
 V=V1(:,1:r);
 A=X1*V*pinv(S)*U';

 Atilde=U'*A*U;
 
  [psiDMD lambdaDMD]=eig(Atilde);
 
 
 dmdbasis=U*psiDMD;
 mu=diag(lambdaDMD);
 omega=log(mu)/dt;
 y0=dmdbasis\X(:,1);

 
end
