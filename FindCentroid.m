function [Cx,Cy]=FindCentroid(DEM)
	% Function to find centroid of drainage basin
	% Input is DEM GRID object 
	% Output is two 1x1 matrices, Cx and Cy with x and y positions of centroid in same coordinates as DEM
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Fall 2015 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


I=isnan(DEM);
I=~I;
Imat=GRIDobj2mat(I);
Imat=double(Imat);

[X,Y]=getcoordinates(DEM);

Ix=bsxfun(@times,Imat,X);
Iy=bsxfun(@times,Imat,Y);

Ixs=nonzeros(Ix);
Iys=nonzeros(Iy);

Cxi=mean(Ixs);
Cyi=mean(Iys);

x_dist=abs(X-Cxi);
y_dist=abs(Y-Cyi);

[~,Ixx]=min(x_dist);
[~,Iyy]=min(y_dist);

Cx=X(Ixx);
Cy=Y(Iyy);

%Function End
end 