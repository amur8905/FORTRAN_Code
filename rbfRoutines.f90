! constructs the coupling matrix H
! Xa = real(:,3),aerodynamic mesh points
! Xs = real(:,3),structural points
! degP = integer,degree of polynomial term (1 is usually fine)
! rbf = character(:), rbf type (see case statement in function CalcRBF)
! r = real, radius
function CalcRBFMatrix(Xa,Xs,degP,rbf,r) result(H)
	implicit none
	
	real(8),			intent(in)	:: Xa(:,:),Xs(:,:),r
	integer,			intent(in)	:: degP
	character(len=*), 	intent(in) 	:: rbf
	integer 						:: Ns,Na,pDim,i,j,k,ind
	real(8), 			allocatable	:: Ps(:,:),Pa(:,:),Pst(:,:),Aas(:,:),F(:,:),Mp(:,:)
	
	real(8), dimension(size(Xa,1),size(Xs,1)) :: right,H
	real(8), dimension(size(Xs,1),size(Xs,1)) :: M,Mi	
			
	! count the number of possible polynomials terms (tetrahedral numbers)
	pDim = (degP+1)*(degP+2)*(degP+3)/6; 			
	Ns = size(Xs,1)
	Na = size(Xa,1)
	
	allocate(Ps(pDim,Ns))	
	allocate(Pa(Na,pDim))
	allocate(Pst(Ns,pDim))
	allocate(Aas(Na,pDim+Ns))
	allocate(F(pDim+Ns,Ns))
	allocate(Mp(pDim,pDim))
	
	! build the components of the Css matrix
	ind = 1
	do k = 0,degP
		do j = 0,degP
			do i = 0,degP
				if (i+j+k <= degP) then					
					Ps(ind,:) = Xs(:,1)**i*Xs(:,2)**j*Xs(:,3)**k										
					ind = ind + 1					
				end if
			end do
		end do
	end do
	
	
	do i = 1,Ns
		do j = 1,Ns
			M(i,j) = CalcRBF(rbf,r,Xs(i,:),Xs(j,:))
		end do
	end do
	
	call inverse(M,Mi,Ns)
	Pst = transpose(Ps)	
	call inverse(matmul(Ps,matmul(Mi,Pst)),Mp,pDim)
	
	! build the Aas matrix	
	ind = 1
	do k = 0,degP
		do j = 0,degP
			do i = 0,degP
				if (i+j+k <= degP) then
					Pa(:,ind) = Xa(:,1)**i*Xa(:,2)**j*Xa(:,3)**k
					ind = ind + 1
				end if
			end do
		end do 
	end do
	
	do i = 1,Na
		do j = 1,Ns
			right(i,j) = CalcRBF(rbf,r,Xa(i,:),Xs(j,:))
		end do
	end do
	
	Aas(:,1:pDim) = Pa
	Aas(:,pDim+1:pDim+Ns) = right
	
	! calculate coupling matrix H
	if (degP == 0) then
		Aas = right
		H = matmul(Aas,Mi)
	else 
		F(1:pDim,:) = matmul(Mp,matmul(Ps,Mi))
		F(pDim+1:pDim+Ns,:) = Mi - matmul(Mi,matmul(Pst,matmul(Mp,matmul(Ps,Mi))))
		H = matmul(Aas,F)
	end if
	
	deallocate(Ps)
	deallocate(Pa)
	deallocate(Pst)
	deallocate(Aas)
	deallocate(F)
	deallocate(Mp)
	return
end function CalcRBFMatrix

! calculates the value of the given RBF at X with center Xc and radius r
! rbf = character(:), rbf type (see case statement in function CalcRBF)
! r = real, radius, available RBFs are given in the case statement below
! X = real(3), point to evaluate rbf at
! Xc = real(3), center of rbf
function CalcRBF(rbf,r,X,Xc) result(phi)
	implicit none
	real(8), 			parameter	:: pi = 3.141592653589793
	character(len=*),	intent(in) 	:: rbf
	real(8), 			intent(in) 	:: r
	real(8),			intent(in)	:: X(3),Xc(3)	
	real(8)							:: Xnorm,Xr,phi	
	
	Xnorm = sqrt((X(1) - Xc(1))**2 + (X(2) - Xc(2))**2 + (X(3) - Xc(3))**2)
	Xr = Xnorm/r

	select case(rbf)
		case("C0")
			phi = (1-Xr)**2;
			if (Xnorm > r) then
				phi = 0
			end if
		case("C2")
			phi = (1 - Xr)**4*(4*Xr + 1);
			if (Xnorm > r) then
				phi = 0
			end if
		case("C4")
			phi = (1-Xr)**6*(35*Xr**2+18*Xr+3)/3;
			if (Xnorm > r) then
				phi = 0;
			end if
		case("C6")
			phi = (1-Xr)**8*(32*Xr**3+25*Xr**2+8*Xr+1);
			if (Xnorm > r) then
				phi = 0
			end if
		case("Euclid")
			phi = pi*((1/12*Xr**3)-0.5**2*Xr+4/3*0.5**3)/(pi*(-1*0.5**2*0+4/3*0.5**3))			
			if (Xnorm > r) then
				phi = 0
			end if
		case("Multiquadric")
			phi = sqrt(1+Xr**2);
		case("InverseMulti")
			phi = 1/sqrt(1+Xr**2);
		case("TPS")
			phi = Xr**2*log(Xr);
		case("Gaussian")
			phi = exp(-1*Xr**2);
		case default
			print *, "Unknown RBF, using Wendland C2"
			phi = (1 - Xr)**4*(4*Xr + 1);
			if (Xnorm > r) then
				phi = 0
			end if
	end select	
	return
end function CalcRBF 
