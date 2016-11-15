
program hello
    implicit none
	integer :: i,j
	real(8), dimension(1:3) 	:: X,Xc
	real(8), dimension(6,1:3)	:: Xa,Xs
	real(8)					:: rad,phi
	real(8), dimension(size(Xa,1),size(Xs,1)) :: H
	
	!X(1:3) = (/ 1,2,3 /)
	!Xc(1:3) = (/ 4,5,6 /)
	Xs(1,1:3) = (/ 1,2,3 /)
	Xs(2,1:3) = (/ 4,5,6 /)
	Xs(3,1:3) = (/ 2,3,6 /)
	Xs(4,1:3) = (/ 3,8,1 /)
	Xs(5,1:3) = (/ 0,1,1 /)
	Xs(6,1:3) = (/ 0,9,1 /)
	Xa = 0.5*Xs
	rad = 200.0
	
	!print *, "X + Xc = ",X(:,1) + Xc(:,1)
	
	! print *, CalcRBF("C0",rad,X,Xc)
	! print *, CalcRBF("C2",rad,X,Xc)
	! print *, CalcRBF("C4",rad,X,Xc)
	! print *, CalcRBF("C6",rad,X,Xc)
	! print *, CalcRBF("Euclid",rad,X,Xc)
	! print *, CalcRBF("Multiquadric",rad,X,Xc)
	! print *, CalcRBF("InverseMulti",rad,X,Xc)
	! print *, CalcRBF("TPS",rad,X,Xc)
	! print *, CalcRBF("Gaussian",rad,X,Xc)
	
	H = CalcRBFMatrix(Xa,Xs,1,"C2",rad)
	call printM(H)

contains 

! prints a matrix to the stream
subroutine printM(M)
	real(8), intent(in) :: M(:,:)
	integer :: i,j
	do, i=1,size(M(:,1))
		write(*,*) ( M(i,j), j=1,size(M(1,:)) )
	enddo
end subroutine printM

! helper function because it's late and i can't combinatorics
! calculates the number of terms in polynomial of degree degP in 3 variables
! conjecture: pDim = degP*(degP+1)*(degP+2)/6, i.e. tetrahedral numbers		
! function CalcPDim(degP) result(pDim)
	! implicit none
	! integer, intent(in) :: degP
	! integer				:: pDim,i,j,k
	
	! pDim = 0
	! do i = 0,degP
		! do j = 0,degP
			! do k = 0,degP
				! if (i+j+k <= degP) then
					! pDim = pDim + 1					
				! end if
			! end do
		! end do
	! end do
	! if (degP == 0) then
		! pDim = 0
	! end if
	! return
! end function

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


  subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k



! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse

end program hello 
