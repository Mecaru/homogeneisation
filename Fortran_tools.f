        subroutine fast_tensor_rotation(A , R, B)
        real, dimension(81),intent(in) :: A
        real, dimension(3,3),intent(in) :: R
        real ::  product,bijkl,term
        real, dimension(81), intent(out) :: B

        do i = 0,2
        do j = 0,2
        do k = 0,2
        do l = 0,2
            B(27*i + 9*j + 3*k + l + 1) = 0
        end do
        end do
        end do
        end do

        do i=0,2
        do j = 0,i
        do k = 0,2
        do l = 0,k    
         do m = 0,2
         do n = 0,2
         do ll = 0,2
         do kk = 0,2
          product = R(i+1,m+1)* R(j+1,n+1)*R(k+1,ll+1)*R(l+1,kk+1)
          bijkl =  B(27*i + 9*j + 3*k + l + 1)
          term = bijkl + A(27*m + 9*n +3*ll + kk +1)*product
          B(27*i + 9*j + 3*k + l + 1) =  term
         end do
         end do
         end do
         end do
         B(27*j + 9*i + 3*k + l + 1) = B(27*i + 9*j + 3*k + l + 1)
         B(27*i + 9*j + 3*l + k + 1) = B(27*i + 9*j + 3*k + l + 1)
         B(27*j + 9*i + 3*l + k + 1) = B(27*i + 9*j + 3*k + l + 1)
        end do
        end do
        end do
        end do
        end subroutine fast_tensor_rotation
        
        
        

        subroutine inversionT(A,S)
        real, dimension(3,3),intent(in) :: A
        real, dimension(3,3),intent(out) :: S
        real :: det,a1,a2,a3,d1,d2,d3

        
        a1 = A(1,1)*A(2,2)*A(3,3)
        a2 = A(1,2)*A(2,3)*A(3,1)
        a3 = A(1,3)*A(2,1)*A(3,2) 
        d1 =  A(3,2)*A(2,3)*A(1,1)
        d2 =  A(3,3)*A(2,1)*A(1,2)
        d3 =  A(3,1)*A(2,2)*A(1,3)
        det = a1+a2+a3 - d1-d2-d3 

        S(1,1) = (A(2,2)*A(3,3)-A(3,2)*A(2,3))/det
        S(2,1) = -(A(1,2)*A(3,3)-A(3,2)*A(1,3))/det
        S(3,1) = (A(1,2)*A(2,3)-A(2,2)*A(1,3))/det

        S(1,2) = -(A(2,1)*A(3,3)-A(3,1)*A(2,3))/det
        S(2,2) = (A(1,1)*A(3,3)-A(3,1)*A(1,3))/det
        S(3,2) = -(A(1,1)*A(2,3)-A(2,1)*A(1,3))/det

        S(1,3) = (A(2,1)*A(3,2)-A(3,1)*A(2,2))/det
        S(2,3) = -(A(1,1)*A(3,2)-A(3,1)*A(1,2))/det
        S(3,3) = (A(1,1)*A(2,2)-A(2,1)*A(1,2))/det


        end subroutine inversionT