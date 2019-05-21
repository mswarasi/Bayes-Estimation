      subroutine multigibbs64(N,Info,maxinfo,GI,a,Se,Sp,p,W,U)

      integer N, maxinfo, GI, a
      integer Info(N,maxinfo), W(N,4)
      real*8  p(4), Se(2), Sp(2), U(N,GI)
      integer g, i, j, Z1, Z2, Y1, Y2, cj, Ysum1, Ysum2, id, stat(4)
      real*8 p1, p2, p3, p4, p11, p01, p10, p00, psum

	p11=p(4)
	p01=p(3)
	p10=p(2)
	p00=p(1)

	DO 101 g=1,GI
	  DO 102 i=1,N
        Info(i,1)=0
	  Info(i,2)=0
	  cj=Info(i,7)
	  Ysum1=0
	  Ysum2=0
	    DO 103 j=1,cj
	    id=Info(i,(7+j))
          Ysum1=Ysum1+Info(id,1)
		Ysum2=Ysum2+Info(id,2)  
103       continue
        Zsum=Info(i,3)+Info(i,4)
	  IF(Ysum1 .gt. 0) THEN
          Ysum1=1
	  ELSE
          Ysum1=0
        END IF
	  IF(Ysum2 .gt. 0) THEN
          Ysum2=1
	  ELSE
          Ysum2=0
        END IF
	  IF(Zsum .gt. 0) THEN
          Zsum=1
	  ELSE
          Zsum=0
        END IF
        Z1=Info(i,3)
	  Z2=Info(i,4)
	  Y1=Info(i,5)
	  Y2=Info(i,6)
        p1=p00*(Se(1)**Z1*(1-Se(1))**(1-Z1))**(Ysum1)
	  p1=p1 *(Sp(1)**(1-Z1)*(1-Sp(1))**(Z1))**(1-Ysum1)
	  p1=p1 *(Se(2)**Z2*(1-Se(2))**(1-Z2))**(Ysum2)
	  p1=p1 *(Sp(2)**(1-Z2)*(1-Sp(2))**(Z2))**(1-Ysum2)
	  p1=p1 *(Sp(1)**(1-Y1)*(1-Sp(1))**(Y1))**Zsum
	  p1=p1 *(Sp(2)**(1-Y2)*(1-Sp(2))**(Y2))**Zsum
      
	  p2=p10*Se(1)**Z1*(1-Se(1))**(1-Z1)
	  p2=p2*(Se(2)**Z2*(1-Se(2))**(1-Z2))**(Ysum2)
	  p2=p2*(Sp(2)**(1-Z2)*(1-Sp(2))**Z2)**(1-Ysum2)
	  p2=p2*(Se(1)**Y1*(1-Se(1))**(1-Y1))**Zsum
	  p2=p2*(Sp(2)**(1-Y2)*(1-Sp(2))**(Y2))**Zsum
       
	  p3=p01*Se(2)**Z2*(1-Se(2))**(1-Z2)
	  p3=p3*(Se(1)**Z1*(1-Se(1))**(1-Z1))**(Ysum1)
	  p3=p3*(Sp(1)**(1-Z1)*(1-Sp(1))**(Z1))**(1-Ysum1)
	  p3=p3*(Se(2)**Y2*(1-Se(2))**(1-Y2))**Zsum
	  p3=p3*(Sp(1)**(1-Y1)*(1-Sp(1))**(Y1))**Zsum
       
	  p4=p11*Se(1)**Z1*(1-Se(1))**(1-Z1)
	  p4=p4 *Se(2)**Z2*(1-Se(2))**(1-Z2)
	  p4=p4*(Se(1)**Y1*(1-Se(1))**(1-Y1))**Zsum
	  p4=p4*(Se(2)**Y2*(1-Se(2))**(1-Y2))**Zsum
        
	  psum=p1+p2+p3+p4
	  p1=p1/psum
	  p2=p2/psum+p1
	  p3=p3/psum+p2

	  IF(U(i,g) .gt. p3) THEN
        stat(1)=0
	  stat(2)=0
	  stat(3)=0
	  stat(4)=1
	  Info(i,1)=1
	  Info(i,2)=1
        ELSE
	    IF(U(i,g) .gt. p2) THEN
            stat(1)=0
	      stat(2)=0
	      stat(3)=1
	      stat(4)=0
	      Info(i,1)=0
	      Info(i,2)=1
 		ELSE
		  IF(U(i,g) .gt. p1) THEN
              stat(1)=0
	        stat(2)=1
	        stat(3)=0
	        stat(4)=0
	        Info(i,1)=1
	        Info(i,2)=0
 		  ELSE
              stat(1)=1
	        stat(2)=0
	        stat(3)=0
	        stat(4)=0
	        Info(i,1)=0
	        Info(i,2)=0
 		  END IF  
          END IF
        END IF
	  IF(g .gt. a) THEN
        W(i,1)=W(i,1)+stat(1)
	  W(i,2)=W(i,2)+stat(2)
	  W(i,3)=W(i,3)+stat(3)
	  W(i,4)=W(i,4)+stat(4)
	  END IF 
102     continue
101   continue 


      return
	end


