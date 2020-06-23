ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DOUBLE PRECISION FUNCTION DGAUSS(F,A,B,EPS)                   

C     ADAPTIVE DOUBLE PRECISION GAUSSIAN QUADRATURE.                   
C                                                                      
C     DGAUSS IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF  
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER  
C     EPS.                                                              
C                                                                      

      DOUBLE PRECISION F,A,B,EPS                                      
      DOUBLE PRECISION W(12),X(12),AA,BB,C1,C2,U,S8,S16,CONST         
      LOGICAL MFLAG,RFLAG                                             
      EXTERNAL F                                                     

      DATA W / 0.10122 85362 90376 259D0,                             
     1         0.22238 10344 53374 471D0,                              
     2         0.31370 66458 77887 287D0,                           
     3         0.36268 37833 78361 983D0,                           
     4         0.27152 45941 17540 949D-1,                        
     5         0.62253 52393 86478 929D-1,                         
     6         0.95158 51168 24927 848D-1,                         
     7         0.12462 89712 55533 872D0,                     
     8         0.14959 59888 16576 732D0,                              
     9         0.16915 65193 95002 538D0,                              
     A         0.18260 34150 44923 589D0,                           
     B         0.18945 06104 55068 496D0/                          
                                                                       
      DATA X / 0.96028 98564 97536 232D0,                             
     1         0.79666 64774 13626 740D0,                             
     2         0.52553 24099 16328 986D0,                       
     3         0.18343 46424 95649 805D0,                        
     4         0.98940 09349 91649 933D0,                          
     5         0.94457 50230 73232 576D0,                         
     6         0.86563 12023 87831 744D0,                        
     7         0.75540 44083 55003 034D0,                          
     8         0.61787 62444 02643 748D0,   
     9         0.45801 67776 57227 386D0,   
     A         0.28160 35507 79258 913D0,  
     B         0.95012 50983 76374 402D-1/ 
C                                      
C     ******************************************************************
C                                                               
C  START.                                            
      DGAUSS=0.0D0                                    
      IF(B.EQ.A) RETURN 
      CONST=0.005D0/(B-A) 
      BB=A  
C          
C  COMPUTATIONAL LOOP.   
    1 AA=BB  
      BB=B   
    2    C1=0.5D0*(BB+AA)  
         C2=0.5D0*(BB-AA)  
         S8=0.0D0          
         DO 3 I=1,4   
            U=C2*X(I) 
            S8=S8+W(I)*(F(C1+U)+F(C1-U)) 
    3    CONTINUE 
         S8=C2*S8 
         S16=0.0D0
         DO 4 I=5,12  
            U=C2*X(I)  
            S16=S16+W(I)*(F(C1+U)+F(C1-U))  
    4    CONTINUE  
         S16=C2*S16  
         IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5 
         BB=C1   
         IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2  
      DGAUSS=0.0D0 
      CALL KERMTR('D103.1',LGFILE,MFLAG,RFLAG) 
      IF(MFLAG) THEN  
         IF(LGFILE.EQ.0) THEN   
c            WRITE(*,6)   
         ELSE  
c            WRITE(LGFILE,6) 
         ENDIF             
      ENDIF                               
      IF(.NOT. RFLAG) CALL ABEND  
      RETURN             
    5 DGAUSS=DGAUSS+S16   
      IF(BB.NE.B) GO TO 1  
      RETURN          
C                     
    6 FORMAT( 4X, 'FUNCTION DGAUSS ... TOO HIGH ACCURACY REQUIRED') 
      END      




      DOUBLE PRECISION FUNCTION DGAUSS3a(F,A,B,EPS)                   

C     ADAPTIVE DOUBLE PRECISION GAUSSIAN QUADRATURE.                   
C                                                                      
C     DGAUSS IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF  
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER  
C     EPS.                                                              
C                                                                      

      DOUBLE PRECISION F,A,B,EPS                                      
      DOUBLE PRECISION W(12),X(12),AA,BB,C1,C2,U,S8,S16,CONST         
      LOGICAL MFLAG,RFLAG                                             
      EXTERNAL F                                                     

      DATA W / 0.10122 85362 90376 259D0,                             
     1         0.22238 10344 53374 471D0,                              
     2         0.31370 66458 77887 287D0,                           
     3         0.36268 37833 78361 983D0,                           
     4         0.27152 45941 17540 949D-1,                        
     5         0.62253 52393 86478 929D-1,                         
     6         0.95158 51168 24927 848D-1,                         
     7         0.12462 89712 55533 872D0,                     
     8         0.14959 59888 16576 732D0,                              
     9         0.16915 65193 95002 538D0,                              
     A         0.18260 34150 44923 589D0,                           
     B         0.18945 06104 55068 496D0/                          
                                                                       
      DATA X / 0.96028 98564 97536 232D0,                             
     1         0.79666 64774 13626 740D0,                             
     2         0.52553 24099 16328 986D0,                       
     3         0.18343 46424 95649 805D0,                        
     4         0.98940 09349 91649 933D0,                          
     5         0.94457 50230 73232 576D0,                         
     6         0.86563 12023 87831 744D0,                        
     7         0.75540 44083 55003 034D0,                          
     8         0.61787 62444 02643 748D0,   
     9         0.45801 67776 57227 386D0,   
     A         0.28160 35507 79258 913D0,  
     B         0.95012 50983 76374 402D-1/ 
C                                      
C     ******************************************************************
C                                                               
C  START.                                            
      DGAUSS3a=0.0D0                                    
      IF(B.EQ.A) RETURN 
      CONST=0.005D0/(B-A) 
      BB=A  
C          
C  COMPUTATIONAL LOOP.   
    1 AA=BB  
      BB=B   
    2    C1=0.5D0*(BB+AA)  
         C2=0.5D0*(BB-AA)  
         S8=0.0D0          
         DO 3 I=1,4   
            U=C2*X(I) 
            S8=S8+W(I)*(F(C1+U)+F(C1-U)) 
    3    CONTINUE 
         S8=C2*S8 
         S16=0.0D0
         DO 4 I=5,12  
            U=C2*X(I)  
            S16=S16+W(I)*(F(C1+U)+F(C1-U))  
    4    CONTINUE  
         S16=C2*S16  
         IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5 
         BB=C1   
         IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2  
      DGAUSS3a=0.0D0 
      CALL KERMTR('D103.1',LGFILE,MFLAG,RFLAG) 
      IF(MFLAG) THEN  
         IF(LGFILE.EQ.0) THEN   
c            WRITE(*,6)   
         ELSE  
c            WRITE(LGFILE,6) 
         ENDIF             
      ENDIF                               
      IF(.NOT. RFLAG) CALL ABEND  
      RETURN             
    5 DGAUSS3a=DGAUSS3a+S16   
      IF(BB.NE.B) GO TO 1  
      RETURN          
C                     
    6 FORMAT( 4X, 'FUNCTION DGAUSS3a ... TOO HIGH ACCURACY REQUIRED') 
      END      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      DOUBLE PRECISION FUNCTION DGAUSS1(F,A,B,EPS)                   

C     ADAPTIVE DOUBLE PRECISION GAUSSIAN QUADRATURE.                   
C                                                                      
C     DGAUSS1 IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF  
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER  
C     EPS.                                                              
C                                                                      

      DOUBLE PRECISION F,A,B,EPS                                      
      DOUBLE PRECISION W(12),X(12),AA,BB,C1,C2,U,S8,S16,CONST         
      LOGICAL MFLAG,RFLAG                                             
      EXTERNAL F                                                     

      DATA W / 0.10122 85362 90376 259D0,                             
     1         0.22238 10344 53374 471D0,                              
     2         0.31370 66458 77887 287D0,                           
     3         0.36268 37833 78361 983D0,                           
     4         0.27152 45941 17540 949D-1,                        
     5         0.62253 52393 86478 929D-1,                         
     6         0.95158 51168 24927 848D-1,                         
     7         0.12462 89712 55533 872D0,                     
     8         0.14959 59888 16576 732D0,                              
     9         0.16915 65193 95002 538D0,                              
     A         0.18260 34150 44923 589D0,                           
     B         0.18945 06104 55068 496D0/                          
                                                                       
      DATA X / 0.96028 98564 97536 232D0,                             
     1         0.79666 64774 13626 740D0,                             
     2         0.52553 24099 16328 986D0,                       
     3         0.18343 46424 95649 805D0,                        
     4         0.98940 09349 91649 933D0,                          
     5         0.94457 50230 73232 576D0,                         
     6         0.86563 12023 87831 744D0,                        
     7         0.75540 44083 55003 034D0,                          
     8         0.61787 62444 02643 748D0,   
     9         0.45801 67776 57227 386D0,   
     A         0.28160 35507 79258 913D0,  
     B         0.95012 50983 76374 402D-1/ 
C                                      
C     ******************************************************************
C                                                               
C  START.                                            
      DGAUSS1=0.0D0                                    
      IF(B.EQ.A) RETURN 
      CONST=0.005D0/(B-A) 
      BB=A  
C          
C  COMPUTATIONAL LOOP.   
    1 AA=BB  
      BB=B   
    2    C1=0.5D0*(BB+AA)  
         C2=0.5D0*(BB-AA)  
         S8=0.0D0          
         DO 3 I=1,4   
            U=C2*X(I) 
            S8=S8+W(I)*(F(C1+U)+F(C1-U)) 
    3    CONTINUE 
         S8=C2*S8 
         S16=0.0D0
         DO 4 I=5,12  
            U=C2*X(I)  
            S16=S16+W(I)*(F(C1+U)+F(C1-U))  
    4    CONTINUE  
         S16=C2*S16  
         IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5 
         BB=C1   
         IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2  
      DGAUSS1=0.0D0 
      CALL KERMTR('D103.1',LGFILE,MFLAG,RFLAG) 
      IF(MFLAG) THEN  
         IF(LGFILE.EQ.0) THEN   
c            WRITE(*,6)   
         ELSE  
c            WRITE(LGFILE,6) 
         ENDIF             
      ENDIF                               
      IF(.NOT. RFLAG) CALL ABEND  
      RETURN             
    5 DGAUSS1=DGAUSS1+S16   
      IF(BB.NE.B) GO TO 1  
      RETURN          
C                     
    6 FORMAT( 4X, 'FUNCTION DGAUSS1 ... TOO HIGH ACCURACY REQUIRED') 
      END      
***********************************************************************
*
      DOUBLE PRECISION FUNCTION DGAUSS2(F,A,B,EPS)                   

C     ADAPTIVE DOUBLE PRECISION GAUSSIAN QUADRATURE.                   
C                                                                      
C     DGAUSS2 IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF  
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER  
C     EPS.                                                              
C                                                                      

      DOUBLE PRECISION F,A,B,EPS                                      
      DOUBLE PRECISION W(12),X(12),AA,BB,C1,C2,U,S8,S16,CONST         
      LOGICAL MFLAG,RFLAG                                             
      EXTERNAL F                                                     

      DATA W / 0.10122 85362 90376 259D0,                             
     1         0.22238 10344 53374 471D0,                              
     2         0.31370 66458 77887 287D0,                           
     3         0.36268 37833 78361 983D0,                           
     4         0.27152 45941 17540 949D-1,                        
     5         0.62253 52393 86478 929D-1,                         
     6         0.95158 51168 24927 848D-1,                         
     7         0.12462 89712 55533 872D0,                     
     8         0.14959 59888 16576 732D0,                              
     9         0.16915 65193 95002 538D0,                              
     A         0.18260 34150 44923 589D0,                           
     B         0.18945 06104 55068 496D0/                          
                                                                       
      DATA X / 0.96028 98564 97536 232D0,                             
     1         0.79666 64774 13626 740D0,                             
     2         0.52553 24099 16328 986D0,                       
     3         0.18343 46424 95649 805D0,                        
     4         0.98940 09349 91649 933D0,                          
     5         0.94457 50230 73232 576D0,                         
     6         0.86563 12023 87831 744D0,                        
     7         0.75540 44083 55003 034D0,                          
     8         0.61787 62444 02643 748D0,   
     9         0.45801 67776 57227 386D0,   
     A         0.28160 35507 79258 913D0,  
     B         0.95012 50983 76374 402D-1/ 
C                                      
C     ******************************************************************
C                                                               
C  START.                                            
      DGAUSS2=0.0D0                                    
      IF(B.EQ.A) RETURN 
      CONST=0.005D0/(B-A) 
      BB=A  
C          
C  COMPUTATIONAL LOOP.   
    1 AA=BB  
      BB=B   
    2    C1=0.5D0*(BB+AA)  
         C2=0.5D0*(BB-AA)  
         S8=0.0D0          
         DO 3 I=1,4   
            U=C2*X(I) 
            S8=S8+W(I)*(F(C1+U)+F(C1-U)) 
    3    CONTINUE 
         S8=C2*S8 
         S16=0.0D0
         DO 4 I=5,12  
            U=C2*X(I)  
            S16=S16+W(I)*(F(C1+U)+F(C1-U))  
    4    CONTINUE  
         S16=C2*S16  
         IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5 
         BB=C1   
         IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2  
      DGAUSS2=0.0D0 
      CALL KERMTR('D103.1',LGFILE,MFLAG,RFLAG) 
      IF(MFLAG) THEN  
         IF(LGFILE.EQ.0) THEN   
            WRITE(*,6)   
         ELSE  
            WRITE(LGFILE,6) 
         ENDIF             
      ENDIF                               
      IF(.NOT. RFLAG) CALL ABEND  
      RETURN             
    5 DGAUSS2=DGAUSS2+S16   
      IF(BB.NE.B) GO TO 1  
      RETURN          
C                     
    6 FORMAT( 4X, 'FUNCTION DGAUSS2 ... TOO HIGH ACCURACY REQUIRED') 
      END      


*
*
      DOUBLE PRECISION FUNCTION DGAUSS3(F,A,B,EPS)                   

C     ADAPTIVE DOUBLE PRECISION GAUSSIAN QUADRATURE.                   
C                                                                      
C     DGAUSS3 IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF  
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER  
C     EPS.                                                              
C                                                                      

      DOUBLE PRECISION F,A,B,EPS                                      
      DOUBLE PRECISION W(12),X(12),AA,BB,C1,C2,U,S8,S16,CONST         
      LOGICAL MFLAG,RFLAG                                             
      EXTERNAL F                                                     

      DATA W / 0.10122 85362 90376 259D0,                             
     1         0.22238 10344 53374 471D0,                              
     2         0.31370 66458 77887 287D0,                           
     3         0.36268 37833 78361 983D0,                           
     4         0.27152 45941 17540 949D-1,                        
     5         0.62253 52393 86478 929D-1,                         
     6         0.95158 51168 24927 848D-1,                         
     7         0.12462 89712 55533 872D0,                     
     8         0.14959 59888 16576 732D0,                              
     9         0.16915 65193 95002 538D0,                              
     A         0.18260 34150 44923 589D0,                           
     B         0.18945 06104 55068 496D0/                          
                                                                       
      DATA X / 0.96028 98564 97536 232D0,                             
     1         0.79666 64774 13626 740D0,                             
     2         0.52553 24099 16328 986D0,                       
     3         0.18343 46424 95649 805D0,                        
     4         0.98940 09349 91649 933D0,                          
     5         0.94457 50230 73232 576D0,                         
     6         0.86563 12023 87831 744D0,                        
     7         0.75540 44083 55003 034D0,                          
     8         0.61787 62444 02643 748D0,   
     9         0.45801 67776 57227 386D0,   
     A         0.28160 35507 79258 913D0,  
     B         0.95012 50983 76374 402D-1/ 
C                                      
C     ******************************************************************
C                                                               
C  START.                                            
      DGAUSS3=0.0D0                                    
      IF(B.EQ.A) RETURN 
      CONST=0.005D0/(B-A) 
      BB=A  
C          
C  COMPUTATIONAL LOOP.   
    1 AA=BB  
      BB=B   
    2    C1=0.5D0*(BB+AA)  
         C2=0.5D0*(BB-AA)  
         S8=0.0D0          
         DO 3 I=1,4   
            U=C2*X(I) 
            S8=S8+W(I)*(F(C1+U)+F(C1-U)) 
    3    CONTINUE 
         S8=C2*S8 
         S16=0.0D0
         DO 4 I=5,12  
            U=C2*X(I)  
            S16=S16+W(I)*(F(C1+U)+F(C1-U))  
    4    CONTINUE  
         S16=C2*S16  
         IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5 
         BB=C1   
         IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2  
      DGAUSS3=0.0D0 
      CALL KERMTR('D103.1',LGFILE,MFLAG,RFLAG) 
      IF(MFLAG) THEN  
         IF(LGFILE.EQ.0) THEN   
            WRITE(*,6)   
         ELSE  
            WRITE(LGFILE,6) 
         ENDIF             
      ENDIF                               
      IF(.NOT. RFLAG) CALL ABEND  
      RETURN             
    5 DGAUSS3=DGAUSS3+S16   
      IF(BB.NE.B) GO TO 1  
      RETURN          
C                     
    6 FORMAT( 4X, 'FUNCTION DGAUSS3 ... TOO HIGH ACCURACY REQUIRED') 
      END      


***********************************************************************
*
      DOUBLE PRECISION FUNCTION DGAUSS4(F,A,B,EPS)                   

C     ADAPTIVE DOUBLE PRECISION GAUSSIAN QUADRATURE.                   
C                                                                      
C     DGAUSS3 IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF  
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER  
C     EPS.                                                              
C                                                                      

      DOUBLE PRECISION F,A,B,EPS                                      
      DOUBLE PRECISION W(12),X(12),AA,BB,C1,C2,U,S8,S16,CONST         
      LOGICAL MFLAG,RFLAG                                             
      EXTERNAL F                                                     

      DATA W / 0.10122 85362 90376 259D0,                             
     1         0.22238 10344 53374 471D0,                              
     2         0.31370 66458 77887 287D0,                           
     3         0.36268 37833 78361 983D0,                           
     4         0.27152 45941 17540 949D-1,                        
     5         0.62253 52393 86478 929D-1,                         
     6         0.95158 51168 24927 848D-1,                         
     7         0.12462 89712 55533 872D0,                     
     8         0.14959 59888 16576 732D0,                              
     9         0.16915 65193 95002 538D0,                              
     A         0.18260 34150 44923 589D0,                           
     B         0.18945 06104 55068 496D0/                          
                                                                       
      DATA X / 0.96028 98564 97536 232D0,                             
     1         0.79666 64774 13626 740D0,                             
     2         0.52553 24099 16328 986D0,                       
     3         0.18343 46424 95649 805D0,                        
     4         0.98940 09349 91649 933D0,                          
     5         0.94457 50230 73232 576D0,                         
     6         0.86563 12023 87831 744D0,                        
     7         0.75540 44083 55003 034D0,                          
     8         0.61787 62444 02643 748D0,   
     9         0.45801 67776 57227 386D0,   
     A         0.28160 35507 79258 913D0,  
     B         0.95012 50983 76374 402D-1/ 
C                                      
C     ******************************************************************
C                                                               
C  START.                                            
      DGAUSS4=0.0D0                                    
      IF(B.EQ.A) RETURN 
      CONST=0.005D0/(B-A) 
      BB=A  
C          
C  COMPUTATIONAL LOOP.   
    1 AA=BB  
      BB=B   
    2    C1=0.5D0*(BB+AA)  
         C2=0.5D0*(BB-AA)  
         S8=0.0D0          
         DO 3 I=1,4   
            U=C2*X(I) 
            S8=S8+W(I)*(F(C1+U)+F(C1-U)) 
    3    CONTINUE 
         S8=C2*S8 
         S16=0.0D0
         DO 4 I=5,12  
            U=C2*X(I)  
            S16=S16+W(I)*(F(C1+U)+F(C1-U))  
    4    CONTINUE  
         S16=C2*S16  
         IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5 
         BB=C1   
         IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2  
      DGAUSS4=0.0D0 
      CALL KERMTR('D103.1',LGFILE,MFLAG,RFLAG) 
      IF(MFLAG) THEN  
         IF(LGFILE.EQ.0) THEN   
            WRITE(*,6)   
         ELSE  
            WRITE(LGFILE,6) 
         ENDIF             
      ENDIF                               
      IF(.NOT. RFLAG) CALL ABEND  
      RETURN             
    5 DGAUSS4=DGAUSS4+S16   
      IF(BB.NE.B) GO TO 1  
      RETURN          
C                     
    6 FORMAT( 4X, 'FUNCTION DGAUSS4 ... TOO HIGH ACCURACY REQUIRED') 
      END      


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DOUBLE PRECISION FUNCTION DGAUSS5(F,A,B,EPS)                   

C     ADAPTIVE DOUBLE PRECISION GAUSSIAN QUADRATURE.                   
C                                                                      
C     DGAUSS IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF  
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER  
C     EPS.                                                              
C                                                                      

      DOUBLE PRECISION F,A,B,EPS                                      
      DOUBLE PRECISION W(12),X(12),AA,BB,C1,C2,U,S8,S16,CONST         
      LOGICAL MFLAG,RFLAG                                             
      EXTERNAL F                                                     

      DATA W / 0.10122 85362 90376 259D0,                             
     1         0.22238 10344 53374 471D0,                              
     2         0.31370 66458 77887 287D0,                           
     3         0.36268 37833 78361 983D0,                           
     4         0.27152 45941 17540 949D-1,                        
     5         0.62253 52393 86478 929D-1,                         
     6         0.95158 51168 24927 848D-1,                         
     7         0.12462 89712 55533 872D0,                     
     8         0.14959 59888 16576 732D0,                              
     9         0.16915 65193 95002 538D0,                              
     A         0.18260 34150 44923 589D0,                           
     B         0.18945 06104 55068 496D0/                          
                                                                       
      DATA X / 0.96028 98564 97536 232D0,                             
     1         0.79666 64774 13626 740D0,                             
     2         0.52553 24099 16328 986D0,                       
     3         0.18343 46424 95649 805D0,                        
     4         0.98940 09349 91649 933D0,                          
     5         0.94457 50230 73232 576D0,                         
     6         0.86563 12023 87831 744D0,                        
     7         0.75540 44083 55003 034D0,                          
     8         0.61787 62444 02643 748D0,   
     9         0.45801 67776 57227 386D0,   
     A         0.28160 35507 79258 913D0,  
     B         0.95012 50983 76374 402D-1/ 
C                                      
C     ******************************************************************
C                                                               
C  START.                                            
      DGAUSS5=0.0D0                                    
      IF(B.EQ.A) RETURN 
      CONST=0.005D0/(B-A) 
      BB=A  
C          
C  COMPUTATIONAL LOOP.   
    1 AA=BB  
      BB=B   
    2    C1=0.5D0*(BB+AA)  
         C2=0.5D0*(BB-AA)  
         S8=0.0D0          
         DO 3 I=1,4   
            U=C2*X(I) 
            S8=S8+W(I)*(F(C1+U)+F(C1-U)) 
    3    CONTINUE 
         S8=C2*S8 
         S16=0.0D0
         DO 4 I=5,12  
            U=C2*X(I)  
            S16=S16+W(I)*(F(C1+U)+F(C1-U))  
    4    CONTINUE  
         S16=C2*S16  
         IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5 
         BB=C1   
         IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2  
      DGAUSS5=0.0D0 
      CALL KERMTR('D103.1',LGFILE,MFLAG,RFLAG) 
      IF(MFLAG) THEN  
         IF(LGFILE.EQ.0) THEN   
            WRITE(*,6)   
         ELSE  
            WRITE(LGFILE,6) 
         ENDIF             
      ENDIF                               
      IF(.NOT. RFLAG) CALL ABEND  
      RETURN             
    5 DGAUSS5=DGAUSS5+S16   
      IF(BB.NE.B) GO TO 1  
      RETURN          
C                     
    6 FORMAT( 4X, 'FUNCTION DGAUSS5 ... TOO HIGH ACCURACY REQUIRED') 
      END      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DOUBLE PRECISION FUNCTION DGAUSS6(F,A,B,EPS)                   

C     ADAPTIVE DOUBLE PRECISION GAUSSIAN QUADRATURE.                   
C                                                                      
C     DGAUSS IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF  
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER  
C     EPS.                                                              
C                                                                      

      DOUBLE PRECISION F,A,B,EPS                                      
      DOUBLE PRECISION W(12),X(12),AA,BB,C1,C2,U,S8,S16,CONST         
      LOGICAL MFLAG,RFLAG                                             
      EXTERNAL F                                                     

      DATA W / 0.10122 85362 90376 259D0,                             
     1         0.22238 10344 53374 471D0,                              
     2         0.31370 66458 77887 287D0,                           
     3         0.36268 37833 78361 983D0,                           
     4         0.27152 45941 17540 949D-1,                        
     5         0.62253 52393 86478 929D-1,                         
     6         0.95158 51168 24927 848D-1,                         
     7         0.12462 89712 55533 872D0,                     
     8         0.14959 59888 16576 732D0,                              
     9         0.16915 65193 95002 538D0,                              
     A         0.18260 34150 44923 589D0,                           
     B         0.18945 06104 55068 496D0/                          
                                                                       
      DATA X / 0.96028 98564 97536 232D0,                             
     1         0.79666 64774 13626 740D0,                             
     2         0.52553 24099 16328 986D0,                       
     3         0.18343 46424 95649 805D0,                        
     4         0.98940 09349 91649 933D0,                          
     5         0.94457 50230 73232 576D0,                         
     6         0.86563 12023 87831 744D0,                        
     7         0.75540 44083 55003 034D0,                          
     8         0.61787 62444 02643 748D0,   
     9         0.45801 67776 57227 386D0,   
     A         0.28160 35507 79258 913D0,  
     B         0.95012 50983 76374 402D-1/ 
C                                      
C     ******************************************************************
C                                                               
C  START.                                            
      DGAUSS6=0.0D0                                    
      IF(B.EQ.A) RETURN 
      CONST=0.005D0/(B-A) 
      BB=A  
C          
C  COMPUTATIONAL LOOP.   
    1 AA=BB  
      BB=B   
    2    C1=0.5D0*(BB+AA)  
         C2=0.5D0*(BB-AA)  
         S8=0.0D0          
         DO 3 I=1,4   
            U=C2*X(I) 
            S8=S8+W(I)*(F(C1+U)+F(C1-U)) 
    3    CONTINUE 
         S8=C2*S8 
         S16=0.0D0
         DO 4 I=5,12  
            U=C2*X(I)  
            S16=S16+W(I)*(F(C1+U)+F(C1-U))  
    4    CONTINUE  
         S16=C2*S16  
         IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5 
         BB=C1   
         IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2  
      DGAUSS6=0.0D0 
      CALL KERMTR('D103.1',LGFILE,MFLAG,RFLAG) 
      IF(MFLAG) THEN  
         IF(LGFILE.EQ.0) THEN   
            WRITE(*,6)   
         ELSE  
            WRITE(LGFILE,6) 
         ENDIF             
      ENDIF                               
      IF(.NOT. RFLAG) CALL ABEND  
      RETURN             
    5 DGAUSS6=DGAUSS6+S16   
      IF(BB.NE.B) GO TO 1  
      RETURN          
C                     
    6 FORMAT( 4X, 'FUNCTION DGAUSS6 ... TOO HIGH ACCURACY REQUIRED') 
      END      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DOUBLE PRECISION FUNCTION DGAUSS7(F,A,B,EPS)                   

C     ADAPTIVE DOUBLE PRECISION GAUSSIAN QUADRATURE.                   
C                                                                      
C     DGAUSS IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF  
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER  
C     EPS.                                                              
C                                                                      

      DOUBLE PRECISION F,A,B,EPS                                      
      DOUBLE PRECISION W(12),X(12),AA,BB,C1,C2,U,S8,S16,CONST         
      LOGICAL MFLAG,RFLAG                                             
      EXTERNAL F                                                     

      DATA W / 0.10122 85362 90376 259D0,                             
     1         0.22238 10344 53374 471D0,                              
     2         0.31370 66458 77887 287D0,                           
     3         0.36268 37833 78361 983D0,                           
     4         0.27152 45941 17540 949D-1,                        
     5         0.62253 52393 86478 929D-1,                         
     6         0.95158 51168 24927 848D-1,                         
     7         0.12462 89712 55533 872D0,                     
     8         0.14959 59888 16576 732D0,                              
     9         0.16915 65193 95002 538D0,                              
     A         0.18260 34150 44923 589D0,                           
     B         0.18945 06104 55068 496D0/                          
                                                                       
      DATA X / 0.96028 98564 97536 232D0,                             
     1         0.79666 64774 13626 740D0,                             
     2         0.52553 24099 16328 986D0,                       
     3         0.18343 46424 95649 805D0,                        
     4         0.98940 09349 91649 933D0,                          
     5         0.94457 50230 73232 576D0,                         
     6         0.86563 12023 87831 744D0,                        
     7         0.75540 44083 55003 034D0,                          
     8         0.61787 62444 02643 748D0,   
     9         0.45801 67776 57227 386D0,   
     A         0.28160 35507 79258 913D0,  
     B         0.95012 50983 76374 402D-1/ 
C                                      
C     ******************************************************************
C                                                               
C  START.                                            
      DGAUSS7=0.0D0                                    
      IF(B.EQ.A) RETURN 
      CONST=0.005D0/(B-A) 
      BB=A  
C          
C  COMPUTATIONAL LOOP.   
    1 AA=BB  
      BB=B   
    2    C1=0.5D0*(BB+AA)  
         C2=0.5D0*(BB-AA)  
         S8=0.0D0          
         DO 3 I=1,4   
            U=C2*X(I) 
            S8=S8+W(I)*(F(C1+U)+F(C1-U)) 
    3    CONTINUE 
         S8=C2*S8 
         S16=0.0D0
         DO 4 I=5,12  
            U=C2*X(I)  
            S16=S16+W(I)*(F(C1+U)+F(C1-U))  
    4    CONTINUE  
         S16=C2*S16  
         IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5 
         BB=C1   
         IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2  
      DGAUSS7=0.0D0 
      CALL KERMTR('D103.1',LGFILE,MFLAG,RFLAG) 
      IF(MFLAG) THEN  
         IF(LGFILE.EQ.0) THEN   
            WRITE(*,6)   
         ELSE  
            WRITE(LGFILE,6) 
         ENDIF             
      ENDIF                               
      IF(.NOT. RFLAG) CALL ABEND  
      RETURN             
    5 DGAUSS7=DGAUSS7+S16   
      IF(BB.NE.B) GO TO 1  
      RETURN          
C                     
    6 FORMAT( 4X, 'FUNCTION DGAUSS7 ... TOO HIGH ACCURACY REQUIRED') 
      END      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DOUBLE PRECISION FUNCTION DGAUSS8(F,A,B,EPS)                   

C     ADAPTIVE DOUBLE PRECISION GAUSSIAN QUADRATURE.                   
C                                                                      
C     DGAUSS IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF  
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER  
C     EPS.                                                              
C                                                                      

      DOUBLE PRECISION F,A,B,EPS                                      
      DOUBLE PRECISION W(12),X(12),AA,BB,C1,C2,U,S8,S16,CONST         
      LOGICAL MFLAG,RFLAG                                             
      EXTERNAL F                                                     

      DATA W / 0.10122 85362 90376 259D0,                             
     1         0.22238 10344 53374 471D0,                              
     2         0.31370 66458 77887 287D0,                           
     3         0.36268 37833 78361 983D0,                           
     4         0.27152 45941 17540 949D-1,                        
     5         0.62253 52393 86478 929D-1,                         
     6         0.95158 51168 24927 848D-1,                         
     7         0.12462 89712 55533 872D0,                     
     8         0.14959 59888 16576 732D0,                              
     9         0.16915 65193 95002 538D0,                              
     A         0.18260 34150 44923 589D0,                           
     B         0.18945 06104 55068 496D0/                          
                                                                       
      DATA X / 0.96028 98564 97536 232D0,                             
     1         0.79666 64774 13626 740D0,                             
     2         0.52553 24099 16328 986D0,                       
     3         0.18343 46424 95649 805D0,                        
     4         0.98940 09349 91649 933D0,                          
     5         0.94457 50230 73232 576D0,                         
     6         0.86563 12023 87831 744D0,                        
     7         0.75540 44083 55003 034D0,                          
     8         0.61787 62444 02643 748D0,   
     9         0.45801 67776 57227 386D0,   
     A         0.28160 35507 79258 913D0,  
     B         0.95012 50983 76374 402D-1/ 
C                                      
C     ******************************************************************
C                                                               
C  START.                                            
      DGAUSS8=0.0D0                                    
      IF(B.EQ.A) RETURN 
      CONST=0.005D0/(B-A) 
      BB=A  
C          
C  COMPUTATIONAL LOOP.   
    1 AA=BB  
      BB=B   
    2    C1=0.5D0*(BB+AA)  
         C2=0.5D0*(BB-AA)  
         S8=0.0D0          
         DO 3 I=1,4   
            U=C2*X(I) 
            S8=S8+W(I)*(F(C1+U)+F(C1-U)) 
    3    CONTINUE 
         S8=C2*S8 
         S16=0.0D0
         DO 4 I=5,12  
            U=C2*X(I)  
            S16=S16+W(I)*(F(C1+U)+F(C1-U))  
    4    CONTINUE  
         S16=C2*S16  
         IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5 
         BB=C1   
         IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2  
      DGAUSS8=0.0D0 
      CALL KERMTR('D103.1',LGFILE,MFLAG,RFLAG) 
      IF(MFLAG) THEN  
         IF(LGFILE.EQ.0) THEN   
            WRITE(*,6)   
         ELSE  
            WRITE(LGFILE,6) 
         ENDIF             
      ENDIF                               
      IF(.NOT. RFLAG) CALL ABEND  
      RETURN             
    5 DGAUSS8=DGAUSS8+S16   
      IF(BB.NE.B) GO TO 1  
      RETURN          
C                     
    6 FORMAT( 4X, 'FUNCTION DGAUSS8 ... TOO HIGH ACCURACY REQUIRED') 
      END      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

***********************************************************************
**********************************************************
*     routines from packlib                              *
*     a) KERSET                                          *
*     b) ABEND                                           *
**********************************************************

*******************************************************************
*
          SUBROUTINE KERSET(ERCODE,LGFILE,LIMITM,LIMITR)
                    PARAMETER(KOUNTE  =  27)         
          CHARACTER*6         ERCODE,   CODE(KOUNTE)  
          LOGICAL             MFLAG,    RFLAG         
          INTEGER             KNTM(KOUNTE),       KNTR(KOUNTE) 
          DATA      LOGF      /  0  /  
          DATA      CODE(1), KNTM(1), KNTR(1)  / 'C204.1', 100, 100 / 
          DATA      CODE(2), KNTM(2), KNTR(2)  / 'C204.2', 100, 100 / 
          DATA      CODE(3), KNTM(3), KNTR(3)  / 'C204.3', 100, 100 / 
          DATA      CODE(4), KNTM(4), KNTR(4)  / 'C205.1', 100, 100 / 
          DATA      CODE(5), KNTM(5), KNTR(5)  / 'C205.2', 100, 100 / 
          DATA      CODE(6), KNTM(6), KNTR(6)  / 'C305.1', 100, 100 / 
          DATA      CODE(7), KNTM(7), KNTR(7)  / 'C308.1', 100, 100 / 
          DATA      CODE(8), KNTM(8), KNTR(8)  / 'C312.1', 100, 100 / 
          DATA      CODE(9), KNTM(9), KNTR(9)  / 'C313.1', 100, 100 / 
          DATA      CODE(10),KNTM(10),KNTR(10) / 'C336.1', 100, 100 / 
          DATA      CODE(11),KNTM(11),KNTR(11) / 'C337.1', 100, 100 / 
          DATA      CODE(12),KNTM(12),KNTR(12) / 'C341.1', 100, 100 / 
          DATA      CODE(13),KNTM(13),KNTR(13) / 'D103.1', 100, 100 / 
          DATA      CODE(14),KNTM(14),KNTR(14) / 'D106.1', 100, 100 / 
          DATA      CODE(15),KNTM(15),KNTR(15) / 'D209.1', 100, 100 / 
          DATA      CODE(16),KNTM(16),KNTR(16) / 'D509.1', 100, 100 / 
          DATA      CODE(17),KNTM(17),KNTR(17) / 'E100.1', 100, 100 / 
          DATA      CODE(18),KNTM(18),KNTR(18) / 'E104.1', 100, 100 / 
          DATA      CODE(19),KNTM(19),KNTR(19) / 'E105.1', 100, 100 / 
          DATA      CODE(20),KNTM(20),KNTR(20) / 'E208.1', 100, 100 / 
          DATA      CODE(21),KNTM(21),KNTR(21) / 'E208.2', 100, 100 / 
          DATA      CODE(22),KNTM(22),KNTR(22) / 'F010.1', 100,   0 / 
          DATA      CODE(23),KNTM(23),KNTR(23) / 'F011.1', 100,   0 / 
          DATA      CODE(24),KNTM(24),KNTR(24) / 'F012.1', 100,   0 / 
          DATA      CODE(25),KNTM(25),KNTR(25) / 'F406.1', 100,   0 / 
          DATA      CODE(26),KNTM(26),KNTR(26) / 'G100.1', 100, 100 / 
          DATA      CODE(27),KNTM(27),KNTR(27) / 'G100.2', 100, 100 /  
          LOGF  =  LGFILE      
          IF(ERCODE .EQ. ' ')  THEN  
             L  =  0  
          ELSE  
             DO 10  L = 1, 6   
                IF(ERCODE(1:L) .EQ. ERCODE)  GOTO 12 
  10            CONTINUE  
  12         CONTINUE     
          ENDIF      
          DO 14     I  =  1, KOUNTE   
             IF(L .EQ. 0)  GOTO 13   
             IF(CODE(I)(1:L) .NE. ERCODE(1:L))  GOTO 14
  13         KNTM(I)  =  LIMITM     
             KNTR(I)  =  LIMITR   
  14         CONTINUE         
          RETURN              
          ENTRY KERMTR(ERCODE,LOG,MFLAG,RFLAG)
          LOG  =  LOGF       
          DO 20     I  =  1, KOUNTE 
             IF(ERCODE .EQ. CODE(I))  GOTO 21 
  20         CONTINUE  
          WRITE(*,1000)  ERCODE  
          CALL ABEND   
          RETURN                                               
  21      RFLAG  =  KNTR(I) .GE. 1  
          IF(RFLAG  .AND.  (KNTR(I) .LT. 100))  KNTR(I)  =  KNTR(I) - 1
          MFLAG  =  KNTM(I) .GE. 1 
          IF(MFLAG  .AND.  (KNTM(I) .LT. 100))  KNTM(I)  =  KNTM(I) - 1
          IF(.NOT. RFLAG)  THEN  
             IF(LOGF .LT. 1)  THEN 
                WRITE(*,1001)  CODE(I)
             ELSE 
                WRITE(LOGF,1001)  CODE(I)  
             ENDIF  
          ENDIF  
          IF(MFLAG .AND. RFLAG)  THEN   
             IF(LOGF .LT. 1)  THEN  
c                WRITE(*,1002)  CODE(I)  
             ELSE  
                WRITE(LOGF,1002)  CODE(I)  
             ENDIF   
          ENDIF     
          RETURN  
1000      FORMAT(' KERNLIB LIBRARY ERROR. ' /    
     +           ' ERROR CODE ',A6,' NOT RECOGNIZED BY KERMTR', 
     +           ' ERROR MONITOR. RUN ABORTED.')  
1001      FORMAT(/' ***** RUN TERMINATED BY CERN LIBRARY ERROR ',  
     +           'CONDITION ',A6)   
1002      FORMAT(/' ***** CERN LIBRARY ERROR CONDITION ',A6) 
          END    


*******************************************************************
*
      SUBROUTINE ABEND 
C                      
C CERN PROGLIB# Z035    ABEND           .VERSION KERNVAX  1.10  811126 
C                      
      STOP '*** ABEND ***' 
      END           

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c modified bessel function first kind K0(x) for positive real x
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         FUNCTION bessk0(x)
         DOUBLE PRECISION bessk0,x
         DOUBLE PRECISION bessi0
         DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,
     &                    q1,q2,q3,q4,q5,q6,q7,y
C
         SAVE p1,p2,p3,p4,p5,p6,p7,
     &        q1,q2,q3,q4,q5,q6,q7
C
         DATA p1,p2,p3,p4,p5,p6,p7/
     &   -0.57721566d0,0.42278420d0,0.23069756d0,
     &   0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
         DATA q1,q2,q3,q4,q5,q6,q7/
     &   1.25331414d0,-0.7832358d-1,0.2189568d-1, 
     &   -0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/
C
         if(x.le.2.0) then
         y = x*x/4.0
         bessk0 = (-log(x/2.0)*bessi0(x)) + 
     &   (p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
         else
         y = 2.0/x
         bessk0 = (exp(-x)/sqrt(x))*
     &   (q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))))
         endif
C
         return
         END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c modified bessel function I0(x) for any postive x, this is used in K0(x)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
         FUNCTION bessi0(x)
         DOUBLE PRECISION bessi0,x
         REAL ax
         DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,
     &                    q1,q2,q3,q4,q5,q6,q7,q8,q9,y
C
         SAVE p1,p2,p3,p4,p5,p6,p7,
     &        q1,q2,q3,q4,q5,q6,q7,q8,q9
C
         DATA p1,p2,p3,p4,p5,p6,p7/
     &   1.0d0,3.5156229d0,3.0899424d0,1.2067492d0,0.2659732d0,
     &   0.360768d-1,0.45813d-2/
         DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/
     &   0.39894228d0,0.1328592d-1,0.225319d-2,-0.157565d-2,
     &   0.916281d-2,-0.2057706d-1,0.2635537d-1,-0.1647633d-1,
     &   0.392377d-2/
C
         if(abs(x).le.3.75) then
         y = (x/3.75)**2
         bessi0 = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
         else
         ax = abs(x)
         y = 3.75/ax
         bessi0 = (exp(ax)/sqrt(ax))*
     &   (q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
         endif
C
         return
         END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c modified bessel function I1(x) for any postive x, this is used in K1(x)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION bessk1(x)
      DOUBLE PRECISION bessk1,x
CU    USES bessi1
      DOUBLE PRECISION bessi1
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,
     *-0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1,
     *
     *0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/
      if (x.le.2.0d0) then
        y=x*x/4.0d0
        bessk1=(log(x/2.0d0)*bessi1(x))+(1.0d0/x)*(p1+y*(p2+y*(p3+y*(p4
     *+y*
     *(p5+y*(p6+y*p7))))))
      else
        y=2.0d0/x
        bessk1=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *q7))))))
      endif
      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      FUNCTION bessi1(x)
      DOUBLE PRECISION bessi1,x
      DOUBLE PRECISION ax
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
     *
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,
     *0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,
     *-0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,-0.2895312d-1,
     *
     *0.1787654d-1,-0.420059d-2/
      if (abs(x).lt.3.75d0) then
        y=(x/3.75d0)**2
        bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
        ax=abs(x)
        y=3.75d0/ax
        bessi1=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *(q7+y*(q8+y*q9))))))))
        if(x.lt.0.d0)bessi1=-bessi1
      endif
      return
      END



cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      function dbesJ0(x)

      implicit none
       
      double precision dbesJ0,x,j,y,jp,yp
      external dbessjy

      call dbessjy(x,0.0d0,j,y,jp,yp)
      dbesJ0 = j

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cII.  Mask to compute J1(x) with same sintax as cernlib
      function dbesJ1(x)

      implicit none
      
      double precision dbesJ1,x,j,y,jp,yp

      call dbessjy(x,1d0,j,y,jp,yp)
      dbesJ1 = j

      return
      end

cccccccccccccccccccccccccccccccccccccccccccc

c II.  Mask to compute J2(x) with same sintax as cernlib

      function dbesJ2(x)

      implicit none
       
      double precision dbesJ2,x,j,y,jp,yp
      external dbessjy

      call dbessjy(x,2.0d0,j,y,jp,yp)
      dbesJ2 = j

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c II.  Mask to compute J2(x) with same sintax as cernlib

      function dbesJ3(x)

      implicit none
       
      double precision dbesJ3,x,j,y,jp,yp
      external dbessjy

      call dbessjy(x,3.0d0,j,y,jp,yp)
      dbesJ3 = j

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

*******************************************************************
cc **   *** double precision bessel functions

      SUBROUTINE dbessjy(x,xnu,rj,ry,rjp,ryp)
      INTEGER MAXIT
      DOUBLE PRECISION rj,rjp,ry,ryp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.e-16,FPMIN=1.e-30,MAXIT=100000,XMIN=2.,
     *PI=3.141592653589793d0)
CU    USES beschb
      INTEGER i,isign,l,nl
      DOUBLE PRECISION a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,
     *f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q,
     *r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1,
     *temp,w,x2,xi,xi2,xmu,xmu2
      if(x.le.0..or.xnu.lt.0.) pause 'bad arguments in bessjy'
      if(x.lt.XMIN)then
        nl=int(xnu+.5d0)
      else
        nl=max(0,int(xnu-x+1.5d0))
      endif
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      w=xi2/PI
      isign=1
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
        b=b+xi2
        d=b-d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b-1.d0/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        del=c*d
        h=del*h
        if(d.lt.0.d0)isign=-isign
        if(abs(del-1.d0).lt.EPS)goto 1
11    continue
      pause 'x too large in bessjy; try asymptotic expansion'
1     continue
      rjl=isign*FPMIN
      rjpl=h*rjl
      rjl1=rjl
      rjp1=rjpl
      fact=xnu*xi
      do 12 l=nl,1,-1
        rjtemp=fact*rjl+rjpl
        fact=fact-xi
        rjpl=fact*rjtemp-rjl
        rjl=rjtemp
12    continue
      if(rjl.eq.0.d0)rjl=EPS
      f=rjpl/rjl
      if(x.lt.XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu).lt.EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call dbeschb(xmu,gam1,gam2,gampl,gammi)
        ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d)
        e=exp(e)
        p=e/(gampl*PI)
        q=1.d0/(e*PI*gammi)
        pimu2=0.5d0*pimu
        if(abs(pimu2).lt.EPS)then
          fact3=1.d0
        else
          fact3=sin(pimu2)/pimu2
        endif
        r=PI*pimu2*fact3*fact3
        c=1.d0
        d=-x2*x2
        sum=ff+r*q
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*(ff+r*q)
          sum=sum+del
          del1=c*p-i*del
          sum1=sum1+del1
          if(abs(del).lt.(1.d0+abs(sum))*EPS)goto 2
13      continue
        pause 'bessy series failed to converge'
2       continue
        rymu=-sum
        ry1=-sum1*xi2
        rymup=xmu*xi*rymu-ry1
        rjmu=w/(rymup-f*rymu)
      else
        a=.25d0-xmu2
        p=-.5d0*xi
        q=1.d0
        br=2.d0*x
        bi=2.d0
        fact=a*xi/(p*p+q*q)
        cr=br+q*fact
        ci=bi+p*fact
        den=br*br+bi*bi
        dr=br/den
        di=-bi/den
        dlr=cr*dr-ci*di
        dli=cr*di+ci*dr
        temp=p*dlr-q*dli
        q=p*dli+q*dlr
        p=temp
        do 14 i=2,MAXIT
          a=a+2*(i-1)
          bi=bi+2.d0
          dr=a*dr+br
          di=a*di+bi
          if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN
          fact=a/(cr*cr+ci*ci)
          cr=br+cr*fact
          ci=bi-ci*fact
          if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN
          den=dr*dr+di*di
          dr=dr/den
          di=-di/den
          dlr=cr*dr-ci*di
          dli=cr*di+ci*dr
          temp=p*dlr-q*dli
          q=p*dli+q*dlr
          p=temp
          if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3
14      continue
        pause 'cf2 failed in bessjy'
3       continue
        gam=(p-f)/q
        rjmu=sqrt(w/((p-f)*gam+q))
        rjmu=sign(rjmu,rjl)
        rymu=rjmu*gam
        rymup=rymu*(p+q/gam)
        ry1=xmu*xi*rymu-rymup
      endif
      fact=rjmu/rjl
      rj=rjl1*fact
      rjp=rjp1*fact
      do 15 i=1,nl
        rytemp=(xmu+i)*xi2*ry1-rymu
        rymu=ry1
        ry1=rytemp
15    continue
      ry=rymu
      ryp=xnu*xi*rymu-ry1
      return
      END


      SUBROUTINE dbeschb(x,gam1,gam2,gampl,gammi)
      INTEGER NUSE1,NUSE2
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x
      PARAMETER (NUSE1=5,NUSE2=5)
CU    USES dchebev
      DOUBLE PRECISION xx,c1(7),c2(8), dchebev
      SAVE c1,c2
      DATA c1/-1.142022680371172d0,6.516511267076d-3,3.08709017308d-4,
     *-3.470626964d-6,6.943764d-9,3.6780d-11,-1.36d-13/
      DATA c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3,
     *-4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/
      xx=8.d0*x*x-1.d0
      gam1=dchebev(-1.d0,1.d0,c1,NUSE1,xx)
      gam2=dchebev(-1.d0,1.d0,c2,NUSE2,xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
      END

      FUNCTION dchebev(a,b,c,m,x)
      INTEGER m
      DOUBLE PRECISION dchebev,a,b,x,c(m)
      INTEGER j
      DOUBLE PRECISION d,dd,sv,y,y2
      if ((x-a)*(x-b).gt.0.) pause 'x not in range in chebev'
      d=0.
      dd=0.
      y=(2.*x-a-b)/(b-a)
      y2=2.*y
      do 11 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv   
11    continue
      dchebev=y*d-dd+0.5*c(1)
      return
      END