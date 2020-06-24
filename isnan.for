
      PROGRAM MAIN
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        I=0
        X=1
        Z=X/I
        M=ISNAN(Z)
        PRINT *,M
        IF (ISNAN(Z)) THEN
          PRINT *,Z,"OK,ITS NAN"
        ELSE 
          PRINT *,Z,"ITS A NUMBER"
        END IF
      END PROGRAM

